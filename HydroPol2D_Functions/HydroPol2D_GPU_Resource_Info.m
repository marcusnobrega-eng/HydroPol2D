%% HydroPol2D — GPU Feasibility + dt/steps + Runtime (ANALYTICAL, Local Inertial)
% This script builds a cfg struct at the top, calls the estimator, and
% defines the estimator as an auxiliary function at the bottom.
%
% ✅ You provide:
%   - GPU memory budget + precision
%   - Grid size (Nx,Ny) OR (Lx,Ly,dx)
%   - Counts of GPU arrays stored (cell floats, qx faces, qy faces, masks, indices)
%   - Temporary arrays as multipliers OR explicit counts
%   - CFL + characteristic U and h + simulation duration
%
% The estimator outputs:
%   - Peak GPU memory (GB) and fit/not-fit
%   - Maximum grid that fits your budget (NxMax, NyMax) and dxMin (if domain size given)
%   - dt from CFL and number of steps
%   - Optional wall-time estimate via an assumed throughput knob

clear; clc;

%% ========================= 1) USER CONFIG =========================
cfg = struct();

% ---------- GPU budget ----------
cfg.gpuBudgetGB     = 128;     % <--- allowed GPU memory you want to use (GB)
cfg.safetyFactor    = 0.85;   % <--- reserve margin (e.g., 0.7–0.9)
cfg.fixedOverheadGB = 0.8;    % <--- CUDA context + fragmentation cushion (GB)
cfg.precision       = 'double';  % 'single' or 'double'

% ---------- Grid (choose ONE option) ----------
% Option A: set Nx,Ny directly
cfg.Nx = 4220;   % e.g., 2000
cfg.Ny = 3569;   % e.g., 1500
cfg.dx_m = 250; % optional if you know it (used for dt)

% Option B: set physical size + dx (Nx,Ny computed automatically)
% cfg.Lx_m = 10e3;    % domain length in x (m)
% cfg.Ly_m = 10e3;    % domain length in y (m)
% cfg.dx_m = 20;      % target resolution (m)

% ---------- Local Inertial dt estimate ----------
cfg.CFL       = 0.5;     % typical 0.3–0.7
cfg.Uchar_mps = 2;     % characteristic velocity (m/s)
cfg.hchar_m   = 5;     % characteristic depth (m)
cfg.Tsim_s    = 365*8600;  % simulation duration (s)
% cfg.dtMin_s   = 0.05;    % optional (s)
% cfg.dtMax_s   = 5.0;     % optional (s)

% ---------- GPU array inventory (PERSISTENT arrays) ----------
% These are COUNTS (not the arrays themselves).
cfg.A = struct();

cfg.A.nCellFloat   = 14;  % # of Nx×Ny floating arrays on GPU (h,wse,z,n,...)
cfg.A.nQxFaces     = 2;   % qx faces: each is (Nx+1)×Ny (you said qx has 2 faces)
cfg.A.nQyFaces     = 2;   % qy faces: each is Nx×(Ny+1) (qy has 2 faces)
cfg.A.nCellLogical = 4;   % masks (wet/dry, nan, boundary flags...) as logical Nx×Ny
cfg.A.nCellInt32   = 0;   % int32 Nx×Ny (receivers, indexing, etc.)
cfg.A.nCellInt64   = 0;   % int64 Nx×Ny (rare)

% ---------- TEMPORARIES (peak extra arrays alive simultaneously) ----------
% Choose ONE approach:

% Approach 1 (recommended): multipliers (simple + robust)
cfg.T = struct();
cfg.T.multiplierCellFloat = 0.8;  % temporaries ~ 80% of persistent cell floats
cfg.T.multiplierFaces     = 0.5;  % temporaries ~ 50% of persistent faces (qx/qy)

% Approach 2 (explicit): uncomment and comment multipliers above
% cfg.T = struct();
% cfg.T.nCellFloat = 12;
% cfg.T.nQxFaces   = 1;
% cfg.T.nQyFaces   = 1;

% ---------- Optional runtime estimate knob ----------
% If you don’t know this yet, leave NaN and you’ll only get dt & steps.
cfg.cellUpdatesPerSecond = NaN;   % e.g., 2e8 (tune later)

% ---------- Optional output cost model (only affects runtime estimate) ----------
cfg.output = struct();
cfg.output.enabled     = true;
cfg.output.writeEvery_s = 300;  % every 5 minutes simulated time
cfg.output.cost_s       = 0.30; % assumed wall seconds per output event

%% ========================= 2) RUN ESTIMATOR =========================
Out = HydroPol2D_GPU_AnalyticalLI(cfg);

%% ========================= 3) OPTIONAL: DISPLAY SUMMARY TABLES =========================
disp(" ");
disp("---- Memory Breakdown (GB) ----");
disp(struct2table(Out.memory, 'AsArray', true))

disp(" ");
disp("---- Capacity ----");
disp(struct2table(Out.capacity, 'AsArray', true))

disp(" ");
disp("---- Timing ----");
disp(struct2table(Out.timing, 'AsArray', true))

%% =====================================================================
%% ======================= AUXILIARY FUNCTION ===========================
%% =====================================================================
function Out = HydroPol2D_GPU_AnalyticalLI(cfg)
%HydroPol2D_GPU_AnalyticalLI
% Analytical estimator for HydroPol2D Local Inertial GPU feasibility + dt/steps + runtime.

%% Defaults
mustHave(cfg, {'gpuBudgetGB','safetyFactor','fixedOverheadGB','precision','CFL','Uchar_mps','hchar_m','Tsim_s','A','T'});
if ~isfield(cfg,'dtMin_s'); cfg.dtMin_s = 0; end
if ~isfield(cfg,'dtMax_s'); cfg.dtMax_s = inf; end
if ~isfield(cfg,'cellUpdatesPerSecond'); cfg.cellUpdatesPerSecond = NaN; end
if ~isfield(cfg,'output'); cfg.output.enabled = false; end
if ~isfield(cfg.output,'enabled'); cfg.output.enabled = false; end
if ~isfield(cfg.output,'writeEvery_s'); cfg.output.writeEvery_s = inf; end
if ~isfield(cfg.output,'cost_s'); cfg.output.cost_s = 0; end

bFloat   = bytesPerFloat(cfg.precision);
bLogical = 1;
bInt32   = 4;
bInt64   = 8;

%% Grid
if isfield(cfg,'Nx') && isfield(cfg,'Ny') && ~isempty(cfg.Nx) && ~isempty(cfg.Ny)
    Nx = cfg.Nx; Ny = cfg.Ny;
    if isfield(cfg,'dx_m') && ~isempty(cfg.dx_m)
        dx = cfg.dx_m;
    elseif isfield(cfg,'Lx_m') && isfield(cfg,'Ly_m') && ~isempty(cfg.Lx_m) && ~isempty(cfg.Ly_m)
        dx = sqrt((cfg.Lx_m*cfg.Ly_m)/(double(Nx)*double(Ny)));
    else
        dx = NaN;
    end
elseif isfield(cfg,'Lx_m') && isfield(cfg,'Ly_m') && isfield(cfg,'dx_m') && ~isempty(cfg.dx_m)
    dx = cfg.dx_m;
    Nx = ceil(cfg.Lx_m/dx);
    Ny = ceil(cfg.Ly_m/dx);
else
    error("Provide either (Nx,Ny) or (Lx_m,Ly_m,dx_m).");
end

%% Temporaries: explicit counts or multipliers
T = cfg.T;
if ~isfield(T,'nCellFloat')
    mustHave(T, {'multiplierCellFloat','multiplierFaces'});
    T.nCellFloat = ceil(T.multiplierCellFloat * cfg.A.nCellFloat);
    T.nQxFaces   = ceil(T.multiplierFaces    * cfg.A.nQxFaces);
    T.nQyFaces   = ceil(T.multiplierFaces    * cfg.A.nQyFaces);
else
    mustHave(T, {'nCellFloat','nQxFaces','nQyFaces'});
end

%% Memory for this grid
nCell  = double(Nx)*double(Ny);
nQx    = double(Nx+1)*double(Ny);    % (Nx+1)×Ny
nQy    = double(Nx)*double(Ny+1);    % Nx×(Ny+1)

persistentBytes = ...
    cfg.A.nCellFloat   * nCell * bFloat + ...
    cfg.A.nQxFaces     * nQx   * bFloat + ...
    cfg.A.nQyFaces     * nQy   * bFloat + ...
    cfg.A.nCellLogical * nCell * bLogical + ...
    cfg.A.nCellInt32   * nCell * bInt32 + ...
    cfg.A.nCellInt64   * nCell * bInt64;

tempBytes = ...
    T.nCellFloat * nCell * bFloat + ...
    T.nQxFaces   * nQx   * bFloat + ...
    T.nQyFaces   * nQy   * bFloat;

peakBytes  = persistentBytes + tempBytes;
fixedBytes = cfg.fixedOverheadGB * 1024^3;

budgetBytes = cfg.gpuBudgetGB * 1024^3;
safeBytes   = budgetBytes * cfg.safetyFactor;

fits = (peakBytes + fixedBytes) <= safeBytes;

%% Invert capacity (max cells for budget) using per-cell approximation
fx = 1 + 1/max(double(Nx),1);
fy = 1 + 1/max(double(Ny),1);

bytesPerCell = ...
    (cfg.A.nCellFloat + T.nCellFloat) * bFloat + ...
    (cfg.A.nQxFaces   + T.nQxFaces)   * bFloat * fx + ...
    (cfg.A.nQyFaces   + T.nQyFaces)   * bFloat * fy + ...
    (cfg.A.nCellLogical)             * bLogical + ...
    (cfg.A.nCellInt32)               * bInt32 + ...
    (cfg.A.nCellInt64)               * bInt64;

avail = max(safeBytes - fixedBytes, 0);
NcellsMax = floor(avail / max(bytesPerCell,1));

% aspect ratio to preserve
if isfield(cfg,'Lx_m') && isfield(cfg,'Ly_m') && ~isempty(cfg.Lx_m) && ~isempty(cfg.Ly_m)
    r = cfg.Lx_m / cfg.Ly_m; % Nx/Ny
else
    r = double(Nx)/double(Ny);
end
NxMax = max(1, floor(sqrt(double(NcellsMax) * r)));
NyMax = max(1, floor(double(NcellsMax) / double(NxMax)));

if isfield(cfg,'Lx_m') && isfield(cfg,'Ly_m') && ~isempty(cfg.Lx_m) && ~isempty(cfg.Ly_m)
    dxMin = sqrt((cfg.Lx_m*cfg.Ly_m) / (double(NxMax)*double(NyMax)));
else
    dxMin = NaN;
end

%% CFL dt estimate (Local Inertial)
g = 9.81;
c = sqrt(g * max(cfg.hchar_m,0));
speed = abs(cfg.Uchar_mps) + c;

dtCFL = cfg.CFL * (dx / max(speed, eps));
dtEst = min(max(dtCFL, cfg.dtMin_s), cfg.dtMax_s);
Nsteps = ceil(cfg.Tsim_s / dtEst);

%% Runtime estimate (pure analytical throughput knob)
if isnan(cfg.cellUpdatesPerSecond) || cfg.cellUpdatesPerSecond <= 0
    runtimeCompute_s = NaN;
else
    runtimeCompute_s = (double(Nsteps) * nCell) / cfg.cellUpdatesPerSecond;
end

% output time (optional)
if cfg.output.enabled && isfinite(cfg.output.writeEvery_s) && cfg.output.writeEvery_s > 0
    Noutputs = floor(cfg.Tsim_s / cfg.output.writeEvery_s) + 1;
else
    Noutputs = 0;
end
runtimeOutput_s = Noutputs * cfg.output.cost_s;

if isnan(runtimeCompute_s)
    runtimeTotal_s = NaN;
else
    runtimeTotal_s = runtimeCompute_s + runtimeOutput_s;
end

%% Pack outputs
Out = struct();
Out.grid = struct('Nx',Nx,'Ny',Ny,'dx_m',dx,'Ncells',nCell);

Out.memory = struct( ...
    'persistentGB', persistentBytes/1024^3, ...
    'temporariesGB', tempBytes/1024^3, ...
    'peakGB', peakBytes/1024^3, ...
    'fixedOverheadGB', cfg.fixedOverheadGB, ...
    'peakPlusFixedGB', (peakBytes+fixedBytes)/1024^3, ...
    'budgetGB', cfg.gpuBudgetGB, ...
    'safeGB', safeBytes/1024^3, ...
    'fits', fits);

Out.capacity = struct('NcellsMax',NcellsMax,'NxMax',NxMax,'NyMax',NyMax,'dxMin_m',dxMin, ...
                      'bytesPerCellApprox', bytesPerCell);

Out.timing = struct('dtCFL_s',dtCFL,'dtEst_s',dtEst,'Nsteps',Nsteps, ...
                    'runtimeCompute_s',runtimeCompute_s, ...
                    'runtimeOutput_s',runtimeOutput_s, ...
                    'runtimeTotal_s',runtimeTotal_s, ...
                    'Noutputs',Noutputs);

%% Console summary
disp("=== HydroPol2D Local Inertial GPU Analytical Estimator ===")
fprintf("Grid: Nx=%d Ny=%d (Ncells=%.3g) dx=%.4g m\n", Nx, Ny, nCell, dx);
fprintf("Precision: %s | Budget: %.2f GB | Safe: %.2f GB | Fixed: %.2f GB\n", ...
    cfg.precision, cfg.gpuBudgetGB, safeBytes/1024^3, cfg.fixedOverheadGB);
fprintf("GPU peak: %.2f GB (arrays) + fixed => %.2f GB | fits=%d\n", ...
    peakBytes/1024^3, (peakBytes+fixedBytes)/1024^3, fits);

fprintf("Capacity ~ Ncells<=%.3g => Nx~%d Ny~%d", NcellsMax, NxMax, NyMax);
if ~isnan(dxMin); fprintf(" | dxMin~%.4g m\n", dxMin); else; fprintf("\n"); end

fprintf("dtCFL=%.4g s | dt=%.4g s | steps=%d for Tsim=%.3g s\n", dtCFL, dtEst, Nsteps, cfg.Tsim_s);

if ~isnan(runtimeTotal_s)
    fprintf("Runtime compute ~ %.2f h (throughput=%.3g cell-updates/s)\n", runtimeCompute_s/3600, cfg.cellUpdatesPerSecond);
    if runtimeOutput_s > 0
        fprintf("Runtime output  ~ %.2f h (Noutputs=%d, %.3gs each)\n", runtimeOutput_s/3600, Noutputs, cfg.output.cost_s);
    end
    fprintf("Runtime total   ~ %.2f h\n", runtimeTotal_s/3600);
else
    fprintf("Runtime: set cfg.cellUpdatesPerSecond to estimate wall time.\n");
end
end

function b = bytesPerFloat(precision)
switch lower(string(precision))
    case "single", b = 4;
    case "double", b = 8;
    otherwise, error("precision must be 'single' or 'double'");
end
end

function mustHave(S, names)
for i = 1:numel(names)
    if ~isfield(S, names{i})
        error("Missing required cfg.%s", names{i});
    end
end
end