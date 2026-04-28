function [SubgridTables, invert_el] = Subgrid_Properties_Lookup(DEM_raster, Reference_raster, coarse_res, varargin)
%--------------------------------------------------------------------------
% Subgrid_Properties_Lookup (PAPER-STYLE LOOKUP TABLE VERSION)
%
% PURPOSE
%   Build subgrid hydraulic lookup tables consistent with the paper:
%
%     1) Cell-centered storage functions:
%           eta  -> wetted volume
%           eta  -> wetted area   (optional but useful)
%
%     2) Shared-face conveyance functions:
%           eta  -> wetted face area
%
%   where eta is the WATER SURFACE ELEVATION (absolute elevation, [m]).
%
% MAIN IDEA
%   - Coarse-cell storage is represented by the fine DEM contained inside
%     each coarse cell.
%   - Shared-face conveyance is represented by a COMMON INTERFACE built from
%     BOTH neighboring coarse cells:
%           z_face = max(edge_current, edge_neighbor)
%   - If one side is NaN, that fine face segment is closed.
%   - If all segments are invalid, the whole face is no-flow.
%
% WHY THIS MATCHES THE PAPER
%   - Runtime local-inertial solver should use:
%         * cell volume <-> eta
%         * face area    <- eta
%   - No single coarse-cell bed elevation is used as the governing storage
%     geometry. Instead, the full within-cell topography defines V(eta).
%
% INPUTS
%   DEM_raster.Z          : fine DEM [m]
%   DEM_raster.cellsize   : fine resolution [m]
%   Reference_raster.Z    : coarse grid shape (dimensions used)
%   coarse_res            : coarse cell size [m]
%
% OPTIONAL NAME-VALUE
%   'dz'       : vertical lookup spacing [m], default = 0.25
%   'maxDepth' : maximum tabulated depth ABOVE local invert [m], default = 10.0
%   'tol_rel'  : tolerance for coarse_res/cellsize non-integer ratio, default = 5e-3
%   'verbose'  : true/false progress printing, default = true
%
% OUTPUTS
%   SubgridTables:
%
%       % Global depth axis (relative, for convenience only)
%       .depth_axis             [1 x nz]
%       .dz                     scalar
%       .maxDepth               scalar
%       .is_uniform             logical = true
%
%       % ---------------- CELL-CENTERED TABLES ----------------
%       .eta_cell               [nyc x nxc x nz]   absolute water surface elevations [m]
%       .area_cell              [nyc x nxc x nz]   wetted plan area [m^2]
%       .volume_cell            [nyc x nxc x nz]   wetted volume [m^3]
%       .invert_el              [nyc x nxc]        minimum fine elevation in cell [m]
%       .eta_top_cell           [nyc x nxc]        highest tabulated eta [m]
%       .Vmax_cell              [nyc x nxc]        highest tabulated volume [m^3]
%
%       % ---------------- SHARED X-FACE TABLES ----------------
%       .eta_x                  [nyc x (nxc-1) x nz]   absolute water surface elevations [m]
%       .area_x                 [nyc x (nxc-1) x nz]   wetted shared-face area [m^2]
%       .invert_x               [nyc x (nxc-1)]        minimum common-face elevation [m]
%       .eta_top_x              [nyc x (nxc-1)]        highest tabulated eta [m]
%
%       % ---------------- SHARED Y-FACE TABLES ----------------
%       .eta_y                  [(nyc-1) x nxc x nz]   absolute water surface elevations [m]
%       .area_y                 [(nyc-1) x nxc x nz]   wetted shared-face area [m^2]
%       .invert_y               [(nyc-1) x nxc]        minimum common-face elevation [m]
%       .eta_top_y              [(nyc-1) x nxc]        highest tabulated eta [m]
%
% NOTES
%   - The paper fundamentally needs:
%         * cell volume vs water surface elevation
%         * shared-face area vs water surface elevation
%   - The inverse relation eta(V) is obtained later by interpolation
%     through the stored monotonic tables volume_cell <-> eta_cell.
%   - This function assumes the fine DEM and coarse grid are aligned.
%--------------------------------------------------------------------------

% ---------------- Parse options ----------------
p = inputParser;
p.addParameter('dz', 0.25, @(x)isscalar(x) && x > 0);
p.addParameter('maxDepth', 10.0, @(x)isscalar(x) && x > 0);
p.addParameter('tol_rel', 5e-3, @(x)isscalar(x) && x > 0);
p.addParameter('verbose', true, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});

dz       = p.Results.dz;
maxDepth = p.Results.maxDepth;
tol_rel  = p.Results.tol_rel;
verbose  = logical(p.Results.verbose);

DEM      = DEM_raster.Z;
cellsize = DEM_raster.cellsize;

% ---------------- Ratio check ----------------
r = coarse_res / cellsize;
r_round = round(r);
rel_err = abs(r - r_round) / max(r, eps);

if rel_err > tol_rel
    warning(['coarse_res/cellsize is not ~integer (r = %.6f). ' ...
             'Shared-face mapping may drift. Strongly recommend resampling DEM ' ...
             'so coarse_res is an exact multiple of cellsize.'], r);
end

[nrows, ncols] = size(DEM);
nrows_coarse = size(Reference_raster.Z,1);
ncols_coarse = size(Reference_raster.Z,2);

% ---------------- Uniform relative depth axis ----------------
% This is only a GENERIC relative axis. Absolute eta axes differ by cell/face
% because each cell/face has its own invert.
depth_axis = 0:dz:maxDepth;
nz = numel(depth_axis);

% ---------------- Preallocate outputs ----------------
SubgridTables = struct();

SubgridTables.depth_axis = depth_axis;
SubgridTables.dz         = dz;
SubgridTables.maxDepth   = maxDepth;
SubgridTables.is_uniform = true;

% ---------------- CELL TABLES ----------------
SubgridTables.eta_cell    = NaN(nrows_coarse, ncols_coarse, nz, 'like', DEM);
SubgridTables.area_cell   = NaN(nrows_coarse, ncols_coarse, nz, 'like', DEM);
SubgridTables.volume_cell = NaN(nrows_coarse, ncols_coarse, nz, 'like', DEM);

invert_el                 = NaN(nrows_coarse, ncols_coarse, 'like', DEM);
SubgridTables.invert_el   = NaN(nrows_coarse, ncols_coarse, 'like', DEM);
SubgridTables.eta_top_cell= NaN(nrows_coarse, ncols_coarse, 'like', DEM);
SubgridTables.Vmax_cell   = NaN(nrows_coarse, ncols_coarse, 'like', DEM);

% ---------------- SHARED X-FACE TABLES ----------------
nxface = max(ncols_coarse-1,1);
SubgridTables.eta_x     = NaN(nrows_coarse, nxface, nz, 'like', DEM);
SubgridTables.area_x    = NaN(nrows_coarse, nxface, nz, 'like', DEM);
SubgridTables.invert_x  = NaN(nrows_coarse, nxface, 'like', DEM);
SubgridTables.eta_top_x = NaN(nrows_coarse, nxface, 'like', DEM);

% ---------------- SHARED Y-FACE TABLES ----------------
nyface = max(nrows_coarse-1,1);
SubgridTables.eta_y     = NaN(nyface, ncols_coarse, nz, 'like', DEM);
SubgridTables.area_y    = NaN(nyface, ncols_coarse, nz, 'like', DEM);
SubgridTables.invert_y  = NaN(nyface, ncols_coarse, 'like', DEM);
SubgridTables.eta_top_y = NaN(nyface, ncols_coarse, 'like', DEM);

% ---------------- Index mapping ----------------
row_idx_start = floor((0:nrows_coarse-1) * coarse_res / cellsize) + 1;
row_idx_end   = floor((1:nrows_coarse)   * coarse_res / cellsize);
col_idx_start = floor((0:ncols_coarse-1) * coarse_res / cellsize) + 1;
col_idx_end   = floor((1:ncols_coarse)   * coarse_res / cellsize);

% Cache coarse-cell DEM patches for face stage
cellPatch = cell(nrows_coarse, ncols_coarse);

% ======================================================================
% 1) CELL-CENTERED STORAGE TABLES: eta -> area, volume
% ======================================================================
total_cells = nrows_coarse * ncols_coarse;
cell_count = 0;

for rowc = 1:nrows_coarse
    for colc = 1:ncols_coarse

        row_idx = row_idx_start(rowc):min(nrows, row_idx_end(rowc));
        col_idx = col_idx_start(colc):min(ncols, col_idx_end(colc));

        sub_DEM = DEM(row_idx, col_idx);
        cellPatch{rowc,colc} = sub_DEM;

        if isempty(sub_DEM) || all(isnan(sub_DEM(:)))
            cell_count = cell_count + 1;
            continue;
        end

        sub_vec = sub_DEM(:);
        valid = ~isnan(sub_vec);

        zmin = min(sub_vec(valid));
        invert_el(rowc,colc)               = zmin;
        SubgridTables.invert_el(rowc,colc) = zmin;

        eta_vec     = zmin + depth_axis;
        area_vec    = zeros(nz,1,'like',DEM);
        volume_vec  = zeros(nz,1,'like',DEM);

        for k = 1:nz
            eta = eta_vec(k);

            d = eta - sub_vec;
            d(~valid) = 0;
            d(d < 0) = 0;

            wet = d > 0;

            area_vec(k)   = sum(wet) * cellsize^2;
            volume_vec(k) = sum(d)   * cellsize^2;
        end

        SubgridTables.eta_cell(rowc,colc,:)    = reshape(eta_vec,    1,1,nz);
        SubgridTables.area_cell(rowc,colc,:)   = reshape(area_vec,   1,1,nz);
        SubgridTables.volume_cell(rowc,colc,:) = reshape(volume_vec, 1,1,nz);

        SubgridTables.eta_top_cell(rowc,colc)  = eta_vec(end);
        SubgridTables.Vmax_cell(rowc,colc)     = volume_vec(end);

        cell_count = cell_count + 1;
        if verbose && (mod(cell_count,100)==0 || cell_count==total_cells)
            fprintf('Subgrid cell preprocessing: %6.2f%% (%d of %d)\n', ...
                100*cell_count/total_cells, cell_count, total_cells);
        end
    end
end

% ======================================================================
% 2) SHARED X-FACE TABLES: eta -> wetted face area
%    Face between cell(rowc,colc) and cell(rowc,colc+1)
% ======================================================================
if ncols_coarse > 1
    total_xfaces = nrows_coarse * (ncols_coarse - 1);
    face_count = 0;

    for rowc = 1:nrows_coarse
        for colc = 1:(ncols_coarse-1)

            left_patch  = cellPatch{rowc,colc};
            right_patch = cellPatch{rowc,colc+1};

            if isempty(left_patch) || isempty(right_patch) || ...
               all(isnan(left_patch(:))) || all(isnan(right_patch(:)))
                face_count = face_count + 1;
                continue;
            end

            % Shared x-face:
            % east edge of left cell and west edge of right cell
            left_edge  = left_patch(:,end);
            right_edge = right_patch(:,1);

            valid_face = ~isnan(left_edge) & ~isnan(right_edge);

            if ~any(valid_face)
                face_count = face_count + 1;
                continue;
            end

            % Conservative common interface as in the paper logic:
            z_face = max(left_edge(valid_face), right_edge(valid_face));

            zf_min = min(z_face);
            eta_vec = zf_min + depth_axis;

            area_face = zeros(nz,1,'like',DEM);

            for k = 1:nz
                eta = eta_vec(k);

                d = eta - z_face;
                d(d < 0) = 0;

                area_face(k) = sum(d) * cellsize;
            end

            SubgridTables.invert_x(rowc,colc)  = zf_min;
            SubgridTables.eta_top_x(rowc,colc) = eta_vec(end);
            SubgridTables.eta_x(rowc,colc,:)   = reshape(eta_vec,   1,1,nz);
            SubgridTables.area_x(rowc,colc,:)  = reshape(area_face, 1,1,nz);

            face_count = face_count + 1;
            if verbose && (mod(face_count,100)==0 || face_count==total_xfaces)
                fprintf('Subgrid x-face preprocessing: %6.2f%% (%d of %d)\n', ...
                    100*face_count/total_xfaces, face_count, total_xfaces);
            end
        end
    end
end

% ======================================================================
% 3) SHARED Y-FACE TABLES: eta -> wetted face area
%    Face between cell(rowc,colc) and cell(rowc+1,colc)
% ======================================================================
if nrows_coarse > 1
    total_yfaces = (nrows_coarse - 1) * ncols_coarse;
    face_count = 0;

    for rowc = 1:(nrows_coarse-1)
        for colc = 1:ncols_coarse

            south_patch = cellPatch{rowc,colc};
            north_patch = cellPatch{rowc+1,colc};

            if isempty(south_patch) || isempty(north_patch) || ...
               all(isnan(south_patch(:))) || all(isnan(north_patch(:)))
                face_count = face_count + 1;
                continue;
            end

            % Shared y-face:
            % north edge of south cell and south edge of north cell
            south_edge = south_patch(1,:);
            north_edge = north_patch(end,:);

            valid_face = ~isnan(south_edge) & ~isnan(north_edge);

            if ~any(valid_face)
                face_count = face_count + 1;
                continue;
            end

            z_face = max(south_edge(valid_face), north_edge(valid_face));

            zf_min = min(z_face);
            eta_vec = zf_min + depth_axis;

            area_face = zeros(nz,1,'like',DEM);

            for k = 1:nz
                eta = eta_vec(k);

                d = eta - z_face;
                d(d < 0) = 0;

                area_face(k) = sum(d) * cellsize;
            end

            SubgridTables.invert_y(rowc,colc)  = zf_min;
            SubgridTables.eta_top_y(rowc,colc) = eta_vec(end);
            SubgridTables.eta_y(rowc,colc,:)   = reshape(eta_vec,   1,1,nz);
            SubgridTables.area_y(rowc,colc,:)  = reshape(area_face, 1,1,nz);

            face_count = face_count + 1;
            if verbose && (mod(face_count,100)==0 || face_count==total_yfaces)
                fprintf('Subgrid y-face preprocessing: %6.2f%% (%d of %d)\n', ...
                    100*face_count/total_yfaces, face_count, total_yfaces);
            end
        end
    end
end

end