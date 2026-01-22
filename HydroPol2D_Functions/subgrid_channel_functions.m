% Sub-grid approach for 1D Channels
% Developer: Marcus Nobrega, Ph.D.
% Goal: define subgrid approach for dealing with channels
%
% Input:
% h: flood depth matrix [m]
% w: river width matrix [m]
% zc: channel elevation matrix [m]
% dx: grid resolution [m]
% nc: channel Manning's roughness coefficient matrix [sm^-1/3]
% Qc_prev: previous outflow for channel [m3/s]
% g: gravity acceleration [m2/s]
% dt: time-step in [sec]
% 
% Output:
% q: outflow [m2/s] with first entry for x and second entry for y
% similarly:
% C_a = river surface area [m2]

function [q,C_a] = subgrid_channel_functions(h, w, zc,dx, nc, Qc_prev,Qci_prev, g, dt, idx_rivers, Subgrid_Properties, SubgridTables)

    % --------- Channel Calculations ----------- %
    yc = h + zc; % Channel head [m]
    yc(isnan(zc)) = nan;

    % Cells that are not river change to nan
    % yc(~idx_rivers) = inf;
    % h(~idx_rivers) = inf;
    % w(~idx_rivers) = inf;
    % nc(~idx_rivers) = inf;
    
    nan_row = nan*ones(1,size(h,2));
    nan_col = nan*ones(size(h,1),1);

    % Water Surface Slope Calculation
    Sc = slope_function(yc,dx,nan_col,nan_row);
  
    % Effective Channel Flow Depth
    hfcflow = hf_function(yc,zc,nan_col,nan_row);

    % Cells in which channel flow occurs
    mask = hfcflow <= 1e-6;    
        
    % Compute flow width for all cells with their respective h_targets
    h(isnan(h)) = 0;
    hfcflow(mask) = 0;
    
    % Polynomial
    % [Ac(:,:,1), Ac(:,:,2), hydraulic_radius(:,:,1), hydraulic_radius(:,:,2)] = Compute_Subgrid_Properties_Polynomialwise(Subgrid_Properties.A_east_spline, ...
    %                                                Subgrid_Properties.A_north_spline, Subgrid_Properties.Rh_east_spline, Subgrid_Properties.Rh_north_spline, ...
    %                                                yc, zc, Subgrid_Properties.invert_el);

    % Lookup Tables
    [Ac(:,:,1), ~, hydraulic_radius(:,:,1), ~] = Compute_Subgrid_Properties_Lookupwise(SubgridTables, yc, zc, 'east');
    [~, Ac(:,:,2), ~, hydraulic_radius(:,:,2)] = Compute_Subgrid_Properties_Lookupwise(SubgridTables, yc, zc, 'north');

    % Channel Effective Flow Width
    wc_flow = width_function(w,nan_col,nan_row);

    % Channel Flow Area
    % Acflow = wc_flow.*hfcflow;
    Acflow = area_function(Ac,nan_col,nan_row);

    % Channel surface area
    % C_a = dx*w;

    C_a = dx*dx*ones(size(zc));
    C_a(isnan(h)) = nan;

    % Channel Hydraulic Radius
    Rc = Rh_function(hydraulic_radius,nan_col,nan_row);

    % Channel Hydraulic Radius
    % Rc = cellfun(@(func, h) func(h) * ~isempty(func), Subgrid_Properties.hydraulic_radius, num2cell(h));

    % Channel Flow
    Qc_prev(:,:,3:4) = Qci_prev * dx;
    Qc = (Qc_prev - g.*Acflow.*dt.*Sc)./(1 + g*dt.*nc.^2.*abs(Qc_prev)./(Rc.^(4/3).*Acflow)); % m3/s
    Qc(mask) = 0;
    Qc(isnan(Qc)) = 0;
    Qc(isinf(Qc)) = 0;

    % Total Flow (Ortogonal Interfaces)
    Q = Qc(:,:,1:2); Q(:,:,3:4) = zeros(size(h,1),size(h,2),2);
    
    q = Q./dx; %  flow per unit width
    q(:,:,3:4) = [];

end


function [S] = slope_function(M,dx,nan_col,nan_row)
    S(:,:,1) = [(M(:,2:end) - M(:,1:(end-1)))/dx,nan_col]; % East
    S(:,:,2) = [nan_row; M(1:(end-1),:) - M(2:end,:)]/dx; % North
    S(2:end,1:(end-1),3) = (M(1:(end-1),2:end) - M(2:end,1:(end-1)))/(sqrt(2)*dx); % NE
    S(1:(end-1),1:(end-1),4) = (M(2:end,2:(end)) - M(1:(end-1),1:(end-1)))/(sqrt(2)*dx); % SE
end

function [hf] = hf_function(y,z,nan_col,nan_row)
    hf(:,:,1) = [max(y(:,2:end), y(:,1:(end-1))) - max(z(:,2:end), z(:,1:(end-1))), nan_col]; % East
    hf(:,:,2) = [nan_row; max(y(1:(end-1),:),y(2:end,:)) - max(z(1:(end-1),:),z(2:end,:))]; % North
    hf(2:end,1:(end-1),3) = max(y(1:(end-1),2:end), y(2:end,1:(end-1))) - max(z(1:(end-1),2:end), z(2:end,1:(end-1))); % NE
    hf(1:(end-1),1:(end-1),4) = max(y(2:end,2:(end)), y(1:(end-1),1:(end-1))) - max(z(2:end,2:(end)), z(1:(end-1),1:(end-1))); % SE
end

function [wc] = width_function(w,nan_col,nan_row)
    wc(:,:,1) = [min(w(:,1:(end-1)),(w(:,2:(end)))), nan_col]; % East
    wc(:,:,2) = [nan_row; min(w(2:(end),:),w(1:(end-1),:))]; % North
    wc(2:end,1:(end-1),3) = [min(w(1:(end-1),2:end),w(2:end,1:(end-1)))]; % NE 
    wc(1:(end-1),1:(end-1),4) = [min(w(2:end,2:(end)),w(1:(end-1),1:(end-1)))]; % SE
    % river_width = repmat(w,[1 1 4]);
    % wc(wc == 0) = river_width(wc == 0);
    % wc(repmat(~idx_rivers,[1 1 4])) = 0;
end

function [wc] = area_function(w,nan_col,nan_row)
    wc(:,:,1) = [min(w(:,1:(end-1),1),(w(:,2:(end),1))), nan_col]; % East
    wc(:,:,2) = [nan_row; min(w(2:(end),:,2),w(1:(end-1),:,2))]; % North
    wc(2:end,1:(end-1),3) = [min(w(1:(end-1),2:end,1),w(2:end,1:(end-1),2))]; % NE 
    wc(1:(end-1),1:(end-1),4) = [min(w(2:end,2:(end),2),w(1:(end-1),1:(end-1),1))]; % SE
    % river_width = repmat(w,[1 1 4]);
    % wc(wc == 0) = river_width(wc == 0);
    % wc(repmat(~idx_rivers,[1 1 4])) = 0;
end

function [wc] = Rh_function(w,nan_col,nan_row)
    wc(:,:,1) = [min(w(:,1:(end-1),1),(w(:,2:(end),1))), nan_col]; % East
    wc(:,:,2) = [nan_row; min(w(2:(end),:,2),w(1:(end-1),:,2))]; % North
    wc(2:end,1:(end-1),3) = [min(w(1:(end-1),2:end,1),w(2:end,1:(end-1),2))]; % NE 
    wc(1:(end-1),1:(end-1),4) = [min(w(2:end,2:(end),2),w(1:(end-1),1:(end-1),1))]; % SE
    % river_width = repmat(w,[1 1 4]);
    % wc(wc == 0) = river_width(wc == 0);
    % wc(repmat(~idx_rivers,[1 1 4])) = 0;
end

function [A, W, Rh] = Compute_Subgrid_Properties_Splinewise(A_splines, W_splines, Rh_splines, H)
% Computes A, W, and Rh using loops to evaluate splines for each depth in H.
% A_splines, W_splines, Rh_splines: [nrows x ncols] cell arrays of spline structures.
% H: [nrows x ncols] depth matrix.

[nrows, ncols] = size(H);
A = NaN(nrows, ncols);
W = NaN(nrows, ncols);
Rh = NaN(nrows, ncols);

% Loop over each cell in the matrix
for rowc = 1:nrows
    for colc = 1:ncols
        % Get the corresponding depth for the current cell
        depth = H(rowc, colc);

        % Check if the spline exists and the depth is valid (not NaN)
        if ~isempty(A_splines{rowc, colc}) && ~isnan(depth)
            % Evaluate the splines for area, width, and hydraulic radius at the given depth
            A(rowc, colc) = feval(A_splines{rowc, colc}, depth);
            W(rowc, colc) = feval(W_splines{rowc, colc}, depth);
            Rh(rowc, colc) = feval(Rh_splines{rowc, colc}, depth);
        end
    end
end

% Ensure physical constraints (no values below the threshold)
A = max(A, 1e-12);
W = max(W, 1e-12);
Rh = max(Rh, 1e-12);

end

function [Area_east, Area_north, Rh_east, Rh_north] = Compute_Subgrid_Properties_Polynomialwise(A_east_coeffs, A_north_coeffs, Rh_east_coeffs, Rh_north_coeffs, yc, zc, invert_el)
    % Computes Area_east, Area_north, Rh_east, Rh_north using stored polynomial coefficients
    % A_east_coeffs, A_north_coeffs, Rh_east_coeffs, Rh_north_coeffs: [nrows x ncols] matrices containing polynomial coefficients.
    % Depths: [nrows x ncols] matrix containing the flow depths at each grid cell.

    % Ensure that the polynomial coefficients are in the correct shape for matrix operations
    [nrows, ncols] = size(yc);

    % Evaluate polynomials across the entire matrix using polyval
    Area_east = zeros(size(yc)); Area_north = zeros(size(yc)); Rh_east = zeros(size(yc)); Rh_north = zeros(size(yc));

    for i = 1:size(A_east_coeffs,3)
        % Area_east = Area_east + A_east_coeffs(:,:,i).*(yc - invert_el).^i;
        % Area_north = Area_north + A_north_coeffs(:,:,i).*(yc - invert_el).^i;
        % Rh_east = Rh_east + Rh_east_coeffs(:,:,i).*(yc - invert_el).^i;
        % Rh_north = Rh_north + Rh_north_coeffs(:,:,i).*(yc - invert_el).^i;

        exp = i - 1;

        Area_east = Area_east + A_east_coeffs(:,:,i).*(yc - zc).^exp;
        Area_north = Area_north + A_north_coeffs(:,:,i).*(yc - zc).^exp;
        Rh_east = Rh_east + Rh_east_coeffs(:,:,i).*(yc - zc).^exp;
        Rh_north = Rh_north + Rh_north_coeffs(:,:,i).*(yc - zc).^exp;
    end
    
    % Area_east = polyval(reshape(A_east_coeffs, nrows * ncols, []), repmat(Depths(:),2)); % Evaluates Area_east for all depths in a vectorized way
    % Area_north = polyval(reshape(A_north_coeffs, nrows * ncols, []), Depths(:));
    % Rh_east = polyval(reshape(Rh_east_coeffs, nrows * ncols, []), Depths(:));
    % Rh_north = polyval(reshape(Rh_north_coeffs, nrows * ncols, []), Depths(:));

    % Reshape back to the original matrix shape
    % Area_east = reshape(Area_east, nrows, ncols);
    % Area_north = reshape(Area_north, nrows, ncols);
    % Rh_east = reshape(Rh_east, nrows, ncols);
    % Rh_north = reshape(Rh_north, nrows, ncols);

    % Ensure physical constraints (no values below the threshold)
    Area_east = max(Area_east, 1e-12);
    Area_north = max(Area_north, 1e-12);
    Rh_east = max(Rh_east, 1e-12);
    Rh_north = max(Rh_north, 1e-12);
end


function [Area_east, Area_north, Rh_east, Rh_north] = Compute_Subgrid_Properties_Lookupwise(SubgridTables, yc, zc, direction)
%--------------------------------------------------------------------------
% Fully matrix-based subgrid interpolation (no for loops, no interp1)
%
% INPUTS:
%   SubgridTables - Struct with .depths, .area_*, .Rh_*, all size [nrows x ncols x nz]
%   yc, zc        - Water surface and bed elevation [nrows x ncols]
%   direction     - 'east' or 'north'
%
% OUTPUTS:
%   Area_east / Area_north - Interpolated flow area [m²]
%   Rh_east   / Rh_north   - Interpolated hydraulic radius [m]
%--------------------------------------------------------------------------

[nrows, ncols, nz] = size(SubgridTables.depths);
N = nrows * ncols;

% Flatten depth
depth = max(yc - zc, 0);         % [nrows x ncols]
depth_flat = reshape(depth, [N 1]);  % [N x 1]

% Choose fields
switch lower(direction)
    case 'east'
        area_table = SubgridTables.area_east;
        rh_table   = SubgridTables.Rh_east;
    case 'north'
        area_table = SubgridTables.area_north;
        rh_table   = SubgridTables.Rh_north;
    otherwise
        error('Invalid direction');
end

% Reshape to [N x nz]
depths_all = reshape(SubgridTables.depths, [N, nz]);      % [N x nz]
areas_all  = reshape(area_table, [N, nz]);
rhs_all    = reshape(rh_table,   [N, nz]);

% Mask invalid rows
valid_mask = ~isnan(depths_all) & ~isnan(areas_all) & ~isnan(rhs_all);

% Initialize output
Area_flat = zeros(N, 1);
Rh_flat   = ones(N, 1) * 1e-12;

% For each cell, find bounding indices
above = depths_all >= depth_flat;     % [N x nz]
below = depths_all <= depth_flat;     % [N x nz]

% First above (min depth >= d)
above_idx = max(above .* (1:nz), [], 2);  % [N x 1]
above_idx(above_idx == 0) = 1;

% Last below (max depth <= d)
below_idx = max(below .* (1:nz), [], 2); % [N x 1]
below_idx(below_idx == 0) = 1;

% Make sure they’re not the same
same_idx = (above_idx == below_idx);
above_idx(same_idx & above_idx < nz) = above_idx(same_idx & above_idx < nz) + 1;

% Linear indices for lookup
idx_upper = sub2ind([N, nz], (1:N)', above_idx);
idx_lower = sub2ind([N, nz], (1:N)', below_idx);

% Fetch values
x1 = depths_all(idx_lower);
x2 = depths_all(idx_upper);
A1 = areas_all(idx_lower);
A2 = areas_all(idx_upper);
R1 = rhs_all(idx_lower);
R2 = rhs_all(idx_upper);

% Avoid divide-by-zero
dx = max(x2 - x1, 1e-6);
alpha = (depth_flat - x1) ./ dx;

% Interpolate
Area_flat = A1 + alpha .* (A2 - A1);
Rh_flat   = R1 + alpha .* (R2 - R1);

% Apply floor
Area_flat = max(Area_flat, 1e-12);
Rh_flat   = max(Rh_flat, 1e-12);

% Reshape to 2D
Area_out = reshape(Area_flat, nrows, ncols);
Rh_out   = reshape(Rh_flat,   nrows, ncols);

% Output assignment
if strcmp(direction, 'east')
    Area_east = Area_out;
    Area_north = [];
    Rh_east = Rh_out;
    Rh_north = [];
else
    Area_north = Area_out;
    Area_east = [];
    Rh_north = Rh_out;
    Rh_east = [];
end

end
