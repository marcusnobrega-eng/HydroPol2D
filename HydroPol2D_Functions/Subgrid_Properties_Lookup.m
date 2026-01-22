function [SubgridTables, invert_el] = Subgrid_Properties_Lookup(DEM_raster, Reference_raster, coarse_res)
%--------------------------------------------------------------------------
% Subgrid_Properties_Lookup
%
% Computes lookup tables of hydraulic properties from a high-resolution DEM
% embedded within a coarse-resolution grid using a HEC-RAS-style subgrid approach.
%
% INPUTS:
%   DEM_raster        - Struct with fields:
%       .Z            - High-resolution DEM (2D matrix) in meters
%       .cellsize     - Resolution of DEM in meters (e.g., 1 m)
%
%   Reference_raster  - Struct with .Z for coarse grid shape (e.g., 50 m cells)
%                       Used to define how many coarse grid cells to compute
%
%   coarse_res        - Coarse cell size in meters (e.g., 50 m)
%
% OUTPUTS:
%   SubgridTables     - Struct containing 3D lookup tables (row, col, depth index):
%       .depths       - Water depths [m] relative to invert
%       .area         - Wetted area of each coarse cell [m²]
%       .volume       - Volume of water in cell [m³]
%       .width_east   - Wetted width along east face [m]
%       .width_north  - Wetted width along north face [m]
%       .Rh_east      - Hydraulic radius on east face [m]
%       .Rh_north     - Hydraulic radius on north face [m]
%
%   invert_el         - Minimum elevation (invert) in each coarse cell [m]
%--------------------------------------------------------------------------
    
    DEM = DEM_raster.Z;                       % High-resolution DEM matrix [m]
    cellsize = DEM_raster.cellsize;           % DEM resolution [m]

    if mod(coarse_res, cellsize) ~= 0
        error('Resample the raster to a multiple resolution of the input raster.');
    end

    [nrows, ncols] = size(DEM);              % Size of the high-resolution DEM
    nrows_coarse = size(Reference_raster.Z, 1); % Number of rows in coarse grid
    ncols_coarse = size(Reference_raster.Z, 2); % Number of cols in coarse grid

    max_depth_points = 200;                  % Max number of depth samples per cell

    % Preallocate 3D lookup tables
    SubgridTables = struct();
    % fields = {'depths', 'area', 'volume', 'area_east', 'area_north', 'width_east', 'width_north', 'Rh_east', 'Rh_north'}; % All Fields
    fields = {'depths', 'area_east', 'area_north', 'width_east', 'width_north', 'Rh_east', 'Rh_north'}; % Required Fields
    for f = fields
        SubgridTables.(f{1}) = NaN(nrows_coarse, ncols_coarse, max_depth_points); % NaN-filled lookup tables
    end
    invert_el = NaN(nrows_coarse, ncols_coarse);  % Minimum elevation per cell [m]

    % Determine fine-resolution row/column indices for each coarse cell
    row_idx_start = floor((0:nrows_coarse-1) * coarse_res / cellsize) + 1;
    row_idx_end = floor((1:nrows_coarse) * coarse_res / cellsize);
    col_idx_start = floor((0:ncols_coarse-1) * coarse_res / cellsize) + 1;
    col_idx_end = floor((1:ncols_coarse) * coarse_res / cellsize);

    % Main loop: for each coarse cell
    for rowc = 1:nrows_coarse
        for colc = 1:ncols_coarse

            % Extract fine-grid indices for this coarse cell
            row_idx = row_idx_start(rowc):min(nrows, row_idx_end(rowc));
            col_idx = col_idx_start(colc):min(ncols, col_idx_end(colc));

            % Get the subgrid DEM patch
            sub_DEM = DEM(row_idx, col_idx); % Local fine-resolution patch

            if all(isnan(sub_DEM(:)))        % Skip empty (all-NaN) regions
                continue;
            end

            sub_DEM_flat = sub_DEM(:);                           % Flattened DEM patch
            invert_el(rowc, colc) = min(sub_DEM_flat);          % Base elevation (invert) for this cell
            elevations = unique(sort(sub_DEM_flat(~isnan(sub_DEM_flat)))); % Unique valid elevations
            min_elev = elevations(1);                           % Minimum elevation [m]
            max_elev = elevations(end);                         % Maximum elevation [m]
            buffer_elev = max_elev + max(0.1, 0.05 * (max_elev - min_elev)); % Buffer elevation for extrapolation

            % Smart sampling: Use unique elevation transitions + buffer
            elev_samples = unique([elevations; linspace(min_elev, max_elev, 15)'; buffer_elev]); 
            depths = elev_samples - min_elev;                   % Water depths relative to invert [m]
            nDepths = length(depths);                           % Number of depth samples

            % Get face indices (east face = last column, north face = first row)
            east_idx = sub_DEM == sub_DEM(:, end); east_idx = east_idx(:); % East face mask (vertical face)
            north_idx = sub_DEM == sub_DEM(1, :); north_idx = north_idx(:); % North face mask (horizontal face)
            sub_DEM_vec = sub_DEM(:)';                          % DEM reshaped for vector operations

            % Allocate output vectors for this coarse cell
            area = zeros(nDepths, 1);                           % Wetted area [m²]
            volume = zeros(nDepths, 1);                         % Volume [m³]
            width_east = zeros(nDepths, 1);                     % Wetted width at east face [m]
            width_north = zeros(nDepths, 1);                    % Wetted width at north face [m]
            area_east = zeros(nDepths, 1);                      % Flow area at east face [m²]
            area_north = zeros(nDepths, 1);                     % Flow area at north face [m²]
            Rh_east = zeros(nDepths, 1);                        % Hydraulic radius at east face [m]
            Rh_north = zeros(nDepths, 1);                       % Hydraulic radius at north face [m]

            % Loop through each sampled depth
            for i = 1:nDepths
                depth = depths(i);                              % Current water depth [m]
                elev = depth + min_elev;                        % Water surface elevation [m]

                depth_matrix = elev - sub_DEM_vec;              % Local depth above terrain
                depth_matrix(depth_matrix < 0) = 0;             % Negative values mean dry cells
                wetted = depth_matrix > 0;                      % Binary mask of wet cells

                area(i) = sum(wetted) * cellsize^2;             % Wetted plan area [m²]
                volume(i) = sum(depth_matrix .* wetted) * cellsize^2; % Volume in the cell [m³]

                width_east(i) = sum(wetted(east_idx)) * cellsize;      % Wetted width along east face [m]
                width_north(i) = sum(wetted(north_idx)) * cellsize;    % Wetted width along north face [m]

                area_east(i) = sum(depth_matrix(east_idx) .* wetted(east_idx)) * cellsize; % Area east face [m²]
                area_north(i) = sum(depth_matrix(north_idx) .* wetted(north_idx)) * cellsize; % Area north face [m²]

                Rh_east(i) = area_east(i) / (width_east(i) + 2 * depth);  % Hydraulic radius east [m]
                Rh_north(i) = area_north(i) / (width_north(i) + 2 * depth); % Hydraulic radius north [m]
            end

            % Compute relative area change between depths
            area_diff = abs([0; diff(area)]);  % Pad to match size
            area_rel_change = area_diff ./ max(area, 1e-6);  % Now same size as area/depths
            important_idx = (area_rel_change > 0.01) | (area_diff > 0.05);  % Logical vector [nDepths x 1]
            important_idx(1) = true;       % Always keep first
            important_idx(end) = true;     % Always keep last

            % Ensure anchor points at regular intervals
            anchor_idx = false(size(depths));
            anchor_spacing = round(length(depths) / 20);  % 20 anchors evenly spread
            anchor_idx(1:anchor_spacing:end) = true;
            
            % Combine and deduplicate
            final_idx = important_idx | anchor_idx;
            final_idx = find(final_idx);
            final_idx = unique(final_idx);  % Sorted indices
            
            % If still too many, sample evenly
            if length(final_idx) > max_depth_points
                final_idx = round(linspace(1, length(depths), max_depth_points));
            end
            
            % Store reduced subset
            d_end = length(final_idx);
            SubgridTables.depths(rowc, colc, 1:d_end)       = depths(final_idx);       % Depth [m]
            SubgridTables.area(rowc, colc, 1:d_end)         = area(final_idx);         % Surface Area [m²]
            SubgridTables.volume(rowc, colc, 1:d_end)       = volume(final_idx);       % Volume [m³]
            SubgridTables.area_east(rowc, colc, 1:d_end)    = area_east(final_idx);    % Area east [m2]
            SubgridTables.area_north(rowc, colc, 1:d_end)    = area_north(final_idx);  % Area north [m2]
            SubgridTables.width_east(rowc, colc, 1:d_end)   = width_east(final_idx);   % Width east [m]
            SubgridTables.width_north(rowc, colc, 1:d_end)  = width_north(final_idx);  % Width north [m]
            SubgridTables.Rh_east(rowc, colc, 1:d_end)      = Rh_east(final_idx);      % Hydraulic radius east [m]
            SubgridTables.Rh_north(rowc, colc, 1:d_end)     = Rh_north(final_idx);     % Hydraulic radius north [m]
            
            % Optional: show progress percentage in command window
            cell_count = (rowc - 1) * ncols_coarse + colc;                   % Current cell index
            total_cells = nrows_coarse * ncols_coarse;                       % Total number of coarse cells
            if mod(cell_count, 100) == 0 || cell_count == total_cells        % Print every 100 or at end
                fprintf('Subgrid Preprocessing Progress: %6.2f%% (%d of %d cells)\n', ...
                    100 * cell_count / total_cells, cell_count, total_cells);
            end


        end
    end
end
