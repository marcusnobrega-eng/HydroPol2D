function [A_coeffs, V_coeffs, Rh_east_coeffs, Rh_north_coeffs, W_east_coeffs, W_north_coeffs, A_east_coeffs, A_north_coeffs, Poly_r2] = Subgrid_Properties_Function(DEM_raster, coarse_res, increment_depth, max_depth)
    % Computes subgrid hydraulic properties: smooth splines for A & V, polynomial fits of order 4 with zero intercept for others.

    DEM = DEM_raster.Z;
    cellsize = DEM_raster.cellsize;
    
    poly_order = 2; % Fixed polynomial order

    if mod(coarse_res, cellsize) ~= 0
        error('Resample the raster to a multiple resolution of the input raster.');
    end

    [nrows, ncols] = size(DEM);
    nrows_coarse = ceil(nrows * cellsize / coarse_res) + 1;
    ncols_coarse = ceil(ncols * cellsize / coarse_res) + 1;

    V_coeffs = NaN(nrows_coarse, ncols_coarse, poly_order);
    A_coeffs = NaN(nrows_coarse, ncols_coarse, poly_order);
    Rh_east_coeffs = NaN(nrows_coarse, ncols_coarse, poly_order);
    Rh_north_coeffs = NaN(nrows_coarse, ncols_coarse, poly_order);
    W_east_coeffs = NaN(nrows_coarse, ncols_coarse, poly_order);
    W_north_coeffs = NaN(nrows_coarse, ncols_coarse, poly_order);
    A_east_coeffs = NaN(nrows_coarse, ncols_coarse, poly_order);
    A_north_coeffs = NaN(nrows_coarse, ncols_coarse, poly_order);

    Poly_r2 = struct('Rh_east', NaN(nrows_coarse, ncols_coarse), ...
                     'Rh_north', NaN(nrows_coarse, ncols_coarse), ...
                     'W_east', NaN(nrows_coarse, ncols_coarse), ...
                     'W_north', NaN(nrows_coarse, ncols_coarse), ...
                     'A_east', NaN(nrows_coarse, ncols_coarse), ...
                     'A_north', NaN(nrows_coarse, ncols_coarse));

    % Precompute row and column indices
    row_idx_start = floor((0:nrows_coarse-1) * coarse_res / cellsize) + 1;
    row_idx_end = floor((1:nrows_coarse) * coarse_res / cellsize);
    col_idx_start = floor((0:ncols_coarse-1) * coarse_res / cellsize) + 1;
    col_idx_end = floor((1:ncols_coarse) * coarse_res / cellsize);

    % Calculate the total number of iterations
    total_iterations = nrows_coarse * ncols_coarse;
    count = 0;

    % Loop for coarse grid
    for rowc = 1:nrows_coarse
        for colc = 1:ncols_coarse
            % Get the row and column indices for the current subgrid
            row_idx = row_idx_start(rowc):min(nrows, row_idx_end(rowc));
            col_idx = col_idx_start(colc):min(ncols, col_idx_end(colc));

            sub_DEM = DEM(row_idx, col_idx);

            % Skip if sub_DEM contains all NaN values
            if sum(sum(isnan(sub_DEM)))  || sum(sum(isempty(sub_DEM)))
                continue;
            end

            east_idx = sub_DEM == sub_DEM(:, end);
            east_idx = east_idx(:);
            north_idx = sub_DEM == sub_DEM(1, :);
            north_idx = north_idx(:);

            min_DEM = min(sub_DEM(:));
            max_Depth = max(sub_DEM(:)) - min_DEM;
            depths = linspace(1e-6, max_Depth + 1e-6, 100)';

            depth_matrix = depths - (sub_DEM(:)' - min_DEM);
            wetted_cells = depth_matrix > 0 & ~isnan(sub_DEM(:)');

            area = sum(wetted_cells * cellsize^2, 2, 'omitnan');
            volume = area .* depths;
            width_east = sum(wetted_cells(:, east_idx), 2, 'omitnan') * cellsize;
            width_north = sum(wetted_cells(:, north_idx), 2, 'omitnan') * cellsize;
            area_east = sum(wetted_cells(:, east_idx) .* depth_matrix(:, east_idx) * cellsize, 2, 'omitnan');
            area_north = sum(wetted_cells(:, north_idx) .* depth_matrix(:, north_idx) * cellsize, 2, 'omitnan');
            Rh_east = area_east ./ (width_east + 2 * depths);
            Rh_north = area_north ./ (width_north + 2 * depths);

            if all(area == 0) || all(isnan(area))
                continue;
            end

            % Polynomial design matrix without intercept (order 4)
            X = depths.^(1:poly_order);

            % Fit using X \ y
            A = X \ area;
            V = X \ volume;
            Rh_east_c = X \ Rh_east;
            Rh_north_c = X \ Rh_north;
            W_east_c = X \ width_east;
            W_north_c = X \ width_north;
            A_east_c = X \ area_east;
            A_north_c = X \ area_north;

            % Save coefficients
            A_coeffs(rowc,colc,:) = A;
            V_coeffs(rowc,colc,:) = V;
            Rh_east_coeffs(rowc, colc, :) = Rh_east_c;
            Rh_north_coeffs(rowc, colc, :) = Rh_north_c;
            W_east_coeffs(rowc, colc, :) = W_east_c;
            W_north_coeffs(rowc, colc, :) = W_north_c;
            A_east_coeffs(rowc, colc, :) = A_east_c;
            A_north_coeffs(rowc, colc, :) = A_north_c;

            % Compute R² for each polynomial fit
            Poly_r2.A_r2(rowc, colc) = compute_r2(X * A, area);
            Poly_r2.A_r2(rowc, colc) = compute_r2(X * V, volume);
            Poly_r2.Rh_east(rowc, colc) = compute_r2(X * Rh_east_c, Rh_east);
            Poly_r2.Rh_north(rowc, colc) = compute_r2(X * Rh_north_c, Rh_north);
            Poly_r2.W_east(rowc, colc) = compute_r2(X * W_east_c, width_east);
            Poly_r2.W_north(rowc, colc) = compute_r2(X * W_north_c, width_north);
            Poly_r2.A_east(rowc, colc) = compute_r2(X * A_east_c, area_east);
            Poly_r2.A_north(rowc, colc) = compute_r2(X * A_north_c, area_north);
            
            count = count + 1;
            progress = 100*count / total_iterations
        end
    end
end

function r2 = compute_r2(y_fit, y_obs)
    r2 = 1 - sum((y_obs - y_fit).^2) / sum((y_obs - mean(y_obs)).^2);
end
