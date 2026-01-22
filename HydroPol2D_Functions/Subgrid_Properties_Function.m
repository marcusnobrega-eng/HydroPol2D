function [A_coeffs, V_coeffs, Rh_east_coeffs, Rh_north_coeffs, W_east_coeffs, W_north_coeffs, A_east_coeffs, A_north_coeffs, Poly_NSE, invert_el] = Subgrid_Properties_Function(DEM_raster, Reference_raster, coarse_res, poly_order)
    % Computes subgrid hydraulic properties: smooth splines for A & V, polynomial fits of order 4 with zero intercept for others.

    DEM = DEM_raster.Z;
    cellsize = DEM_raster.cellsize;
    
    if mod(coarse_res, cellsize) ~= 0
        warning('Resample the raster to a multiple resolution of the input raster.');
    end

    [nrows, ncols] = size(DEM);
    % nrows_coarse = nrows * cellsize / coarse_res; %%%% Probably wrong
    % ncols_coarse = ncols * cellsize / coarse_res;
    nrows_coarse = size(Reference_raster.Z,1);
    ncols_coarse = size(Reference_raster.Z,2);

    V_coeffs = NaN(nrows_coarse, ncols_coarse, poly_order + 1);
    A_coeffs = NaN(nrows_coarse, ncols_coarse, poly_order + 1);
    Rh_east_coeffs = NaN(nrows_coarse, ncols_coarse, poly_order + 1);
    Rh_north_coeffs = NaN(nrows_coarse, ncols_coarse, poly_order + 1);
    W_east_coeffs = NaN(nrows_coarse, ncols_coarse, poly_order + 1);
    W_north_coeffs = NaN(nrows_coarse, ncols_coarse, poly_order + 1);
    A_east_coeffs = NaN(nrows_coarse, ncols_coarse, poly_order + 1);
    A_north_coeffs = NaN(nrows_coarse, ncols_coarse, poly_order + 1);
    invert_el = NaN(nrows_coarse, ncols_coarse);

    Poly_NSE = struct('Rh_east', NaN(nrows_coarse, ncols_coarse), ...
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
            invert_el(rowc,colc) = min_DEM; % Invert elevation [m]
            max_Depth = max(sub_DEM(:)) - min_DEM;
            depths = linspace(1e-6, max_Depth + 1e-6, 100)'; % Depths [m]
            elevations = depths + min_DEM; % Elevations [m]

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

            % Add Intercept
            A = [0; A];
            V =[0; V];
            Rh_east_c = [0; Rh_east_c];
            Rh_north_c = [0; Rh_north_c];
            W_east_c = [0; W_east_c];
            W_north_c = [0; W_north_c];
            A_east_c = [0; A_east_c];
            A_north_c = [0; A_north_c];

            if max_Depth == 0
                A = 0*A; A(1) = area(1);
                V = 0*V; V(1) = volume(1);
                Rh_east_c = 0*Rh_east_c; Rh_east_c(2) = 1; % Constant with depth
                Rh_north_c = 0*Rh_north_c; Rh_north_c(2) = 1; % Constant with depth
                W_east_c = 0*W_east_c; W_east_c(1) = coarse_res;
                W_east_c = 0*W_east_c; W_east_c(1) = coarse_res;
                W_north_c = 0*W_north_c; W_north_c(1) = coarse_res;
                A_east_c = 0*A_east_c; A_east_c(2) = coarse_res; % Linear wih depth
                A_north_c = 0*A_north_c; A_north_c(2) = coarse_res; % Linear wih depth
            end


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
            X_star = [zeros(size(X,1),1), X];
            if max_Depth ~= 0
                Poly_NSE.A_NSE(rowc, colc) = compute_NSE(X_star * A, area);
                Poly_NSE.A_NSE(rowc, colc) = compute_NSE(X_star * V, volume);
                Poly_NSE.Rh_east(rowc, colc) = compute_NSE(X_star * Rh_east_c, Rh_east);
                Poly_NSE.Rh_north(rowc, colc) = compute_NSE(X_star * Rh_north_c, Rh_north);
                Poly_NSE.W_east(rowc, colc) = compute_NSE(X_star * W_east_c, width_east);
                Poly_NSE.W_north(rowc, colc) = compute_NSE(X_star * W_north_c, width_north);
                Poly_NSE.A_east(rowc, colc) = compute_NSE(X_star * A_east_c, area_east);
                Poly_NSE.A_north(rowc, colc) = compute_NSE(X_star * A_north_c, area_north);
            end
            count = count + 1;
            if count == 26119
                ttt = 1;
            end
            progress = 100*count / total_iterations
        end
    end
end

function NSE = compute_NSE(y_fit, y_obs)
    NSE = 1 - sum((y_obs - y_fit).^2) / sum((y_obs - mean(y_obs)).^2);
end
