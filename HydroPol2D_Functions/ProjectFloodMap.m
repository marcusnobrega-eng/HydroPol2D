function high_res_flood_map = ProjectFloodMap(high_res_DEM, coarse_DEM, coarse_flood_map)
    % Projects a coarse-resolution flood map onto a high-resolution DEM
    
    % Get sizes
    [nrows_high, ncols_high] = size(high_res_DEM.Z);
    
    high_res_cellsize = high_res_DEM.cellsize;
    coarse_res_cellsize = coarse_DEM.cellsize;
    
    % Compute approximate scale factor
    scale_factor = coarse_res_cellsize / high_res_cellsize;
    scale_factor_round = round(scale_factor);
    
    % Allow small mismatch
    if abs(scale_factor - scale_factor_round) > 1e-6
        warning('Resolution ratio is not exact. Using rounded scale factor.');
    end
    
    scale_factor = scale_factor_round;
    
    if scale_factor < 1
        error('Coarse resolution must be >= high-resolution cell size.');
    end
    
    % Repeat the coarse flood map
    high_res_flood_surface = repelem(double(coarse_flood_map), scale_factor, scale_factor);
    
    % Make output exactly match DEM size
    temp = zeros(nrows_high, ncols_high);
    nr = min(nrows_high, size(high_res_flood_surface,1));
    nc = min(ncols_high, size(high_res_flood_surface,2));
    temp(1:nr, 1:nc) = high_res_flood_surface(1:nr, 1:nc);
    high_res_flood_surface = temp;
    
    % Compute flood depth
    high_res_flood_map = high_res_flood_surface - double(high_res_DEM.Z);
    
    % Ensure non-negative flood depths
    high_res_flood_map(high_res_flood_map < 0) = 0;
end