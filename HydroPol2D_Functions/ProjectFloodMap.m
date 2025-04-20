function high_res_flood_map = ProjectFloodMap(high_res_DEM, coarse_DEM, coarse_flood_map)
    % Projects a coarse-resolution flood map onto a high-resolution DEM
    
    % Get resolutions
    [nrows_high, ncols_high] = size(high_res_DEM.Z);
    [nrows_coarse, ncols_coarse] = size(coarse_DEM.Z);
    
    high_res_cellsize = high_res_DEM.cellsize;
    coarse_res_cellsize = coarse_DEM.cellsize;
    
    % Compute scale factor
    scale_factor = coarse_res_cellsize / high_res_cellsize;
    if mod(scale_factor, 1) ~= 0
        error('High-resolution cell size must be an integer fraction of coarse resolution.');
    end
    
    % Repeat the coarse flood map to match high resolution grid
    high_res_flood_surface = repelem(coarse_flood_map, scale_factor, scale_factor);
    
    % Ensure dimensions match (crop if necessary)
    high_res_flood_surface = high_res_flood_surface(1:nrows_high, 1:ncols_high);
    
    % Compute flood depth at high resolution
    high_res_flood_map = high_res_flood_surface - double(high_res_DEM.Z);
    
    % Ensure non-negative flood depths
    high_res_flood_map(high_res_flood_map < 0) = 0;
end