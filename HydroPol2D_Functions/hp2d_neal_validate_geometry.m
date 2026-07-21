function hp2d_neal_validate_geometry(River_Width, River_Depth, Resolution)
%HP2D_NEAL_VALIDATE_GEOMETRY Validate channel geometry for Neal-mode subgrid.

has_width = River_Width > 0;
has_depth = River_Depth > 0;

invalid_negative = (River_Width < 0) | (River_Depth < 0);
if any(invalid_negative(:))
    error('HydroPol2D:NealSubgridInvalidGeometry', ...
        'Neal subgrid received negative River_Width or River_Depth values.');
end

invalid_partial = xor(has_width, has_depth);
if any(invalid_partial(:))
    error('HydroPol2D:NealSubgridIncompleteGeometry', ...
        'Neal subgrid requires River_Width and River_Depth together for every active channel cell.');
end

invalid_resolution = has_width & (River_Width > Resolution);
if any(invalid_resolution(:))
    bad_count = nnz(invalid_resolution);
    error('HydroPol2D:NealSubgridWidthTooLarge', ...
        ['Neal subgrid cannot place a channel wider than one cell in a ', ...
         'single centerline cell. %d channel cells have River_Width > ', ...
         'Resolution. Represent these reaches with a multi-cell channel ', ...
         'bathymetry or use a finer grid.'], bad_count);
end

end
