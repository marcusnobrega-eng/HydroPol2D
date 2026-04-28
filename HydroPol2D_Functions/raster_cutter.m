function rainfall_raster = raster_cutter(DEM_raster, rR, rr2, reset_cache)
%RASTER_CUTTER  Fast resample/reproject an input raster (rainfall/ETP/etc) to DEM grid.
%
% Optimized design
% ----------------
% 1) Cache DEM geometry and source-grid mapping
% 2) Rebuild cache only when:
%       - reset_cache == true
%       - cache is empty
%       - DEM/source geometry changes
% 3) Interpolate directly onto DEM domain
% 4) Reuse a cached griddedInterpolant object (only update Values each call)
% 5) Crop source raster to DEM-overlap window before interpolation
% 6) Avoid expensive clip(...)
% 7) Minimize repeated masking / cleaning
%
% Inputs
% ------
% DEM_raster   : DEM raster (projected)
% rR           : raster reference of source raster
% rr2          : source raster values
% reset_cache  : logical flag, usually (k == 1)
%
% Output
% ------
% rainfall_raster : DEM-clone raster with aligned source values

    persistent CACHE

    if nargin < 4 || isempty(reset_cache)
        reset_cache = false;
    end
    reset_cache = logical(reset_cache);

    % ---------------------------------------------------------------------
    % Decide whether to rebuild cache
    % ---------------------------------------------------------------------
    mustRebuild = false;

    if isempty(CACHE) || reset_cache
        mustRebuild = true;
    else
        try
            currentSig = local_build_signature(DEM_raster, rR, size(rr2));
            mustRebuild = ~isequal(CACHE.signature, currentSig);
        catch
            mustRebuild = true;
        end
    end

    if mustRebuild
        CACHE = local_build_cache(DEM_raster, rR, size(rr2));
    end

    % ---------------------------------------------------------------------
    % Read only overlapping source window
    % ---------------------------------------------------------------------
    rr_sub = double(rr2(CACHE.rowIdx, CACHE.colIdx));

    % Make row orientation consistent with ascending ySrcSub
    if CACHE.flipRows
        rr_sub = flipud(rr_sub);
    end

    % ---------------------------------------------------------------------
    % Reuse cached interpolant object: only update Values
    % ---------------------------------------------------------------------
    CACHE.F.Values = rr_sub;
    new_mesh = CACHE.F(CACHE.Yq, CACHE.Xq);

    % ---------------------------------------------------------------------
    % Fast cleanup
    % ---------------------------------------------------------------------
    % Outside DEM domain or outside interpolation support -> 0
    new_mesh(~CACHE.keepMask) = 0;

    % Remove NaN / impossible values
    bad = isnan(new_mesh) | (new_mesh < 0) | (new_mesh > 300);
    if any(bad(:))
        new_mesh(bad) = 0;
    end

    % ---------------------------------------------------------------------
    % Return DEM clone
    % ---------------------------------------------------------------------
    rainfall_raster = DEM_raster;
    rainfall_raster.Z = new_mesh;
    rainfall_raster.name = 'rainfall_resampled';
end


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function CACHE = local_build_cache(DEM_raster, rR, rr2_size)

    CACHE = struct();

    % ---------------------------------------------------------------------
    % DEM geometry
    % ---------------------------------------------------------------------
    demRef = DEM_raster.georef.SpatialRef;

    if ~isprop(demRef, "ProjectedCRS") || isempty(demRef.ProjectedCRS)
        error("DEM raster has no valid ProjectedCRS. DEM must be projected.");
    end

    DEM_CRS = demRef.ProjectedCRS;

    nRowsDEM = DEM_raster.size(1);
    nColsDEM = DEM_raster.size(2);

    xDEM = demRef.XWorldLimits(1) + (0:nColsDEM-1) * demRef.CellExtentInWorldX;

    if isprop(demRef, "RowsStartFrom") && strcmpi(demRef.RowsStartFrom, "north")
        yDEM = demRef.YWorldLimits(2) - (0:nRowsDEM-1) * demRef.CellExtentInWorldY;
    else
        yDEM = demRef.YWorldLimits(1) + (0:nRowsDEM-1) * demRef.CellExtentInWorldY;
    end

    [X_dem, Y_dem] = meshgrid(xDEM, yDEM);

    CACHE.signature = local_build_signature(DEM_raster, rR, rr2_size);
    CACHE.demMask   = ~isnan(DEM_raster.Z);

    % DEM bounds
    demXmin = min(xDEM);
    demXmax = max(xDEM);
    demYmin = min(yDEM);
    demYmax = max(yDEM);

    % ---------------------------------------------------------------------
    % Detect source type
    % ---------------------------------------------------------------------
    isRainGeographic = isa(rR, 'map.rasterref.GeographicCellsReference');
    isRainProjected  = isa(rR, 'map.rasterref.MapCellsReference') || isprop(rR, "ProjectedCRS");

    if ~isRainGeographic && ~isRainProjected
        error('Unsupported raster reference type for source raster.');
    end

    Rain_CRS = [];
    if isRainGeographic
        if isprop(rR,'GeographicCRS') && ~isempty(rR.GeographicCRS)
            Rain_CRS = rR.GeographicCRS;
        else
            Rain_CRS = geocrs(4326);
        end
    else
        if isprop(rR,"ProjectedCRS") && ~isempty(rR.ProjectedCRS)
            Rain_CRS = rR.ProjectedCRS;
        end
    end

    % ---------------------------------------------------------------------
    % Source coordinate vectors
    % ---------------------------------------------------------------------
    nRowsSrc = rr2_size(1);
    nColsSrc = rr2_size(2);

    if isRainGeographic
        xSrc = rR.LongitudeLimits(1) + (0:nColsSrc-1) * rR.CellExtentInLongitude;

        if isprop(rR,"RowsStartFrom") && strcmpi(rR.RowsStartFrom,"north")
            ySrc = rR.LatitudeLimits(2) - (0:nRowsSrc-1) * rR.CellExtentInLatitude;
            flipRows = true;
        else
            ySrc = rR.LatitudeLimits(1) + (0:nRowsSrc-1) * rR.CellExtentInLatitude;
            flipRows = false;
        end
    else
        xSrc = rR.XWorldLimits(1) + (0:nColsSrc-1) * rR.CellExtentInWorldX;

        if isprop(rR,"RowsStartFrom") && strcmpi(rR.RowsStartFrom,"north")
            ySrc = rR.YWorldLimits(2) - (0:nRowsSrc-1) * rR.CellExtentInWorldY;
            flipRows = true;
        else
            ySrc = rR.YWorldLimits(1) + (0:nRowsSrc-1) * rR.CellExtentInWorldY;
            flipRows = false;
        end
    end

    if flipRows
        ySrcAsc = fliplr(ySrc);
    else
        ySrcAsc = ySrc;
    end

    CACHE.flipRows = flipRows;

    % ---------------------------------------------------------------------
    % Build DEM query points in source coordinate system
    % ---------------------------------------------------------------------
    if isRainGeographic
        [lat_q, lon_q] = projinv(DEM_CRS, X_dem, Y_dem);
        Xq = lon_q;
        Yq = lat_q;

    else
        if ~isempty(Rain_CRS) && ~isequal(Rain_CRS, DEM_CRS)
            [lat_q, lon_q] = projinv(DEM_CRS, X_dem, Y_dem);
            [x_q, y_q] = projfwd(Rain_CRS, lat_q, lon_q);
            Xq = x_q;
            Yq = y_q;
        else
            Xq = X_dem;
            Yq = Y_dem;

            if isempty(Rain_CRS)
                srcXmin = min(xSrc);
                srcXmax = max(xSrc);
                srcYmin = min(ySrcAsc);
                srcYmax = max(ySrcAsc);

                overlapX = ~(srcXmax < demXmin || srcXmin > demXmax);
                overlapY = ~(srcYmax < demYmin || srcYmin > demYmax);

                if ~(overlapX && overlapY)
                    error(['Source raster is projected but CRS metadata is missing and extents do not overlap DEM. ' ...
                           'Fix the GeoTIFF CRS metadata or provide compatible input.']);
                end
            end
        end
    end

    CACHE.Xq = Xq;
    CACHE.Yq = Yq;

    % ---------------------------------------------------------------------
    % Crop source grid to DEM overlap
    % ---------------------------------------------------------------------
    xqMin = min(Xq(:), [], 'omitnan');
    xqMax = max(Xq(:), [], 'omitnan');
    yqMin = min(Yq(:), [], 'omitnan');
    yqMax = max(Yq(:), [], 'omitnan');

    dx = local_mean_step(xSrc);
    dy = local_mean_step(ySrcAsc);

    xPad = max(abs(dx), eps);
    yPad = max(abs(dy), eps);

    xMinCrop = xqMin - xPad;
    xMaxCrop = xqMax + xPad;
    yMinCrop = yqMin - yPad;
    yMaxCrop = yqMax + yPad;

    colIdx = find(xSrc    >= xMinCrop & xSrc    <= xMaxCrop);
    rowIdxAsc = find(ySrcAsc >= yMinCrop & ySrcAsc <= yMaxCrop);

    if isempty(colIdx) || isempty(rowIdxAsc)
        colIdx    = 1:nColsSrc;
        rowIdxAsc = 1:nRowsSrc;
    end

    if flipRows
        rowIdx = nRowsSrc - rowIdxAsc + 1;
        rowIdx = sort(rowIdx, 'ascend');
    else
        rowIdx = rowIdxAsc;
    end

    xSrcSub = xSrc(colIdx);
    ySrcSub = ySrcAsc(rowIdxAsc);

    CACHE.rowIdx  = rowIdx;
    CACHE.colIdx  = colIdx;
    CACHE.xSrcSub = xSrcSub;
    CACHE.ySrcSub = ySrcSub;

    % ---------------------------------------------------------------------
    % Prebuild interpolant once and reuse
    % ---------------------------------------------------------------------
    CACHE.F = griddedInterpolant( ...
        {ySrcSub, xSrcSub}, ...
        zeros(numel(ySrcSub), numel(xSrcSub)), ...
        'linear', 'none');

    % ---------------------------------------------------------------------
    % Keep-mask: only valid DEM cells are retained
    % Also reject invalid projected query coordinates
    % ---------------------------------------------------------------------
    finiteQueryMask = isfinite(Xq) & isfinite(Yq);
    CACHE.keepMask = CACHE.demMask & finiteQueryMask;
end


function sig = local_build_signature(DEM_raster, rR, rr2_size)

    demRef = DEM_raster.georef.SpatialRef;

    sig = struct();
    sig.demSize = DEM_raster.size;
    sig.srcSize = rr2_size;

    sig.demXLim = local_try_get(demRef, 'XWorldLimits', []);
    sig.demYLim = local_try_get(demRef, 'YWorldLimits', []);
    sig.demDX   = local_try_get(demRef, 'CellExtentInWorldX', []);
    sig.demDY   = local_try_get(demRef, 'CellExtentInWorldY', []);
    sig.demRSF  = local_try_get(demRef, 'RowsStartFrom', '');

    if isprop(demRef,'ProjectedCRS') && ~isempty(demRef.ProjectedCRS)
        sig.demCRSName = char(string(demRef.ProjectedCRS.Name));
    else
        sig.demCRSName = '';
    end

    sig.srcClass = class(rR);

    if isa(rR,'map.rasterref.GeographicCellsReference')
        sig.srcLonLim = local_try_get(rR, 'LongitudeLimits', []);
        sig.srcLatLim = local_try_get(rR, 'LatitudeLimits', []);
        sig.srcDX     = local_try_get(rR, 'CellExtentInLongitude', []);
        sig.srcDY     = local_try_get(rR, 'CellExtentInLatitude', []);
        sig.srcRSF    = local_try_get(rR, 'RowsStartFrom', '');

        if isprop(rR,'GeographicCRS') && ~isempty(rR.GeographicCRS)
            sig.srcCRSName = char(string(rR.GeographicCRS.Name));
        else
            sig.srcCRSName = 'WGS84_assumed';
        end
    else
        sig.srcXLim = local_try_get(rR, 'XWorldLimits', []);
        sig.srcYLim = local_try_get(rR, 'YWorldLimits', []);
        sig.srcDX   = local_try_get(rR, 'CellExtentInWorldX', []);
        sig.srcDY   = local_try_get(rR, 'CellExtentInWorldY', []);
        sig.srcRSF  = local_try_get(rR, 'RowsStartFrom', '');

        if isprop(rR,'ProjectedCRS') && ~isempty(rR.ProjectedCRS)
            sig.srcCRSName = char(string(rR.ProjectedCRS.Name));
        else
            sig.srcCRSName = '';
        end
    end
end


function out = local_try_get(s, fieldName, defaultVal)
    try
        out = s.(fieldName);
    catch
        out = defaultVal;
    end
end


function d = local_mean_step(v)
    if numel(v) < 2
        d = 1;
    else
        d = mean(diff(v), 'omitnan');
    end
end