function DEM_filled = fillsinks_for_routing(DEM, flags)
%FILLSINKS_FOR_ROUTING Fill DEM sinks using the active routing connectivity.
%
% TopoToolbox fillsinks uses 8-neighbor image reconstruction. That is
% appropriate for D8 routing, but it can leave depressions that are only
% drainable diagonally when the hydraulic solver uses D4 faces. For D4
% routing, fill with 4-neighbor reconstruction so DEM conditioning and
% routing connectivity are consistent.

if nargin < 2 || ~isstruct(flags) || ~isfield(flags, 'flag_D8') || flags.flag_D8 == 1
    DEM_filled = fillsinks(DEM);
    return;
end

DEM_filled = DEM;
dem = DEM.Z;

inan = isnan(dem);
dem(inan) = -inf;

marker = -dem;
interior = false(size(dem));
if size(dem, 1) > 2 && size(dem, 2) > 2
    interior(2:end-1, 2:end-1) = true;
end
marker(interior & ~inan) = -inf;

demfs = -imreconstruct(marker, -dem, 4);
if isfloat(demfs)
    demfs(inan) = nan;
end

DEM_filled.Z = demfs;
end
