function flat = hp2d_remap_intensive(mapping, values, minimum_coverage, fill)
%HP2D_REMAP_INTENSIVE Coverage-normalised mesh-to-raster remap for INTENSIVE fields.
%
%   Divides the conservative remap by each raster cell's mesh coverage, so a
%   constant field maps to that constant everywhere the mesh reaches, and returns
%   FILL where coverage is below MINIMUM_COVERAGE.
%
%   Without this, water surface elevation was the visible casualty. Measured on a
%   90 m Pune mesh over a 5 m raster: 184 raster cells had below 50% coverage and
%   805 interior cells had none, so a 605 m water surface was written as 302 m at
%   the worst partly covered cell and as 0 m where the mesh did not reach. The
%   reported WSE range was 0.00-1030.99 m over a bed of 605.29-1039.73 m.
%
%   Use HydroPol2D_Read_Overlap's mapping.mesh_to_raster directly only for
%   EXTENSIVE quantities that are meant to sum.

if nargin < 3 || isempty(minimum_coverage), minimum_coverage = 0.5; end
if nargin < 4 || isempty(fill), fill = NaN; end
assert(numel(values) == size(mapping.mesh_to_raster,2), ...
    'hp2d_remap_intensive:size', 'Field length must match the overlap mesh.');

flat = mapping.mesh_to_raster * double(values(:));
coverage = mapping.raster_coverage;
positive = coverage > 0;
flat(positive) = flat(positive) ./ coverage(positive);
flat(~positive) = fill;
flat(coverage < minimum_coverage) = fill;
end
