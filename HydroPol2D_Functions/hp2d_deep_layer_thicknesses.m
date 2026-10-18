function layers = hp2d_deep_layer_thicknesses(max_depth_m, first_m, ratio, max_m, explicit_m)
%HP2D_DEEP_LAYER_THICKNESSES Deep-slot thicknesses (m) reaching max_depth_m below the soil base.
%   Ruler: first_m, multiplied by ratio per slot, capped at max_m. An explicit list is
%   used as given and its last value repeated until the sum reaches max_depth_m.
%   Same rule as HydroPol2D-Python modules.subsurface.deep_layer_thicknesses.
if nargin < 2 || isempty(first_m), first_m = 0.25; end
if nargin < 3 || isempty(ratio), ratio = 2.0; end
if nargin < 4 || isempty(max_m), max_m = 4.0; end
if nargin < 5, explicit_m = []; end
if max_depth_m <= 0
    error('HydroPol2D:Subsurface', 'max_depth_m must be positive.');
end
if ~isempty(explicit_m)
    layers = double(explicit_m(:)).';
    if any(layers <= 0)
        error('HydroPol2D:Subsurface', 'Explicit deep-layer thicknesses must be positive.');
    end
    while sum(layers) < max_depth_m - 1e-9
        layers(end+1) = layers(end); %#ok<AGROW>
    end
    return
end
if first_m <= 0 || ratio < 1 || max_m < first_m
    error('HydroPol2D:Subsurface', 'Deep-layer rule needs first > 0, ratio >= 1 and max >= first.');
end
layers = [];
thickness = first_m;
while sum(layers) < max_depth_m - 1e-9
    layers(end+1) = min(thickness, max_m); %#ok<AGROW>
    thickness = thickness*ratio;
end
end
