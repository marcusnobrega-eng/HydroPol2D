function s = hydropol2d_make_parameter_spec(name, target, field, idx, ...
        operation, lower, upper, initial, transform, enabled)
%HYDROPOL2D_MAKE_PARAMETER_SPEC Create a LULC/SOIL calibration parameter.
%
% target:
%   'LULC' or 'SOIL'
%
% operation:
%   'value'      -> write the calibrated value directly.
%   'multiplier' -> multiply the baseline field values for idx.
%   'offset'     -> add the calibrated value to baseline field values.

s = hydropol2d_empty_parameter_spec();
s.name = char(name);
s.target = char(target);
s.field = char(field);
s.index = idx(:);
s.operation = char(operation);
s.lower = double(lower);
s.upper = double(upper);
s.initial = double(initial);
s.transform = char(transform);
s.enabled = logical(enabled);
end
