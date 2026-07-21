function s = hydropol2d_make_inputdata_parameter_spec(name, path, lower, ...
        upper, initial, transform, enabled)
%HYDROPOL2D_MAKE_INPUTDATA_PARAMETER_SPEC Calibrate InputData_Bypass fields.
%
% path is a dot-delimited path inside InputData_Bypass, for example:
%   'general.slope_outlet'
%   'flags.flag_infiltration'

s = hydropol2d_empty_parameter_spec();
s.name = char(name);
s.target = 'InputData';
s.path = char(path);
s.operation = 'value';
s.lower = double(lower);
s.upper = double(upper);
s.initial = double(initial);
s.transform = char(transform);
s.enabled = logical(enabled);
end
