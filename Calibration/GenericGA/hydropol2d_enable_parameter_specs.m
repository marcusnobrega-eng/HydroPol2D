function specs = hydropol2d_enable_parameter_specs(specs, enablePatterns)
%HYDROPOL2D_ENABLE_PARAMETER_SPECS Enable specs by name or wildcard pattern.
%
% Examples
% --------
% Enable exact names:
%   specs = hydropol2d_enable_parameter_specs(specs, [
%       "lulc_roughness_multiplier_all"
%       "soil_theta_i_multiplier_all"]);
%
% Enable all per-class soil theta parameters:
%   specs = hydropol2d_enable_parameter_specs(specs, "soil_*_theta_i");

if ischar(enablePatterns) || isstring(enablePatterns)
    enablePatterns = string(enablePatterns);
else
    enablePatterns = string(enablePatterns(:));
end

for i = 1:numel(specs)
    specs(i).enabled = false;
end

names = string({specs.name});
for p = 1:numel(enablePatterns)
    pattern = enablePatterns(p);
    if contains(pattern, "*")
        expr = "^" + regexptranslate('wildcard', char(pattern)) + "$";
        match = ~cellfun('isempty', regexp(cellstr(names), char(expr), 'once'));
    else
        match = names == pattern;
    end
    for i = find(match)
        specs(i).enabled = true;
    end
end
end
