function specs = hydropol2d_finalize_parameter_specs(specs)
%HYDROPOL2D_FINALIZE_PARAMETER_SPECS Add internal GA coordinates to specs.
%
% Use `transform = 'log'` for strictly positive parameters where genetic
% search should be multiplicative. Use `transform = 'linear'` for bounded
% fractions, switches, offsets, and values that can be zero.

for i = 1:numel(specs)
    validate_spec(specs(i));
    if strcmpi(specs(i).transform, 'log')
        specs(i).lower_internal = log(specs(i).lower);
        specs(i).upper_internal = log(specs(i).upper);
        specs(i).initial_internal = log(specs(i).initial);
    else
        specs(i).lower_internal = specs(i).lower;
        specs(i).upper_internal = specs(i).upper;
        specs(i).initial_internal = specs(i).initial;
    end
end
end

function validate_spec(s)
    assert(~isempty(s.name), 'Every calibration parameter needs a name.');
    assert(isfinite(s.lower) && isfinite(s.upper) && s.upper > s.lower, ...
        'Invalid bounds for calibration parameter "%s".', s.name);
    assert(isfinite(s.initial), ...
        'Invalid initial value for calibration parameter "%s".', s.name);
    assert(s.initial >= s.lower && s.initial <= s.upper, ...
        'Initial value for "%s" must be inside [lower, upper].', s.name);
    if strcmpi(s.transform, 'log')
        assert(s.lower > 0 && s.upper > 0 && s.initial > 0, ...
            'Log-transformed parameter "%s" must be strictly positive.', s.name);
    end
end
