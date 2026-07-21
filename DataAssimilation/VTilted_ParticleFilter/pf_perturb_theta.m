function theta = pf_perturb_theta(theta, specs, scale)
%PF_PERTURB_THETA Apply bounded Gaussian post-resampling parameter noise.

if nargin < 3 || isempty(scale)
    scale = 1;
end

for i = 1:numel(specs)
    spec = specs(i);
    if ~spec.enabled
        continue;
    end
    if strcmpi(char(spec.target_type), 'IC')
        continue; % ponytail: initial conditions are sampled at t=0, not jittered mid-run.
    end

    value = get_theta_value(theta, spec);
    sigma = scale * spec.perturb_std;
    if sigma <= 0
        continue;
    end

    if spec.transform == "log"
        value = exp(log(max(value, eps)) + sigma * randn());
    else
        value = value + sigma * randn();
    end

    value = min(max(value, spec.lower), spec.upper);
    theta = set_theta_value(theta, spec, value);
end
end

function value = get_theta_value(theta, spec)
target = lower(char(spec.target_type));
property = char(spec.property);
class_id = spec.class_id;

switch target
    case 'lulc'
        value = theta.lulc(class_id).(property);
    case 'soil'
        value = theta.soil(class_id).(property);
    case 'gw'
        value = theta.gw(class_id).(property);
    case 'ic'
        value = theta.ic(class_id).(property);
    otherwise
        error('pf_perturb_theta:target', 'Unknown target_type: %s', spec.target_type);
end
end

function theta = set_theta_value(theta, spec, value)
target = lower(char(spec.target_type));
property = char(spec.property);
class_id = spec.class_id;

switch target
    case 'lulc'
        theta.lulc(class_id).(property) = value;
    case 'soil'
        theta.soil(class_id).(property) = value;
    case 'gw'
        theta.gw(class_id).(property) = value;
    case 'ic'
        theta.ic(class_id).(property) = value;
    otherwise
        error('pf_perturb_theta:target', 'Unknown target_type: %s', spec.target_type);
end
end
