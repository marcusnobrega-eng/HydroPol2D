function [weights, objective, log_likelihood, Info] = pf_paper_weights(PF, Domain, obs_table, sim_matrix)
%PF_PAPER_WEIGHTS Particle weights following the WRR paper pseudocode.
%
% w_m = sum_k beta_k * N(z_k - z_{p,m,k}; 0, sigma_k)
% followed by w_m <- w_m / sum_m w_m.
%
% Discharge is converted from m3/s to catchment-normalized mm/h before
% weighting, matching the paper's seepage flux units.

N = size(sim_matrix, 1);
types = ["discharge", "groundwater_depth", "soil_moisture"];
group_likelihood = zeros(N, numel(types));
group_log_likelihood = -inf(N, numel(types));
group_beta = zeros(1, numel(types));

for g = 1:numel(types)
    type = types(g);
    beta = observation_group_weight(PF, type);
    group_beta(g) = beta;
    if beta <= 0
        continue;
    end

    idx = obs_table.type == type;
    if ~any(idx)
        continue;
    end

    obs_values = obs_table.value(idx)';
    sim_values = sim_matrix(:, idx);
    [obs_w, sim_w, sigma] = convert_to_weighting_units(PF, Domain, type, obs_values, sim_values);
    residual = sim_w - obs_w;
    log_pdf = -0.5 .* (residual ./ sigma) .^ 2 - log(sigma .* sqrt(2 * pi));
    group_log_likelihood(:, g) = logmeanexp_local(log_pdf, 2);
    group_likelihood(:, g) = exp(group_log_likelihood(:, g) - max(group_log_likelihood(:, g), [], 'omitnan'));
end

active = group_beta > 0 & any(isfinite(group_log_likelihood), 1);
if ~any(active)
    weights = ones(N, 1) ./ N;
    objective = ones(N, 1);
    log_likelihood = log(weights);
else
    beta = group_beta(active);
    beta = beta ./ sum(beta);
    log_terms = group_log_likelihood(:, active) + log(beta(:))';
    log_likelihood = logsumexp_local(log_terms, 2);
    relative_likelihood = exp(log_likelihood - max(log_likelihood, [], 'omitnan'));
    weights = relative_likelihood ./ sum(relative_likelihood);
    objective = 1 - relative_likelihood ./ max(max(relative_likelihood), realmin);
end

Info = struct();
Info.method = "paper_beta_normal_pdf";
Info.types = types;
Info.beta = group_beta;
Info.group_likelihood = group_likelihood;
Info.total_likelihood = exp(log_likelihood - max(log_likelihood, [], 'omitnan'));
end

function y = logmeanexp_local(x, dim)
m = max(x, [], dim, 'omitnan');
y = m + log(mean(exp(x - m), dim, 'omitnan'));
y(~isfinite(m)) = -inf;
end

function y = logsumexp_local(x, dim)
m = max(x, [], dim, 'omitnan');
y = m + log(sum(exp(x - m), dim, 'omitnan'));
y(~isfinite(m)) = -inf;
end

function beta = observation_group_weight(PF, type)
beta = 0;
if ~isfield(PF, 'weighting') || ~isfield(PF.weighting, 'group_weights')
    return;
end
field = char(type);
if isfield(PF.weighting.group_weights, field)
    beta = PF.weighting.group_weights.(field);
end
end

function [obs_w, sim_w, sigma] = convert_to_weighting_units(PF, Domain, type, obs_values, sim_values)
switch char(type)
    case 'discharge'
        area_m2 = sum(Domain.valid(:)) * Domain.cell_area_m2;
        factor = 1000 * 3600 / area_m2; % m3/s -> mm/h
        obs_w = obs_values .* factor;
        sim_w = sim_values .* factor;
        sigma = PF.weighting.sigma.discharge_mm_h;
    case 'groundwater_depth'
        obs_w = obs_values;
        sim_w = sim_values;
        sigma = PF.weighting.sigma.groundwater_depth_m;
    case 'soil_moisture'
        obs_w = obs_values;
        sim_w = sim_values;
        sigma = PF.weighting.sigma.soil_moisture_m3m3;
    otherwise
        error('pf_paper_weights:type', 'Unknown observation type: %s', type);
end
sigma = max(sigma, realmin);
end
