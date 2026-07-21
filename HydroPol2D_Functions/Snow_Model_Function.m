function [SWE_t, H_snow_t, M_snow, P_snow, P_rain, rho_snow, E_s, mass_balance_error_t] = ...
    Snow_Model_Function(SWE_prev, H_snow_prev, rho_prev, T_air, T_min, P, wind, ...
    lat, DOY, alpha, epsilon, C_e, DDF, T_snow_all, T_rain_all, ...
    rho_snow_init, rho_max, k_t, k_swe, k_D, dt_s)
%SNOW_MODEL_FUNCTION Advance the HydroPol2D snow store by one model step.
%
% All water terms are in mm water equivalent. Snow depth is in mm physical
% snow depth and density is in kg m^-3. Rate parameters are expressed per
% day and are scaled here by the actual routing time step.

if nargin ~= 21
    error('Snow_Model_Function expects 21 inputs, including rho_prev and dt_s.');
end
if ~isscalar(dt_s) || ~isfinite(dt_s) || dt_s <= 0
    error('Snow model time step dt_s must be a positive scalar in seconds.');
end

dt_day = dt_s / 86400;
rho_water = 1000; % kg m^-3, equal to mm water per m2

% A linear mixed-phase transition is configured per land-cover class.
valid_thresholds = isfinite(T_rain_all) & isfinite(T_snow_all);
if any(T_rain_all(valid_thresholds) < T_snow_all(valid_thresholds))
    error('Snow_Model_Function requires T_rain_all >= T_snow_all in every valid cell.');
end
snow_fraction = precipitation_snow_fraction(T_air, T_snow_all, T_rain_all);
P_snow = max(P, 0) .* snow_fraction;
P_rain = max(P, 0) - P_snow;

% Existing and new snow are mixed before compaction. An empty snowpack uses
% the LULC-specific fresh-snow density as its dormant state value.
SWE_after_snow = max(SWE_prev, 0) + P_snow;
rho_prev = max(rho_prev, rho_snow_init);
rho_prev = min(rho_prev, rho_max);
rho_mixed = rho_snow_init;
has_snow = SWE_after_snow > 0;
rho_mixed(has_snow) = (rho_prev(has_snow) .* max(SWE_prev(has_snow), 0) + ...
    rho_snow_init(has_snow) .* P_snow(has_snow)) ./ SWE_after_snow(has_snow);

% The radiation term is an energy flux [W m^-2]. Multiplying by seconds and
% dividing by latent heat [J kg^-1] gives kg m^-2, numerically equal to mm.
[S_down, L_down] = compute_radiation(T_air, lat, DOY);
Q_net = compute_net_radiation(S_down, L_down, T_air, alpha, epsilon);
latent_heat_fusion = 3.34e5; % J kg^-1
melt_degree_day = DDF .* max(T_air, 0) .* dt_day;
melt_radiation = max(Q_net, 0) .* dt_s / latent_heat_fusion;
M_snow = min(max(melt_degree_day + melt_radiation, 0), SWE_after_snow);

SWE_after_melt = SWE_after_snow - M_snow;
E_s = compute_sublimation(wind, T_air, T_min, SWE_after_melt, C_e, dt_day);
SWE_t = max(SWE_after_melt - E_s, 0);

% Compaction is a rate process. Use the carried density rather than resetting
% density to its initial value at every routing step.
H_before_compaction = zeros(size(SWE_after_snow), 'like', SWE_after_snow);
H_before_compaction(has_snow) = SWE_after_snow(has_snow) .* rho_water ./ rho_mixed(has_snow);
rho_snow = rho_mixed + dt_day .* ( ...
    k_t .* max(T_air, 0) + ...
    k_swe .* SWE_after_snow + ...
    k_D .* H_before_compaction);
rho_snow = min(max(rho_snow, rho_snow_init), rho_max);
rho_snow(SWE_t <= 0) = rho_snow_init(SWE_t <= 0);

H_snow_t = zeros(size(SWE_t), 'like', SWE_t);
wet_snow = SWE_t > 0;
H_snow_t(wet_snow) = SWE_t(wet_snow) .* rho_water ./ rho_snow(wet_snow);

% Preserve domain nodata consistently across all outputs.
invalid = ~isfinite(SWE_prev) | ~isfinite(H_snow_prev) | ~isfinite(rho_prev) | ...
    ~isfinite(T_air) | ~isfinite(T_min) | ~isfinite(P) | ~isfinite(wind) | ...
    ~isfinite(lat) | ~isfinite(alpha) | ~isfinite(epsilon) | ~isfinite(C_e) | ...
    ~isfinite(DDF) | ~isfinite(T_snow_all) | ~isfinite(T_rain_all) | ...
    ~isfinite(rho_snow_init) | ~isfinite(rho_max) | ~isfinite(k_t) | ...
    ~isfinite(k_swe) | ~isfinite(k_D);
SWE_t(invalid) = NaN;
H_snow_t(invalid) = NaN;
M_snow(invalid) = NaN;
P_snow(invalid) = NaN;
P_rain(invalid) = NaN;
rho_snow(invalid) = NaN;
E_s(invalid) = NaN;

mass_balance_error_t = nansum((SWE_prev(:) + P_snow(:)) - ...
    (M_snow(:) + E_s(:) + SWE_t(:)));
end

function snow_fraction = precipitation_snow_fraction(T_air, T_snow_all, T_rain_all)
snow_fraction = zeros(size(T_air), 'like', T_air);
step_threshold = abs(T_rain_all - T_snow_all) <= eps(max(abs(T_rain_all), 1));
snow_fraction(step_threshold) = T_air(step_threshold) <= T_snow_all(step_threshold);

transition = ~step_threshold;
snow_fraction(transition) = (T_rain_all(transition) - T_air(transition)) ./ ...
    (T_rain_all(transition) - T_snow_all(transition));
snow_fraction(T_air <= T_snow_all) = 1;
snow_fraction(T_air >= T_rain_all) = 0;
snow_fraction = min(max(snow_fraction, 0), 1);
end

function [S_down, L_down] = compute_radiation(T_air, lat, DOY)
sigma = 5.67e-8; % W m^-2 K^-4
T_kelvin = T_air + 273.15;
decl = 23.45 * sind(360 * (DOY - 81) / 365);
lat_rad = deg2rad(lat);
decl_rad = deg2rad(decl);
cos_theta = sin(lat_rad) .* sin(decl_rad) + cos(lat_rad) .* cos(decl_rad);
theta_z = acos(min(max(cos_theta, -1), 1));
S_down = 0.75 * 1367 .* max(cos(theta_z), 0);
epsilon_a = 0.642 + 0.035 * sqrt(max(T_kelvin - 273.15, 0));
L_down = epsilon_a .* sigma .* T_kelvin.^4;
end

function Q_net = compute_net_radiation(S_down, L_down, T_air, alpha, epsilon)
sigma = 5.67e-8;
T_snow_kelvin = 273.15;
L_up = epsilon .* sigma .* T_snow_kelvin^4;
Q_net = S_down .* (1 - alpha) + L_down - L_up - ...
    epsilon .* sigma .* (T_snow_kelvin^4 - (T_air + 273.15).^4);
end

function E_s = compute_sublimation(wind, T_air, T_min, SWE_available, C_e, dt_day)
e_s = 6.112 .* exp((17.67 .* T_air) ./ (T_air + 243.5));
e_a = 10 .* (0.61 .* exp((17.27 .* T_min) ./ (T_min + 237.3)));
q_s = 0.622 .* e_s ./ (1013.25 - e_s);
q_a = 0.622 .* e_a ./ (1013.25 - e_a);
E_s = C_e .* max(wind, 0) .* max(q_s - q_a, 0) .* max(SWE_available, 0) .* dt_day;
E_s = min(E_s, SWE_available);
end
