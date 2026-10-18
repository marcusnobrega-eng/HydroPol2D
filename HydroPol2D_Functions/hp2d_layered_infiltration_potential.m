function rate_mm_h = hp2d_layered_infiltration_potential(V, available_mm, theta_r, theta_s, alpha, n_vg, ksat_mm_h, dt_s, l_vg, ltop_m, dh_max_m, ksurf_weight, valid)
%HP2D_LAYERED_INFILTRATION_POTENTIAL Bucket-scale Darcy intake rate (mm/h) on a soil view.
%   Same formulation as Infiltration_Module (and HydroPol2D-Python infiltration_potential):
%   C = Ksoil*[min(h0 - h(theta_bucket), dh_max)/Ltop + 1], Ksoil the weighted interface
%   conductivity between Ks and K(theta_top); limited by supply and remaining storage.
available = max(available_mm(:), 0);
alpha = max(alpha, 1e-12);
m = 1 - 1./n_vg;
ks = max(ksat_mm_h(:), 0)/1000/3600;
ltop = max(ltop_m + 0*available, 1e-6);
ltop = max(min(ltop, max(V.near_surface_thickness_m, 1e-6)), 1e-6);
dh_max = max(dh_max_m, 0);
wk = min(max(ksurf_weight, 0), 1);
dt_h = dt_s/3600;
supply = available/dt_h;
remaining_mm = max(V.total_capacity_mm - sum(V.water_mm, 2), 0);
storage_limit = remaining_mm/dt_h;
pore = max(V.pervious_fraction, 0);
near_thick = max(V.near_surface_thickness_m, 1e-6);
near_m = max(V.water_mm(:, 1), 0)/1000;
theta_bucket = min(max(theta_r + near_m./max(pore.*near_thick, 1e-9), theta_r), theta_s);
se = min(max((theta_bucket - theta_r)./max(theta_s - theta_r, 1e-12), 1e-6), 1);
psi = -(1./alpha).*((se.^(-1./m) - 1).^(1./n_vg));
psi(se >= 0.999999) = 0;
available_m = available/1000;
near_deficit_m = max(V.capacity_mm(:, 1) - V.water_mm(:, 1), 0)/1000;
wetting = min(min(available_m, remaining_mm/1000), near_deficit_m);
theta_top = min(max(theta_bucket + wetting./max(pore.*ltop, 1e-12), theta_r), theta_s);
se_top = min(max((theta_top - theta_r)./max(theta_s - theta_r, 1e-12), 1e-6), 1);
term = min(max(1 - se_top.^(1./m), 0), 1);
kr = min(max(se_top.^l_vg.*(1 - term.^m).^2, 0), 1);
k_top = ks.*kr;
wet = supply > 0 | available_m > 0;
conductivity = k_top;
blend = wk.*ks + (1 - wk).*k_top;
conductivity(wet) = blend(wet);
conductivity = min(max(conductivity, 0), ks);
gradient = min(max(available_m - psi, 0), dh_max)./ltop + 1;
capacity = conductivity.*gradient*1000*3600.*pore;
rate_mm_h = min(min(supply, capacity), storage_limit);
rate_mm_h(pore <= 1e-9 | V.vadose_depth_m <= 0 | V.total_capacity_mm <= 0) = 0;
rate_mm_h(~isfinite(rate_mm_h)) = 0;
rate_mm_h = max(rate_mm_h, 0);
rate_mm_h(~valid) = NaN;
end
