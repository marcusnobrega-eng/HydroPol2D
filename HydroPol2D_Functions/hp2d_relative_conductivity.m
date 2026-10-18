function kr = hp2d_relative_conductivity(column, theta, l_vg)
%HP2D_RELATIVE_CONDUCTIVITY Mualem-van Genuchten relative conductivity.
if nargin < 3, l_vg = 0.5; end
m = 1 - 1./column.n_vg;
se = (theta - column.theta_r)./max(column.theta_s - column.theta_r, 1e-12);
se = min(max(se, 1e-9), 1);
inner = min(max(1 - se.^(1./m), 0), 1);
kr = min(max(se.^l_vg.*(1 - inner.^m).^2, 0), 1);
end
