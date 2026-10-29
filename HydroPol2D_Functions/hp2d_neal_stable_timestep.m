function dt_s = hp2d_neal_stable_timestep(h, Qc, Qf, River_Width, River_Depth, Resolution, courant, max_dt_s, idx_nan)
%HP2D_NEAL_STABLE_TIMESTEP CFL time step for the Neal channel-subgrid scheme [s].
%   h is the representative depth above the channel bed [m]. Qc and Qf are
%   the channel and floodplain discharges on the east (:,:,1) and north
%   (:,:,2) faces [m3/s]. The cell discharge is the mean of its two faces in
%   each direction. The velocity divides it by the wetted section: channel
%   width x depth in bank, cell width x depth over the banks. Then
%       dt = courant * dx / max(|v| + sqrt(g h))
%   over the valid cells, capped at max_dt_s (estimate_neal_subgrid_dt_s in
%   HydroPol2D-Python).

g = 9.81;
[ny,nx] = size(h);
h = max(h, 0);
h(~isfinite(h)) = 0;

Q = Qc(:,:,1:2) + Qf(:,:,1:2);
Q(~isfinite(Q)) = 0;
q_x = 0.5 * (Q(:,:,1) + [zeros(ny,1,'like',h), Q(:,1:end-1,1)]);
q_y = 0.5 * (Q(:,:,2) + [Q(2:end,:,2); zeros(1,nx,'like',h)]);

area = hp2d_neal_cell_area(h, River_Width, River_Depth, Resolution);
velocity = hypot(q_x, q_y) .* Resolution ./ max(area .* h, 1e-12);
speed = velocity + sqrt(g .* h);
speed(idx_nan | ~isfinite(speed)) = 0;

max_speed = max(speed(:));
if max_speed <= 1e-12
    dt_s = max_dt_s;
else
    dt_s = min(max_dt_s, courant * Resolution / max_speed);
end

end
