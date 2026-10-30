function speed_m_s = hp2d_d8_cell_velocity(volume_m3, discharge_m3_s, D8)
%HP2D_D8_CELL_VELOCITY Fastest D8 link velocity at each cell [m/s].
%   A link's velocity is |Q| / A with its current flow area; each cell takes
%   the maximum over the links that start or end in it (0 elsewhere). This
%   is d8_river_cell_velocity_jax in HydroPol2D-Python. Returns ny x nx.

stage = hp2d_d8_cell_stage(D8.Cells, max(volume_m3(:), 0));
area = hp2d_d8_link_hydraulics(D8.Faces, D8.Links, stage);
speed = abs(discharge_m3_s(:)) ./ max(area, 1e-30);
speed_m_s = zeros(D8.N, 1, 'like', stage);
if ~isempty(speed)
    speed_m_s = max(accumarray(D8.Links.source, speed, [D8.N 1], @max), ...
        accumarray(D8.Links.target, speed, [D8.N 1], @max));
end
speed_m_s = reshape(speed_m_s, D8.ny, D8.nx);

end
