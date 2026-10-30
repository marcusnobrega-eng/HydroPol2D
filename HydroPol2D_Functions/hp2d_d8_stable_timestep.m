function dt_s = hp2d_d8_stable_timestep(volume_m3, D8, courant, max_dt_s)
%HP2D_D8_STABLE_TIMESTEP Link-length CFL time step of the D8 links.
%   dt = courant * min over wet links of L / sqrt(g h), with L the link
%   length and h the link depth (higher cell stage minus the link datum);
%   max_dt_s when no link is wet. This is the Wflow/Bates wave-celerity
%   limit, estimate_d8_river_timestep_jax in HydroPol2D-Python.

stage = hp2d_d8_cell_stage(D8.Cells, max(volume_m3(:), 0));
edge_depth = max(max(stage(D8.Links.source), stage(D8.Links.target)) - D8.Faces.datum, 0);
wet = edge_depth > 0;
if ~any(wet)
    dt_s = max_dt_s;
    return
end
candidate = courant * D8.Links.length(wet) ./ sqrt(9.81 * edge_depth(wet));
dt_s = min(min(candidate), max_dt_s);

end
