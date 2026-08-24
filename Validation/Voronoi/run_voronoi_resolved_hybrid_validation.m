function summary = run_voronoi_resolved_hybrid_validation(fine_file, hybrid_file)
%RUN_VORONOI_RESOLVED_HYBRID_VALIDATION Fine resolved vs coarse Neal flow.

fine = HydroPol2D_Read_UGRID(fine_file); hybrid = HydroPol2D_Read_UGRID(hybrid_file);
assert(fine.channel.n_nodes == 0 && hybrid.channel.n_nodes > 0);
river = isfinite(fine.cell_roughness);
depth = zeros(fine.n_cells,1); depth(river) = 1.5;
slope = 1e-3; width = 30; roughness = 0.035;
area = width * 1.5; radius = area / (width + 3);
Q = area / roughness * radius^(2/3) * sqrt(slope);

internal = fine.edge_neighbor > 0;
river_edge = internal & river(fine.edge_owner) & river(max(fine.edge_neighbor,1));
initial_q = zeros(fine.n_edges,1);
initial_q(river_edge) = (Q/width) .* fine.edge_normal_x(river_edge);
boundary = find(fine.edge_neighbor == 0);
inlet = boundary(abs(fine.edge_normal_x(boundary)+1) < 1e-8 & river(fine.edge_owner(boundary)));
outlet = boundary(abs(fine.edge_normal_x(boundary)-1) < 1e-8 & river(fine.edge_owner(boundary)));
assert(~isempty(inlet) && ~isempty(outlet));
inlet_share = fine.edge_length(inlet) / sum(fine.edge_length(inlet));
fine_forcing.surface_boundary = struct('edge_id',[inlet;outlet], ...
    'type',[repmat("inflow",numel(inlet),1);repmat("normal_flow",numel(outlet),1)], ...
    'value',[Q*inlet_share;repmat(slope,numel(outlet),1)]);
fine_config = struct('duration_s',30,'initial_surface_depth_m',depth, ...
    'initial_surface_discharge_per_width_m2_s',initial_q,'max_dt_s',0.25, ...
    'min_dt_s',1e-5,'output_interval_s',10,'surface_roughness',roughness);
fine_result = HydroPol2D_Voronoi_Run(fine_file,fine_config,fine_forcing);
section_x = 500;
crossing = internal & river(fine.edge_owner) & river(max(fine.edge_neighbor,1)) & ...
    (fine.cell_x(fine.edge_owner)-section_x).*(fine.cell_x(max(fine.edge_neighbor,1))-section_x) <= 0;
fine_outlet_q = sum(fine_result.edge_discharge_per_width_m2_s(crossing).* ...
    fine.edge_length(crossing).*fine.edge_normal_x(crossing));

channel = hybrid.channel;
channel_inlet = setdiff((1:channel.n_nodes)',channel.link_down);
channel_outlet = setdiff((1:channel.n_nodes)',channel.link_up);
hybrid_forcing.channel_boundary = struct('node_id',[channel_inlet;channel_outlet], ...
    'type',["inflow";"normal_flow"],'value',[Q;slope], ...
    'width_m',[width;width],'roughness',[roughness;roughness]);
hybrid_config = struct('duration_s',30,'initial_channel_depth_m',1.5, ...
    'initial_channel_discharge_m3_s',Q,'max_dt_s',0.25,'min_dt_s',1e-5,'output_interval_s',10);
hybrid_result = HydroPol2D_Voronoi_Run(hybrid_file,hybrid_config,hybrid_forcing);
hybrid_link_midpoint = 0.5*(hybrid.cell_x(channel.host_cell(channel.link_up)) + hybrid.cell_x(channel.host_cell(channel.link_down)));
[~,middle_link] = min(abs(hybrid_link_midpoint-section_x));
hybrid_outlet_q = hybrid_result.channel_discharge_m3_s(middle_link);

[~,mid] = min(abs(hybrid.cell_x(channel.host_cell)-500));
comparison_x = hybrid.cell_x(channel.host_cell(mid));
fine_mid = river & abs(fine.cell_x-comparison_x) < 15;
fine_stage = mean(fine.surface_bed(fine_mid)+fine_result.final_surface_volume_m3(fine_mid)./fine.surface_area(fine_mid));
hybrid_stage = channel.bed(mid)+hybrid_result.final_channel_volume_m3(mid)/channel.plan_area(mid);
discharge_error = abs(hybrid_outlet_q-fine_outlet_q)/max(abs(fine_outlet_q),eps);
level_error_fraction_bankfull = abs(hybrid_stage-fine_stage)/2;
fine_volume = sum(fine_result.final_surface_volume_m3);
hybrid_volume = sum(hybrid_result.final_surface_volume_m3) + sum(hybrid_result.final_channel_volume_m3);
volume_error = abs(hybrid_volume-fine_volume)/fine_volume;
fine_inundated_area = sum(fine.cell_area(fine_result.final_surface_volume_m3./fine.surface_area > 1e-6));
hybrid_inundated_area = sum(channel.plan_area(hybrid_result.final_channel_volume_m3./channel.plan_area > 1e-6)) + ...
    sum(hybrid.surface_area(hybrid_result.final_surface_volume_m3./hybrid.surface_area > 1e-6));
inundated_area_error = abs(hybrid_inundated_area-fine_inundated_area)/fine_inundated_area;
fprintf('fine Q %.12g, hybrid Q %.12g, fine stage %.12g, hybrid stage %.12g\n', ...
    fine_outlet_q,hybrid_outlet_q,fine_stage,hybrid_stage);
assert(discharge_error <= 0.05,'Hybrid discharge differs from fine resolved reference by %.3g.',discharge_error);
assert(level_error_fraction_bankfull <= 0.05,'Hybrid level differs by %.3g bankfull depths.',level_error_fraction_bankfull);
assert(volume_error <= 0.02,'Hybrid volume differs from fine resolved reference by %.3g.',volume_error);
assert(inundated_area_error <= 0.05,'Hybrid inundated area differs by %.3g.',inundated_area_error);
summary = struct('fine_outlet_discharge_m3_s',fine_outlet_q, ...
    'hybrid_outlet_discharge_m3_s',hybrid_outlet_q,'relative_discharge_error',discharge_error, ...
    'water_level_error_fraction_bankfull',level_error_fraction_bankfull, ...
    'relative_volume_error',volume_error,'relative_inundated_area_error',inundated_area_error);
end
