function summary = run_voronoi_bidirectional_transition_validation(node_to_cell_file, cell_to_node_file)
%RUN_VORONOI_BIDIRECTIONAL_TRANSITION_VALIDATION Both orientations and reversal.

files = {node_to_cell_file, cell_to_node_file};
errors = zeros(2,2); discharges = zeros(2,2);
for fixture = 1:2
    mesh = HydroPol2D_Read_UGRID(files{fixture});
    assert(mesh.channel.n_transitions >= 1);
    expected_orientation = fixture == 1;
    assert(all(mesh.channel.transition_positive_node_to_cell == expected_orientation));
    node = mesh.channel.transition_node(:); cell_id = mesh.channel.transition_cell(:);
    for direction = 1:2
        surface = zeros(mesh.n_cells,1); channel = zeros(mesh.channel.n_nodes,1);
        if (direction == 1) == expected_orientation
            channel(unique(node)) = max(20 - mesh.channel.bed(unique(node)), 0) .* mesh.channel.plan_area(unique(node));
        else
            surface(cell_id) = max(20 - mesh.surface_bed(cell_id), 0) .* mesh.surface_area(cell_id);
        end
        initial = sum(surface) + sum(channel);
        [surface, channel, q, ~] = Voronoi_Channel_Transition_Step(mesh, surface, channel, zeros(mesh.channel.n_transitions,1), 0.01);
        final = sum(surface) + sum(channel);
        errors(fixture,direction) = abs(final - initial) / max(initial, eps);
        discharges(fixture,direction) = sum(q);
        assert(errors(fixture,direction) <= 1e-12);
        if direction == 1, assert(all(q > 0)); else, assert(all(q < 0)); end
    end
end
summary = struct('relative_mass_errors', errors, 'signed_transition_discharges_m3_s', discharges);
end
