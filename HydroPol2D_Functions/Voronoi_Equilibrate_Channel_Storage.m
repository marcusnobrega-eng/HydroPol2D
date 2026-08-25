function [surface_volume, channel_volume, channel_stage, transfer_to_channel_m3] = Voronoi_Equilibrate_Channel_Storage(mesh, surface_volume, channel_volume)
%VORONOI_EQUILIBRATE_CHANNEL_STORAGE Share stage without double-counting volume.

channel = mesh.channel;
channel_stage = zeros(channel.n_nodes, 1);
channel_before = channel_volume;
for node = 1:channel.n_nodes
    cell_id = channel.host_cell(node);
    total = surface_volume(cell_id) + channel_volume(node);
    channel_area = channel.plan_area(node);
    floodplain_area = mesh.surface_area(cell_id);
    bank_depth = max(channel.bank(node) - channel.bed(node), 0);
    bank_volume = channel_area * bank_depth;
    if total <= bank_volume
        stage = channel.bed(node) + total / channel_area;
        channel_volume(node) = total;
        surface_volume(cell_id) = 0;
    else
        stage = channel.bank(node) + (total - bank_volume) / (channel_area + floodplain_area);
        channel_volume(node) = bank_volume + channel_area * (stage - channel.bank(node));
        surface_volume(cell_id) = floodplain_area * (stage - channel.bank(node));
    end
    channel_stage(node) = stage;
end
transfer_to_channel_m3 = sum(channel_volume - channel_before);
end
