function speed = hp2d_cell_centre_velocity(flow_rate, depth_mm, resolution, idx_nan)
%HP2D_CELL_CENTRE_VELOCITY Cell-centre speed from the D4 face discharges [m/s].
%   The discharge per width at a cell's centre is the mean of the net
%   discharges through its two opposite faces in each direction, divided by
%   the face width; the speed is its magnitude over the cell depth, in cells
%   deeper than 1 mm (0 elsewhere). This is _cell_velocity_m_s in
%   HydroPol2D-Python, which both human-risk modes use.
%   flow_rate.qout_left_t/right_t/up_t/down_t are the cell outflows [mm/h]
%   through each face; depth_mm is the cell depth [mm].

to_m3_s = resolution^2 / 1000 / 3600;
left = flow_rate.qout_left_t * to_m3_s; right = flow_rate.qout_right_t * to_m3_s;
up = flow_rate.qout_up_t * to_m3_s; down = flow_rate.qout_down_t * to_m3_s;
left(~isfinite(left)) = 0; right(~isfinite(right)) = 0;
up(~isfinite(up)) = 0; down(~isfinite(down)) = 0;
% Net eastward discharge through every vertical face, and net southward
% discharge through every horizontal face (row 1 is the north).
qx = [-left(:,1), right(:,1:end-1) - left(:,2:end), right(:,end)];
qy = [-up(1,:); down(1:end-1,:) - up(2:end,:); down(end,:)];
hu = 0.5 * (qx(:,1:end-1) + qx(:,2:end)) / resolution;
hv = 0.5 * (qy(1:end-1,:) + qy(2:end,:)) / resolution;
depth = depth_mm / 1000;
speed = zeros(size(depth), 'like', depth);
wet = depth > 1e-3;
speed(wet) = hypot(hu(wet), hv(wet)) ./ depth(wet);
if nargin >= 4 && ~isempty(idx_nan), speed(idx_nan) = 0; end
end
