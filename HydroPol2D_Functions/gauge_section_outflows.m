function q = gauge_section_outflows(flow_rate, flags)
%GAUGE_SECTION_OUTFLOWS Per-direction cell outflows (mm/h) of the last routing step.
%   For hp2d_gauge_section_discharge; the diagonal ones only with D8 routing.
%   signed_faces: the D4 non-CA solvers (Local_Inertial_Model_D4, Routing_Wave_Model_D4,
%   Full_Momentum_Model_D4) return signed face fluxes; CA_Routing, CA_Routing_8D and
%   Bates_Inertial_8D return one-sided outflows (the solver choice of HydroPol2D_Main_While).
q = struct('left', [], 'right', [], 'up', [], 'down', []);
d8 = isfield(flags, 'flag_D8') && flags.flag_D8 == 1;
ca = isfield(flags, 'flag_CA') && flags.flag_CA == 1 && ~(isfield(flags, 'flag_inertial') && flags.flag_inertial == 1);
q.signed_faces = ~d8 && ~ca;
if isfield(flow_rate, 'qout_left_t'), q.left = flow_rate.qout_left_t; end
if isfield(flow_rate, 'qout_right_t'), q.right = flow_rate.qout_right_t; end
if isfield(flow_rate, 'qout_up_t'), q.up = flow_rate.qout_up_t; end
if isfield(flow_rate, 'qout_down_t'), q.down = flow_rate.qout_down_t; end
if d8
    for name = {'ne', 'se', 'sw', 'nw'}
        field = ['qout_' name{1} '_t'];
        if isfield(flow_rate, field), q.(name{1}) = flow_rate.(field); end
    end
end
end
