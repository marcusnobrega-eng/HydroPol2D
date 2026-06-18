run('/Users/mngomes/Downloads/input_data_bypass_script.m');

if ~isfield(InputData_Bypass.flags, 'flag_full_momentum')
    InputData_Bypass.flags.flag_full_momentum = 0;
end

routing_mode = lower(strtrim(string(getenv('HYDROPOL2D_VTILT_ROUTING'))));
if strlength(routing_mode) == 0
    routing_mode = "local_inertial";
end

InputData_Bypass.flags.flag_full_momentum = 0;
InputData_Bypass.flags.flag_inertial = 0;
InputData_Bypass.flags.flag_CA = 0;
InputData_Bypass.flags.flag_kinematic = 0;
InputData_Bypass.flags.flag_diffusive = 0;

switch routing_mode
    case {"full_momentum", "fullmomentum", "fm"}
        InputData_Bypass.flags.flag_full_momentum = 1;
    case {"local_inertial", "localinertial", "li"}
        InputData_Bypass.flags.flag_inertial = 1;
    case {"cellular_automata", "ca"}
        InputData_Bypass.flags.flag_CA = 1;
    case {"kinematic", "kinematic_wave", "kw"}
        InputData_Bypass.flags.flag_kinematic = 1;
    case {"diffusive", "diffusive_wave", "dw"}
        InputData_Bypass.flags.flag_diffusive = 1;
    otherwise
        error('Unsupported HYDROPOL2D_VTILT_ROUTING mode: %s', routing_mode);
end
