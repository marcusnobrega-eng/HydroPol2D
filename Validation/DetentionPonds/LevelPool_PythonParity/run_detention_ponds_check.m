% VAL-POND-001: detention ponds (level pools) against HydroPol2D-Python, and their
% water balance in full HydroPol2D runs.
%
% 1. Kernel parity (VAL-POND-001K). The case of make_python_reference.py
%    (kernel_case.json: two ponds, DTM and power law, with an orifice, a spillway, a
%    pump and a disabled orifice) runs through the MATLAB kernels for 600 steps of
%    10 s and is compared, step by step, with the HydroPol2D-Python kernel
%    (python_reference_steps.csv): storage, stage, rain, evaporation, seepage and
%    recharge, critical inflow, device discharge and pump switches.
% 2. Full model (VAL-POND-001G and VAL-POND-001N). run_detention_ponds_full_model.m
%    runs the normal preprocessing and main loop in bypass mode with two ponds,
%    with groundwater on (seepage recharges the aquifer) and off (seepage leaves the
%    model). The event mass ledger, with the pond storage and pond losses, must close;
%    so must each pond's own balance.
%
% Writes Outputs/Validation/Pass_Fail.csv (case_id, status, report_ready, ...),
% Kernel_Parity.csv and FullModel_Summary.csv.

clear; clc;

% 1. Kernel parity with HydroPol2D-Python (writes Kernel_Parity.csv).
case_dir = fileparts(mfilename('fullpath'));
run(fullfile(case_dir, 'run_detention_ponds_kernel_check.m'));
case_dir = fileparts(mfilename('fullpath'));
summary_dir = fullfile(case_dir, 'Outputs', 'Validation');

% 2. Full-model runs (preprocessing clears the workspace except listed variables,
% among them Cases, icase, ValidationCase, case_dir and summary_dir).
full_summary = fullfile(summary_dir, 'FullModel_Summary.csv');
if isfile(full_summary); delete(full_summary); end
Cases = {'groundwater', 'no_groundwater'};
for icase = 1:numel(Cases)
    ValidationCase = struct('variant', Cases{icase});
    run(fullfile(case_dir, 'run_detention_ponds_full_model.m'));
end

% 3. Pass/Fail.
Kernel = readtable(fullfile(summary_dir, 'Kernel_Parity.csv'), 'TextType', 'string');
Full = readtable(fullfile(summary_dir, 'FullModel_Summary.csv'), 'TextType', 'string');
ids = ["VAL-POND-001K"; "VAL-POND-001G"; "VAL-POND-001N"];
passed = [logical(Kernel.passed(1)); ...
    logical(Full.passed(Full.variant == "groundwater")); ...
    logical(Full.passed(Full.variant == "no_groundwater"))];
metric = [Kernel.max_volume_error_m3(1); ...
    Full.system_residual_m3(Full.variant == "groundwater"); ...
    Full.system_residual_m3(Full.variant == "no_groundwater")];
notes = [ ...
    sprintf("MATLAB vs Python kernel (commit %s): max volume error %.2e m3, stage %.2e m, discharge %.2e m3/s, %d pump-switch mismatches.", ...
        extractBefore(Kernel.python_commit(1), 9), Kernel.max_volume_error_m3(1), ...
        Kernel.max_stage_error_m(1), Kernel.max_discharge_error_m3_s(1), Kernel.pump_switch_mismatches(1)); ...
    sprintf("Full run, groundwater on: ledger residual %.2e m3 (%.1f m3 without the pond terms); pond seepage %.3f m3 to the aquifer.", ...
        Full.system_residual_m3(Full.variant == "groundwater"), ...
        Full.residual_without_pond_terms_m3(Full.variant == "groundwater"), ...
        Full.pond_seepage_m3(Full.variant == "groundwater")); ...
    sprintf("Full run, groundwater off: ledger residual %.2e m3 (%.1f m3 without the pond terms); pond seepage %.3f m3 out of the model.", ...
        Full.system_residual_m3(Full.variant == "no_groundwater"), ...
        Full.residual_without_pond_terms_m3(Full.variant == "no_groundwater"), ...
        Full.pond_seepage_out_m3(Full.variant == "no_groundwater"))];
status = repmat("FAIL", 3, 1);
status(passed) = "PASS";
Pass = table(ids, status, double(passed), double(passed), metric, notes, ...
    'VariableNames', {'case_id', 'status', 'passed', 'report_ready', 'metric', 'notes'});
writetable(Pass, fullfile(summary_dir, 'Pass_Fail.csv'));
disp(Pass);
