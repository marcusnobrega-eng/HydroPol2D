% Run all Phase 1 snow checks, including the normal HydroPol2D V-tilted run.
%
% The full-model script is intentionally last because regular preprocessing
% clears the script workspace as part of the standard model workflow.

case_dir = fileparts(mfilename('fullpath'));
run(fullfile(case_dir, 'run_snow_model.m'));
run(fullfile(case_dir, 'run_snow_full_model_validation.m'));
