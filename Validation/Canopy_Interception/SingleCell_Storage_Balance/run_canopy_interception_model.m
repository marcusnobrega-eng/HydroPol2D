% P1-CANOPY-001: run HydroPol2D interceptionModel against the Phase 1 reference.
%
% Run from the repository root or from this case folder in MATLAB:
%   run('HydroPol2D_Model/Validation/Canopy_Interception/SingleCell_Storage_Balance/run_canopy_interception_model.m')

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Model', 'HydroPol2D_Functions');
reference_file = fullfile(repo_root, 'HydroPol2D_Model', 'Validation', ...
    'Reference_Outputs', 'Phase1', 'canopy_bucket', ...
    'P1-CANOPY-001_reference.csv');
output_dir = fullfile(case_dir, 'Outputs', 'Validation');

if ~exist(functions_dir, 'dir')
    error('HydroPol2D functions directory not found: %s', functions_dir);
end
if ~exist(reference_file, 'file')
    error('Reference file not found. Generate it with phase1_reference_solutions.py: %s', reference_file);
end
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

addpath(functions_dir);

ref = readtable(reference_file, 'TextType', 'string');
model_rows = table();

scenarios = unique(ref.scenario, 'stable');
for s = 1:numel(scenarios)
    scenario = scenarios(s);
    idx = ref.scenario == scenario;
    ref_s = ref(idx, :);

    storage_prev = 0;
    for i = 1:height(ref_s)
        P = ref_s.gross_rainfall_mm(i);
        Ep = ref_s.potential_evaporation_mm(i);
        LAI = ref_s.lai(i);
        C = ref_s.coefficient_mm_per_lai(i);

        [S, T, E, F, mb_error] = interceptionModel(P, Ep, LAI, storage_prev, C, 0);

        next_row = table( ...
            scenario, ...
            ref_s.time_s(i), ...
            P, ...
            Ep, ...
            LAI, ...
            C * LAI, ...
            S, ...
            T, ...
            E, ...
            F, ...
            mb_error, ...
            'VariableNames', {'scenario','time_s','gross_rainfall_mm', ...
            'potential_evaporation_mm','lai','smax_mm','canopy_storage_mm', ...
            'throughfall_mm','evaporation_mm','stemflow_mm','mass_residual_mm'});
        model_rows = [model_rows; next_row]; %#ok<AGROW>

        storage_prev = S;
    end
end

writetable(model_rows, fullfile(output_dir, 'HydroPol2D_Canopy_Model.csv'));

disp('Wrote HydroPol2D canopy model output to:');
disp(fullfile(output_dir, 'HydroPol2D_Canopy_Model.csv'));
disp('Run compare_canopy_interception.py to generate metrics and pass/fail diagnostics.');
