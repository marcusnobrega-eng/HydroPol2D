% HP2D_Wrapper.m
% HydroPol2D Sherlock batch wrapper
%
% Assumptions:
%   - this wrapper lives in the model root folder
%   - the main script is always HydroPol2D_V115.m
%   - results are always written to ./Modeling_Results

clc;
clear;

fprintf('============================================================\n');
fprintf('HydroPol2D batch execution started\n');
fprintf('Start time        : %s\n', datestr(now));
fprintf('============================================================\n');

%% Fixed paths based on current folder
project_root     = pwd;
main_script      = 'HydroPol2D_V115.m';
main_script_path = fullfile(project_root, main_script);
results_dir      = fullfile(project_root, 'Modeling_Results');

%% Check main script exists
if ~isfile(main_script_path)
    error('Main MATLAB script not found:\n%s', main_script_path);
end

%% Make sure results folder exists
if ~isfolder(results_dir)
    mkdir(results_dir);
end

%% Set working folder and paths
cd(project_root);
addpath(genpath(project_root));

%% Headless-safe plotting
set(0, 'DefaultFigureVisible', 'off');

%% Optional CUDA forward compatibility
try
    parallel.gpu.enableCUDAForwardCompatibility(true);
    fprintf('CUDA forward compatibility enabled.\n');
catch ME_fc
    fprintf('Could not enable CUDA forward compatibility: %s\n', ME_fc.message);
end

%% Diagnostics
fprintf('Working directory : %s\n', pwd);
fprintf('Project root      : %s\n', project_root);
fprintf('Main script       : %s\n', main_script);
fprintf('Results directory : %s\n', results_dir);
fprintf('Host name         : %s\n', getenv('HOSTNAME'));
fprintf('SLURM Job ID      : %s\n', getenv('SLURM_JOB_ID'));
fprintf('CPUs allocated    : %s\n', getenv('SLURM_CPUS_PER_TASK'));
fprintf('GPUs visible      : %s\n', getenv('CUDA_VISIBLE_DEVICES'));
fprintf('MATLAB version    : %s\n', version);
fprintf('============================================================\n');

%% Put results directory in base workspace so HydroPol2D_V115.m can use it
assignin('base', 'exportDirectory', results_dir);
assignin('base', 'project_root', project_root);

%% Run main file
try
    fprintf('Running main script: %s\n', main_script);

    run(main_script_path);

    fprintf('============================================================\n');
    fprintf('HydroPol2D finished successfully\n');
    fprintf('Finish time       : %s\n', datestr(now));
    fprintf('============================================================\n');

    exit(0);

catch ME
    fprintf(2, '============================================================\n');
    fprintf(2, 'HydroPol2D execution FAILED\n');
    fprintf(2, 'Finish time       : %s\n', datestr(now));
    fprintf(2, 'Error message     : %s\n', ME.message);
    fprintf(2, '------------------------------------------------------------\n');
    fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    fprintf(2, '============================================================\n');

    exit(1);
end