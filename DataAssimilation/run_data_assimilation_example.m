function Results = run_data_assimilation_example(varargin)
%RUN_DATA_ASSIMILATION_EXAMPLE Entry point for the isolated DA prototype.

da_dir = fullfile(fileparts(mfilename('fullpath')), 'VTilted_ParticleFilter');
addpath(da_dir);
Results = run_vtilted_augmented_pf_full_hydropol2d(varargin{:});
end
