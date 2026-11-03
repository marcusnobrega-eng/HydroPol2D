%% HydroPol2D | former name of the main script
% The main script is now run_HydroPol2D.m. This file runs it, so existing cases,
% cluster jobs and Case Builder exports that call HydroPol2D_V115.m keep working.
run(fullfile(fileparts(mfilename('fullpath')), 'run_HydroPol2D.m'));
