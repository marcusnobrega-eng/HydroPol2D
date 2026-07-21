function run_snisb_hydropol2d_batch(varargin)
%RUN_SNISB_HYDROPOL2D_BATCH Loop HydroPol2D over prepared SNISB dam folders.
%
% Examples:
%   run_snisb_hydropol2d_batch
%   run_snisb_hydropol2d_batch('MaxDams', 2)
%   run_snisb_hydropol2d_batch('DamCodes', {'002151','020275'})
%   run_snisb_hydropol2d_batch('DryRun', true)
%   run_snisb_hydropol2d_batch('DamCodes', {'002151','020275'}, 'SimulationMinutes', 15)

p = inputParser;
addParameter(p, 'OutputsRoot', fullfile(fileparts(fileparts(mfilename('fullpath'))), 'outputs', 'dams'), ...
    @(x) ischar(x) || isstring(x));
addParameter(p, 'DamCodes', {}, @(x) iscell(x) || isstring(x) || isnumeric(x));
addParameter(p, 'MaxDams', inf, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'CleanOutputFolder', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'RunPostprocessing', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'EnableLogging', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'DryRun', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'SimulationMinutes', 2 * 24 * 60, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'RecordTimeMapsMinutes', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'RoutingModel', 'local_inertial', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});

outputs_root = char(p.Results.OutputsRoot);
dam_codes = normalize_dam_codes(p.Results.DamCodes);
max_dams = p.Results.MaxDams;

if exist(outputs_root, 'dir') ~= 7
    error('SNISB outputs root not found: %s', outputs_root);
end

dam_dirs = discover_dam_dirs(outputs_root, dam_codes);
if isempty(dam_dirs)
    error('No prepared dam folders found under: %s', outputs_root);
end

if isfinite(max_dams)
    dam_dirs = dam_dirs(1:min(numel(dam_dirs), max_dams));
end

fprintf('\nSNISB HydroPol2D batch will run %d dam(s).\n', numel(dam_dirs));
for i = 1:numel(dam_dirs)
    fprintf('  %2d. %s\n', i, dam_dirs{i});
end

failures = {};
for i = 1:numel(dam_dirs)
    dam_dir = dam_dirs{i};
    fprintf('\n============================================================\n');
    fprintf('SNISB HydroPol2D batch %d/%d\n', i, numel(dam_dirs));
    fprintf('============================================================\n');
    try
        run_snisb_hydropol2d_case( ...
            dam_dir, ...
            'CleanOutputFolder', logical(p.Results.CleanOutputFolder), ...
            'RunPostprocessing', logical(p.Results.RunPostprocessing), ...
            'EnableLogging', logical(p.Results.EnableLogging), ...
            'DryRun', logical(p.Results.DryRun), ...
            'SimulationMinutes', double(p.Results.SimulationMinutes), ...
            'RecordTimeMapsMinutes', double(p.Results.RecordTimeMapsMinutes), ...
            'RoutingModel', char(p.Results.RoutingModel));
    catch ME
        failures{end+1,1} = dam_dir; %#ok<AGROW>
        failures{end,2} = ME.message;
        warning('HydroPol2D failed for %s\n%s', dam_dir, ME.message);
    end
end

if ~isempty(failures)
    fprintf('\nHydroPol2D batch finished with %d failure(s):\n', size(failures,1));
    for i = 1:size(failures,1)
        fprintf('  - %s\n    %s\n', failures{i,1}, failures{i,2});
    end
    error('One or more HydroPol2D cases failed.');
end

fprintf('\nHydroPol2D batch completed successfully.\n');
end

function dam_dirs = discover_dam_dirs(outputs_root, dam_codes)
    entries = dir(fullfile(outputs_root, 'snisb_*'));
    dam_dirs = {};
    for i = 1:numel(entries)
        if ~entries(i).isdir
            continue;
        end
        dam_dir = fullfile(entries(i).folder, entries(i).name);
        if ~has_required_case_files(dam_dir)
            continue;
        end
        if ~isempty(dam_codes) && ~matches_any_code(entries(i).name, dam_codes)
            continue;
        end
        dam_dirs{end+1,1} = dam_dir; %#ok<AGROW>
    end
    dam_dirs = sort(dam_dirs);
end

function ok = has_required_case_files(dam_dir)
    ok = exist(fullfile(dam_dir, 'rasters', 'domain_clipped', 'dem_fabdem_30m_domain.tif'), 'file') == 2 && ...
         exist(fullfile(dam_dir, 'rasters', 'domain_clipped', 'lulc_mapbiomas_30m_domain.tif'), 'file') == 2 && ...
         exist(fullfile(dam_dir, 'rasters', 'domain_clipped', 'soil_texture_usda_30m_domain.tif'), 'file') == 2 && ...
         exist(fullfile(dam_dir, 'hydrograph', 'inlet_cells_hydrograph.csv'), 'file') == 2;
end

function tf = matches_any_code(folder_name, dam_codes)
    tf = false;
    for i = 1:numel(dam_codes)
        code = dam_codes{i};
        if ~isempty(strfind(folder_name, code)) %#ok<STREMP>
            tf = true;
            return;
        end
    end
end

function codes = normalize_dam_codes(raw)
    if isempty(raw)
        codes = {};
        return;
    end
    if isnumeric(raw)
        raw = arrayfun(@(x) sprintf('%06d', x), raw(:), 'UniformOutput', false);
    elseif isstring(raw)
        raw = cellstr(raw(:));
    end
    codes = cell(size(raw));
    for i = 1:numel(raw)
        s = char(raw{i});
        digits_only = regexp(s, '\d+', 'match', 'once');
        if isempty(digits_only)
            codes{i} = s;
        else
            codes{i} = sprintf('%06d', str2double(digits_only));
        end
    end
end
