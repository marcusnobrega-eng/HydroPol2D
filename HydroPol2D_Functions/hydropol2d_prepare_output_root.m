function exportRootDir = hydropol2d_prepare_output_root(exportRootDir, cleanOutputFolder)
%HYDROPOL2D_PREPARE_OUTPUT_ROOT Prepare an output directory without aborting
% a run when Windows keeps an old result file locked.

arguments
    exportRootDir (1,:) char
    cleanOutputFolder (1,1) logical = true
end

if ~isfolder(exportRootDir)
    create_directory(exportRootDir);
    return;
end

if is_dir_empty(exportRootDir)
    return;
end

if ~cleanOutputFolder
    warning('HydroPol2D:ExportRootNotEmpty', ...
        ['Export root already contains files and clean_output_folder is false. ' ...
         'Existing contents will be preserved:\n  %s'], exportRootDir);
    return;
end

warning('HydroPol2D:ExportRootNotEmpty', ...
    ['Export root already contains files. HydroPol2D will remove the old ' ...
     'contents before exporting new results:\n  %s'], exportRootDir);

[cleaned, failureSummary] = delete_dir_contents(exportRootDir);
if cleaned
    return;
end

requestedRoot = exportRootDir;
exportRootDir = create_fallback_root(requestedRoot);
warning('HydroPol2D:ExportCleanupFailed', ...
    ['Some previous outputs could not be removed, usually because a file is ' ...
     'open in MATLAB, Excel, File Explorer, or another program. The current ' ...
     'run will continue in a new directory.\nRequested root:\n  %s\n' ...
     'Current run root:\n  %s\nCleanup details:\n  %s'], ...
    requestedRoot, exportRootDir, failureSummary);
end

function [cleaned, failureSummary] = delete_dir_contents(folder)
entries = dir(folder);
failures = strings(0,1);

for index = 1:numel(entries)
    name = entries(index).name;
    if strcmp(name, '.') || strcmp(name, '..')
        continue;
    end

    target = fullfile(folder, name);
    try
        make_writable(target);
        if entries(index).isdir
            [removed, message, messageId] = rmdir(target, 's');
            if ~removed
                error(normalize_message_id(messageId), '%s', message);
            end
        else
            delete(target);
            if isfile(target)
                error('HydroPol2D:OutputDeleteFailed', ...
                    'The file still exists after deletion was requested.');
            end
        end
    catch ME
        failures(end + 1, 1) = sprintf('%s: %s', target, ME.message); %#ok<AGROW>
    end
end

cleaned = is_dir_empty(folder);
if cleaned
    failureSummary = '';
elseif isempty(failures)
    failureSummary = 'One or more output entries remain locked.';
else
    failureSummary = strjoin(failures, newline);
end
end

function make_writable(target)
try
    if ispc
        fileattrib(target, '+w -h -s', '', 's');
    else
        fileattrib(target, '+w', 'u', 's');
    end
catch
    % Deletion below reports the actionable path if attributes cannot change.
end
end

function exportRootDir = create_fallback_root(requestedRoot)
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
baseName = ['Run_' timestamp];
lastMessage = 'All fallback directory names already exist.';

for suffix = 0:999
    if suffix == 0
        candidate = fullfile(requestedRoot, baseName);
    else
        candidate = fullfile(requestedRoot, sprintf('%s_%03d', baseName, suffix));
    end
    if ~isfolder(candidate) && ~isfile(candidate)
        [created, message] = mkdir(candidate);
        if created
            exportRootDir = candidate;
            return;
        end
        lastMessage = message;
    end
end

error('HydroPol2D:OutputFallbackFailed', ...
    'Could not create a fallback output folder under %s. Last error: %s', ...
    requestedRoot, lastMessage);
end

function create_directory(folder)
[created, message, messageId] = mkdir(folder);
if ~created
    error(normalize_message_id(messageId), ...
        'Could not create output directory %s: %s', folder, message);
end
end

function tf = is_dir_empty(folder)
entries = dir(folder);
tf = all(ismember({entries.name}, {'.', '..'}));
end

function messageId = normalize_message_id(messageId)
if isempty(messageId)
    messageId = 'HydroPol2D:OutputCleanupFailed';
end
end
