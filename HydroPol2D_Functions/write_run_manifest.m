function manifest_paths = write_run_manifest(output_dir, manifest, run_name)
%WRITE_RUN_MANIFEST Save a debug manifest for HydroPol2D runs.
%
% Saves:
%   - MAT file
%   - JSON file
%   - TXT summary with nested fields expanded

    if nargin < 3 || strlength(string(run_name)) == 0
        run_name = "run";
    end

    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    ts = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    base_name = sprintf('%s_manifest_%s', char(run_name), char(ts));

    mat_path  = fullfile(output_dir, [base_name '.mat']);
    json_path = fullfile(output_dir, [base_name '.json']);
    txt_path  = fullfile(output_dir, [base_name '.txt']);

    manifest.generated_at = char(datetime('now', ...
        'Format', 'yyyy-MM-dd HH:mm:ss'));

    % Save MAT
    save(mat_path, 'manifest', '-v7.3');

    % Save JSON
    try
        json_text = jsonencode(manifest, 'PrettyPrint', true);
    catch
        json_text = jsonencode(manifest);
    end

    fid = fopen(json_path, 'w');
    if fid == -1
        error('Could not open JSON file for writing: %s', json_path);
    end
    fwrite(fid, json_text, 'char');
    fclose(fid);

    % Save TXT with nested expansion
    fid = fopen(txt_path, 'w');
    if fid == -1
        error('Could not open TXT file for writing: %s', txt_path);
    end

    fprintf(fid, 'HydroPol2D preprocessing manifest\n');
    fprintf(fid, '================================\n\n');

    write_struct_recursive(fid, manifest, '');

    fclose(fid);

    manifest_paths = struct();
    manifest_paths.mat  = mat_path;
    manifest_paths.json = json_path;
    manifest_paths.txt  = txt_path;

    fprintf('Manifest written:\n');
    fprintf('  %s\n', mat_path);
    fprintf('  %s\n', json_path);
    fprintf('  %s\n', txt_path);
end

function write_struct_recursive(fid, s, prefix)
    fields = fieldnames(s);

    for i = 1:numel(fields)
        name = fields{i};
        value = s.(name);

        if isempty(prefix)
            full_name = name;
        else
            full_name = [prefix '.' name];
        end

        if isstruct(value)
            fprintf(fid, '%s:\n', full_name);
            write_struct_recursive(fid, value, full_name);
        else
            fprintf(fid, '%s = %s\n', full_name, value_to_string(value));
        end
    end
end

function out = value_to_string(v)
    if isempty(v)
        out = '';
    elseif ischar(v)
        out = v;
    elseif isstring(v)
        out = char(join(v, ", "));
    elseif islogical(v)
        if isscalar(v)
            out = mat2str(v);
        else
            out = mat2str(double(v));
        end
    elseif isnumeric(v)
        if isscalar(v)
            out = num2str(v);
        else
            sz = size(v);
            out = ['[' strjoin(string(sz), 'x') ' ' class(v) ']'];
        end
    elseif iscell(v)
        out = ['{' num2str(numel(v)) ' cell}'];
    elseif isdatetime(v)
        if isscalar(v)
            out = char(v);
        else
            out = ['[' num2str(numel(v)) ' datetime]'];
        end
    else
        out = ['[' class(v) ']'];
    end
end