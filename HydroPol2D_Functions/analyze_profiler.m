function analyze_profiler(topN_funcs, topN_lines)
% ANALYZE_PROFILER Analyze MATLAB profiler results
% topN_funcs: number of top functions to display
% topN_lines: number of top lines per function

    % Get profiler info
    p = profile('info');

    if isempty(p.FunctionTable)
        error('No profiling data found. Run profile on first.');
    end

    % Extract total time per function
    funcTimes = [p.FunctionTable.TotalTime];

    % Sort functions by total time (descending)
    [~, idx] = sort(funcTimes, 'descend');

    fprintf('\n===== Top %d Functions by Total Time =====\n', topN_funcs);

    % Loop over top N functions
    for i = 1:min(topN_funcs, length(idx))
        f = p.FunctionTable(idx(i));

        fprintf('\n--- Function #%d ---\n', i);
        fprintf('Name       : %s\n', f.FunctionName);
        fprintf('File       : %s\n', f.FileName);
        fprintf('Total Time : %.6f sec\n', f.TotalTime);
        fprintf('Num Calls  : %d\n', f.NumCalls);

        % If no line info available, skip
        if isempty(f.ExecutedLines)
            fprintf('No line-level profiling data available.\n');
            continue;
        end

        % ExecutedLines format: [lineNumber, hits, time]
        lines = f.ExecutedLines;

        % Sort lines by time spent (3rd column)
        [~, lineIdx] = sort(lines(:,3), 'descend');

        fprintf('Top %d lines in this function:\n', topN_lines);

        % Read file content
        if exist(f.FileName, 'file')
            fileText = fileread(f.FileName);
            fileLines = splitlines(fileText);
        else
            fileLines = {};
        end

        % Loop through top lines
        for j = 1:min(topN_lines, size(lines,1))
            ln = lines(lineIdx(j), 1);
            hits = lines(lineIdx(j), 2);
            t = lines(lineIdx(j), 3);

            % Get code line if available
            if ~isempty(fileLines) && ln <= length(fileLines)
                codeLine = strtrim(fileLines{ln});
            else
                codeLine = '[Code not available]';
            end

            fprintf('  Line %d | Time: %.6f sec | Hits: %d\n', ln, t, hits);
            fprintf('    -> %s\n', codeLine);
        end
    end

    fprintf('\n===== Analysis Complete =====\n');
end