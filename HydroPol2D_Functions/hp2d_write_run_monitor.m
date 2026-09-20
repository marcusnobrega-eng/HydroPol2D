function hp2d_write_run_monitor(progress_path, metrics_path, values, append_metrics)
%HP2D_WRITE_RUN_MONITOR Write lightweight, atomic run telemetry.

if nargin < 4
    append_metrics = false;
end
if isempty(progress_path)
    return
end

folder = fileparts(progress_path);
if ~isempty(folder) && ~isfolder(folder)
    mkdir(folder);
end
values.updated_utc = char(datetime('now', 'TimeZone', 'UTC', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ss.SSSXXX'));

temporary = [progress_path '.tmp'];
fid = fopen(temporary, 'w');
if fid < 0
    warning('HydroPol2D:RunMonitorWriteFailed', ...
        'Could not open run-progress file: %s', temporary);
    return
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(values));
clear cleanup
movefile(temporary, progress_path, 'f');

if ~append_metrics || isempty(metrics_path)
    return
end

new_file = ~isfile(metrics_path);
fid = fopen(metrics_path, 'a');
if fid < 0
    warning('HydroPol2D:RunMetricsWriteFailed', ...
        'Could not open live-metrics file: %s', metrics_path);
    return
end
cleanup = onCleanup(@() fclose(fid));
if new_file
    fprintf(fid, ['updated_utc,simulated_minutes,total_minutes,percent,', ...
        'time_step_seconds,elapsed_seconds,eta_seconds,rainfall_mean_mm_h,', ...
        'rainfall_max_mm_h,max_surface_depth_m,outlet_discharge_m3_s,', ...
        'total_storage_m3,mass_balance_error_m3,mean_groundwater_depth_m\n']);
end
fprintf(fid, ['%s,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,', ...
    '%.9g,%.9g,%.9g,%.9g,%.9g\n'], values.updated_utc, ...
    values.simulated_minutes, values.total_minutes, values.percent, ...
    values.time_step_seconds, values.elapsed_seconds, values.eta_seconds, ...
    values.rainfall_mean_mm_h, values.rainfall_max_mm_h, ...
    values.max_surface_depth_m, values.outlet_discharge_m3_s, ...
    values.total_storage_m3, values.mass_balance_error_m3, ...
    values.mean_groundwater_depth_m);
end
