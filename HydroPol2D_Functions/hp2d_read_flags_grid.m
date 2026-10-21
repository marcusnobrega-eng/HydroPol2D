function F = hp2d_read_flags_grid(workbookPath, settingsGrid)
%HP2D_READ_FLAGS_GRID The grid to scan for flag_* name/value pairs.
%   Older workbooks keep every flag on a Flags sheet. The single case workbook has no
%   Flags sheet: each switch sits in the block it turns on, so the flags are read from
%   the settings tabs (settingsGrid, from hp2d_read_settings_grid, read here if not given).
if hp2d_workbook_has_sheet(workbookPath, 'Flags')
    F = readcell(workbookPath, 'Sheet', 'Flags');
elseif nargin >= 2 && ~isempty(settingsGrid)
    F = settingsGrid;
else
    F = hp2d_read_settings_grid(workbookPath);
end
end
