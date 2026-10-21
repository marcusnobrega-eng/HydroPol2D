function T = hp2d_read_case_table(workbookPath, tabName, legacyFile, varargin)
%HP2D_READ_CASE_TABLE A class table or time series: its tab in the case workbook, if any.
%   The single case workbook holds each table on its own tab, cell for cell where the old
%   separate workbook had it, so readtable returns the same table either way. Without
%   that tab the old workbook is read exactly as before (legacyFile, found by name on the
%   MATLAB path). varargin: readtable options (e.g. 'Range', 'A2').
%   Tabs: Soil classes, Land use classes, Human instability, Rain series, Rain gauges,
%   Climate stations, Inflow hydrographs, Stage hydrographs, Control structures.
if hp2d_workbook_has_sheet(workbookPath, tabName)
    T = readtable(workbookPath, 'Sheet', tabName, varargin{:});
else
    T = readtable(legacyFile, varargin{:});
end
end
