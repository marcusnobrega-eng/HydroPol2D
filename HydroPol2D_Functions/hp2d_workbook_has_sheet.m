function tf = hp2d_workbook_has_sheet(workbookPath, sheetName)
%HP2D_WORKBOOK_HAS_SHEET True when the workbook has a sheet of this name (any case).
%   False for an empty path, a missing file or an unreadable workbook.
tf = false;
if isempty(workbookPath) || ~isfile(workbookPath)
    return
end
try
    names = cellstr(string(sheetnames(workbookPath)));
catch
    return
end
tf = any(strcmpi(names, char(sheetName)));
end
