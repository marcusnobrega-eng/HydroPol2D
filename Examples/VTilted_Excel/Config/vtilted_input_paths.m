function InputPaths = vtilted_input_paths(model_root, Overrides)
%VTILTED_INPUT_PATHS Portable paths for the bundled V-Tilted example.

if nargin < 2 || isempty(Overrides)
    Overrides = struct();
end
config_root = fileparts(mfilename('fullpath'));
case_root = fileparts(config_root);
static_root = fullfile(case_root, 'Static');
forcing_root = fullfile(case_root, 'Forcing');

InputPaths = struct();
InputPaths.case_root = case_root;
InputPaths.topo_path = fullfile(model_root, 'third_party', 'topotoolbox_lite');
InputPaths.hydropol2d_tools = fullfile(model_root, 'HydroPol2D_Functions');
InputPaths.DEM_path = fullfile(static_root, 'DEM.tif');
InputPaths.LULC_path = fullfile(static_root, 'LULC.tif');
InputPaths.SOIL_path = fullfile(static_root, 'SOIL.tif');
InputPaths.DTB_path = fullfile(static_root, 'DTB.tif');
InputPaths.LAI_path = fullfile(static_root, 'LAI.tif');
InputPaths.Albedo_path = fullfile(static_root, 'Albedo.tif');
InputPaths.Initial_Soil_Moisture_path = fullfile(static_root, 'Initial_SM.tif');
InputPaths.Outlet_Cells_CSV = fullfile(forcing_root, 'Outlet', 'outlet_cells.csv');
InputPaths.Rainfall_Timeseries_File = fullfile( ...
    case_root, 'Input_Data_Sheets', 'Rainfall_Intensity_Data.xlsx');
InputPaths.Rainfall_Rasters_Folder = fullfile(forcing_root, 'Rainfall');
InputPaths.Transpiration_Rasters_Folder = fullfile(forcing_root, 'Transpiration');
InputPaths.Evaporation_Rasters_Folder = fullfile(forcing_root, 'Evaporation');
InputPaths.Rainfall_Raster_Files = {};
InputPaths.Transpiration_Raster_Files = {};
InputPaths.Evaporation_Raster_Files = {};
InputPaths.Inflow_Hydrograph_CSV = fullfile(forcing_root, 'Inflow', 'inflow_hydrograph.csv');
InputPaths.Stage_Hydrograph_CSV = fullfile(forcing_root, 'Stage', 'stage_hydrograph.csv');
InputPaths.Observed_Gauges_CSV = fullfile(forcing_root, 'Observed_Gauges', 'observed_gauges.csv');
InputPaths.ETP_input_spreadsheet = fullfile(forcing_root, 'Evapotranspiration', 'ETP_input_data.xlsx');

names = fieldnames(Overrides);
for index = 1:numel(names)
    InputPaths.(names{index}) = Overrides.(names{index});
end
end
