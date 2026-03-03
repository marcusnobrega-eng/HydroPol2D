%% HydroPol2D Model
% Developer: Marcus Nobrega, Ph.D.
% Main Script
% ---------- Version 15.0 -------------
% Last Update - 7/18/2024

%% Pre-Processing
clear; clc;

% -------------------------------------------------------------------------
% Main input spreadsheet (General_Data)
% -------------------------------------------------------------------------
model_folder = '/Users/mngomes/Documents/GitHub/HydroPol2D/Input_Data_Sheets/General_Data_HydroPol2D.xlsx';

% Optional: make sure folder containing spreadsheets is visible
addpath('Input_Data_Sheets');

% Read General_Data as a raw cell grid (label/value format)
GD = readcell(model_folder, 'Sheet', 'General_Data');

% -------------------------------------------------------------------------
% Add required toolboxes / model functions (name-based)
% -------------------------------------------------------------------------
hydropol2d_tools = string(xlget(GD,'hydropol2d_tools'));
topo_path        = string(xlget(GD,'topo_path'));   % if you want TopoToolbox available globally

addpath(genpath(char(hydropol2d_tools)));
addpath(genpath(char(topo_path)));

% -------------------------------------------------------------------------
% Run preprocessing 
% -------------------------------------------------------------------------
HydroPol2D_preprocessing;   %

%% Main While Loop (solver)
HydroPol2D_Main_While;

%% Post-Processing Results
close all
post_processing;


% ========================================================================
% Local helper (same as you already have in preprocessing)
% ========================================================================
function v = xlget(GD, key)
%XLGET Find a label (key) anywhere in a readcell() grid and return the cell to the right.
% Case-insensitive. Errors if key not found.

    S = strings(size(GD));
    for r = 1:size(GD,1)
        for c = 1:size(GD,2)
            if ischar(GD{r,c}) || isstring(GD{r,c})
                S(r,c) = string(GD{r,c});
            end
        end
    end

    [rr, cc] = find(strcmpi(strtrim(S), key), 1, 'first');
    if isempty(rr)
        error("General_Data: key '%s' not found.", key);
    end
    if cc == size(GD,2)
        error("General_Data: key '%s' found in last column; no value to the right.", key);
    end

    v = GD{rr, cc+1};
end