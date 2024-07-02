%% HydroPol2D Model
% Developer: Marcus Nobrega, Ph.D.
% Main Script
% % % % % % % % % % % % % Model Status % % % % % % % % % % % % %
% ---------- Version 14.0 -------------
% Last Update - 5/10/2024

% % % % % % % % % % % % % All Rights Reserved % % % % % % % % % % % % %
%% Pre-Processing

% Clearing All Previous Data
clear all; clc;

% Adding Paths
addpath 'Input_Data_Sheets'
model_folder = 'Input_Data_Sheets\General_Data_HydroPol2D.xlsx';
input_table = readtable(model_folder);

% Load Model Functions
HydroPol2D_tools = char(table2cell(input_table(9,31)));
addpath(genpath(char(HydroPol2D_tools)));

HydroPol2D_preprocessing
% Pre_Processing_Data
%% Sensitivity Analysis
% Activate this code below to run a sensitivity analysis
% ---- Variation Range ---- %
% var_range = [0.5 0.75 1 1.25 1.5]'; % variable/baseline, such that 0.5
% means the parameters will be changed in 50%.
% Sensitivity_Analysis_HydroPol2D;

%% Main While Loop
% This code runs the HydroPol2D solver
% If you want to assign outlet normal flow boundary conditions to all
% domain perimeter, activate the following code
% outlet_index = Wshed_Properties.perimeter;

outlet_index = Wshed_Properties.perimeter;

% Main While to Run the Model
HydroPol2D_Main_While; % Runs the Main While Loop of the Model

%% Post-Processing Results
close all
post_processing
