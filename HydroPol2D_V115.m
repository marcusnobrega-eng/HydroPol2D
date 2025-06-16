%% HydroPol2D Model
% Developer: Marcus Nobrega, Ph.D.
% Main Script
% % % % % % % % % % % % % Model Status % % % % % % % % % % % % %
% ---------- Version 15.0 -------------
% Last Update - 7/18/2024

% % % % % % % % % % % % % All Rights Reserved % % % % % % % % % % % % %
%% Pre-Processing

% Clearing All Previous Data
clear all; clc;

% Adding Paths
addpath 'Input_Data_Sheets' % Folder where the excel files are available
model_folder = 'Input_Data_Sheets\General_Data_HydroPol2D.xlsx'; % Name of the main general data
input_table = readtable(model_folder);

% Load Model Functions
HydroPol2D_tools = char(table2cell(input_table(19,27))); % File path of the model functions
addpath(genpath(char(HydroPol2D_tools)));

HydroPol2D_preprocessing % Preprocessing code

%% Sensitivity Analysis
% Activate this code below to run a sensitivity analysis
% clear all
% Adding Paths
% addpath 'Input_Data_Sheets'
% addpath HydroPol2D_Functions\
% model_folder = 'Input_Data_Sheets\General_Data_HydroPol2D.xlsx';
% input_table = readtable(model_folder);
% ---- Variation Range ---- %
% var_range = [0.25 0.5 0.75 1 1.25 1.5 1.75]'; % variable/baseline, such that 0.5
% means the parameters will be changed in 50%.
% Sensitivity_Analysis_HydroPol2D;

%% Main While Loop
% This code runs the HydroPol2D solver
% Main While to Run the Model
HydroPol2D_Main_While; % Runs the Main While Loop of the Model

%% Post-Processing Results
close all
post_processing
