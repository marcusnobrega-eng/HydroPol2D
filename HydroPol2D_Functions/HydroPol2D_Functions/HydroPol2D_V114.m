
%% HydroPol2D Model
% Developer: Marcus Nobrega
% Main Script
% % % % % % % % % % % % % Model Status % % % % % % % % % % % % %
% ---------- Version 12.0 -------------
% Last Update - 6/22/2023
% Current Update - Added Huff and Alternated Blocks


%% Pre-Processing

clear all; clc;
model_folder = 'General_Data_HydroPol2D.xlsx';
input_table = readtable(model_folder);

% Load Model Functions
HydroPol2D_tools = char(table2cell(input_table(9,31)));
addpath(genpath(char(HydroPol2D_tools)));

HydroPol2D_preprocessing
% Pre_Processing_Data
%% Sensitivity Analysis
% Sensitivity_Analysis_HydroPol2D;

%% Main While Loop
outlet_index = Wshed_Properties.perimeter;
HydroPol2D_Main_While; % Runs the Main While Loop of the Model

%% Post-Processing Results
close all
post_processing