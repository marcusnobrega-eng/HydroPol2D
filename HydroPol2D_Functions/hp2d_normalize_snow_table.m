function SnowConfig = hp2d_normalize_snow_table(source)
%HP2D_NORMALIZE_SNOW_TABLE Read the public per-LULC snow parameter table.

if istable(source)
    T = source;
elseif isstruct(source) && isfield(source, 'table') && istable(source.table)
    T = source.table;
else
    error('Snow parameters must be supplied as a table or as a struct with a table field.');
end

fields = { ...
    'LULC_Index', 'alpha'; ...
    'Snow_Albedo', 'alpha'; ...
    'Snow_Emissivity', 'epsilon'; ...
    'Sublimation_Coefficient_d_1', 'C_e'; ...
    'Degree_Day_Factor_mm_C_day', 'DDF'; ...
    'T_Snow_All_C', 'T_snow_all'; ...
    'T_Rain_All_C', 'T_rain_all'; ...
    'Rho_Snow_Init_kg_m3', 'rho_snow_init'; ...
    'Rho_Snow_Max_kg_m3', 'rho_max'; ...
    'Compaction_Temperature_kg_m3_C_day', 'k_t'; ...
    'Compaction_SWE_kg_m3_mm_day', 'k_swe'; ...
    'Compaction_Depth_kg_m3_mm_day', 'k_D'};

public_names = string(fields(2:end, 1))';
internal_names = string(fields(2:end, 2))';
% The public LULC table uses `Index` for its class code.  Retain
% `LULC_Index` as an accepted alias so older snow-only tables still work.
index_column = find_column(T.Properties.VariableNames, fields{1,1});
if isempty(index_column)
    index_column = find_column(T.Properties.VariableNames, 'Index');
end
if isempty(index_column)
    error('Snow table is missing required column "%s" (or the LULC table column "Index").', fields{1,1});
end

class_index = numeric_column(T, index_column, fields{1,1});
values = nan(height(T), numel(public_names));
for i = 1:numel(public_names)
    column = find_column(T.Properties.VariableNames, public_names(i));
    if isempty(column)
        error('Snow table is missing required column "%s".', public_names(i));
    end
    values(:,i) = numeric_column(T, column, public_names(i));
end

valid = isfinite(class_index) & all(isfinite(values), 2);
class_index = class_index(valid);
values = values(valid,:);
if isempty(class_index) || any(abs(class_index - round(class_index)) > 1e-9) || ...
        numel(unique(class_index)) ~= numel(class_index)
    error('Snow table must contain unique finite integer LULC_Index values.');
end

if any(values(:,1) < 0 | values(:,1) > 1) || any(values(:,2) < 0 | values(:,2) > 1)
    error('Snow albedo and emissivity must lie between 0 and 1.');
end
if any(values(:,3) < 0 | values(:,4) < 0 | values(:,8) < values(:,7) | ...
        values(:,7) <= 0 | values(:,9) < 0 | values(:,10) < 0 | values(:,11) < 0)
    error('Snow rate and density parameters must be nonnegative and rho_max must exceed rho_snow_init.');
end
if any(values(:,6) < values(:,5))
    error('T_Rain_All_C must be greater than or equal to T_Snow_All_C.');
end

SnowConfig = struct();
SnowConfig.class_index = class_index;
SnowConfig.parameter_names = internal_names;
SnowConfig.parameter_values = values;
end

function idx = find_column(variable_names, required_name)
target = normalize_name(required_name);
names = cellfun(@normalize_name, cellstr(variable_names), 'UniformOutput', false);
idx = find(strcmp(names, target), 1);
end

function value = numeric_column(T, column, label)
value = T{:, column};
if iscell(value)
    value = str2double(string(value));
elseif isstring(value) || ischar(value)
    value = str2double(string(value));
else
    value = double(value);
end
if ~isnumeric(value)
    error('Snow table column "%s" must be numeric.', label);
end
value = double(value(:));
end

function normalized = normalize_name(value)
value = regexprep(char(string(value)), '\s*\[.*\]\s*$', '');
normalized = lower(regexprep(value, '[^a-zA-Z0-9]', ''));
end
