function specs = hydropol2d_build_parameter_specs(baseOverrides, Catalog)
%HYDROPOL2D_BUILD_PARAMETER_SPECS Build generic GA parameter specs.
%
% Inputs
% ------
% baseOverrides:
%   Struct with optional LULC and SOIL override tables:
%     baseOverrides.LULC.Index
%     baseOverrides.LULC.<field>
%     baseOverrides.SOIL.Index
%     baseOverrides.SOIL.<field>
%
% Catalog:
%   Struct from hydropol2d_default_parameter_catalog or a user-edited copy.
%
% Output
% ------
% specs:
%   Struct array consumed by run_hydropol2d_ga_calibration.

if nargin < 1 || isempty(baseOverrides)
    baseOverrides = struct();
end
if nargin < 2 || isempty(Catalog)
    Catalog = hydropol2d_default_parameter_catalog(baseOverrides);
end

specs = repmat(hydropol2d_empty_parameter_spec(), 0, 1);
specs = [specs; build_class_group_specs(baseOverrides, Catalog, 'LULC')];
specs = [specs; build_class_group_specs(baseOverrides, Catalog, 'SOIL')];
specs = [specs; build_inputdata_specs(Catalog)];
end

function specs = build_class_group_specs(baseOverrides, Catalog, group)
    specs = repmat(hydropol2d_empty_parameter_spec(), 0, 1);
    if ~isfield(Catalog, group)
        return;
    end
    G = Catalog.(group);
    idx = [];
    if isfield(G, 'Index') && ~isempty(G.Index)
        idx = G.Index(:);
    elseif isfield(baseOverrides, group) && isfield(baseOverrides.(group), 'Index')
        idx = baseOverrides.(group).Index(:);
    end
    if isempty(idx)
        return;
    end
    fields = G.Fields;
    for f = 1:numel(fields)
        rule = fields(f);
        if isfield(rule, 'make_global_multiplier') && rule.make_global_multiplier
            globalName = rule.global_name;
            if isempty(globalName)
                globalName = sprintf('%s_%s_multiplier_all', lower(group), rule.field);
            end
            specs(end + 1, 1) = hydropol2d_make_parameter_spec( ...
                globalName, group, rule.field, idx, 'multiplier', ...
                rule.multiplier_lower, rule.multiplier_upper, ...
                rule.multiplier_initial, rule.transform, false); %#ok<AGROW>
        end

        if ~isfield(rule, 'make_class_values') || rule.make_class_values
            baseVals = get_base_values(baseOverrides, group, rule.field, idx, ...
                rule.initial_default);
            for i = 1:numel(idx)
                specs(end + 1, 1) = hydropol2d_make_parameter_spec( ...
                    sprintf('%s_%g_%s', lower(group), idx(i), rule.field), ...
                    group, rule.field, idx(i), 'value', ...
                    rule.lower, rule.upper, baseVals(i), rule.transform, false); %#ok<AGROW>
            end
        end
    end
end

function specs = build_inputdata_specs(Catalog)
    specs = repmat(hydropol2d_empty_parameter_spec(), 0, 1);
    if ~isfield(Catalog, 'InputData') || isempty(Catalog.InputData)
        return;
    end
    for i = 1:numel(Catalog.InputData)
        rule = Catalog.InputData(i);
        specs(end + 1, 1) = hydropol2d_make_inputdata_parameter_spec( ...
            rule.name, rule.path, rule.lower, rule.upper, rule.initial, ...
            rule.transform, rule.enabled); %#ok<AGROW>
    end
end

function values = get_base_values(baseOverrides, group, field, idx, defaultValue)
    values = repmat(defaultValue, numel(idx), 1);
    if ~isfield(baseOverrides, group) || ~isfield(baseOverrides.(group), 'Index') || ...
            ~isfield(baseOverrides.(group), field)
        return;
    end
    baseIdx = baseOverrides.(group).Index(:);
    baseVals = baseOverrides.(group).(field);
    for i = 1:numel(idx)
        row = find(baseIdx == idx(i), 1);
        if isempty(row)
            continue;
        end
        if numel(baseVals) == 1
            values(i) = baseVals;
        else
            values(i) = baseVals(row);
        end
    end
end
