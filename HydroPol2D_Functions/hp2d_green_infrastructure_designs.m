function designs = hp2d_green_infrastructure_designs(designs_in)
%HP2D_GREEN_INFRASTRUCTURE_DESIGNS Green-infrastructure design library with defaults.
%   designs = hp2d_green_infrastructure_designs() returns the three starting
%   templates of HydroPol2D-Python (default_green_infrastructure_designs):
%   'bioretention_default', 'infiltration_trench_default' and
%   'permeable_pavement_default'.
%   designs = hp2d_green_infrastructure_designs(designs_in) fills the fields a
%   design leaves out with the defaults of GreenInfrastructureDesign
%   (config/models.py) and checks the same limits. designs_in is a struct
%   array, a cell array of structs or a table, one design per row, with the
%   Python field names and units:
%     identifier, name, practice_type ('bioretention', 'infiltration_trench'
%     or 'permeable_pavement'), surface_storage_depth_m (0.15),
%     surface_manning_n (0.10), media_depth_m (0.60), media_porosity (0.35),
%     media_residual_water_content (0.05), media_initial_water_content (0.15),
%     media_ksat_mm_h (100), media_alpha_vg_m_inverse (3.5), media_n_vg (1.60),
%     gravel_depth_m (0.30), gravel_void_fraction (0.40),
%     gravel_ksat_mm_h (1000), et_crop_coefficient (1.0),
%     surface_exchange_coefficient (0.60), inlet_width_m (1.0),
%     overflow_crest_depth_m (0.15), overflow_coefficient (1.70),
%     underdrain_coefficient_m3_s_m_exponent (0), underdrain_exponent (0.5),
%     underdrain_sill_depth_m (0).

defaults = struct( ...
    'identifier', '', 'name', '', 'practice_type', 'bioretention', ...
    'surface_storage_depth_m', 0.15, 'surface_manning_n', 0.10, ...
    'media_depth_m', 0.60, 'media_porosity', 0.35, ...
    'media_residual_water_content', 0.05, 'media_initial_water_content', 0.15, ...
    'media_ksat_mm_h', 100.0, 'media_alpha_vg_m_inverse', 3.5, 'media_n_vg', 1.60, ...
    'gravel_depth_m', 0.30, 'gravel_void_fraction', 0.40, 'gravel_ksat_mm_h', 1000.0, ...
    'et_crop_coefficient', 1.0, 'surface_exchange_coefficient', 0.60, ...
    'inlet_width_m', 1.0, 'overflow_crest_depth_m', 0.15, 'overflow_coefficient', 1.70, ...
    'underdrain_coefficient_m3_s_m_exponent', 0.0, 'underdrain_exponent', 0.5, ...
    'underdrain_sill_depth_m', 0.0);

if nargin < 1 || isempty(designs_in)
    designs_in = { ...
        struct('identifier', 'bioretention_default', 'name', 'Bioretention / rain garden', ...
            'practice_type', 'bioretention', 'surface_storage_depth_m', 0.20, ...
            'surface_manning_n', 0.12, 'media_depth_m', 0.60, 'media_ksat_mm_h', 100.0, ...
            'gravel_depth_m', 0.30, 'underdrain_coefficient_m3_s_m_exponent', 0.02), ...
        struct('identifier', 'infiltration_trench_default', 'name', 'Infiltration trench', ...
            'practice_type', 'infiltration_trench', 'surface_storage_depth_m', 0.10, ...
            'surface_manning_n', 0.08, 'media_depth_m', 0.30, 'media_ksat_mm_h', 200.0, ...
            'gravel_depth_m', 0.90, 'gravel_void_fraction', 0.40, ...
            'overflow_crest_depth_m', 0.10), ...
        struct('identifier', 'permeable_pavement_default', 'name', 'Permeable pavement', ...
            'practice_type', 'permeable_pavement', 'surface_storage_depth_m', 0.03, ...
            'surface_manning_n', 0.02, 'media_depth_m', 0.20, 'media_ksat_mm_h', 250.0, ...
            'gravel_depth_m', 0.35, 'et_crop_coefficient', 0.0, ...
            'overflow_crest_depth_m', 0.03)};
end

items = as_cell(designs_in);
names = fieldnames(defaults);
designs = repmat(defaults, numel(items), 1);
for i = 1:numel(items)
    item = items{i};
    unknown = setdiff(fieldnames(item), names);
    if ~isempty(unknown)
        error('Green-infrastructure design %d has unknown fields: %s.', i, strjoin(unknown, ', '));
    end
    d = defaults;
    given = fieldnames(item);
    for j = 1:numel(given)
        value = item.(given{j});
        if isstring(value), value = char(value); end
        if iscell(value) && isscalar(value), value = value{1}; end
        if ~(ischar(value) && isempty(value)) && ~(isnumeric(value) && isempty(value))
            d.(given{j}) = value;
        end
    end
    d.identifier = char(d.identifier);
    if isempty(strtrim(d.identifier))
        error('Green-infrastructure design %d needs an identifier.', i);
    end
    if isempty(d.name), d.name = d.identifier; end
    d.name = char(d.name);
    d.practice_type = char(d.practice_type);
    if ~any(strcmp(d.practice_type, {'bioretention', 'infiltration_trench', 'permeable_pavement'}))
        error('GI design %s: practice_type must be bioretention, infiltration_trench or permeable_pavement.', ...
            d.identifier);
    end
    for j = 4:numel(names)
        d.(names{j}) = double(d.(names{j}));
        if ~isscalar(d.(names{j})) || ~isfinite(d.(names{j}))
            error('GI design %s: %s must be a finite number.', d.identifier, names{j});
        end
    end
    check(d, d.surface_storage_depth_m >= 0, 'surface_storage_depth_m >= 0');
    check(d, d.surface_manning_n > 0, 'surface_manning_n > 0');
    check(d, d.media_depth_m > 0, 'media_depth_m > 0');
    check(d, d.media_porosity > 0 && d.media_porosity <= 1, '0 < media_porosity <= 1');
    check(d, d.media_residual_water_content >= 0 && d.media_residual_water_content < 1, ...
        '0 <= media_residual_water_content < 1');
    check(d, d.media_initial_water_content >= 0 && d.media_initial_water_content <= 1, ...
        '0 <= media_initial_water_content <= 1');
    check(d, d.media_ksat_mm_h >= 0, 'media_ksat_mm_h >= 0');
    check(d, d.media_alpha_vg_m_inverse > 0, 'media_alpha_vg_m_inverse > 0');
    check(d, d.media_n_vg > 1, 'media_n_vg > 1');
    check(d, d.gravel_depth_m > 0, 'gravel_depth_m > 0');
    check(d, d.gravel_void_fraction > 0 && d.gravel_void_fraction <= 1, '0 < gravel_void_fraction <= 1');
    check(d, d.gravel_ksat_mm_h >= 0, 'gravel_ksat_mm_h >= 0');
    check(d, d.et_crop_coefficient >= 0, 'et_crop_coefficient >= 0');
    check(d, d.surface_exchange_coefficient >= 0, 'surface_exchange_coefficient >= 0');
    check(d, d.inlet_width_m > 0, 'inlet_width_m > 0');
    check(d, d.overflow_crest_depth_m >= 0, 'overflow_crest_depth_m >= 0');
    check(d, d.overflow_coefficient >= 0, 'overflow_coefficient >= 0');
    check(d, d.underdrain_coefficient_m3_s_m_exponent >= 0, 'underdrain_coefficient_m3_s_m_exponent >= 0');
    check(d, d.underdrain_exponent > 0, 'underdrain_exponent > 0');
    check(d, d.underdrain_sill_depth_m >= 0, 'underdrain_sill_depth_m >= 0');
    check(d, d.media_residual_water_content < d.media_porosity, ...
        'media residual water content below porosity');
    check(d, d.media_residual_water_content <= d.media_initial_water_content && ...
        d.media_initial_water_content <= d.media_porosity, ...
        'media initial water content between residual content and porosity');
    check(d, d.overflow_crest_depth_m <= d.surface_storage_depth_m, ...
        'overflow crest not above the surface storage depth');
    check(d, d.underdrain_sill_depth_m <= d.gravel_depth_m, ...
        'underdrain sill not above the gravel depth');
    designs(i) = d;
end
ids = lower({designs.identifier});
if numel(unique(ids)) ~= numel(ids)
    error('Green-infrastructure design identifiers must be unique.');
end
end

function check(d, ok, rule)
if ~ok
    error('GI design %s breaks the rule: %s.', d.identifier, rule);
end
end

function items = as_cell(value)
if istable(value)
    items = cell(height(value), 1);
    for i = 1:height(value)
        items{i} = table2struct(value(i, :));
    end
elseif iscell(value)
    items = value(:);
elseif isstruct(value)
    items = num2cell(value(:));
else
    error('Green-infrastructure designs must be a struct array, a cell array of structs or a table.');
end
end
