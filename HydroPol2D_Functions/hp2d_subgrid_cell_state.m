function state = hp2d_subgrid_cell_state(SubgridTables, query, query_type, Resolution)
%HP2D_SUBGRID_CELL_STATE Query lookup-subgrid cell storage state.
%
% query_type:
%   'depth'  representative depth above cell invert [m]
%   'eta'    absolute water-surface elevation [m]
%   'volume' stored volume [m3]

if nargin < 4 || isempty(Resolution)
    Resolution = [];
end

switch lower(string(query_type))
    case "depth"
        drep = max(query, 0);
    case "eta"
        drep = max(query - SubgridTables.invert_el, 0);
    case "volume"
        if isempty(Resolution)
            error('hp2d_subgrid_cell_state:missingResolution', ...
                'Resolution is required when querying cell state by volume.');
        end
        drep = hp2d_subgrid_inverse_volume( ...
            SubgridTables.volume_cell, query, ...
            SubgridTables.dz, SubgridTables.maxDepth, Resolution);
    otherwise
        error('hp2d_subgrid_cell_state:badQueryType', ...
            'query_type must be depth, eta, or volume.');
end

state = struct();
state.drep = drep;
state.eta = SubgridTables.invert_el + drep;
state.volume = hp2d_subgrid_lookup_depth( ...
    SubgridTables.volume_cell, drep, SubgridTables.dz, SubgridTables.maxDepth);
state.area = hp2d_subgrid_lookup_depth( ...
    SubgridTables.area_cell, drep, SubgridTables.dz, SubgridTables.maxDepth);
state.area(~isfinite(state.area)) = 0;
state.volume(~isfinite(state.volume)) = 0;
end
