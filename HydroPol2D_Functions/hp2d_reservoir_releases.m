function [exchange_m3, external_m3, release_m3] = hp2d_reservoir_releases( ...
    V_m3, reservoir_x, reservoir_y, k1, h1, k2, k3, h2, k4, yds1, xds1, yds2, xds2, dt_s, geometry)
%HP2D_RESERVOIR_RELEASES Rating-curve releases of the reservoir controls.
%   V_m3 is each cell's water volume after this step's face fluxes [m3].
%   Each control drains its source cell through up to two outlets, in order:
%       Q = k * max(h - h0, 0)^b                [m3/s]
%       release = min(Q * dt, volume above h0)  [m3]
%   h is the source depth after the face fluxes and the earlier outlets; h0
%   is the sill above the cell invert. The release leaves the source and
%   enters the destination cell, so the transfer conserves mass. A
%   destination index of 0 in both directions (workbook coordinates 0, 0)
%   sends the release out of the model; NaN switches the outlet off. This is
%   reservoir_outlet_transfers in HydroPol2D-Python.
%
%   geometry says how depth and volume relate in a cell:
%       struct('type','area','area',A)   depth = V / A (A scalar or matrix)
%       struct('type','neal','River_Width',W,'River_Depth',D,'Resolution',dx)
%           depth above the channel bed of a Neal (2012) cell
%       struct('type','handles','depth',@(V,i) ...,'volume',@(h,i) ...)
%
%   exchange_m3: change of each cell's volume [m3], negative at the sources.
%   external_m3: volume that left the model, at its source cell [m3].
%   release_m3:  total release of each control [m3].

exchange_m3 = zeros(size(V_m3), 'like', V_m3);
external_m3 = zeros(size(V_m3), 'like', V_m3);
n_controls = numel(reservoir_y);
release_m3 = zeros(n_controls, 1);

for ii = 1:n_controls
    source = sub2ind(size(V_m3), reservoir_y(ii), reservoir_x(ii));
    for outlet = 1:2
        if outlet == 1
            k = k1(ii); sill = h1(ii); b = k2(ii); row = yds1(ii); col = xds1(ii);
        else
            k = k3(ii); sill = h2(ii); b = k4(ii); row = yds2(ii); col = xds2(ii);
        end
        if isnan(row) || isnan(col)
            continue
        end

        V_source = max(V_m3(source) + exchange_m3(source), 0);
        head = max(cell_depth(V_source, source, geometry) - sill, 0);
        available = max(V_source - cell_volume(sill, source, geometry), 0);
        dV = min(k * head^b * dt_s, available);
        if ~(dV > 0)
            continue
        end

        exchange_m3(source) = exchange_m3(source) - dV;
        if row == 0 && col == 0
            external_m3(source) = external_m3(source) + dV;
        else
            destination = sub2ind(size(V_m3), row, col);
            exchange_m3(destination) = exchange_m3(destination) + dV;
        end
        release_m3(ii) = release_m3(ii) + dV;
    end
end

end

function h = cell_depth(V, idx, geometry)
switch geometry.type
    case 'area'
        h = V / pick(geometry.area, idx);
    case 'neal'
        h = hp2d_neal_depth_from_volume(V, pick(geometry.River_Width, idx), ...
            pick(geometry.River_Depth, idx), geometry.Resolution);
    case 'handles'
        h = geometry.depth(V, idx);
    otherwise
        error('hp2d_reservoir_releases: unknown geometry type %s.', geometry.type);
end
end

function V = cell_volume(h, idx, geometry)
switch geometry.type
    case 'area'
        V = h * pick(geometry.area, idx);
    case 'neal'
        V = hp2d_neal_cell_volume(h, pick(geometry.River_Width, idx), ...
            pick(geometry.River_Depth, idx), geometry.Resolution);
    case 'handles'
        V = geometry.volume(h, idx);
    otherwise
        error('hp2d_reservoir_releases: unknown geometry type %s.', geometry.type);
end
end

function value = pick(field, idx)
if isscalar(field)
    value = field;
else
    value = field(idx);
end
end
