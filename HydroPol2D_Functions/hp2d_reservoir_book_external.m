function outlet_flow = hp2d_reservoir_book_external(outlet_flow, external_m3, cell_area, dt_s)
%HP2D_RESERVOIR_BOOK_EXTERNAL Add reservoir releases that leave the model to
%   the outlet flow [mm/h over the coarse cell] at their source cells, so the
%   outlet record and the mass ledger count them. external_m3 comes from
%   hp2d_reservoir_releases and is already removed from the depth.

idx = external_m3 > 0;
if ~any(idx(:))
    return
end

if isscalar(cell_area)
    area = cell_area;
else
    area = cell_area(idx);
end
outlet_flow(idx) = outlet_flow(idx) + 1000 * 3600 .* external_m3(idx) ./ (area .* dt_s);

end
