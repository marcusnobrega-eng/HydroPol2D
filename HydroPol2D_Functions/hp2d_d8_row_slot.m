function slot = hp2d_d8_row_slot(keys, query, count, shared_keys)
%HP2D_D8_ROW_SLOT Bracketing breakpoint of each table row (1-based).
%   keys is rows x points, each row non-decreasing (rows shorter than the
%   table width repeat their last value). slot(r) is the last breakpoint of
%   row r that is <= query(r), kept within [1, count(r)]. This is
%   searchsorted(row, query, side='right') - 1, clipped to [0, count - 1],
%   in HydroPol2D-Python (_row_slot_jax), plus one.
%   shared_keys (optional): one strictly increasing row used by every table
%   row; the lookup then uses discretize instead of comparing every column.

query = query(:);
if nargin > 3 && ~isempty(shared_keys)
    if isgpuarray(query)
        slot = sum(shared_keys(:).' <= query, 2);
    else
        slot = discretize(query, [shared_keys(:).', inf]);
        slot(isnan(slot)) = 1;
    end
else
    slot = sum(keys <= query, 2);
end
slot = min(max(slot, 1), count(:));

end
