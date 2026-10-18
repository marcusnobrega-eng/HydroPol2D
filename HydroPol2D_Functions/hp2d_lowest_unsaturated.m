function lowest = hp2d_lowest_unsaturated(unsat)
%HP2D_LOWEST_UNSATURATED Mask [N x L] of the deepest slot that has an unsaturated part.
has = unsat > 0;
from_bottom = fliplr(cumsum(fliplr(double(has)), 2));
lowest = has & (from_bottom - double(has) == 0);
end
