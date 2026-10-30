function stage_m = hp2d_detention_pond_stage(storage_m3, P)
%HP2D_DETENTION_POND_STAGE Level-pool water level of every detention pond [m].
%   storage_m3: storage of each pond [m3]; P: hp2d_prepare_detention_ponds.
%   A DTM pond fills its footprint from the lowest cell up: with the n lowest
%   cells under water, stage = (V/A + sum of their elevations)/n, the candidate
%   that lies between the n-th and (n+1)-th elevations; an empty pond sits at its
%   lowest cell. A power-law pond inverts V = C (stage - invert)^b. This is
%   _stage_from_storage in HydroPol2D-Python (routing/detention.py).

storage = max(double(storage_m3(:)), 0);
stage_m = zeros(P.n_ponds, 1);
for p = 1:P.n_ponds
    V = storage(p);
    if P.dtm(p)
        z = P.sorted_elevation_m{p};
        if ~(V > 0)
            stage_m(p) = z(1);
            continue
        end
        candidates = (V / P.cell_area_m2 + P.cumulative_elevation_m{p}) ./ (1:numel(z))';
        valid = candidates >= z & candidates <= [z(2:end); inf];
        if any(valid)
            stage_m(p) = min(candidates(valid));
        else
            % Round-off between two equal elevations: the last bracket that holds.
            stage_m(p) = candidates(find(candidates >= z, 1, 'last'));
        end
    else
        stage_m(p) = P.invert_elevation_m(p) + ...
            (V / P.storage_coefficient(p))^(1 / P.storage_exponent(p));
    end
end
end
