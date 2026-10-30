function F = hp2d_detention_pond_forcing(storage_m3, rainfall_mm_h, evaporation_mm_day, P, dt_s)
%HP2D_DETENTION_POND_FORCING Rain, open-water evaporation and bottom seepage of the ponds.
%   storage_m3: pond storage [m3] (n_ponds x 1); rainfall_mm_h and
%   evaporation_mm_day: rates at the pond cells P.cells (same order); P:
%   hp2d_prepare_detention_ponds; dt_s: step [s].
%   In order, as pond_forcing_step in HydroPol2D-Python (routing/detention.py):
%     1. rain on the whole footprint enters the pond;
%     2. open water evaporates: a DTM pond from its wet cells (cell elevation below
%        the stage), a power-law pond at the footprint-mean rate over its curve area
%        dV/dh = C b (h - invert)^(b-1); capped at the storage;
%     3. an unlined pond seeps at its bottom Ksat over the wetted area, capped at
%        the storage.
%   F.storage_m3 (after the three), F.rainfall_m3, F.evaporation_m3, F.seepage_m3
%   (per pond) and F.recharge_rate_m_s at P.cells: the seepage as a rate over the
%   wet cells of a DTM pond or uniformly over a power-law footprint. The caller
%   sends it to the groundwater under the pond, or out of the model when the
%   groundwater is off (as Python does).

if ~(dt_s > 0)
    error('HydroPol2D:DetentionPonds', 'Pond forcing dt_s must be positive.');
end
n = P.n_ponds;
A = P.cell_area_m2;
idx = P.cell_pond;
storage = max(double(storage_m3(:)), 0);
rain = max(double(rainfall_mm_h(:)), 0);
evaporation_rate = max(double(evaporation_mm_day(:)), 0);

rainfall = accumarray(idx, rain, [n 1]) / 1000 / 3600 * A * dt_s;
current = storage + rainfall;

[wet, wet_count, area] = wetted(current, P);
dtm_evaporation = accumarray(idx, evaporation_rate .* wet, [n 1]) * A;
mean_evaporation = accumarray(idx, evaporation_rate, [n 1]) ./ max(P.footprint_count, 1);
evaporation = mean_evaporation .* area;
evaporation(P.dtm) = dtm_evaporation(P.dtm);
evaporation = evaporation / 1000 / 86400 * dt_s;
evaporation = min(current, max(evaporation, 0));
current = current - evaporation;

[wet, wet_count, area] = wetted(current, P);
seepage = min(current, P.bottom_ksat_m_s .* area * dt_s);
current = current - seepage;

% Seepage as a recharge rate under the pond [m/s].
nominal = wet_count .* P.bottom_ksat_m_s * A * dt_s;
dtm_scale = zeros(n, 1);
has = nominal > 0;
dtm_scale(has) = seepage(has) ./ nominal(has);
footprint_area = P.footprint_count * A;
uniform_rate = zeros(n, 1);
has = footprint_area > 0;
uniform_rate(has) = seepage(has) ./ footprint_area(has) / dt_s;
recharge = uniform_rate(idx);
on_dtm = P.dtm(idx);
recharge(on_dtm) = wet(on_dtm) .* P.bottom_ksat_m_s(idx(on_dtm)) .* dtm_scale(idx(on_dtm));

F = struct('storage_m3', current, 'rainfall_m3', rainfall, 'evaporation_m3', evaporation, ...
    'seepage_m3', seepage, 'recharge_rate_m_s', recharge);
end

function [wet, wet_count, area] = wetted(storage, P)
% Wet cells (strictly below the stage) and the water-surface area of each pond.
stage = hp2d_detention_pond_stage(storage, P);
wet = P.cell_elevation_m < stage(P.cell_pond);
wet_count = accumarray(P.cell_pond, double(wet), [P.n_ponds 1]);
area = wet_count * P.cell_area_m2;
curve = ~P.dtm;
depth = max(stage(curve) - P.invert_elevation_m(curve), 0);
area(curve) = P.storage_coefficient(curve) .* P.storage_exponent(curve) .* ...
    depth .^ (P.storage_exponent(curve) - 1);
end
