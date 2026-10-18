function column = hp2d_slot_column(z_top, z_bot, P)
%HP2D_SLOT_COLUMN Assemble the layered column struct (slot fields [N x L], cell fields [N x 1]).
%   P fields ([N x 1] or scalar): soil_theta_s, soil_theta_r, soil_alpha_m_inverse,
%   soil_n_vg, soil_kv_m_s ([N x 3] or [1 x 3]), soil_lateral_anisotropy,
%   regolith_theta_s, regolith_theta_r, regolith_alpha_m_inverse, regolith_n_vg,
%   regolith_kv_m_s, regolith_kh_m_s, retention_head_m, ground_m, column_base_depth_m,
%   valid, and optional regolith_specific_yield (NaN = derive).
%   A falling table leaves theta_retained, so Sy = theta_s - theta_retained. It is theta
%   at field capacity (retention_head_m) unless a finite rock Sy is given for deep slots.
%   Same rule as HydroPol2D-Python build_slot_column.
[N, L] = size(z_top);
deep = (1:L) > 3;
    function m = stack(soil_value, deep_value)
        soil_col = soil_value(:) + zeros(N, 1, 'like', z_top);
        deep_col = deep_value(:) + zeros(N, 1, 'like', z_top);
        m = repmat(soil_col, 1, L);
        m(:, deep) = repmat(deep_col, 1, L - 3);
    end
column.z_top_m = z_top;
column.z_bot_m = z_bot;
column.theta_s = stack(P.soil_theta_s, P.regolith_theta_s);
column.theta_r = stack(P.soil_theta_r, P.regolith_theta_r);
column.alpha_m_inverse = stack(P.soil_alpha_m_inverse, P.regolith_alpha_m_inverse);
column.n_vg = stack(P.soil_n_vg, P.regolith_n_vg);
soil_kv = P.soil_kv_m_s;
if size(soil_kv, 1) == 1, soil_kv = repmat(soil_kv, N, 1); end
kv = zeros(N, L, 'like', z_top);
kv(:, 1:3) = soil_kv;
kv(:, deep) = repmat(P.regolith_kv_m_s(:) + zeros(N, 1, 'like', z_top), 1, L - 3);
kh = zeros(N, L, 'like', z_top);
kh(:, 1:3) = soil_kv*P.soil_lateral_anisotropy;
kh(:, deep) = repmat(P.regolith_kh_m_s(:) + zeros(N, 1, 'like', z_top), 1, L - 3);
column.kv_m_s = kv;
column.kh_m_s = kh;
retained = hp2d_theta_from_matric_head(-abs(P.retention_head_m), column.theta_r, ...
    column.theta_s, column.alpha_m_inverse, column.n_vg);
if isfield(P, 'regolith_specific_yield') && ~isempty(P.regolith_specific_yield)
    sy = P.regolith_specific_yield(:) + zeros(N, 1, 'like', z_top);
    override = repmat(isfinite(sy), 1, L) & repmat(deep, N, 1);
    replacement = column.theta_s - sy;
    retained(override) = replacement(override);
end
column.theta_retained = retained;
column.ground_m = P.ground_m(:) + zeros(N, 1, 'like', z_top);
column.base_depth_m = P.column_base_depth_m(:) + zeros(N, 1, 'like', z_top);
column.valid = logical(P.valid(:)) & true(N, 1);
end
