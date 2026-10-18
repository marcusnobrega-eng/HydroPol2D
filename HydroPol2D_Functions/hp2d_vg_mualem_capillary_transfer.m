function moved = hp2d_vg_mualem_capillary_transfer( ...
        ksat_m_s,alpha_1_m,n_vg,l_vg,gap_m,field_capacity_head_m, ...
        rooted_capacity_m,rooted_storage_m,available_groundwater_m,dt_s)
%HP2D_VG_MUALEM_CAPILLARY_TRANSFER Storage-limited capillary rise.
m=1-1./max(n_vg,1+eps);
se_gap=(1+(max(alpha_1_m,eps).*max(gap_m,0)).^max(n_vg,1+eps)).^(-m);
term=min(max(1-se_gap.^(1./m),0),1);
conductivity_factor=se_gap.^l_vg.*(1-term.^m).^2;
potential=max(ksat_m_s,0).*conductivity_factor.*dt_s;
se_fc=(1+(max(alpha_1_m,eps).*field_capacity_head_m).^max(n_vg,1+eps)).^(-m);
deficit=max(se_fc.*max(rooted_capacity_m,0)-max(rooted_storage_m,0),0);
moved=min(min(deficit,potential),max(available_groundwater_m,0));
moved(~isfinite(moved))=0;
end
