function head_m = hp2d_matric_head_from_theta(theta, theta_r, theta_s, alpha_m_inverse, n_vg)
%HP2D_MATRIC_HEAD_FROM_THETA Inverse van Genuchten matric head (m, <= 0); saturation in [1e-9, 1].
m = 1 - 1./n_vg;
saturation = (theta - theta_r)./max(theta_s - theta_r, 1e-12);
saturation = min(max(saturation, 1e-9), 1);
head_m = min(-((saturation.^(-1./m) - 1).^(1./n_vg))./alpha_m_inverse, 0);
end
