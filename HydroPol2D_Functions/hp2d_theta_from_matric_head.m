function theta = hp2d_theta_from_matric_head(head_m, theta_r, theta_s, alpha_m_inverse, n_vg)
%HP2D_THETA_FROM_MATRIC_HEAD Van Genuchten water content at a matric head (m, <= 0).
m = 1 - 1./n_vg;
suction = max(-head_m, 0);
saturation = (1 + (alpha_m_inverse.*suction).^n_vg).^(-m);
theta = theta_r + (theta_s - theta_r).*saturation;
end
