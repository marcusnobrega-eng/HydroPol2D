function alpha = hp2d_feddes_alpha(head_m, feddes, rate_mm_day)
%HP2D_FEDDES_ALPHA Feddes et al. (1978) reduction factor in [0, 1] for matric head (m).
%   0 wetter than h1, linear to 1 at h2, 1 down to h3, linear to 0 at the wilting head h4.
%   feddes: struct with h1_m, h2_m, h3_high_m, h3_low_m, h4_m.
if nargin < 3, rate_mm_day = 0; end
h1 = feddes.h1_m; h2 = feddes.h2_m; h4 = feddes.h4_m;
h3 = hp2d_feddes_h3(feddes, rate_mm_day);
wet_ramp = (h1 - head_m)/(h1 - h2);
dry_ramp = (head_m - h4)./max(h3 - h4, 1e-12);
alpha = ones(size(head_m + h3), 'like', head_m + h3);
alpha(head_m + 0*h3 > h1) = 0;
wet = head_m + 0*h3 <= h1 & head_m + 0*h3 > h2;
alpha(wet) = wet_ramp(wet);
dry = head_m + 0*h3 < h3;
dry_values = dry_ramp + 0*alpha;
alpha(dry) = dry_values(dry);
alpha = min(max(alpha, 0), 1);
end
