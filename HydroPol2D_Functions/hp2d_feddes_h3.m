function h3 = hp2d_feddes_h3(feddes, rate_mm_day)
%HP2D_FEDDES_H3 Demand-dependent h3 (SWAP): h3_high at >= 5 mm/d, h3_low at <= 1 mm/d.
weight = min(max((rate_mm_day - 1)/(5 - 1), 0), 1);
h3 = feddes.h3_low_m + weight.*(feddes.h3_high_m - feddes.h3_low_m);
end
