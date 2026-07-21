Keep the supplied workbook names and existing row layout. In
LULC_parameters.xlsx, Kc [-] scales internally calculated reference ET to the
potential soil-ET demand of each land-cover class. It is not applied when
evapotranspiration maps are prescribed, and it does not scale ponded open-water
evaporation, which uses Ep directly. Enter snow parameters in the snow columns
of LULC_parameters.xlsx; HydroPol2D no longer uses a separate
Snow_Parameters.xlsx workbook.

The supplied LULC rows are generic starting values, not a calibration. Their
literature basis, recommended application, and parameters that must be defined
for a site-specific water-quality or snow study are documented in
LULC_Defaults.md.
