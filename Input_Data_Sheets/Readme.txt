As supplied, these workbooks run a small example on a laptop: the San
Francisquito Creek case of Examples/Stanford on its 30 m grid (100 mm/h of rain
for 1 h, 3 h simulated, on the CPU, 15 to 20 min). The maps are read from
Examples/Stanford/Static; the dates (2000-01-01 00:00 to 03:00) are only labels.
Run it in MATLAB from the HydroPol2D folder with

    setenv('HYDROPOL_RUN_MODE', 'excel'); run('run_HydroPol2D.m')

and expect the same results as Examples/Stanford/run_excel.m. For your own case,
copy HydroPol2D_Inputs.xlsx and replace the maps, dates and settings.

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
