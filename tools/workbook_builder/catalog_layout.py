"""Layout of the HydroPol2D input workbook: tabs, blocks and the inputs in each.

Each input row is (key, label, kind, unit, default, used_when). kind is one of
number, text, path, date, switch (0/1) or a tuple of allowed integer options for
multi-choice switches. default None means required; "" means optional with no value.
The first input of a block may be its switch; the rest of the block greys out when it
is 0. A block may name another block's switch in "controlled_by". Notes (hover text)
come from the existing workbook comments, merged by make_catalog.py.
"""

REQ = None  # required: no default

TABS = [
    {
        "name": "Run control",
        "purpose": "Simulation period, time step, performance and what the model records.",
        "blocks": [
            {"title": "Simulation period", "inputs": [
                ("date_begin", "Start date and time", "date", "date", REQ, "always"),
                ("date_end", "End date and time", "date", "date", REQ, "always"),
            ]},
            {"title": "Surface time step", "inputs": [
                ("time_step_model", "Initial time step", "number", "s", REQ, "always"),
                ("min_time_step", "Minimum time step", "number", "s", REQ, "always"),
                ("max_time_step", "Maximum time step", "number", "s", REQ, "always"),
                ("alfa_min", "Courant factor, lower", "number", "-", REQ, "always"),
                ("flag_adaptive_timestepping", "Adaptive time stepping", "switch", "0 / 1", REQ, "always"),
            ]},
            {"title": "Computing", "inputs": [
                ("flag_GPU", "Run on the GPU", "switch", "0 / 1", REQ, "always"),
                ("flag_single", "Single precision", "switch", "0 / 1", REQ, "always"),
            ]},
            {"title": "Recording intervals", "inputs": [
                ("record_time_maps", "Save maps every", "number", "min", REQ, "always"),
                ("record_time_hydrographs", "Save hydrographs every", "number", "min", REQ, "always"),
                ("record_time_spatial_rainfall", "Save rainfall maps every", "number", "min", REQ, "spatial rainfall on"),
                ("record_time_spatial_ETP", "Save ET maps every", "number", "min", REQ, "ET on"),
            ]},
            {"title": "Outputs and display", "inputs": [
                ("flag_export_maps", "Export maps", "switch", "0 / 1", REQ, "always"),
                ("flag_export_groundwater_maps", "Export groundwater maps", "switch", "0 / 1", 0, "groundwater on"),
                ("flag_wse", "Water-surface elevation maps", "switch", "0 / 1", REQ, "always"),
                ("depth_wse", "Minimum depth for those maps", "number", "m", REQ, "water-surface maps on"),
                ("flag_waterbalance", "Water-balance report", "switch", "0 / 1", REQ, "always"),
                ("flag_dashboard", "Live dashboard", "switch", "0 / 1", REQ, "always"),
                ("flag_elapsed_time", "Show elapsed time", "switch", "0 / 1", REQ, "always"),
            ]},
            {"title": "Warm start", "inputs": [
                ("flag_warmup", "Start from a warm-up state", "switch", "0 / 1", REQ, "always"),
            ]},
            {"title": "Advanced (keep the defaults unless you know why)", "style": "advanced", "inputs": [
                ("alfa_max", "Courant factor, upper", "number", "-", REQ, "always"),
                ("time_step_change", "Time-step refresh increment", "number", "s", REQ, "always"),
                ("flag_timestep", "Full time-step stability check", "switch", "0 / 1", REQ, "always"),
            ]},
        ],
    },
    {
        "name": "Input maps",
        "purpose": "Paths to the GeoTIFF maps. Leave optional maps empty to use the class tables or scalars.",
        "wide_value": True,
        "blocks": [
            {"title": "Required maps", "inputs": [
                ("DEM_path", "Terrain elevation (DEM)", "path", "m", REQ, "always"),
                ("LULC_path", "Land use / land cover classes", "path", "class", REQ, "always"),
                ("SOIL_path", "Soil classes", "path", "class", REQ, "always"),
            ]},
            {"title": "Land surface", "inputs": [
                ("Albedo_path", "Albedo", "path", "0 to 1", "", "spatial albedo on"),
                ("LAI_path", "Leaf area index", "path", "m2/m2", "", "abstraction on"),
                ("ImperviousFraction_path", "Impervious fraction", "path", "0 to 1", "", "optional"),
            ]},
            {"title": "Rivers and terrain", "inputs": [
                ("RiverWidths_path", "River width", "path", "m", "", "channel sub-grid with river maps"),
                ("RiverDepths_path", "River depth", "path", "m", "", "channel sub-grid with river maps"),
            ]},
            {"title": "Initial states", "inputs": [
                ("Warmup_Depth_path", "Initial surface water depth", "path", "m", "", "warm start on"),
                ("Initial_Soil_Moisture_path", "Initial soil water", "path", "see note", "", "warm start on"),
                ("GW_table_path", "Initial water-table elevation", "path", "m", "", "groundwater on"),
                ("Initial_SWE_path", "Initial snow water equivalent", "path", "mm", "", "snow on"),
                ("Initial_Snow_Depth_path", "Initial snow depth", "path", "mm", "", "snow on"),
            ]},
            {"title": "Soil-property maps (override the SOIL classes)", "inputs": [
                ("ThetaS_path", "Soil theta_s", "path", "m3/m3", "", "optional"),
                ("ThetaR_path", "Soil theta_r", "path", "m3/m3", "", "optional"),
                ("AlphaVG_path", "Soil van Genuchten alpha", "path", "1/m", "", "optional"),
                ("NVG_path", "Soil van Genuchten n", "path", "-", "", "optional"),
                ("KsatNearSurface_path", "Ksat, near-surface layer", "path", "mm/h", "", "optional"),
                ("KsatRootZone_path", "Ksat, root-zone layer", "path", "mm/h", "", "optional"),
                ("KsatTransmission_path", "Ksat, transmission layer", "path", "mm/h", "", "optional"),
                ("GroundwaterKsat_path", "Aquifer Ksat (kgsat)", "path", "mm/h", "", "optional"),
                ("SpecificYield_path", "Aquifer specific yield", "path", "-", "", "legacy groundwater"),
            ]},
            {"title": "Water-quality maps", "inputs": [
                ("Initial_Buildup_path", "Initial pollutant mass", "path", "kg/cell", "", "initial build-up on"),
                ("B1_path", "Build-up asymptote C1", "path", "see note", "", "water quality with maps"),
                ("B2_path", "Build-up rate C2", "path", "see note", "", "water quality with maps"),
                ("W1_path", "Wash-off coefficient C3", "path", "see note", "", "water quality with maps"),
                ("W2_path", "Wash-off exponent C4", "path", "see note", "", "water quality with maps"),
            ]},
        ],
    },
    {
        "name": "Soil and groundwater",
        "purpose": "Infiltration, evapotranspiration, groundwater, the layered subsurface and the deep aquifer. Each map sits under the value it overrides.",
        "wide_value": True,
        "blocks": [
            {"title": "Infiltration and evapotranspiration", "inputs": [
                ("flag_infiltration", "Infiltration", "switch", "0 / 1", REQ, "always"),
                ("flag_ETP", "Evapotranspiration", "switch", "0 / 1", REQ, "always"),
                ("flag_abstraction", "Canopy interception", "switch", "0 / 1", REQ, "always"),
                ("flag_capillary_rise", "Capillary rise", "switch", "0 / 1", 1, "groundwater on"),
                ("flag_snow_modeling", "Snow", "switch", "0 / 1", REQ, "always"),
                ("flag_spatial_albedo", "Albedo from a map", "switch", "0 / 1", REQ, "snow or ET on"),
                ("flag_neglect_infiltration_river", "No infiltration in river cells", "switch", "0 / 1", REQ, "infiltration on"),
            ]},
            {"title": "Groundwater", "switch": "flag_groundwater_modeling", "inputs": [
                ("flag_groundwater_modeling", "Groundwater", "switch", "0 / 1", REQ, "switch for this block"),
                ("flag_baseflow", "Lateral groundwater flow", "switch", "0 / 1", REQ, "groundwater on"),
                ("DTB_path", "Depth to bedrock map", "path", "m", "", "legacy: required; layered: optional"),
                ("GW_Dirichlet_Head_path", "Fixed groundwater head map", "path", "m", "", "groundwater on"),
                ("flag_prescribed_recharge", "Prescribed recharge series", "switch", "0 / 1", 0, "groundwater on"),
                ("Prescribed_Recharge_Timeseries_File", "Prescribed recharge file", "path", "file", "", "prescribed recharge on"),
                ("groundwater_courant", "Groundwater Courant factor", "number", "-", 0.25, "lateral flow on"),
                ("groundwater_relative_tolerance", "Lateral clock relative tolerance", "number", "-", 1e-3, "lateral flow on"),
                ("groundwater_head_tolerance_m", "Lateral clock head tolerance", "number", "m", 1e-3, "lateral flow on"),
                ("riverbed_thickness_m", "Riverbed thickness", "number", "m", 0.5, "groundwater on"),
            ]},
            {"title": "Layered subsurface: column", "switch": "flag_layered_subsurface", "inputs": [
                ("flag_layered_subsurface", "Layered subsurface", "switch", "0 / 1", 0, "switch for these blocks"),
                ("subsurface_soil_base_depth_m", "Soil thickness", "number", "m", 2.0, "layered on"),
                ("SoilBaseDepth_path", "Soil thickness map", "path", "m", "", "layered on (wins where it has data)"),
                ("subsurface_depth_to_bedrock_m", "Depth to bedrock", "number", "m", "", "layered, where the DTB map is empty"),
                ("subsurface_max_depth_m", "Deepest column", "number", "m", 50.0, "layered on"),
                ("subsurface_deep_first_thickness_m", "First rock layer thickness", "number", "m", 0.25, "layered on"),
                ("subsurface_deep_growth_ratio", "Rock layer growth ratio", "number", "-", 2.0, "layered on"),
                ("subsurface_deep_max_thickness_m", "Thickest rock layer", "number", "m", 4.0, "layered on"),
            ]},
            {"title": "Layered subsurface: rock properties", "controlled_by": "flag_layered_subsurface", "inputs": [
                ("regolith_theta_s", "Rock theta_s", "number", "m3/m3", 0.35, "layered on"),
                ("RegolithThetaS_path", "Rock theta_s map", "path", "m3/m3", "", "layered on (wins where it has data)"),
                ("regolith_theta_r", "Rock theta_r", "number", "m3/m3", 0.05, "layered on"),
                ("RegolithThetaR_path", "Rock theta_r map", "path", "m3/m3", "", "layered on (wins where it has data)"),
                ("regolith_alpha_vg_m_inverse", "Rock van Genuchten alpha", "number", "1/m", 2.0, "layered on"),
                ("RegolithAlphaVG_path", "Rock van Genuchten alpha map", "path", "1/m", "", "layered on (wins where it has data)"),
                ("regolith_n_vg", "Rock van Genuchten n", "number", "-", 1.5, "layered on"),
                ("RegolithNVG_path", "Rock van Genuchten n map", "path", "-", "", "layered on (wins where it has data)"),
                ("regolith_ksat_v_mm_h", "Rock vertical Ksat", "number", "mm/h", 5.0, "layered on"),
                ("RegolithKsatV_path", "Rock vertical Ksat map", "path", "mm/h", "", "layered on (wins where it has data)"),
                ("regolith_ksat_h_mm_h", "Rock horizontal Ksat", "number", "mm/h", "", "layered on (empty: kgsat)"),
                ("RegolithKsatH_path", "Rock horizontal Ksat map", "path", "mm/h", "", "layered on (wins where it has data)"),
                ("regolith_specific_yield", "Rock specific yield", "number", "-", "", "layered on (empty: from curve)"),
                ("RegolithSpecificYield_path", "Rock specific yield map", "path", "-", "", "layered on (wins where it has data)"),
            ]},
            {"title": "Layered subsurface: soil and numerics", "controlled_by": "flag_layered_subsurface", "inputs": [
                ("soil_lateral_anisotropy", "Soil Kh / Kv", "number", "-", 1.0, "layered on"),
                ("capillary_field_capacity_head_m", "Field-capacity suction", "number", "m", 3.3, "layered on"),
                ("subsurface_drainage_substeps", "Drainage sub-steps", "number", "-", 1, "layered on"),
                ("groundwater_target_dt_s", "Lateral update interval", "number", "s", 3600, "layered on"),
                ("groundwater_min_dt_s", "Shortest lateral interval", "number", "s", 60, "layered on"),
                ("flag_implicit_water_table_flux", "Implicit water-table flux", "switch", "0 / 1", 1, "layered on"),
            ]},
            {"title": "Lateral conductivity multipliers", "controlled_by": "flag_layered_subsurface", "inputs": [
                ("soil_lateral_conductivity_multiplier", "Soil Kh multiplier", "number", "-", 1.0, "layered on"),
                ("regolith_lateral_conductivity_multiplier", "Rock Kh multiplier", "number", "-", 1.0, "layered on"),
            ]},
            {"title": "Feddes root-water stress", "switch": "flag_feddes", "inputs": [
                ("flag_feddes", "Feddes stress", "switch", "0 / 1", 0, "switch for this block"),
                ("feddes_h1_m", "h1, anaerobiosis", "number", "m", -0.1, "Feddes on"),
                ("feddes_h2_m", "h2, optimal from", "number", "m", -0.25, "Feddes on"),
                ("feddes_h3_m", "h3, stress from", "number", "m", -4.0, "Feddes on"),
                ("feddes_h3_low_m", "h3 at low demand", "number", "m", "", "Feddes on (empty: h3)"),
                ("feddes_h3_high_m", "h3 at high demand", "number", "m", "", "Feddes on (empty: h3)"),
                ("feddes_h4_m", "h4, wilting point", "number", "m", -80.0, "Feddes on"),
            ]},
            {"title": "Deep aquifer", "switch": "flag_deep_aquifer", "inputs": [
                ("flag_deep_aquifer", "Deep aquifer", "switch", "0 / 1", 0, "switch for this block"),
                ("DeepAquiferSubbasins_path", "Deep-aquifer subbasins map", "path", "ID", "", "deep aquifer on: required"),
                ("DeepAquiferChannelWeight_path", "Deep-aquifer return weights map", "path", "-", "", "deep aquifer on: required"),
                ("deep_recharge_fraction", "Share of recharge to the deep aquifer", "number", "-", 0.05, "deep aquifer on"),
                ("deep_aquifer_tau_days", "Response time", "number", "days", 180, "deep aquifer on"),
            ]},
            {"title": "Deep aquifer: starting storage (give one, or none to start empty)", "controlled_by": "flag_deep_aquifer", "inputs": [
                ("deep_aquifer_initial_baseflow_mm_day", "Starting deep baseflow (from a gauge)", "number", "mm/day", "", "option 1"),
                ("deep_aquifer_mean_recharge_mm_year", "Long-term mean recharge (no gauge)", "number", "mm/year", "", "option 2: steady state"),
                ("DeepAquiferInitialState_path", "Calibrated starting storage (JSON)", "path", "file", "", "option 3: from the Builder"),
            ]},
        ],
    },
    {
        "name": "Surface flow",
        "purpose": "Routing scheme, outlets, river channels, terrain treatment and structures.",
        "blocks": [
            {"title": "Routing scheme", "inputs": [
                ("flag_inertial", "Local-inertial routing", "switch", "0 / 1", REQ, "always"),
                ("flag_numerical_scheme", "Inertial discretization", (1, 2, 3), "1 / 2 / 3", REQ, "inertial on"),
                ("flag_CA", "Cellular-automata routing", "switch", "0 / 1", REQ, "always"),
                ("flag_kinematic", "Kinematic wave", "switch", "0 / 1", REQ, "always"),
                ("flag_diffusive", "Diffusive wave", "switch", "0 / 1", REQ, "always"),
                ("flag_D8", "D8 flow directions", "switch", "0 / 1", REQ, "always"),
            ]},
            {"title": "Outlets and boundaries", "inputs": [
                ("flag_boundary", "Open boundary treatment", "switch", "0 / 1", REQ, "always"),
                ("flag_outlet_type", "Outlet relation", (1, 2), "1 / 2", REQ, "always"),
                ("flag_critical", "Critical-flow outlet", "switch", "0 / 1", REQ, "always"),
                ("slope_outlet", "Outlet energy slope", "number", "m/m", REQ, "normal-flow outlet"),
                ("n_outlets_data", "Number of outlet cells", "number", "-", 1, "always"),
            ]},
            {"title": "River channel (sub-grid)", "inputs": [
                ("flag_neal_channel", "Neal (2012) channel sub-grid", "switch", "0 / 1", REQ, "always"),
                ("flag_river_rasters", "River width and depth from maps", "switch", "0 / 1", REQ, "channel sub-grid on"),
                ("Manning", "Channel Manning n", "number", "s m^-1/3", REQ, "always"),
                ("alfa_1", "Depth law coefficient", "number", "-", REQ, "channel from drainage area"),
                ("alfa_2", "Depth law exponent", "number", "-", REQ, "channel from drainage area"),
                ("beta_1", "Width law coefficient", "number", "-", REQ, "channel from drainage area"),
                ("beta_2", "Width law exponent", "number", "-", REQ, "channel from drainage area"),
            ]},
            {"title": "Terrain treatment", "inputs": [
                ("flag_resample", "Resample the grid", "switch", "0 / 1", REQ, "always"),
                ("resolution_resample", "Model resolution", "number", "m", REQ, "resample on"),
                ("flag_fill_DEM", "Fill sinks", "switch", "0 / 1", REQ, "always"),
                ("flag_smoothening", "Smooth river profiles", "switch", "0 / 1", REQ, "always"),
                ("min_area", "Stream threshold area", "number", "km2", REQ, "smoothing on"),
                ("K_value", "Smoothing strength", "number", "0 to 20", REQ, "smoothing on"),
                ("tau", "Smoothing regularization", "number", "0 to 1", REQ, "smoothing on"),
                ("sl", "Minimum imposed slope", "number", "m/m", REQ, "smoothing on"),
                ("flag_smooth_cells", "Smooth cells", "switch", "0 / 1", REQ, "always"),
                ("flag_trunk", "Trunk streams only", "switch", "0 / 1", REQ, "always"),
                ("flag_reduce_DEM", "Reduce the DEM", "switch", "0 / 1", REQ, "always"),
                ("flag_DTM", "DTM filter", "switch", "0 / 1", REQ, "always"),
                ("slope_DTM", "DTM filter slope", "number", "%", REQ, "DTM filter on"),
            ]},
            {"title": "Unstructured mesh", "inputs": [
                ("flag_voronoi", "Voronoi mesh", "switch", "0 / 1", 0, "always"),
                ("flag_voronoi_subgrid", "Voronoi sub-grid tables", "switch", "0 / 1", REQ, "Voronoi on"),
            ]},
            {"title": "Inflows and structures", "inputs": [
                ("flag_inflow", "Inflow hydrographs", "switch", "0 / 1", REQ, "always"),
                ("flag_stage_hydrograph", "Stage hydrographs", "switch", "0 / 1", REQ, "always"),
                ("flag_reservoir", "Reservoirs", "switch", "0 / 1", REQ, "always"),
            ]},
            {"title": "Paused (keep at 0: setting 1 stops the run)", "style": "paused", "inputs": [
                ("flag_structured_lookup_subgrid", "Lookup-table sub-grid (paused)", "switch", "0 / 1", 0, "paused"),
                ("Subgrid_DEM_path", "High-resolution DEM for that sub-grid (paused)", "path", "m", "", "paused"),
            ]},
        ],
    },
    {
        "name": "Rain and ET",
        "purpose": "Where rainfall and evapotranspiration come from, and the lists of map files.",
        "wide_value": True,
        "blocks": [
            {"title": "Rainfall source", "inputs": [
                ("flag_rainfall", "Rainfall", "switch", "0 / 1", REQ, "always"),
                ("flag_spatial_rainfall", "Spatially varied rainfall", "switch", "0 / 1", REQ, "rainfall on"),
                ("flag_input_rainfall_map", "Rainfall from the map list below", "switch", "0 / 1", REQ, "rainfall on"),
                ("flag_satellite_rainfall", "Satellite rainfall", "switch", "0 / 1", REQ, "rainfall on"),
                ("flag_real_time_satellite_rainfall", "Real-time satellite rainfall", "switch", "0 / 1", REQ, "rainfall on"),
            ]},
            {"title": "Evapotranspiration source", "inputs": [
                ("flag_input_ETP_map", "ET from the map lists below", "switch", "0 / 1", REQ, "ET on"),
            ]},
            {"title": "Advanced (keep the defaults unless you know why)", "style": "advanced", "inputs": [
                ("flag_rainfall_multiple_runs", "Multiple rainfall runs (legacy workflow)", "switch", "0 / 1", REQ, "rainfall on"),
            ]},
        ],
        "lists": [
            {"title": "Rainfall maps", "anchor": "Satellite or Radar Rainfall",
             "columns": ["Time [min]", "Raster Directory with values in mm/h"],
             "note": "One row per rainfall map: minutes since the start date and the GeoTIFF path (mm/h). Read when flag_input_rainfall_map or a satellite rainfall switch is 1."},
            {"title": "Transpiration maps", "anchor": "Satellite transpiration",
             "columns": ["Time [day]", "Raster Directory with values in mm/day"],
             "note": "One row per transpiration map: days since the start date and the GeoTIFF path (mm/day). Read when flag_input_ETP_map = 1; same times as the evaporation list."},
            {"title": "Evaporation maps", "anchor": "Satellite Evaporation",
             "columns": ["Time [day]", "Raster Directory with values in mm/day"],
             "note": "One row per evaporation map: days since the start date and the GeoTIFF path (mm/day). Read when flag_input_ETP_map = 1; same times as the transpiration list."},
        ],
    },
    {
        "name": "Design storm",
        "purpose": "Synthetic storms from a local IDF relation i = K RP^a / (b + duration)^c.",
        "blocks": [
            {"title": "Hyetograph", "inputs": [
                ("flag_alternated_blocks", "Alternated blocks", "switch", "0 / 1", REQ, "design storm"),
                ("flag_huff", "Huff curves", "switch", "0 / 1", REQ, "design storm"),
                ("RP", "Return period", "number", "years", REQ, "design storm"),
                ("Rainfall Duration", "Storm duration", "number", "min", REQ, "design storm"),
                ("dt_design", "Hyetograph time step", "number", "min", REQ, "design storm"),
            ]},
            {"title": "IDF relation", "inputs": [
                ("K", "IDF scale K", "number", "see note", REQ, "design storm"),
                ("a", "IDF return-period exponent a", "number", "-", REQ, "design storm"),
                ("b", "IDF duration offset b", "number", "min", REQ, "design storm"),
                ("c", "IDF duration exponent c", "number", "-", REQ, "design storm"),
            ]},
        ],
    },
    {
        "name": "Water quality",
        "purpose": "Pollutant build-up and wash-off.",
        "blocks": [
            {"title": "Water quality", "switch": "flag_waterquality", "inputs": [
                ("flag_waterquality", "Water quality", "switch", "0 / 1", REQ, "switch for this block"),
                ("flag_WQ_Rasters", "Parameters from maps", "switch", "0 / 1", REQ, "water quality on"),
                ("flag_initial_buildup", "Initial build-up from a map", "switch", "0 / 1", REQ, "water quality on"),
                ("ADD", "Antecedent dry days", "number", "days", REQ, "water quality on"),
                ("min_Bt", "Minimum surface mass", "number", "g/m2", REQ, "water quality on"),
                ("Bmin", "Build-up lower bound", "number", "g/m2", REQ, "water quality on"),
                ("Bmax", "Build-up upper bound", "number", "g/m2", REQ, "water quality on"),
                ("Pol_min", "Map display threshold", "number", "mg/L", REQ, "water quality on"),
            ]},
            {"title": "Advanced (keep the defaults unless you know why)", "style": "advanced", "controlled_by": "flag_waterquality", "inputs": [
                ("flag_wq_model", "Formulation (0 = mass-based build-up / wash-off)", "number", "0 / ...", REQ, "water quality on"),
            ]},
        ],
    },
    {
        "name": "Risk and gauges",
        "purpose": "Human-instability risk and observation points.",
        "blocks": [
            {"title": "Human instability", "inputs": [
                ("flag_human_instability", "Human-instability risk", (0, 1, 3), "0 / 1 / 3", REQ, "always"),
            ]},
            {"title": "Observation points", "inputs": [
                ("flag_obs_gauges", "Record at the gauges below", "switch", "0 / 1", REQ, "always"),
            ]},
        ],
        "lists": [
            {"title": "Gauges", "anchor": "Gauge",
             "columns": ["Gauge", "Easting (m)", "Northing (m)", "Label Name"],
             "optional_columns": [
                 {"name": "Cross-section width (cells)", "default": 5,
                  "note": "Discharge cross-section: number of cells in a line across the flow, centred on the gauge cell. Odd numbers keep the gauge in the middle (an even one is raised by 1). Empty: 5. Depth and water level stay at the gauge cell."},
                 {"name": "Downstream direction", "default": "auto", "options": "auto,N,S,E,W",
                  "note": "Which way the water leaves the section. auto: downhill from the terrain around the gauge. N or S gives an east-west line; E or W a north-south line. Set it where the terrain is flat. Empty: auto."},
             ],
             "header_is_anchor": True,
             "note": "One row per observation point: a number, its coordinates in the DEM's projection (m) and a label. Read when flag_obs_gauges = 1. The discharge is summed across a cross-section centred on the gauge (last two columns)."},
        ],
    },
]

# Class tables and time series, one tab each. Every tab is a cell-for-cell copy of the old
# separate workbook (same rows and columns), so readtable gives the model the same table;
# the builder only styles it and adds notes. A series is copied in full only when the
# case's switches read it ("read_when": all "on" switches = 1, no "off" switch = 1);
# otherwise only its header rows are kept, as a template. "files": where the old
# workbook may sit, relative to the source workbook's folder (then the repository's
# Input_Data_Sheets as a template). "header_rows": rows styled as column headers.
TABLE_TABS = [
    {"name": "Soil classes", "kind": "class", "files": ["SOIL_parameters.xlsx"], "header_rows": [2],
     "purpose": "One row per soil class of the SOIL map: infiltration, van Genuchten and groundwater parameters.",
     "note": "Soil class table (was SOIL_parameters.xlsx). Keep the rows and the column order: the model reads the columns by position. Add a class as a new row."},
    {"name": "Land use classes", "kind": "class", "files": ["LULC_parameters.xlsx"], "header_rows": [2],
     "purpose": "One row per land-use class of the LULC map: roughness, losses, build-up and wash-off, roots, crop and snow parameters.",
     "note": "Land-use class table (was LULC_parameters.xlsx). Keep the rows and the column order: the model reads it from row 2 by position."},
    {"name": "Human instability", "kind": "class", "files": ["human_risk.xlsx", "Human_risk.xlsx", "Outdated/Human_risk.xlsx", "Outdated/human_risk.xlsx"],
     "header_rows": [2],
     "purpose": "Body and water parameters of the three human-instability methods.",
     "note": "Human-instability parameters (was human_risk.xlsx). Read when flag_human_instability = 1 (method 1) or 3 (method 3). Keep every cell in place: the model reads the values by position."},
    {"name": "Mesh", "kind": "class", "files": ["Voronoi_Settings.xlsx"], "header_rows": [4],
     "purpose": "Voronoi and quadtree mesh settings.",
     "note": "Mesh settings (was Voronoi_Settings.xlsx). Read when flag_voronoi = 1. Each value sits right of its Parameter name."},
    {"name": "Control structures", "kind": "class", "files": ["BC_Control_structure.xlsx", "BC_Control_Structure.xlsx"], "header_rows": [2, 3], "template_rows": 3,
     "purpose": "Internal boundary conditions: one row per reservoir or structure with its two outflow laws.",
     "note": "Control structures (was BC_Control_structure.xlsx). Read when flag_reservoir = 1. Outflow Q = k1 (h - h01)^k2 to cell (xd1, yd1) and Q = k3 (h - h02)^k4 to cell (xd2, yd2); coordinates in the DEM projection (m).",
     "column_notes": {"Index": "Structure number.", "x(m)": "Easting of the structure cell (m, DEM projection).", "y(m)": "Northing of the structure cell (m).",
                      "k1": "Coefficient of the first outflow law.", "h01 [m]": "Threshold depth of the first law (m).", "k2": "Exponent of the first law.",
                      "xd1(m)": "Easting of the cell receiving the first outflow (m).", "yd1(m)": "Northing of that cell (m).",
                      "k3": "Coefficient of the second outflow law.", "h02 [m]": "Threshold depth of the second law (m).", "k4": "Exponent of the second law.",
                      "xd2(m)": "Easting of the cell receiving the second outflow (m).", "yd2(m)": "Northing of that cell (m)."}},
    {"name": "Rain series", "kind": "series", "files": ["Rainfall_Intensity_Data.xlsx"], "header_rows": [2], "template_rows": 2,
     "read_when": {"on": ["flag_rainfall"], "off": ["flag_alternated_blocks", "flag_huff", "flag_input_rainfall_map", "flag_spatial_rainfall",
                                                    "flag_real_time_satellite_rainfall", "flag_satellite_rainfall", "flag_rainfall_multiple_runs"]},
     "purpose": "One rain gauge: rainfall intensity over time.",
     "note": "Rain series (was Rainfall_Intensity_Data.xlsx). Read when rainfall is on and no design storm, map, satellite or spatial rainfall is chosen.",
     "column_notes": {"Time(min)": "Minutes from the start of the run, one row per step. The step is the difference of the first two rows.",
                      "Rainfall Intensity": "Rainfall intensity over the step (mm/h)."}},
    {"name": "Rain gauges", "kind": "series", "files": ["Rainfall_Spatial_Input.xlsx"], "header_rows": [3, 4, 5, 6, 7], "template_rows": 7,
     "read_when": {"on": ["flag_rainfall", "flag_spatial_rainfall"], "off": ["flag_real_time_satellite_rainfall", "flag_satellite_rainfall", "flag_input_rainfall_map"]},
     "purpose": "Several rain gauges, interpolated over the domain.",
     "note": "Spatial rain gauges (was Rainfall_Spatial_Input.xlsx). Read when flag_spatial_rainfall = 1 and no rainfall maps are given. One column per gauge: its index, name and coordinates (m, DEM projection), then its intensity (mm/h) at each time (min) of column B."},
    {"name": "Climate stations", "kind": "series", "files": ["ETP_input_data.xlsx"], "header_rows": [2, 3], "template_rows": 3,
     "read_when": {"on": ["flag_ETP"], "off": ["flag_input_ETP_map"]},
     "purpose": "Daily weather at one or more stations, for the reference evapotranspiration.",
     "note": "Climate stations (was ETP_input_data.xlsx). Read when flag_ETP = 1 and no ET maps are given. Row 2 names each station block (index, easting and northing in m); row 3 heads its columns; one row per day from row 4."},
    {"name": "Inflow hydrographs", "kind": "series", "files": ["Inflow_Hydrograph.xlsx"], "header_rows": [2, 3], "template_rows": 3,
     "read_when": {"on": ["flag_inflow"], "off": []},
     "purpose": "Inflow hydrographs entering the domain, up to 50 gauges.",
     "note": "Inflow hydrographs (was Inflow_Hydrograph.xlsx). Read when flag_inflow = 1. Five columns per gauge, 50 blocks: time (min), flow (m3/s), and the easting and northing (m) of each cell the flow enters, one cell per row. Keep all 50 blocks, even empty."},
    {"name": "Stage hydrographs", "kind": "series", "files": ["Stage_Hydrograph.xlsx"], "header_rows": [2, 3], "template_rows": 3,
     "read_when": {"on": ["flag_stage_hydrograph"], "off": []},
     "purpose": "Water depths imposed at boundary cells, up to 50 gauges.",
     "note": "Stage hydrographs (was Stage_Hydrograph.xlsx). Read when flag_stage_hydrograph = 1. Five columns per gauge, 50 blocks: time (min), depth (m), and the easting and northing (m) of each cell, one per row. Keep all 50 blocks, even empty."},
]

# Workbook keys left out on purpose: General_Data cells no reader asks for, and keys the
# model derives itself.
EXCLUDED = {
    "ro_person": "human_risk.xlsx holds the human-instability parameters",
    "weight_person": "human_risk.xlsx", "height_person": "human_risk.xlsx",
    "width1_person": "human_risk.xlsx", "width2_person": "human_risk.xlsx", "ro_water": "human_risk.xlsx",
    "flag_subgrid": "derived from the three channel sub-grid switches",
    "flag_overbanks": "derived from flag_neal_channel",
    "terrain_runtime": "only in a code comment", "NDVI_path": "its reader line is commented out",
    "deep_aquifer_response_time_days": "alias of deep_aquifer_tau_days",
    "River_K_coeff": "read, then deactivated in input_data_script",
    "time_save_ETP": "read, then not used by a model run",
    "flag_river_heigth_compensation": "used only by old code a run does not call",
    "flag_dam_break": "used only by old main loops (ensemble, synthetic, workspace)",
    "flag_data_source": "used only by old code a run does not call",
}

# Reader changes the new workbook needs (input_data_script.m, HydroPol2D_preprocessing.m).
# Reader changes the workbook needed; all done on 2026-10-20 (hp2d_read_settings_grid,
# hp2d_read_flags_grid, hp2d_read_case_table, requiredFlags, xlblock_2col).
READER_FOLLOWUPS = []
