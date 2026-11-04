<p align="center">
  <a href="https://marcusnobrega-eng.github.io/HydroPol2D-docs/">
    <img src="https://marcusnobrega-eng.github.io/HydroPol2D-docs/img/social_card.png" width="720" alt="HydroPol2D: a distributed hydrologic-hydrodynamic framework for coupled surface-subsurface processes" />
  </a>
</p>

<p align="center">
  <a href="https://github.com/marcusnobrega-eng/HydroPol2D/stargazers"><img src="https://img.shields.io/github/stars/marcusnobrega-eng/HydroPol2D?style=social" alt="GitHub stars" /></a>
  <a href="https://github.com/marcusnobrega-eng/HydroPol2D/forks"><img src="https://img.shields.io/github/forks/marcusnobrega-eng/HydroPol2D?style=social" alt="GitHub forks" /></a>
  <a href="https://github.com/marcusnobrega-eng/HydroPol2D/watchers"><img src="https://img.shields.io/github/watchers/marcusnobrega-eng/HydroPol2D?style=social" alt="GitHub watchers" /></a>
</p>

<p align="center">
  <a href="https://github.com/marcusnobrega-eng/HydroPol2D/tags"><img src="https://img.shields.io/github/v/tag/marcusnobrega-eng/HydroPol2D?sort=semver&include_prereleases&label=version&color=0a6e4f" alt="Latest version" /></a>
  <a href="https://github.com/marcusnobrega-eng/HydroPol2D/commits/main"><img src="https://img.shields.io/github/last-commit/marcusnobrega-eng/HydroPol2D/main?label=last%20commit" alt="Last commit" /></a>
  <a href="https://github.com/marcusnobrega-eng/HydroPol2D/commits/main"><img src="https://img.shields.io/github/commit-activity/m/marcusnobrega-eng/HydroPol2D?label=commits" alt="Commit activity" /></a>
  <a href="https://github.com/marcusnobrega-eng/HydroPol2D/issues"><img src="https://img.shields.io/github/issues/marcusnobrega-eng/HydroPol2D" alt="Open issues" /></a>
  <a href="LICENSE"><img src="https://img.shields.io/github/license/marcusnobrega-eng/HydroPol2D" alt="License: GPL-3.0" /></a>
  <img src="https://img.shields.io/badge/MATLAB-R2020a%2B-e16737" alt="MATLAB R2020a or newer" />
  <a href="https://marcusnobrega-eng.github.io/HydroPol2D-docs/"><img src="https://img.shields.io/badge/docs-website-2f6fb5" alt="Documentation" /></a>
  <a href="https://github.com/marcusnobrega-eng/HydroPol2D-Studio-releases/releases"><img src="https://img.shields.io/badge/Studio-desktop%20app-6b3fa0" alt="HydroPol2D Studio" /></a>
  <a href="https://doi.org/10.1016/j.jhydrol.2023.129982"><img src="https://img.shields.io/badge/paper-J.%20Hydrology%202023-555555" alt="HydroPol2D paper" /></a>
</p>

<p align="center">
  <a href="https://marcusnobrega-eng.github.io/HydroPol2D-docs/"><b>Documentation</b></a> ·
  <a href="#quick-start"><b>Quick start</b></a> ·
  <a href="#examples-you-can-run"><b>Examples</b></a> ·
  <a href="#gallery"><b>Gallery</b></a> ·
  <a href="#the-hydropol2d-family"><b>Studio and Case Builder</b></a> ·
  <a href="#licenses"><b>Licenses</b></a> ·
  <a href="#cite"><b>Cite</b></a>
</p>

HydroPol2D is a computer model of floods and water quality. It follows rain across a landscape,
cell by cell. It computes how much water plants and snow hold, how much soaks into the soil and
reaches the groundwater, and where the rest flows, how deep it gets and what pollutants it
carries. It works for cities, suburbs and rural catchments, including places with few
measurements. This repository is the MATLAB version of the model. It is free and open source.

<p align="center">
  <a href="https://marcusnobrega-eng.github.io/profile/files/stanford-campus-coupled-states.mp4">
    <img src=".github/assets/stanford-campus-10m-coupled-states.gif" width="900" alt="Animated maps of eight synchronized HydroPol2D outputs for a 10 m simulation of Stanford campus" />
  </a>
  <br />
  <sub>Stanford campus, CA, simulated on a 10 m grid. Top: flood depth, flow speed, flood hazard and water level. Bottom: infiltration, rainfall, groundwater depth and groundwater recharge. Click for the full video.</sub>
</p>

## What it does

<p align="center">
  <img src="https://marcusnobrega-eng.github.io/HydroPol2D-docs/img/getting-started/questions-modules-outputs.svg" width="760" alt="Six questions, the part of the model that answers each and the outputs that show it: surface flow, rain to runoff, snow, groundwater, people at risk and water quality" />
</p>

| Area | What HydroPol2D computes |
|---|---|
| **Flood flow** | Water flowing over the land and in rivers. Five ways to compute it are available. The default is local inertial. |
| **Rain to runoff** | Rain caught by leaves, snow that builds up and melts, evaporation and plant water use, and water soaking into three soil layers. |
| **Groundwater** | The water table, groundwater flow between cells, a deep aquifer and water exchange with rivers. |
| **Structures** | Reservoirs and gates, dam breaks, storm drains, detention ponds and green infrastructure such as rain gardens. |
| **Water quality and risk** | Pollutants that build up on surfaces and wash off with runoff, and the risk for people standing in floodwater. |
| **Grids and computers** | A regular grid or a flexible mesh. It runs on a normal computer or on an NVIDIA graphics card, and its time step adapts so that water is conserved. |
| **Outputs** | Maps (GeoTIFF), animations (MP4), hydrographs and water balances (CSV), and figures (PDF and PNG). |

## Quick start

```bash
git clone https://github.com/marcusnobrega-eng/HydroPol2D.git
```

Open MATLAB in the `HydroPol2D` folder and run an example:

```matlab
run('Examples/VTilted_Excel/run_excel.m')   % a simple synthetic catchment, about 3 min
run('Examples/Stanford/run_excel.m')        % San Francisquito Creek, CA, about 14 min
```

**You need** MATLAB R2020a or newer with the Mapping Toolbox and the Image Processing Toolbox.
Runs on a graphics card also need the Parallel Computing Toolbox and an NVIDIA card. Terrain
smoothing needs the Optimization Toolbox.

**Your own case** is one Excel workbook. Copy
[`Input_Data_Sheets/HydroPol2D_Inputs.xlsx`](Input_Data_Sheets/HydroPol2D_Inputs.xlsx), fill in
its tabs and run [`run_HydroPol2D.m`](run_HydroPol2D.m). As supplied, the workbook runs the
30 m San Francisquito Creek example on the CPU
([`Input_Data_Sheets/Readme.txt`](Input_Data_Sheets/Readme.txt)). The website explains every tab, with a
picture of each one, in
[Case Workbook Settings](https://marcusnobrega-eng.github.io/HydroPol2D-docs/manual/configuration/general_data)
and [Run from Excel](https://marcusnobrega-eng.github.io/HydroPol2D-docs/manual/getting-started/excel-workflow).
For automated studies, a script can replace the workbook
([`Config/input_data_bypass_script.m`](Config/input_data_bypass_script.m)).

## Examples you can run

| Example | Run in MATLAB | What it shows | Time | Expected results |
|---|---|---|---|---|
| [V-tilted catchment](Examples/VTilted_Excel/README.md) (synthetic) | `run('Examples/VTilted_Excel/run_excel.m')` | 10.8 mm/h of rain for 90 min, 180 min simulated, on a regular grid or a mesh | ~3 min | rain 16.20 mm, peak outflow 1.525 m³/s, runoff coefficient 0.94, maximum depth 0.15 m |
| [San Francisquito Creek](Examples/Stanford/README.md), CA (30 m grid) | `run('Examples/Stanford/run_excel.m')` | 100 mm/h of rain for 1 h on a 157 km² catchment, 180 min simulated | ~14 min | rain 99.97 mm, peak outflow across the domain edge 402.8 m³/s, runoff coefficient 0.16, maximum depth 5.87 m |
| [One check per process](Validation/README.md) | each case's run script | about 80 small synthetic cases, each with a pass or fail rule | 1-10 min each | [Validation_Cases.csv](Validation/Validation_Cases.csv) |

Each example writes its results to its own `Outputs/` folder. A 24 h version of the Stanford case,
with maps and a video, is on the
[website](https://marcusnobrega-eng.github.io/HydroPol2D-docs/manual/Examples/RainfallRunoff).

## Gallery

Results from studies. Their inputs are not in this repository.

<table>
  <tr>
    <td align="center" width="50%">
      <img src="https://marcusnobrega-eng.github.io/profile//files/Rain_on_the_grid.gif" width="300" alt="Animated flood simulation of an urban catchment" /><br />
      <sub>A 50-year storm over an urban catchment with storm drains, São Paulo, Brazil</sub>
    </td>
    <td align="center" width="50%">
      <img src="https://github.com/user-attachments/assets/adfcdfe4-66e7-42aa-9335-de0e6d56a729" width="300" alt="Dam-break flood over a city" /><br />
      <sub>A full dam break over a city in Pernambuco, Brazil</sub>
    </td>
  </tr>
</table>

**A rainfall-runoff event in a catchment near Palo Alto, CA**

https://github.com/user-attachments/assets/739e39eb-18e8-4ec3-bee6-a7a647793774

**A 22-day simulation with rainfall that varies in space (MSWEP), Pune, India**

https://github.com/user-attachments/assets/3ad69bb3-3e24-4df3-b5ff-584aafcc6296

## The HydroPol2D family

| Product | What it is | License |
|---|---|---|
| **HydroPol2D (MATLAB)** | This repository. The model's source code in MATLAB. You prepare a case in an Excel workbook and run it in MATLAB. | Free and open source ([GPL-3.0](LICENSE)) |
| **[HydroPol2D Studio](https://marcusnobrega-eng.github.io/HydroPol2D-docs/manual/studio/overview)** | A desktop app for Mac and Windows. You build a model on a map, run it, calibrate it and view the results in one window, with no MATLAB and no coding. It runs the Python version of the model. | Proprietary, needs a license code |
| **[HydroPol2D Case Builder](https://marcusnobrega-eng.github.io/HydroPol2D-docs/manual/builder/overview)** | A tool that prepares a ready-to-run case for a city or a catchment from public global data: the model area, terrain, rivers, land use, soil, rain and weather. It writes cases for this MATLAB version and for the Python version. | Proprietary, needs a license code |

<p align="center">
  <img src="https://marcusnobrega-eng.github.io/HydroPol2D-docs/img/getting-started/data-to-results.svg" width="720" alt="Maps, parameter tables, forcing and settings make one case, which HydroPol2D Studio or HydroPol2D (MATLAB) runs into maps, time series, animations and a results package" />
</p>

## Licenses

- **HydroPol2D (MATLAB)** is free under the [GNU General Public License v3.0](LICENSE). You may use,
  study, change and share it. Software that you distribute and that is built on it must stay open
  source under the same license. The bundled TopoToolbox code keeps its own GPL-3.0 license
  ([`third_party/NOTICE.md`](third_party/NOTICE.md)).
- **HydroPol2D Studio** and **HydroPol2D Case Builder** are proprietary. To get a license, write to
  Marcus N. Gomes Jr. at [marcusnobrega.engcivil@gmail.com](mailto:marcusnobrega.engcivil@gmail.com).
  Studio downloads are on its
  [release page](https://github.com/marcusnobrega-eng/HydroPol2D-Studio-releases/releases).

## Repository layout

| Path | Contents |
|---|---|
| [`run_HydroPol2D.m`](run_HydroPol2D.m) | Main script. It reads the case workbook and runs the model. `HydroPol2D_V115.m`, its former name, still works. |
| [`HydroPol2D_Functions/`](HydroPol2D_Functions) | Model code |
| [`Input_Data_Sheets/`](Input_Data_Sheets) | Case workbook template |
| [`Config/`](Config) | Script settings for automated studies |
| [`Examples/`](Examples) | Examples you can run |
| [`Validation/`](Validation/README.md) | Process checks with pass or fail rules |
| [`Applications/`](Applications/README.md) | Cluster (SLURM) template and dam-case runners |
| [`tools/`](tools) | Workbook builder, which converts older `General_Data.xlsx` cases, and other tools |
| [`third_party/`](third_party/NOTICE.md) | Bundled TopoToolbox terrain code |

## Documentation

The [HydroPol2D website](https://marcusnobrega-eng.github.io/HydroPol2D-docs/), in English,
Portuguese and Spanish, has the full manual:
[getting started](https://marcusnobrega-eng.github.io/HydroPol2D-docs/manual/getting-started/overview),
[inputs](https://marcusnobrega-eng.github.io/HydroPol2D-docs/manual/inputs/overview),
[case workbook](https://marcusnobrega-eng.github.io/HydroPol2D-docs/manual/configuration/general_data),
[theory](https://marcusnobrega-eng.github.io/HydroPol2D-docs/manual/theory/overview),
[validation](https://marcusnobrega-eng.github.io/HydroPol2D-docs/manual/validation/overview),
[Studio](https://marcusnobrega-eng.github.io/HydroPol2D-docs/manual/studio/overview),
[Case Builder](https://marcusnobrega-eng.github.io/HydroPol2D-docs/manual/builder/overview) and the
[changelog](https://marcusnobrega-eng.github.io/HydroPol2D-docs/manual/validation/changelog).

## Versions

This repository is at **v1.19.0-rc2** ([release notes](RELEASE_NOTES_v1.19.0-rc2.md)), a release
candidate. It matches HydroPol2D Studio 0.9.0rc2. Version `v1.17.0` stays the stable reference
until the last hardware and packaging checks pass.

## Data sources of the examples

- **San Francisquito Creek, terrain:** USGS 3D Elevation Program (3DEP), 10 m averaged to 30 m
  (`Examples/Stanford/Static/DEM.tif`) and 10 m (`Examples/Stanford/Static_10m_LiDAR`), U.S. public domain.
- **San Francisquito Creek, land use:** ESA WorldCover 2020 (v100), CC BY 4.0. © ESA WorldCover
  project 2020 / Contains modified Copernicus Sentinel data (2020) processed by ESA WorldCover
  consortium. The class codes 10 to 90 and their names in the land-use table are WorldCover's.
- **San Francisquito Creek, other maps:** single values set for the example (soil class 8, depth
  to bedrock 1.5 m, LAI 1, albedo 0.25, initial soil moisture 0.1104).
- **V-tilted catchment:** all inputs are generated by
  [`prepare_standard_case.m`](Examples/VTilted_Excel/prepare_standard_case.m).

## Cite

If you use HydroPol2D, please cite:

- Gomes Jr, M.N., et al., 2023. [HydroPol2D: Distributed hydrodynamic and water quality model: Challenges and opportunities in poorly-gauged catchments](https://doi.org/10.1016/j.jhydrol.2023.129982). *Journal of Hydrology*, 625, 129982.
- Gomes Jr, M.N., Castro, M.A., Castillo, L.M., Sánchez, M.H., Giacomoni, M.H., de Paiva, R.C. and Bates, P.D., 2025. [Spatio-temporal performance of 2D local inertial hydrodynamic models for urban drainage and dam-break applications](https://doi.org/10.1016/j.jhydrol.2025.134661). *Journal of Hydrology*.

<details>
<summary><b>More publications using HydroPol2D</b></summary>

- Gomes Jr, M.N., Jalihal, V., Castro, M. and Mendiondo, E.M., 2025. [Exploring the impact of rainfall temporal distribution and critical durations on flood hazard modeling](https://doi.org/10.1007/s11069-025-07186-3). *Natural Hazards*.
- Gomes Jr, M.N., et al., 2024. [Global optimization-based calibration algorithm for a 2D distributed hydrologic-hydrodynamic and water quality model](https://doi.org/10.1016/j.envsoft.2024.106128). *Environmental Modelling & Software*, 179, 106128.
- Rápalo, L.M., Gomes Jr, M.N. and Mendiondo, E.M., 2024. [Developing an open-source flood forecasting system adapted to data-scarce regions: A digital twin coupled with hydrologic-hydrodynamic simulations](https://doi.org/10.1016/j.jhydrol.2024.131929). *Journal of Hydrology*, 644, 131929.
- Rápalo, L.M., Gomes Jr, M.N. and Mendiondo, E.M., 2025. [Multiple levels of human instability due to urban overland flow within the 21st century: An urban catchment study case in Brazil](https://doi.org/10.1016/j.ijdrr.2025.105931). *International Journal of Disaster Risk Reduction*.
- Castro, M., Gomes Jr, M., Rápalo, L. and Mendiondo, E., 2025. [Performance of low-complexity hydrodynamic models for dam-break flood mapping: Trade-offs between full momentum, diffusive-wave, and local-inertial models](https://doi.org/10.22541/essoar.175242171.16927854/v1) (under review).
- Castro, M.D.A.R.A., Gomes Jr, M.N. and Mendiondo, E.M., 2025. *Probabilistic dam-break flood mapping via Monte Carlo simulations using a 2D local-inertial model* (DOI not available).
- Sanchez, M.H., Gomes Jr, M.N., Rápalo, L., Castro, M.D.A.R.A. and Mendiondo, E.M., 2026. [Assessment of urban vulnerability to floods under current and future demographic and climate scenarios](https://doi.org/10.2139/ssrn.6381670). SSRN (under review).
- Sousa, M.R.D., Mendiondo, E.M. and Gomes, M.N., 2026. [Hydrological-hydrodynamic modeling of climate-induced urban flooding of design storms using HydroPol2D: a case study in São Carlos, Brazil](https://www.scielo.br/j/rbrh/a/KFQjyK4vxHPscnyCCRrBGjH/?format=pdf&lang=en). *RBRH*.
- Sánchez, M.H., 2025. [A novel framework for flood vulnerability assessment and fuzzy-controlled reservoir optimization in data-scarce urban watersheds](https://teses.usp.br/teses/disponiveis/18/18138/tde-30072025-105642/publico/CORRIGIDO_Mateo_Hernandez_Sanchez.pdf). Master's thesis, University of São Paulo.
- Castro, M.D.A.R.A., 2024. [Métodos determinísticos e probabilísticos para a avaliação do impacto de rompimentos de barragens via modelagem hidrodinâmica](https://teses.usp.br/teses/disponiveis/18/18138/tde-21032025-085649/publico/Dissertacao_versao_corrigida.pdf). Master's thesis, University of São Paulo.
- Rápalo, L.M.C., 2024. [Open-source tools for flood risk assessment for multiple spatio-temporal scales under scenarios of change](https://teses.usp.br/teses/disponiveis/18/18138/tde-07102024-094249/publico/ThesisCastilloRapaloLuisMiguelCorrected.pdf). Doctoral dissertation, University of São Paulo.
- Gomes Jr, M.N., 2023. [Advances in open source hydroinformatics for flood modeling and disaster education](https://teses.usp.br/teses/disponiveis/18/18138/tde-18042024-142515/publico/TeseMarcusNobregaGomesJuniorVersaoCorrigidaCompressed.pdf). Doctoral dissertation, University of São Paulo.

</details>

## Contact

**Dr. Marcus N. Gomes Jr.**, [mngomes@stanford.edu](mailto:mngomes@stanford.edu). Please report bugs
and ask questions as [GitHub issues](https://github.com/marcusnobrega-eng/HydroPol2D/issues).
For licenses and for a short course (2-hour classes with slides), write to
[marcusnobrega.engcivil@gmail.com](mailto:marcusnobrega.engcivil@gmail.com).

**Acknowledgments:** University of São Paulo, São Carlos School of Engineering. University of
Texas at San Antonio, Civil and Environmental Engineering. Stanford University, Stanford Doerr
School of Sustainability, Department of Earth System Science.
