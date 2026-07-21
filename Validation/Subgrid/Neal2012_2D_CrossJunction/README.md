# Two-Dimensional Cross-Junction Benchmark

This case extends the steady straight-channel test to an orthogonal river
junction. Constant flows enter from the west, north, and south. The three
branches meet at the center of the domain and discharge through the eastern
boundary.

## Configuration

- Domain: 990 m by 990 m
- Fine reference resolution: 10 m
- Coarse resolution: 30 m
- Channel width: 10 m
- Bankfull depth: 1 m
- Longitudinal slope: 0.005 m/m
- Lateral hillslope: 0.005 m/m toward the channel
- Outlet: normal flow over the complete eastern boundary
- Outlet slope: 0.005 m/m
- Channel and floodplain Manning coefficient: 0.035
- Inflow: 100 m3/s at each of three inlets
- Total steady inflow: 300 m3/s
- Simulation duration: 180 min

The driver compares a fine explicit local-inertial model, an ordinary coarse
local-inertial model, and a coarse Neal (2012) channel-subgrid model. All
three simulations receive the same water volume. Every cell on the eastern
boundary uses the normal-flow outlet condition.

## Run

```matlab
cd HydroPol2D_Model/Validation/Subgrid/Neal2012_2D_CrossJunction
run_neal_2d_cross_junction_validation
run_neal_2d_cross_junction_validation("inbank")
```

Outputs include inlet-branch, downstream, and outlet hydrographs; storage and
wet-area time series; maximum-depth maps; event snapshots; and mass-balance
and pass/fail tables.

## Overbank results

The Neal (2012) simulation reproduces the fine-grid transient and steady
responses more closely than the ordinary 30 m model. Outlet NSE increases
from 0.937 to 0.997, and outlet RMSE decreases from 17.06 to 3.40 m3/s. NSE
values at the west, north, and south branch gauges are 0.996, 0.996, and
0.997, respectively. Main-channel and downstream-gauge NSE values are 0.997
and 0.997.

Maximum-depth RMSE decreases from 0.123 m for the ordinary 30 m model to
0.029 m for the Neal model. Wet-area CSI at the 0.10 m threshold increases
from 0.677 to 0.949, and the Neal wet-area error is -0.31%. Storage RMSE
decreases from 30,652 to 6,649 m3.

During the final 20 minutes, the Neal outlet discharge averages 300.00 m3/s,
equal to the imposed total flow to the reported precision. The largest
branch-flow error is 0.209%. Water balance closes to numerical precision in
all three simulations.

## In-bank results

The 50 m3/s test uses the same planform, slopes, roughness, and boundary
conditions, with a 3 m bankfull depth. The Manning uniform-flow capacity is
92.16 m3/s, so the imposed flow is 54.3% of bankfull capacity. Each inlet
supplies 16.67 m3/s.

The fine, ordinary coarse, and Neal simulations remain confined to the
channel, with no water deeper than 0.01 m on the floodplain. The Neal outlet
NSE is 0.976 and its outlet RMSE is 1.61 m3/s, compared with 0.643 and
6.17 m3/s for the ordinary 30 m model. Maximum-depth RMSE decreases from
0.091 to 0.023 m, and storage RMSE decreases from 10,908 to 3,238 m3. The
steady outlet flow is 50.00 m3/s and water balance closes to numerical
precision in every simulation.

The channel centerline is rasterized separately at 10 and 30 m so that it
reaches the eastern boundary cell in both grids. This prevents an artificial
bank at the fine-grid outlet and ensures that the entire eastern boundary
applies the same normal-flow condition.
