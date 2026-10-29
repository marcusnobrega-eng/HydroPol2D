# ET after infiltration (VAL-ET-ORDER-001)

Purpose: check the order of the legacy vertical step in `Hydrological_Model.m`.
Snow, canopy, infiltration, open-water evaporation, then soil ET, as in
HydroPol2D-Python and the MATLAB layered path. A cell counts as wet only if
water is left on the surface after infiltration.

One step of 15 min on a 3 x 3 grid with 2 mm of rain, ETP 4 mm/day, Kc 0.8 and
Ep 5 mm/day:

- Row 1, the soil takes all the rain: the cell ends dry, has no open-water
  evaporation and loses soil ET = ETP x Kc.
- Row 2, the soil is full: the rain ponds and loses Ep; no soil ET.
- Row 3, no rain: soil ET = ETP x Kc.
- Ledger: rain = change in soil + change in surface + evaporation + soil ET.

Before the fix, ET ran before infiltration. The rained-on cell in row 1 then
counted as open water: it lost Ep from the rain and got no soil ET.

## Run

```matlab
run('Validation/Evapotranspiration/StepOrder_RainyCell/run_et_infiltration_order_validation.m')
```

## Outputs

- `Outputs/Validation/ET_Order_Cells.csv`
- `Outputs/Validation/ET_Order_Checks.csv`
- `Outputs/Validation/ET_Order_Ledger.csv`
- `Outputs/Validation/Pass_Fail.csv`

Acceptance: every cell value within 1e-12 of the expected value; ledger
residual below 1e-10 m3.
