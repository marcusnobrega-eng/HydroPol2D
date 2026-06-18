# ROUTE-001: Kinematic/Diffusive Routing Audit

Purpose: audit and correct the kinematic and diffusive routing implementations before validation claims are made.

Audit checklist:
- Confirm governing equations and discretization for kinematic routing.
- Confirm pressure-gradient/diffusive terms and numerical stability limits for diffusive routing.
- Verify sign conventions, slope handling, boundary fluxes, and wet/dry behavior.
- Compare kinematic, diffusive, and CA/D8 routing on the same controlled synthetic planes only after the audit passes.

Required diagnostics after correction:
- implementation audit notes
- `Outputs/Validation/Mass_Balance.csv`
- `Outputs/Validation/Metric_Summary.csv`
- `Outputs/Validation/Pass_Fail.csv`
- hydrograph comparison figures for all routing alternatives

Acceptance threshold: no validation claim until the audit passes; after correction, domain water-balance error `< 0.1%` and hydrograph metrics are reported.
