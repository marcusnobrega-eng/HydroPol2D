# Validation Scripts

Run the Phase 1 suite from the repository root:

```bash
python3 Validation/scripts/run_phase1_release.py \
  --matlab /Applications/MATLAB_R2025b.app/bin/matlab
```

`run_phase1_release.py` starts each driver from a clean MATLAB path and
writes `Validation/Release/Outputs/Phase1_Release_Summary.csv`. It fails
when a report-ready result fails; documented diagnostic rows remain in the
summary without blocking the run.

After a completed run, check that every registry row has its source folder,
driver, metadata, and pass/fail output:

```bash
python3 Validation/scripts/audit_validation_registry.py
```

`phase1_reference_solutions.py` produces independent reference series for
the analytical Phase 1 cases. Its outputs are written under
`Validation/Reference_Outputs/Phase1`.
