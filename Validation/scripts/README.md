# Validation Scripts

Run validation suite from the repository root:

```bash
python3 Validation/scripts/run_validation_release.py \
  --matlab /Applications/MATLAB_R2025b.app/bin/matlab
```

`run_validation_release.py` starts each driver from a clean MATLAB path and
writes `Validation/Release/Outputs/Validation_Release_Summary.csv`. It fails
when a report-ready result fails; documented diagnostic rows remain in the
summary without blocking the run.

After a completed run, check that every registry row has its source folder,
driver, metadata, and pass/fail output:

```bash
python3 Validation/scripts/audit_validation_registry.py
```

`reference_solutions.py` produces independent reference series for
the analytical validation cases. Its outputs are written under
`Validation/Reference_Outputs/Validation`.
