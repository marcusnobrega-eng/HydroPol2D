# Validation Scripts

`audit_validation_registry.py` checks that each registered validation case has a folder, acceptance threshold, metrics, and report references. Use it before compiling the validation report and before marking a case report-ready.

Example:

```bash
python3 HydroPol2D_Model/Validation/scripts/audit_validation_registry.py HydroPol2D_Model/Validation/Validation_Cases.csv
```

Use `--strict` when the campaign is close to submission and warnings should fail the audit.

`phase1_reference_solutions.py` generates analytical/reference CSVs for the Phase 1 formula and dynamics verification cases. These outputs live under `HydroPol2D_Model/Validation/Reference_Outputs/Phase1` and should be compared against HydroPol2D/module outputs.
