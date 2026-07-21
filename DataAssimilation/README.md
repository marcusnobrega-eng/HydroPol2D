# HydroPol2D Data Assimilation

This folder contains optional data-assimilation workflows. The production
HydroPol2D model source remains in `HydroPol2D_Functions`.

Current runnable prototype:

```matlab
addpath('/Users/mngomes/Documents/HydroPol2D_CodexUpdate/HydroPol2D_Model/DataAssimilation')
run_data_assimilation_example('SmokeTest', true, 'NParticles', 1, 'NWindows', 1)
```

Generated DA files stay inside each prototype folder under `Runs`, `Outputs`,
and `Figures`.
