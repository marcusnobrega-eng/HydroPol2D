# SLURM Batch Template

`HP2D_HPC_Config.sbatch` submits the standard HydroPol2D launcher to a SLURM cluster. `HP2D_Wrapper.m` locates the repository root, registers the bundled runtime, disables figure display, and runs `HydroPol2D_V115.m`.

Submit the template from the repository root:

```bash
sbatch Applications/HPC/HP2D_HPC_Config.sbatch
```

Before submitting, set the partition, time limit, CPU, memory, GPU, and module commands for the target cluster. The supplied template requests one GPU and loads MATLAB; it should be adapted to the local scheduler environment.

The model configuration remains in `Config/` or the selected input spreadsheet. This folder only provides the scheduler entry point.
