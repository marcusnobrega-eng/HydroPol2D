# VAL-GW-LINRES-001: Linear Reservoir Recession

Purpose: verify groundwater/baseflow recession dynamics against an exponential linear-reservoir solution.

Truth source: `S(t) = S0 exp(-k t)` and `Q(t) = k S(t)`.

This is also a simple inverse-modeling case: generate synthetic discharge from known `k`, then recover `k` from the log-linear recession.

Status: planned formula/dynamics verification.
