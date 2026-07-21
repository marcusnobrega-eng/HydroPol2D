# Constant-flow composite-channel test

This test compares three local-inertial simulations of the same straight
compound channel:

- a 10 m grid that resolves the 10 m-wide channel explicitly;
- an ordinary 30 m grid obtained by averaging the 10 m terrain;
- a 30 m Neal grid with a 10 m-wide, 1 m-deep embedded channel.

The reach is represented by 990 m so that both grids share identical
downstream boundaries. The cross section contains a 10 m channel, 90 m
floodplain benches on each side, and 10 m outer containment berms. The
longitudinal slope is 0.005 m/m. Manning n is 0.035 for the channel and
floodplain. A constant 100 m3/s inflow is applied for 120 minutes from a dry
initial condition.

A second case uses a Nash-type hydrograph with shape parameter 4, a peak
of 100 m3/s at 30 minutes, and a 180-minute simulation. Run it with:

`Q(t) = Qp (t/tp)^(n-1) exp[(n-1)(1-t/tp)]`

where `Qp = 100 m3/s`, `tp = 30 min`, and `n = 4`.

```matlab
run_neal_composite_channel_nash100
```

Transient outputs are written under `Outputs/Validation/NashTransient`,
and transient figures are written under `Figures/NashTransient`.

Run in MATLAB with:

```matlab
run_neal_composite_channel_100cms
```

Outputs are written under `Outputs/Validation`, and figures are written to
`Figures`.
