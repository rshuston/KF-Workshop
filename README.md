# KF Workshop

A MATLAB workspace to explore various Kalman filtering techniques.

The following filter types are explored:
- Kalman Filter, Transformed Range and Bearing Measurements
- Extended Kalman Filter, Range and Bearing Measurements
- Extended Kalman Filter, Range and Direction Cosine Measurements
- Unscented Kalman Filter (measurement relation only), Range and Bearing Measurements
- Unscented Kalman Filter (measurement relation ony), Range and Direction Cosine Measurements
- Unscented Kalman Filter (classic Julier and Uhlmann algorithm), Range and Direction Cosine Measurements

Two model types are explored:
- Constant velocity with acceleration disturbances
- Constant turn with acceleration disturbances

Six flight scenarios are simulated:
1. circle
2. curve
3. line
4. s-curve
5. square
6. wiggle

Run the "setup_env" script in Octave to add the necessary subdirectories to
the function search path.

To generate plots, run "fly scenario" in Octave, where "scenario" is one of
the above scenario names. If just "fly" is specified, the "line" scenario is
run.

An "engineer's derivation" of the Kalman filter and its various forms can be
found in the `doc` directory. Also, a writeup of the various target motion
models is given.
