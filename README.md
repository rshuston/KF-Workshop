# KF Workshop

A workspace to explore various Kalman filtering techniques.

The following filter types are explored:
- Kalman Filter, Transformed Range and Bearing Measurements
- Extended Kalman Filter, Range and Bearing Measurements
- Extended Kalman Filter, Range and Direction Cosine Measurements
- Unscented Kalman Filter (measurement relation only), Range and Bearing Measurements
- Unscented Kalman Filter (measurement relation ony), Range and Direction Cosine Measurements
- Unscented Kalman Filter (classic Julier and Uhlmann algorithm), Range and Direction Cosine Measurements

The filters are implemented in Octave, although things should port easily to
MATLAB.

Five flight scenarios are simulated:
1. line
2. curve
3. circle
4. wiggle
5. square

To generate plots, run "fly scenario" in Octave, where "scenario" is one of
the above scenario names. If just "fly" is specified, the "line" scenario is
run.

An "engineer's derivation" of the Kalman filter and its various forms can be
found in the `doc` directory.
