# KF Workshop

A workspace to explore various Kalman Filtering techniques.

The following filter types are explored:
- Kalman Filter, Transformed Range and Bearing Measurements
- Extended Kalman Filter, Range and Bearing Measurements
- Extended Kalman Filter, Range and Direction Cosine Measurements
- Unscented Kalman Filter, Range and Bearing Measurements
- Unscented Kalman Filter, Range and Direction Cosine Measurements

The filters are implemented in Octave, although things should port
easily to MATLAB.

Five flight scenarios are simulated:
1. line
2. curve
3. circle
4. wiggle
5. square

To generate plots, run "fly scenario" in Octave, where "scenario"
is one of the above scenario names. If just "fly" is specified,
the "line" scenario is run.
