# Exact-fast validated snapshot

This directory is an isolated optimization laboratory for BSMPT v3.1.8.  The
repository root remains the unmodified upstream implementation.  The validated
implementation is stored as a small source overlay under `upstream/`, together
with wrappers, test inputs, exact comparison outputs, and the optimization
report.

The committed overlay consists of these files:

- `include/BSMPT/ThermalFunctions/thermalcoefficientcalculator.h`
- `include/BSMPT/bounce_solution/action_calculation.h`
- `include/BSMPT/bounce_solution/calcgw_profiler.h`
- `include/BSMPT/models/ClassPotentialOrigin.h`
- `include/BSMPT/utility/NumericalDerivatives.h`
- `include/BSMPT/utility/spline/spline.h`
- `src/ThermalFunctions/ThermalFunctions.cpp`
- `src/ThermalFunctions/thermalcoefficientcalculator.cpp`
- `src/bounce_solution/CMakeLists.txt`
- `src/bounce_solution/action_calculation.cpp`
- `src/bounce_solution/bounce_solution.cpp`
- `src/bounce_solution/calcgw_profiler.cpp`
- `src/models/ClassPotentialOrigin.cpp`
- `src/models/ClassPotentialR2HDM.cpp`
- `src/transition_tracer/transition_tracer.cpp`
- `src/utility/CMakeLists.txt`
- `src/utility/NumericalDerivatives.cpp`
- `src/utility/spline/spline.cpp`

To reconstruct the lab source tree in a fresh checkout, copy the repository
root source tree into a new directory under `bsmpt_speed_lab/` while excluding
`bsmpt_speed_lab` itself, then overwrite it with the committed files from the
`upstream/` overlay.  Keep the reconstruction and every build directory inside
`bsmpt_speed_lab/`; do not modify the repository-root source tree.

The production entry point is `run_calcgw_exact_fast.sh`.  Correctness was
checked against `run_calcgw_dense_control.sh` with `compare_outputs.py --rtol 0
--atol 0`; only the `runtime` column is ignored.  See
`OPTIMIZATION_REPORT_ZH.md` for the validation matrix and rejected candidates.
