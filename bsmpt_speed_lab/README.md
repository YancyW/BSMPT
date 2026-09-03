# BSMPT speed laboratory

> 当前完整中文总结、最终推荐配置、全部接受/拒绝实验和复现方法，请以
> [`OPTIMIZATION_REPORT_ZH.md`](OPTIMIZATION_REPORT_ZH.md) 为准。本文后续内容
> 保留了开发过程记录，其中部分早期 benchmark 已被更新结果取代。

This directory is isolated from the existing BSMPT installation. Nothing here
is installed, linked, or copied over `install/bin/CalcGW`.

## CalcGW code-level optimization result

Phase instrumentation on the `mH < 125`, SNR 11.55 reference point split the
original 201.9 s internal runtime into approximately 20.2 s of vacuum/phase
tracing, 167.1 s of bounce construction, 13.6 s of percolation/completion, and
0.4 s of GW-parameter evaluation. The bounce solver is therefore the dominant
target.

The isolated source contains three default optimizations:

1. `BounceActionInt::SetPath` no longer builds a 1001-point `dV/dl` grid which
   `Solve1DBounce` always overwrote before its first use.
2. Numerical gradient/Hessian routines reuse temporary vectors, cache the
   repeated `V(phi)` value, reserve fixed grid storage, and avoid repeated
   spline evaluation in `Calc_d2Vdl2`.
3. `BackwardsPropagation` now initializes `Initial_lmin` before the exact
   threshold scan uses it. The true-vacuum Hessian, which is constant for that
   path, is evaluated once rather than about 1001 times inside the scan. This
   removes undefined-value sensitivity and tens of thousands of repeated
   `VEff`/mass-eigensolver calls per action.

The third change intentionally fixes undefined behavior in the upstream order
of operations, so its stable result need not reproduce the accidental value
of an older binary bit-for-bit. Repeated runs of the corrected build agreed in
all non-runtime fields bit-for-bit in both mass regions:

| Reference point | Earlier optimized runtime | Corrected/cached runtime | Repeated corrected output |
|---|---:|---:|---|
| `mH < 125`, SNR about 11.53 | 181.56 s | 145.47-146.17 s | bit-for-bit stable |
| `mH > 125`, very small SNR | 77.62 s | 55.93-60.71 s | bit-for-bit stable |

Further exact common-subexpression optimizations now:

- reuse the numerical gradient when `dV/dl` and `d2V/dl2` are requested at the
  same threshold point;
- reuse adjacent Simpson endpoints in the final potential-action integral;
- optionally skip strictly zero R2HDM field components while assembling the
  Higgs mass matrix (`BSMPT_SKIP_ZERO_HIGGS_FIELDS=1`).

Two additional matrix experiments were validated on both sides of the
`mH=125` split:

- `BSMPT_USE_R2HDM_HIGGS_PAIR=1` builds the common field-dependent R2HDM
  Higgs matrix once for the Arnold--Espinosa `M(v,T)`/`M(v,0)` pair, while
  retaining two independent self-adjoint eigensolver calls.  Every non-runtime
  output field was bit-for-bit identical on the tested `mH<125` high-SNR and
  `mH>125` ultraweak points.  On the latter, a same-binary A/B reduced internal
  runtime from 101.726 s to 96.209 s (5.4%) under the then-current machine
  load.
- `BSMPT_USE_ANALYTIC_GAUGE_MASSES=1` uses the exact R2HDM 4x4 arrowhead
  spectrum for `diff==0`, after checking the analytic matrix against the full
  assembled matrix and falling back on any mismatch.  It was also bit-for-bit
  identical in both mass regions, but its measured benefit was only about
  2.3% on the `mH>125` point and noise-level on `mH<125`, so it remains
  experimental.

The two switches together were bit-for-bit identical in both regions and gave
about 6.0% on the `mH>125` A/B (101.726 s to 95.624 s).  They remain opt-in
until a broader classified-point validation is complete.

The subsequently added `BSMPT_USE_PATH_GEOMETRY_JET=1` evaluates the path
geometry (`phi`, `dphi/dl`, and `d2phi/dl2`) from one shared spline inversion.
Together with the Higgs-pair path it reduced a same-load `mH>125` run from
61.765 s to 53.072 s (14.1%) and an `mH<125` run from 157.066 s to 142.068 s
(9.5%), with bit-for-bit identical physical output.

`BSMPT_SKIP_ZERO_R2HDM_HIGGS_TERMS=1` separately skips exactly-zero R2HDM
L3/L4 tensor contributions without changing the order of nonzero additions.
After correcting an experimental control-flow bug (a zero L3 term must not
skip a potentially nonzero L4 term), it is bit-for-bit identical in both mass
regions.  It gave a further 7.2% on the tested `mH>125` point and noise-level
change on the `mH<125` point.

The current recommended exact validation combination is therefore:

```bash
BSMPT_USE_PATH_GEOMETRY_JET=1 \
BSMPT_USE_R2HDM_HIGGS_PAIR=1 \
BSMPT_SKIP_ZERO_R2HDM_HIGGS_TERMS=1
```

`BSMPT_USE_R2HDM_HIGGS_PAIR_DIFF=1` extends the same local T/0 matrix
assembly reuse to the `diff=1..8` Arnold--Espinosa paths.  It retains two
independent eigenvalue-derivative solves and only shares the common mass and
derivative matrix construction.  It was bit-for-bit identical on both main
mass-region references and on all four additional classified points.  Its
incremental gain was about 1.1% on the ultraweak `mH>125` reference, 6.5% on
the `mH<125` high-SNR reference, and 2.3% summed over the four-point set
(204.338 s to 199.666 s).  Add this switch to the recommended exact
combination.

With this addition, the `mH<125` high-SNR reference is 132.934 s versus
203.37 s for the original installed binary (34.6% less time, 1.53x
throughput).  The `mH>125` reference is 48.691 s versus 88.66 s originally
(45.1% less time, 1.82x throughput).  Machine load varies, so paired internal
runtime comparisons remain the acceptance criterion.

## Current exact opt-in set

Further strict implementations replace repeated zero tests with per-model
index lists, pair the R2HDM gauge matrix assembly, and retain the previous
spline interval during continuous queries.  The current fully validated set is:

```bash
BSMPT_USE_PATH_GEOMETRY_JET=1 \
BSMPT_USE_SPLINE_INTERVAL_HINT=1 \
BSMPT_USE_COMBINED_NUMERICAL_DERIVATIVES=1 \
BSMPT_USE_R2HDM_HIGGS_PAIR=1 \
BSMPT_USE_R2HDM_HIGGS_PAIR_DIFF=1 \
BSMPT_USE_R2HDM_HIGGS_INDEX_CACHE=1 \
BSMPT_USE_R2HDM_GAUGE_PAIR=1 \
BSMPT_USE_R2HDM_GAUGE_INDEX_CACHE=1 \
BSMPT_USE_R2HDM_QUARK_INDEX_CACHE=1 \
BSMPT_USE_R2HDM_LEPTON_INDEX_CACHE=1 \
BSMPT_USE_R2HDM_QUARK_FIXED12_DIFF0=1 \
BSMPT_USE_R2HDM_LEPTON_FIXED9_DIFF0=1 \
BSMPT_USE_V1LOOP_THERMAL_CONTEXT=1 \
BSMPT_USE_R2HDM_VTREE_INDEX_CACHE=1 \
BSMPT_USE_R2HDM_COUNTERTERM_INDEX_CACHE=1
```

The same validated set is wrapped by `run_calcgw_exact_fast.sh`; normal
CalcGW arguments pass through unchanged:

```bash
./bsmpt_speed_lab/run_calcgw_exact_fast.sh \
  --model=r2hdm --input=points.tsv --output=gw.tsv \
  --firstline=2 --lastline=100
```

The Higgs and gauge index caches store only ordered nonzero tensor indices in
each model instance; they never share parameter values between points.  The
gauge cache has both the `(i,j)` mass-contraction order and a per-derivative-
field `j` list, so it covers `diff==0` and `diff>0` without changing the
contraction order.  Gauge pairing shares matrix construction but retains
independent T and zero-temperature eigensolver calls.  All switches
automatically fall back outside their R2HDM and Arnold--Espinosa scopes.

This set remains bit-for-bit identical on both main references and all four
additional classified points.  The four-point corrected-baseline total is now
248.635 s to 126.270 s (49.2% less time).  The main references are:

| Region | Original installed | Current exact | Time reduction | Throughput |
|---|---:|---:|---:|---:|
| `mH < 125` | 203.37 s | 95.490 s | 53.0% | 2.13x |
| `mH > 125` | 88.66 s | 30.535 s | 65.6% | 2.90x |

The latest set includes sparse ordered caches for the tree potential,
counterterms, quark matrix, and lepton matrix.  The combined numerical
derivative path preserves the original four-point gradient and full Hessian
stencils but reuses the identical `phi +/- 2 eps e_i` potential samples.  No
finite-difference order, step, or threshold grid is reduced.

The fixed-size R2HDM fermion paths now directly assemble the complete 12x12
quark and 9x9 lepton matrices and leave their dynamic fallback matrices
unallocated unless that fallback is actually selected.  This removes dead
heap allocation without changing matrix entries or eigensolver arithmetic.
They are now enabled by the generic wrapper after both main references, all
four classified validation points, and the independent `mH<125`
SNR-about-150 point remained bit-for-bit identical in every non-runtime field
(seven points total).  The four-point total changed from 126.270 s to
121.957 s in the validation run.

An active-direction Hessian experiment found all four R2HDM path directions
strictly nonzero on the `mH>125` reference (102488 active dimensions over
25622 combined evaluations), so it could not skip any matrix entries and is
not recommended.  A local mass-contribution hash cache and a fixed-storage
RK5 experiment were also slower and remain disabled.

The earlier `BSMPT_SKIP_ZERO_R2HDM_HIGGS_TERMS` branch remains useful as a
fallback experiment, but the ordered per-instance index cache supersedes it.

On four additional classified points not used while developing the changes
(two on each side of `mH=125`, SNR approximately `3e-27` through `2.6e-10`),
all non-runtime fields agreed bit-for-bit with the corrected exact baseline.
Their summed internal runtime fell from 248.635 s to 204.338 s, a 17.8%
reduction.  The input and both outputs are retained as
`classified_validation_4*.tsv`.

Two further exact-output experiments were not performance wins:
`BSMPT_SAFE_R2HDM_HIGGS_DERIVATIVE_UPPER_ONLY=1` increased the tested runtime
from 49.246 s to 52.121 s, while `BSMPT_SKIP_ZERO_R2HDM_GAUGE_TERMS=1`
increased it to 50.337 s.  They remain disabled.

The mathematically valid `BSMPT_USE_R2HDM_QUARK_BLOCK6=1` reduction was
rejected for production use: on the ultraweak `mH>125` reference it was slower
and shifted beta/H by about 0.31%, despite preserving the exact block spectrum
in real arithmetic.  This illustrates why matrix-size reductions are not
promoted without full CalcGW validation.

The first two remain bit-for-bit identical on both mass-region references. The
current exact `mH > 125` reference runs in about 51.3 s before the optional
zero-field matrix shortcut (49.7 s with it), compared with 88.7 s for the
original installed binary. The exact `mH < 125` high-SNR reference is about
137-140 s with the tested exact shortcuts, compared with 203.4 s originally.

`BSMPT_ADAPTIVE_THRESHOLD=1` is still experimental. It reduced threshold
samples by 68% and produced excellent agreement on the SNR~11.5 point, but on
an independent SNR~150 point it shifted beta/H by about 3% and SNR by about
1.7% relative to the corrected dense scan. It is therefore not an exact-mode
default. A dense 2000-step convergence run also showed that the historical
1000-step threshold has non-negligible discretization error of its own.

Measured end-to-end results with identical CalcGW options:

| Reference point | Control | Optimized | Speedup | Physical output |
|---|---:|---:|---:|---|
| `mH < 125`, SNR 11.55 | 203.37 s | 181.63 s | 1.120x | bit-for-bit identical |
| `mH > 125`, very small SNR | 88.66 s | 77.69 s | 1.141x | bit-for-bit identical |

The comparison ignores only the `runtime` column and uses zero numerical
tolerance for every other status and floating-point output field.

## Experimental fast-screen mode

Setting `BSMPT_USE_CENTRAL2_GRADIENT=1` changes the bounce solver's four-point
gradient stencil to a two-point central stencil. This is a screening mode, not
an exact replacement for the default calculation.

| Region / reference | Original runtime | Fast runtime | Dominant SNR relative error |
|---|---:|---:|---:|
| `mH < 125`, SNR 11.55 | 181.56 s | 129.35 s | -1.07e-5 |
| `mH < 125`, SNR 150.18 | 304.54 s | 206.78 s | -2.93e-5 |
| `mH > 125`, SNR 54.69 | 619.50 s | 323.66 s | -1.08e-3 |

All dominant-transition statuses agreed in these high-SNR checks. On the last
point, however, a secondary transition with SNR around `7e-31` changed its
nucleation status. Consequently, use the default exact mode whenever the goal
is to classify arbitrarily weak `SNR > 0` signals without boundary ambiguity.

The analytic-gradient, reduced-raster, directional-derivative, LTO, native
vectorization, and mixed-gradient experiments were rejected: they either
changed boundary statuses/important GW values, provided too little speedup, or
were slower/unstable. Their environment-controlled code paths remain disabled
by default for reproducible research only.

## Source provenance

The official upstream Git object already present in the parent repository was
exported from commit `04cb17d1233522f3c423cbd957a8922be037241e` (BSMPT 3.3.1).
Direct cloning was blocked by the configured GitHub proxy, so `upstream/` has
no `.git` directory and cannot modify the parent repository.

## Validated control build

`build-conda-control/bin/CalcGW` was compiled with the same Conda compiler and
the same non-vectorized options as the current installation. On
`benchmark_input.tsv`, all 98 physical/status output fields were bit-for-bit
identical to `install/bin/CalcGW`; only `runtime` differed.

The official `BSMPTUseVectorization=ON` option is not safe with the current
binary dependency set: the resulting executable aborts with heap corruption.
Do not use the `build-native*` experiments.

## Opt-in R2HDM Higgs T/zero-temperature matrix-pair experiment

The isolated source now contains a disabled-by-default experiment controlled by
`BSMPT_USE_R2HDM_HIGGS_PAIR=1`. It applies only when the model ID is R2HDM,
`diff == 0`, and the Arnold--Espinosa path is active. Inside one `V1Loop`, it
assembles the common zero-temperature Higgs field polynomial once, copies it,
adds the Debye term to the thermal copy, and still runs the original
`SelfAdjointEigenSolver<MatrixXd>` independently on both matrices. Other models,
derivatives, dimensions, and Parwani mode automatically use the original path.

This is an unvalidated opt-in source change awaiting a control-vs-experiment
build and exact TSV comparison. Do not enable it for production scans until the
parent agent has checked both `mH < 125` and `mH > 125` points with
`compare_outputs.py --rtol 0 --atol 0`.

`BSMPT_USE_R2HDM_HIGGS_PAIR_DIFF=1` is a separate, more restrictive experiment
for `diff=1..8` field derivatives. It reuses the R2HDM Higgs field-derivative
matrix assembly between the thermal and zero-temperature evaluations, while
retaining two independent `FirstDerivativeOfEigenvalues` calls. It is also
disabled unless R2HDM, Arnold--Espinosa, and the expected `8` Higgs / `4` VEV
dimensions are detected. This derivative mode has not yet passed end-to-end
TSV validation.

`BSMPT_USE_R2HDM_GAUGE_PAIR=1` is a further opt-in for the R2HDM gauge spectra
in `V1Loop`, covering `diff==0` and `diff=1..8` in Arnold--Espinosa mode. It
shares the field-dependent gauge matrix assembly between T and T=0 while
retaining separate eigenvalue/eigenvalue-derivative solves. It preserves the
`BSMPT_SAFE_GAUGE_UPPER_ONLY` triangle policy and automatically falls back for
`diff==-1`, Parwani, non-R2HDM, or unexpected dimensions. It also awaits exact
TSV validation.

`BSMPT_USE_R2HDM_GAUGE_INDEX_CACHE=1` enables an additional per-instance cache
inside that gauge-pair path. It skips exactly zero entries of
`Curvature_Gauge_G2H2` while traversing the original `i,j` order for the mass
matrix and the original `j` order for each `diff` field. The cache is rebuilt
after every `SetCurvatureArrays`/`initVectors` cycle, so changing model
parameters cannot leave stale tensor indices. It is read-only during `V1Loop`,
with no global or cross-point state. If the cache is unavailable, the
gauge-pair branch uses its original dense loops.

`BSMPT_USE_R2HDM_FIXED_EIGENSOLVER=1` is an additional `diff==0` experiment
for the paired R2HDM path. It copies the 8-by-8 Higgs and 4-by-4 gauge mass
matrices into fixed-size Eigen types and runs two independent
`SelfAdjointEigenSolver` instances for T and T=0. The dynamic-size solver is
used automatically when this switch is absent or outside the R2HDM,
Arnold--Espinosa, and expected-dimension guards. Because fixed-size Eigen may
take a different floating-point path, exact output and runtime must be
checked separately before combining it with the validated switches.

`BSMPT_USE_R2HDM_PAIR_STACK_MATRICES=1` is the recommended follow-up
allocation experiment. When both R2HDM Higgs and Gauge pair switches are
active, it assembles the 8-by-8 and 4-by-4 `diff==0` matrices in fixed-size
stack storage, then reuses one dynamic `MatrixXd` workspace for the two
unchanged dynamic eigensolver calls. It is disabled for `v.size()!=8` and all
non-paired, derivative, Parwani, non-R2HDM, or unexpected-dimension paths.
Unlike the fixed-eigensolver experiment, it does not change the eigensolver
algorithm; exact output should nevertheless be checked because the stack
assembly is a separate implementation path.

`BSMPT_USE_R2HDM_COUNTERTERM_INDEX_CACHE=1` enables an independent R2HDM
counterterm cache for `CounterTerm(diff==0)` and `CounterTerm(diff>0)`. It
stores only nonzero L1--L4 indices in an `i -> j -> k -> l` hierarchy, so the
nonzero additions retain the original nested order and coefficients. The
cache is built after `set_CT_Pot_Par` (when the counterterms are populated),
cleared by `initVectors`, and rebuilt for a new parameter point. It is
instance-local and read-only during evaluation; other models and cache-off
calls use the original dense loops.

`BSMPT_USE_V1LOOP_MASS_CONTRIBUTION_CACHE=1` enables a local, per-`V1Loop`
cache for repeated `boson`/`fermion` return values. Keys use the exact bit
representation of every floating-point argument (`MassSquared`, `Temp`, the
boson `cb` and derivative-mass argument) plus the derivative selector. Cache
hits still execute every original `res +=` statement independently, so
degenerate states are not combined and the accumulation order is unchanged.
The cache is discarded at the end of each call and is disabled by default.
The underlying helpers are read-only/pure for these arguments, so caching does
not suppress model state changes.

Internal minimizer multithreading reduced wall time from about 203.3 s to
196.6 s on the benchmark point, but introduced small numerical variation
(about 2.65e-5 relative in total SNR and 3.06e-4 in beta/H). It is therefore
not the default recommendation for boundary-sensitive scans.

## Exact point-level parallelism

For a TSV containing many independent R2HDM points:

```bash
python3 parallel_calcgw.py \
  --binary ../install/bin/CalcGW \
  --input points.tsv \
  --output results.tsv \
  --jobs 10
```

Each worker runs the unchanged serial CalcGW algorithm. Rows are merged back
in input order. Keep `jobs` below the physical-core count and below the memory
limit; start with 8-10 on the i7-13700K.

The two-point validation used two copies of the high-SNR reference point. A
serial estimate is about 406.6 s, while `--jobs 2` completed in 204.5 s
(1.99x throughput). Both output rows matched the serial reference in every
field except `runtime`, including an identical total SNR of
`11.551909703375555`.

Compare an output against a reference with:

```bash
python3 compare_outputs.py reference.tsv candidate.tsv
```

## Isolated GCC PGO experiment

`build-pgo-gcc` is a separate build tree and does not modify
`build-conda-control` or `install`.  It uses the existing Conan toolchain and
the Conda GCC compiler, with the following effective compile flags:

```text
-O3 -g -fprofile-generate -fno-fast-math -fno-tree-vectorize -fPIC
-fstack-protector-strong -fno-plt -ffunction-sections -pipe -m64
```

`BSMPTUseVectorization=OFF` is set explicitly.  In particular, neither
`-march=native` nor fast-math is enabled.  The profile-generating configure
command is:

```bash
cmake -S upstream -B build-pgo-gcc -DCMAKE_BUILD_TYPE=Release \
  -DBSMPTUseVectorization=OFF \
  -DCMAKE_C_FLAGS='-O3 -g -fprofile-generate -fno-fast-math -fno-tree-vectorize -fPIC -fstack-protector-strong -fno-plt -ffunction-sections -pipe -m64' \
  -DCMAKE_CXX_FLAGS='-O3 -g -fprofile-generate -fno-fast-math -fno-tree-vectorize -fPIC -fstack-protector-strong -fno-plt -ffunction-sections -pipe -m64' \
  -DCMAKE_EXE_LINKER_FLAGS='-fprofile-generate -m64 -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,/home/yancy/Package_Management/anaconda3/envs/root63604py311/lib -Wl,-rpath-link,/home/yancy/Package_Management/anaconda3/envs/root63604py311/lib -L/home/yancy/Package_Management/anaconda3/envs/root63604py311/lib -Wl,--defsym,__TMC_END__=0 -L/home/yancy/Package_Management/anaconda3/envs/root63604py311/lib/gcc/x86_64-conda-linux-gnu/14.3.0'
```

Train with the trusted `mH>125` input (the legend is line 1 and the point is
line 2), using the validated opt-in environment block above:

```bash
cmake --build build-pgo-gcc --target CalcGW -j 1
<validated-opt-in-environment> build-pgo-gcc/bin/CalcGW R2HDM \
  benchmark_input_mh_gt_125.tsv pgo_train_mh_gt.tsv 2 2
```

After training, reconfigure the same tree with `-fprofile-use -fprofile-correction`
in place of `-fprofile-generate` (retain `-O3`, `-fno-fast-math`,
`-fno-tree-vectorize`, and `BSMPTUseVectorization=OFF`), then build `CalcGW`
again.  The `.gcda` files must remain beside their generating objects until
the profile-use build is complete; any source/control-flow change requires a
fresh generate-and-train pass.

## Safety notes

## Current bounce hot-path audit (2026-09-01)

Phase-only profiling of the exact-fast configuration on the `mH>125`
reference point measured 8.316 s in vacuum/phase tracing and 26.894 s in
bounce construction.  The later nucleation, percolation, completion and GW
parameter stages together used less than 0.14 s.  Aggregating the per-action
instrumentation further split bounce work as follows:

```text
solve_1d          24 calls   21.523 s
path_check        12 calls    1.352 s
action_integral   12 calls    0.664 s
path_deformation  12 calls    0.442 s
```

Thus the remaining exact bottleneck is the dense exact-solution threshold
scan and its numerical gradient/Hessian potential evaluations, not the final
GW integration.

Three allocation/control-flow experiments were bit-for-bit exact but rejected
after paired testing.  Quark eigensolver workspace reuse averaged 32.72 s
versus 32.43 s control; RK5 vector workspace reuse measured 32.684 s versus
31.623 s; and skipping the unused saturated logistic branch measured 32.870 s.
Their hot-path branches were removed again so failed experiments do not bloat
the production benchmark binary.

The VEff component profile on the `mH>125` reference point attributes the
largest remaining costs to the complete quark, Higgs, and lepton spectra.
After direct fixed-matrix construction was added, the bit-for-bit flat
CounterTerm stream also improved the tested `mH<125` workloads and is now in
the generic wrapper.

The generic wrapper enables full fixed-size Eigen solves for
the R2HDM `diff==0` fermion spectra: a complete complex 12x12 quark matrix and
a complete complex 9x9 lepton matrix.  These are not block reductions and do
not duplicate analytically degenerate eigenvalues.  On `mH>125`, the two
fixed-size solvers reduce a same-build control from 31.48 s to 30.50 s; with
the flat CounterTerm stream the specialized wrapper measures 29.90 s.  Every
non-runtime output field is bit-for-bit identical.  With direct fixed-matrix
construction, current `mH<125` tests also show a benefit: the main point was
93.511 s with fixed9+fixed12 and the four-point set total was 121.957 s.

`BSMPT_USE_R2HDM_LEPTON_ANALYTIC_DIFF0=1` is deliberately excluded from all
exact wrappers.  It reduces runtime by 13.8--18.2%, but validation on the
`mH<125`, SNR~150 point changed total SNR by -2.74% and beta/H by -0.757%.
It is therefore an approximate research switch, not a production scan mode.

An isolated build disabling stack protectors in Models and BounceSolution was
bit-for-bit exact but took 38.65 s on the `mH>125` reference and was rejected.

Do not invoke the upstream `Build.py` in this environment without isolating
`CONAN_HOME`: its setup step attempts to remove and recreate the existing
`~/.conan2/profiles/BSMPT` directory.

The exact wrapper also enables an explicitly unrolled order-four low-temperature
fermion series.  It preserves the original `pow` calls and accumulation order.
Seven validation points are bit-for-bit identical outside `runtime`; paired
measurements showed a modest 0.3--2.0% improvement, including 73.395 s versus
74.449 s on the SNR~150 point.

The paired Higgs/gauge paths now reserve their small result vectors before
appending eigenvalues.  This allocation-only change passed the same seven-point
zero-tolerance gate and measured 30.527 s versus 30.821 s on `mH>125`, and
95.433 s versus 96.280 s on `mH<125`.

When detailed VEff profiling is disabled, the wrapper also caches its
process-wide timing gate inside V1Loop instead of querying it for every
particle contribution.  Profiling-on behavior is retained.  The seven-point
exact gate passed; paired main-point timings improved by 0.31--1.75%.

Three additional R2HDM points, not used by the earlier CalcGW acceptance gate,
confirm that the latest three micro-optimizations generalize: their paired
batch improved by 2.86% and remained bit-for-bit identical outside `runtime`.
A fermion-low4 ablation attributed 2.12% to that specialization, while the
result-reserve plus profiler-gate pair contributed 1.26% in a separate run.
An analogous boson-low3 specialization was rejected: despite positive timing
on five paired points, the SNR~150 gate changed 53 non-runtime fields.

A second independent validation batch adds eight previously unused R2HDM
points.  Seven complete the GW calculation and one retains a `nan` GW boundary
status.  All eight match the pre-micro exact control bit-for-bit outside
`runtime`.  Four representative points were also compared with every
exact-fast opt-in disabled: all fields still matched exactly, while aggregate
runtime fell from 392.366 s to 170.037 s (56.66%, or 2.308x).

An NLO-boundary stress test interpolates between one fully valid point and one
`no_nlo_stability` point.  Five valid points immediately around alternating
valid/invalid regions match the dense control bit-for-bit and improve from
383.043 s to 183.053 s (52.21%, 2.093x).  Four adjacent invalid points also
match exactly and remain rejected as `no_nlo_stability`, confirming that the
optimized path does not move the tested acceptance boundary.

Two further ten-point stratified groups broaden the parameter coverage.  Group
A mixes low/intermediate soft masses, higher `tan(beta)`, both signs of `L5`,
and varied couplings; group B stresses high soft masses, low `tan(beta)`, and
long-tail paths.  Both groups match the dense control bit-for-bit in every
non-runtime field (20/20 rows).  Aggregate runtime improved by 50.61% in A
and 49.32% in B, while covering successful GW, `nan`, `not_set`, `non_bfb`,
and `no_nlo_stability` outcomes.

A third ten-point broad group stresses negative `L3`, large positive/negative
quartics, and additional mass branches.  All 10 rows match the dense control
bit-for-bit outside `runtime`, including five complete GW results, four `nan`
GW outcomes, one NLO rejection, and a 313 s optimized long tail.  Aggregate
runtime improved from 1055.075 s to 556.091 s (47.29%, 1.897x).

Yukawa coverage now includes types 2, 3, and 4.  A nine-row scalar-point/type
matrix matches exactly and improves from 810.131 s to 405.369 s (49.96%,
1.999x), directly exercising the complete fixed 12x12 quark and 9x9 lepton
Eigen paths.  A separate three-row high-SNR scalar point also matches exactly:
type 3 retains SNR 150.868 while types 2 and 4 enter weak-signal branches;
runtime improves from 430.514 s to 220.162 s (48.86%).

The high-SNR type-3 point was additionally checked under every CalcGW
`multistepmode`: `0`, `1`, `2`, and `auto`.  All four optimized/dense pairs are
bit-for-bit identical outside `runtime`.  Their aggregate runtime is 358.033 s
versus 702.271 s (49.02% lower); mode 2 exercises the heavier global-minimum
tracing path.  Across the expanded boundary, broad, Yukawa, and mode matrices,
at least 73 new row-level dense comparisons now pass the zero-tolerance gate.

An experimental, non-strict profile is available through
`run_calcgw_approx_safe.sh`.  It uses the analytic R2HDM lepton spectrum only
for a first pass and automatically reruns a single risky point with the
validated exact-fast profile.  Risk includes failed/non-finite/incomplete
statuses, very small total SNR, and configurable margins around SNR cuts.  The
current 42-row study found two counterexamples to unguarded approximation; both
are caught by this policy.  Among the six rows accepted by the default guard,
the maximum observed SNR-component error is 4.36%.  This is empirical evidence,
not a whole-parameter-space guarantee; see `APPROX_SAFE_REPORT_ZH.md`.

The default guarded first pass now combines the analytic lepton spectrum with
the central two-point gradient.  On the 42-row matrix its raw runtime is 26.52%
lower than exact-fast.  The guard accepts six rows (maximum observed SNR
component error 3.89%) and sends 36 risky/weak rows to exact-fast.  In
particular, it catches a central2 false-positive GW transition because its
signal lies below the conservative SNR floor.  The strict wrapper and binary
are unchanged.

A four-point relative-perturbation cloud around that false-positive boundary
was also checked.  Central2 flipped the GW outcome on all four points in both
directions, but every case is rejected by the guard: three false successes
have total SNR below `1e-20`, and the fourth reports a failed/non-finite GW.
This boundary region has essentially no raw speed benefit and becomes slower
after exact fallback; the guarded profile targets stable signal regions, not
uniform acceleration of every point.

The current default first pass additionally enables the 64-point adaptive
exact-solution threshold search.  Across the same 42-row matrix it reduces raw
runtime by 49.11% (2173.759 s to 1106.233 s), while the six guard-accepted rows
have a maximum observed SNR-component error of 3.72%.  Four CalcGW multistep
modes and the four-point false-positive perturbation cloud also pass the guard
policy.  Analytic-only and central2-only wrappers remain available for
ablation; none of these changes rebuilds or modifies the strict binary.

The validated guarded default now also uses a 500-interval bounce raster.  On
the 42-row matrix its raw first-pass runtime is 940.724 s versus 2173.759 s for
exact-fast (56.72% lower), while the six accepted rows remain below 3.89%
maximum observed SNR-component error.  NLO boundary, broad A/B/C, high-SNR
Yukawa, and all four multistep modes were rerun for this configuration.

Further raster reductions to 250 and 100 intervals were not promoted.  Their
incremental gains over raster500 shrink to roughly 3--6% on the high-SNR and
NLO probes, while the lower interpolation density expands untested numerical
risk.  The committed wrappers remain experimental ablations; raster500 stays
the guarded default.

Two further ablations were rejected.  Replacing central2 with the analytic
gradient exceeded the 10% component-SNR limit on the five-point NLO boundary
set (13.19%) and was slower on the accepted high-SNR probe.  Reducing the
adaptive grid from 64 to 32 produced identical non-runtime output on the
high-SNR and NLO probes but was slower in both measurements.  Neither changes
the guarded default.

An intermediate 400-interval raster was also rejected.  It preserved all
status/history fields on the NLO, high-SNR, and broad-C probes, but improved
the three-row high-SNR runtime by only 0.45% and was slower than raster500 on
the ten-row broad-C set.  Its known A/B boundary divergences remained caught
by the guard.

For a chronological, commit-linked audit trail of every accepted and rejected
research round, see `RESEARCH_JOURNAL_ZH.md`.  New experiments must update that
journal with their hypothesis, sample scope, correctness/runtime result,
decision, artifacts, and commit/push state.

The safe runner also performs a conservative input prefilter using
`approx_exact_direct_anchors.tsv`.  It sends only explicitly validated risk
anchors/neighborhoods directly to exact-fast, avoiding an approximate pass
that is already known to require fallback.  A non-match never bypasses the
existing output guard.  New anchors must not be added without saved
exact/approx neighborhood evidence.

The current guarded first pass additionally uses the isolated
`build-approx-thermal/bin/CalcGW` and `BSMPT_USE_THERMAL_FAST_POWERS=1` through
`run_calcgw_approx_c2_adaptive_r500_thermal_fast.sh`.  The switch replaces
only integer/half-integer powers in the existing low-temperature expansions;
it remains off in the strict binary.  Across the 42-row matrix this reduced
the recorded raw first-pass total from 940.724 s to about 908.155 s, with the
same six accepted rows and at most 3.80% accepted SNR-component error.  The
strict fallback is launched after the first-pass subprocess exits and retains
the unchanged exact-fast binary.

Remaining work is planned before execution in
`REMAINING_OPTIMIZATION_ROADMAP_ZH.md`.  It fixes candidate priority, staged
validation and stop conditions, the two-CalcGW memory limit, and a one-at-a-
time read-only Luna audit workflow to control compute and model usage.
