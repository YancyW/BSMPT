# R2HDM gauge-mass analytic reduction (read-only design)

This note is a proof and implementation design only.  No model source was
modified.  The statements below refer to the isolated source in
`upstream/`, in particular `ClassPotentialR2HDM::SetCurvatureArrays`,
`Class_Potential_Origin::GaugeMassesSquared`, and
`Class_Potential_R2HDM::CalculateDebyeGaugeSimplified`.

## 1. Exact matrix structure for all eight fields

The R2HDM interaction-basis field order is

\[
 v=(\rho_1,\eta_1,\rho_2,\eta_2,\zeta_1,\psi_1,\zeta_2,\psi_2).
\]

Define

\[
\begin{split}
 R^2 &= \rho_1^2+\eta_1^2+\rho_2^2+\eta_2^2+
        \zeta_1^2+\psi_1^2+\zeta_2^2+\psi_2^2,\\
 q &= \rho_1\zeta_1+\eta_1\psi_1+\rho_2\zeta_2+\eta_2\psi_2,\\
 p &= \rho_1\psi_1-\eta_1\zeta_1+\rho_2\psi_2-eta_2\zeta_2,\\
 d &= \rho_1^2+\eta_1^2+\rho_2^2+\eta_2^2-
      \zeta_1^2-\psi_1^2-\zeta_2^2-\psi_2^2.
\end{split}
\]

Writing \(g=C_g\), \(g'=C_{gs}\), and \(\kappa=gg'/2\), direct contraction
of the nonzero entries of `Curvature_Gauge_G2H2` with

\[
 M_{ab}^{\rm field}=\frac12\sum_{ij}G_{abij}v_i v_j
\]

gives, for *every* point in the eight-dimensional field space,

\[
 M_{\rm field}^2=
 \begin{pmatrix}
 A&0&0&u\\
 0&A&0&w_1\\
 0&0&A&w_2\\
 u&w_1&w_2&B
 \end{pmatrix},
\quad
 A=\frac{g^2R^2}{4},\quad B=\frac{g'^2R^2}{4},
\]

with

\[
 (u,w_1,w_2)=\kappa(q,p,d).
\]

The first three gauge indices are `(W1,W2,W3)` and the fourth is `B`.  The
factor \(1/2\) in the code is important: each off-diagonal field entry occurs
twice, so for example `M(0,3)=gg' q/2`, not `gg' q/4`.

This is an arrowhead matrix, not a block-diagonal matrix in a generic point.
It remains valid when charged fields are nonzero (charge-breaking, CB) and
when CP-odd fields are nonzero (CP-violating).  Therefore a proposed fast
path must not assume the usual neutral CP-conserving vacuum.

## 2. Debye terms and the exact eigenvalue formula

For R2HDM, `CalculateDebyeGaugeSimplified` sets

\[
 D=\operatorname{diag}(2g^2,2g^2,2g^2,2g'^2).
\]

`GaugeMassesSquared` adds \(T^2D\) whenever `Temp != 0`, so define

\[
 A_T=\frac{g^2R^2}{4}+2g^2T^2,\qquad
 B_T=\frac{g'^2R^2}{4}+2g'^2T^2.
\]

The Debye correction changes only the two diagonal entries and does not alter
the arrowhead vector \(t=(u,w_1,w_2)\).  At `Temp == 0`, the terms vanish
exactly in the current code.  With \(t^2=u^2+w_1^2+w_2^2\), the exact
eigenvalues are

\[
 \lambda_A=A_T\quad\text{(multiplicity two)},
\]

and

\[
 \lambda_\pm=\frac{A_T+B_T\pm
 \sqrt{(A_T-B_T)^2+4t^2}}{2}.
\]

Proof: any vector `(x,0)` with \(x\cdot t=0\) is an eigenvector with
eigenvalue \(A_T\); the orthogonal complement is spanned by the normalized
`t` direction and `B`, where the remaining matrix is
\(\left(\begin{smallmatrix}A_T&|t|\\|t|&B_T\end{smallmatrix}\right)\).

Useful special cases:

* Neutral CP-conserving configurations have charged fields and `psi` fields
  zero.  Thus `q=p=0`, `d=-R^2`; only the `(W3,B)` 2x2 sector mixes and the
  `T=0` lower eigenvalue is the photon zero mode.
* Neutral CP-violating configurations still have `q=p=0` when charged fields
  vanish; the same gauge-mass form follows without imposing `psi=0`.
* CB configurations generally have nonzero `q` and/or `p`, so all three W
  directions can couple to `B`.  The two-fold `A_T` eigenvalue is still exact.

## 3. Eigen ordering and zero-mass conventions

For `diff == 0`, Eigen's `SelfAdjointEigenSolver` returns eigenvalues in
nondecreasing order.  `GaugeMassesSquared` preserves that order and then
replaces each value satisfying

\[
 |\lambda|<10^{-5}
\]

by exactly zero.  This threshold is on mass-squared in the code (nominally
GeV\(^2\)), not on the mass.  An analytic implementation must apply the same
strict `< 1e-5` clipping *after* sorting; do not clip before comparing the
analytic values.

For `diff > 0`, and for `diff == -1`, the current implementation calls
`FirstDerivativeOfEigenvalues`.  Its first `NGauge` entries are the sorted
eigenvalues and its next `NGauge` entries are the corresponding derivatives.
It uses `EVThres=1e-6`: eigenvalues and derivatives with absolute value below
that threshold are set to zero.  A replacement must keep this two-block output
layout and threshold behavior.

## 4. Derivatives and degeneracies

Away from degeneracy, for a field coordinate `x`, let dots denote derivatives
with respect to `x`.  Then

\[
 \dot\lambda_\pm=\frac12\left[\dot A_T+\dot B_T
 \pm\frac{(A_T-B_T)(\dot A_T-\dot B_T)+4t\cdot\dot t}
 {\sqrt{(A_T-B_T)^2+4t^2}}\right],
\]

where Debye terms have no field derivative,
\(\dot A_T=g^2x/2\), \(\dot B_T=g'^2x/2\), and `dot t` follows by differentiating
`q,p,d`.  The two `A_T` branches have derivative `dot A_T` when they are
isolated from `lambda_+` and `lambda_-`.

The formula above must **not** be used when the denominator is numerically
degenerate.  The existing `FirstDerivativeOfEigenvalues` handles this by
grouping eigenvalues within `1e-6`, projecting `MDiff` onto the degenerate
eigenspace, and diagonalizing that projected matrix.  An analytic opt-in must
fall back to this exact existing routine when

\[
 \Delta=\sqrt{(A_T-B_T)^2+4t^2}\le 1e-6,
\]

or when an `A_T` branch is within `1e-6` of a \(\lambda_\pm\) branch.  At
`t=0`, the W block is a three-fold `A_T` degeneracy (unless it meets `B_T`);
its first-order projected derivatives are three copies of `dot A_T`.  If also
`A_T=B_T`, all four states are degenerate and the derivative matrix must be
diagonalized as a full 4x4 projected block; off-diagonal `dot t` can then split
the first-order derivatives.  These cases are why an unconditional closed-form
derivative replacement is unsafe.

The two-fold `A_T` degeneracy is not a numerical accident.  A generic Eigen
eigensolver may choose any orthonormal basis in that subspace, but the
projected-derivative prescription is basis independent and is the convention
that must be reproduced.

## 5. Floating-point implementation details

Use

```text
delta = hypot(A_T - B_T, 2 * hypot(u, hypot(w1, w2)))
```

instead of explicitly squaring large numbers.  The lower root can suffer
cancellation when `A_T+B_T` and `delta` are close.  Prefer the determinant
identity

\[
 \lambda_-=(A_TB_T-t^2)/\lambda_+
\]

when `lambda_plus` is safely away from zero, with the direct expression as a
fallback near zero.  The result must still be clipped with the current
`ZeroMass=1e-5` rule.  Keep `double`, matching Eigen's `MatrixXd`; using long
double only in the fast path would produce small last-bit differences in the
thermal potential and make exact output comparisons less meaningful.

The `diff == 0` analytic path can be exact up to ordinary floating-point
rounding.  It is a good candidate for speed optimization because it avoids
constructing a dynamic 4x4 matrix and invoking a general eigensolver.  The
`diff > 0`/`-1` path should initially retain Eigen for ordinary and degenerate
points until a component-wise comparison has been completed.

## 6. Safe opt-in plan

The safe sequence is:

1. Add a non-default environment switch, for example
   `BSMPT_USE_ANALYTIC_GAUGE_MASSES=1`, in the isolated build only.
2. In the R2HDM gauge routine, use the arrowhead formula only for `diff == 0`;
   retain the existing implementation for derivatives.
3. Apply Eigen ordering and the exact `1e-5` clipping after the formula.
4. For every test point, compare the four returned values against the current
   Eigen result using both absolute and relative errors.  Include random
   eight-dimensional points, neutral CP-conserving points, neutral CP-odd
   points, CB points, `T=0`, and finite `T`.
5. Require extra tests near `t=0`, `A_T=B_T`, the photon zero mode, and values
   around `1e-5`.  Disable the analytic path and fall back to Eigen whenever
   the degeneracy guard triggers.

Only after this comparison is stable should a derivative (`diff > 0` or
`diff == -1`) reduction be attempted.  This preserves the current treatment
of repeated eigenvalues and avoids silently changing the ordering used by the
thermal one-loop potential.
