# R2HDM Type-I: exact quark-matrix block reduction design

This note is a read-only analysis of the isolated upstream checkout.  No
upstream source file is changed by this note.  The target is the
`QuarkMassesSquared` call used by `CalcGW`, for the CP-conserving R2HDM with
Yukawa Type-I, including arbitrary real CB/CP field components in the bounce
configuration.

## Conclusion

For R2HDM Type-I, the 12 by 12 matrix constructed by
`Class_Potential_Origin::QuarkMassMatrix` has the exact chiral form

\[
 M=\begin{pmatrix}0_{6}&K\\K^T&0_{6}\end{pmatrix},
\]

where `T` is an ordinary transpose, not a Hermitian adjoint.  Therefore the
matrix used by the current code,

\[
 M^*M=\begin{pmatrix}K^*K^T&0\\0&K^\dagger K\end{pmatrix},
\]

is exactly two six-dimensional Hermitian blocks with identical real spectra:

\[
 K^*K^T=(K^\dagger K)^T.
\]

Consequently one 6 by 6 self-adjoint eigensolve is sufficient and each of its
six eigenvalues must be returned twice.  The same reduction applies to the
field derivative used by `QuarkMassesSquared(v,diff)`: if `L=\partial K/\partial
v_m`, the lower block derivative is

\[
 \partial(K^\dagger K)=L^\dagger K+K^\dagger L.
\]

The derivative eigenvalues are also duplicated.  This is an exact algebraic
reduction, not a zero-entry approximation.

## Source-level index mapping

`ClassPotentialR2HDM.h` documents the 12-component basis as

```text
0  u_R       1  c_R       2  t_R
3  d_R       4  s_R       5  b_R
6  bar(u)_L  7  bar(c)_L  8  bar(t)_L
9  bar(d)_L  10 bar(s)_L  11 bar(b)_L
```

Let the first six indices be `R(r)`, `r=0..5`, and the last six be
`L(l)=6+l`, `l=0..5`.  Then

```text
M[R(r),R(s)] = 0
M[L(l),L(q)] = 0
M[R(r),L(l)] = K[r,l]
M[L(l),R(r)] = K[r,l]
```

The last equality is the source-level symmetry `M=M^T`; it is not the
Hermiticity condition `M=M^dagger`.

## Explicit Type-I R2HDM form

In `Class_Potential_R2HDM::SetCurvatureArrays`, `NQuarks=12` and the VEV
components are

```text
rho1=0, eta1=1, rho2=2, eta2=3,
zeta1=4, psi1=5, zeta2=6, psi2=7.
```

For Type-I the quark Yukawa entries on fields `0,1,4,5` are zero.  Thus CB
components can be active in the bounce, but they do not enter the quark mass
matrix in this model.  Define

\[
 x=v_2^{(\mathrm{field})},\quad y=v_3^{(\mathrm{field})},\quad
 z=v_6^{(\mathrm{field})},\quad w=v_7^{(\mathrm{field})},
\]

and

\[
 p=x-i y,\qquad q=z-i w,\qquad r=z+i w.
\]

With `U=diag(m_u,m_c,m_t)/v_2`, `D=diag(m_d,m_s,m_b)/v_D`, and the CKM
matrix `V` (whose row index is up-type and column index is down-type), the
nonzero 3 by 3 Yukawa blocks read directly from the source:

\[
 A=-U V^*,\qquad B=U,\qquad C=D V^T,\qquad E=D.
\]

Here `v_D=v_2` for Type-I.  In the reduced ordering

```text
K rows:    [u_R,c_R,t_R | d_R,s_R,b_R]
K columns: [bar(u)_L,bar(c)_L,bar(t)_L | bar(d)_L,bar(s)_L,bar(b)_L]
```

the exact matrix is

\[
 K=\begin{pmatrix}
      q B & p A\\
      p C & r E
    \end{pmatrix}.
\]

This formula remains valid for arbitrary real values of all eight scalar
components: the four components absent from this formula simply have zero
Yukawa coefficient in Type-I.  It also shows why there is generally no second
exact reduction to independent 3 by 3 matrices: physical non-diagonal CKM
mixing makes the 6 by 6 `K` matrix coupled between the two 3-generation
sectors.

The source constructs the individual Yukawa matrices, then explicitly copies
the lower triangle from the upper triangle in `SetCurvatureArrays`.  In
R2HDM, `Curvature_Quark_F2` is initialized to zero, so there is no constant
term that spoils this form.  `QuarkMassMatrix` then adds the real field
components to those symmetric Yukawa tensors.  These facts establish
`M=M^T` for every real CB/CP configuration, not only for a neutral path.

## Why `.conjugate()` still permits this proof

The implementation uses

```cpp
MassMatrix = MIJ.conjugate() * MIJ;
```

`conjugate()` is element-wise complex conjugation.  It is not generally an
adjoint.  However, for this R2HDM Type-I matrix `MIJ=MIJ.transpose()` exactly,
so element-wise conjugation happens to equal the adjoint:

\[
 M^* = (M^T)^* = M^\dagger.
\]

The reduction above does not silently replace the requested operation: it
starts from the literal product `M.conjugate()*M`, uses the chiral block form,
and obtains the two blocks `K^*K^T` and `K^dagger K`.  An implementation must
retain a symmetry check or a safe full-matrix fallback if it is moved to a
generic model.

## Derivative form and ordering

For `diff>0`, the source sets `m=diff-1` and constructs

\[
 D_M=M^*M_m+M_m^*M,
 \qquad M_m=\partial M/\partial v_m.
\]

The derivative has the same chiral form with

\[
 L=\partial K/\partial v_m,
\quad
 D_M=\operatorname{diag}\left(L^*K^T+K^*L^T,
                               L^\dagger K+K^\dagger L\right).
\]

For Type-I, `L=0` for `m=0,1,4,5`; for `m=2,3,6,7` it is obtained by
replacing the corresponding scalar factor in `K` by its derivative.  The
lower six-dimensional block is therefore sufficient for
`FirstDerivativeOfEigenvalues`.

The current function returns first all eigenvalues and then all derivatives.
The reduced implementation must preserve this layout:

```text
reduced spectrum:      lambda[0..5]
full-compatible output: lambda[0],lambda[0], ..., lambda[5],lambda[5],
                         d[0],d[0], ..., d[5],d[5]
```

`SelfAdjointEigenSolver` sorts eigenvalues in ascending order.  Exact duplicate
copying preserves the physical sorted spectrum.  The full 12 by 12 solver can
choose either member of an exactly degenerate pair first, and floating-point
roundoff can split the two copies by tiny amounts; code must not depend on the
identity of a duplicate.  `V1Loop` pairs each mass with its derivative by
index, but duplicated masses give the same scalar contribution, so duplicate
copying is the mathematically correct stable convention.

The existing `FirstDerivativeOfEigenvalues` groups eigenvalues whose absolute
difference is below `EVThres=1e-6`, then diagonalizes the perturbation in the
degenerate subspace.  The reduced implementation should use the same function
and the same threshold.  In particular, do not use a new relative threshold
or drop small matrix entries: zero-mass postprocessing in
`QuarkMassesSquared` uses a separate `ZeroMass=1e-10` rule that must remain
unchanged.

## Safe opt-in implementation plan

The recommended experiment is an opt-in branch inside the isolated checkout,
for example controlled by `BSMPT_USE_R2HDM_QUARK_BLOCK6`.  The default path
must remain the current 12 by 12 implementation until validation is complete.

1. Build `K` and `L` from the existing `MIJ` and `Curvature_Quark_F2H1`
   entries using the index map above.  For maximum semantic safety, initially
   form the current full `MIJ` and verify that its four blocks reconstruct `K`
   and that `max(abs(M-M.transpose()))` is below a scale-aware tolerance.
2. If the model is not exactly R2HDM Type-I, dimensions are not 12, the
   symmetry check fails, an eigensolver reports failure, or a non-finite value
   appears, fall back to the original full path.
3. Compute only
   `H=K.adjoint()*K` (the lower block of the literal product) and solve the
   6 by 6 self-adjoint problem.  Return each value twice, applying the original
   `ZeroMass=1e-10` cleanup after duplication.
4. For `diff>0`, compute
   `dH=L.adjoint()*K+K.adjoint()*L`, call the existing derivative helper on
   the 6 by 6 pair `(H,dH)`, and duplicate both the six values and six
   derivatives.  Keep `diff<=0` and `diff>NHiggs` behavior identical to the
   current function.
5. During validation, at representative points from both `mH<125` and
   `mH>125`, compare the complete returned vectors from the full and reduced
   paths, including all eight field derivatives.  Require finite values,
   ascending mass order, duplicate-pair agreement, and agreement in the final
   `VEff`/`CalcGW` output within a predeclared numerical tolerance.  Also test
   the symmetric point where all Yukawa masses vanish and points with active
   CB/CP components.

A debug-only stronger check can compare the reduced `H` eigenvalues with the
two six-dimensional spectra obtained by explicitly partitioning the full
`M.conjugate()*M`; this isolates eigensolver ordering from any later thermal
function differences.  No threshold-based sparsification should be enabled in
the production opt-in path.

## Expected benefit and limitations

The eigensolver and derivative-subspace work changes from 12 by 12 to 6 by 6;
for cubic dense eigensolver work this is approximately an eightfold reduction
in that portion, in addition to halving the matrix storage and multiplication
work.  It does not remove the repeated calls to `QuarkMassMatrix` or thermal
function interpolation, so the total `CalcGW` speedup must be benchmarked rather
than inferred from the dimension ratio.  The optimization is nevertheless
exact for the stated R2HDM Type-I structure and is suitable as a low-risk
opt-in after the full-vector and end-to-end comparisons pass.
