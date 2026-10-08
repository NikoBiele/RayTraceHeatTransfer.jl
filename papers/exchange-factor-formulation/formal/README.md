# Formal Verification of the Graph Equilibrium Radiative Transfer (GERT) Exchange Factor Formulation

Machine-checked proofs (Lean 4) of all theorems in:

> N. M. Bielefeld, *A Radiation Exchange Factor Formulation with Proven
> Non-Negativity and Unconditional Energy Conservation*, arXiv:2512.22157v2.

Every theorem in Sections 3–4 of the paper is formalized and verified.
The development also produced two corrections and several sharpenings,
which will appear in v3 of the paper; they are summarized below and
recorded in detail in the file headers.

## Layout

Two independent tracks, per theorem:

- **Core track** (`ThmXX_basic.lean`): self-contained files with **zero
  dependencies** — no Mathlib, no imports. Each states its algebraic
  axioms explicitly (all satisfied by ℝ; a ℤ instance certifies
  consistency) and carries `#eval` smoke tests. Checkable with a bare
  `lean` binary. These cover the entrywise/algebraic layer.
- **Mathlib track** (`ThmXX_mathlib.lean`): full-strength statements over
  ℝ and ℂ — complex eigenvalues, determinants, genuine invertibility.
  Require a Mathlib project (see *Building*).

## Theorem ↔ file map

| Paper | Content | Core | Mathlib | Status |
|---|---|---|---|---|
| Thm 3.1 | A + R = F; Perron vector 𝟙 | `Thm31_basic.lean` | `Thm31_mathlib.lean` | ✅ verified as published |
| Thm 3.2 | ρ(R) ≤ b_max; D invertible | `Thm32_basic.lean` | `Thm32_mathlib.lean` | ✅ verified as published |
| Thm 3.3 | C singular M-matrix | `Thm33_basic.lean` | `Thm33_mathlib.lean` | ✅ verified as published (incl. det C = 0, Re μ ≥ 0 ∀ eigenvalues) |
| Thm 3.4 | D non-singular M-matrix | `Thm34_basic.lean` | `Thm34_mathlib.lean` | ✅ verified as published (incl. Re μ ≥ 1 − b_max) |
| Thm 3.5 | Mixed-boundary non-singularity | — | `Thm35multi_mathlib.lean` | ✅ verified, **corrected** (see Findings 1–2) |
| Thm 3.5 (single-row analysis) | — | — | `Thm35single_mathlib.lean` | superseded by the multi-row file; retained as the record of the original proof's gap and its reciprocity-based repair |
| Thm 4.1 | Non-negative radiation | — | `Thm41_mathlib.lean` | ✅ verified via an **alternative self-contained proof** (see Finding 4) |
| Thm 4.2 | Unconditional energy conservation | `Thm42_basic.lean` | `Thm42_mathlib.lean` | ✅ verified as published |

Files are standalone by design (definitions repeated per file); each can
be checked in isolation.

## Findings (to appear in paper v3)

1. **Thm 3.5 proof gap (repaired).** The published Step 3 infers
   nullity(Mᵀ) < nullity(Cᵀ) from 𝟙 ∉ null(Mᵀ); this presupposes
   null(Mᵀ) ⊆ null(Cᵀ), which does not hold in general. The case
   x_j ≠ 0 requires the strict positivity of the right null vector of C —
   supplied either by Perron–Frobenius or, more directly, by the
   reciprocity relation already asserted in the paper (eq. 7), whose
   capacity vector w is an explicit positive null vector.
   (`Thm35single_mathlib.lean`.)
2. **Thm 3.5 statement correction.** The theorem requires the additional
   hypothesis ∃ j₀ ∈ S with b_{j₀} < 1 (S the set of replaced rows): if
   b ≡ 1 on S, every replaced row of D equals the row of C it replaces
   and M = C is singular. This is the abstract's standing assumption
   "maximum reflection-scattering coefficient strictly less than unity,"
   localized to S. (`Thm35multi_mathlib.lean`.)
3. **Thm 3.5, simplified proof, zero classical citations.** The verified
   multi-row proof is a discrete maximum principle on the masked
   reflection operator F∘ḡ (g = b on S, g = 1 off S) and eliminates the
   Perron–Frobenius input entirely; the Horn & Johnson citation becomes
   background. Theorem 3.2(c) is recovered as the special case S = all
   rows, unifying the two non-singularity results.
4. **Thm 4.1 alternative route.** The published proof cites the classical
   inverse-positivity of non-singular M-matrices (not presently in
   Mathlib). The verified proof instead uses a minimum principle in
   reciprocity-scaled variables y = j/w — the exact dual of Finding 3's
   maximum principle. Consequence: uniqueness (3.5) and non-negativity
   (4.1) are the same argument, the discrete elliptic maximum principle.
   The classical M-matrix route remains valid but is not machine-checked.
5. **Hypothesis accounting.** Every theorem consumes exact
   row-stochasticity (conservation); Theorem 4.1's self-contained proof
   additionally consumes reciprocity. Conservation is a physical law;
   reciprocity is a kernel symmetry contingent on the isotropic-emission
   model, and it is what removes all classical citations. Raw Monte Carlo
   estimates satisfy neither exactly; the smoothing procedure of the
   companion paper enforces the joint structure {conservation ∧
   reciprocity} to machine precision — a genuine constraint-intersection
   problem (a single row normalization restores conservation alone but
   breaks symmetry, and vice versa). If the framework is extended to
   anisotropic distributions, reciprocity fails and the machine-checked
   positivity guarantee no longer applies as formalized; positivity then
   rests on the classical M-matrix theorem alone.
6. **Minor sharpenings** recorded in file headers: Thm 3.1 Step 1 needs
   no finiteness; Thm 3.4's diagonal needs only b ≤ 1 (not 0 ≤ b), and
   D = C + Aᵀ is hypothesis-free; in-edge existence needs no
   non-negativity; irreducibility is formalized as strong connectivity of
   the positive-entry relation (`StronglyConnected`), verbatim §2.1.

## Building

**Core files** — no setup:
```
lean Thm31_basic.lean    # any Lean >= 4.15; developed and checked on 4.33.0 and 4.34.0-rc1
```

**Mathlib files** — inside a Mathlib project:
```
lake new gertformal math
cd gertformal && lake exe cache get
# copy the *_mathlib.lean files into the project, then
lake env lean Thm35multi_mathlib.lean
```
The `lean-toolchain` in this directory pins the verified compiler
version. Mathlib evolves; lemma-name drift may require occasional
one-line fixes in future snapshots.

## Provenance and scope of the verification claim

The Lean formalizations — both the formal statements and the proofs —
were formulated and drafted by Anthropic's Claude (Fable 5) in
collaboration with the author. The author compiled and machine-checked
every file (the Mathlib track against a local Mathlib build; the core
track additionally checked in-session), iterating with Claude until all
files verified without errors or `sorry`.

**Scope of the claim.** The author is not proficient in Lean and does
not independently certify that the formal statements are faithful
renderings of the paper's theorems. The precise claim is therefore
conditional: *if* the formal statements below correspond to the paper's
theorem statements, *then* every theorem of the paper is machine
verified. This is the standard trusted base of any formalization (the
"specification gap"): proof checkers verify proofs against statements;
statement faithfulness is a human-audit question.

**Auditing faithfulness.** The audit surface is deliberately small:
only the theorem *statements* need human review, not the proofs. Each
statement is a few lines, cross-referenced in its docstring to the
paper's equation numbers, and each core-track file carries `#eval`
smoke tests on concrete matrices. Readers proficient in Lean are
invited to audit the statements against the paper; corrections and
issue reports are welcome.