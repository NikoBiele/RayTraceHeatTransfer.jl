# Formal Verification of the GERT Exchange Factor Formulation

Machine-checked proofs (Lean 4, Mathlib) of the main theorems in:

> N. M. Bielefeld, *Non-Negativity, Non-Singularity and Energy Conservation of
> the Zonal Method with Non-Reciprocal Exchange Factors*, arXiv:2512.22157v3.

The proofs of the paper's third version rest on a single discrete minimum
principle (Lemma 3.4), from which the non-singularity of the mixed-boundary
system (Theorem 3.5) and the non-negativity of its solution (Theorem 4.1)
both follow. This development machine-checks that lemma and both theorems,
together with the identities of Theorems 3.1 and 4.2 and the spectral bound of
Theorem 3.2. No reciprocity of the exchange factors is assumed anywhere, as
in the paper.

## Layout

The verification is a Lake project. All modules share one set of definitions.

- `GERT/Defs.lean` — the objects of the paper's Section 2 and the masked form
  of the mixed-boundary matrix (eqs. 23–25): `Amat`, `Rmat`, `Cmat`, `Dmat`,
  `MmatS` (the mixed matrix for a prescribed set S), `bmask` (the paper's b̃),
  `StronglyConnected` (irreducibility, §2.1), and the entrywise reduction
  Mᵀx = c ⟹ xₖ = Σₗ Fₖₗ b̃ₗ xₗ + cₖ.
- `GERT/MinPrinciple.lean` — Lemma 3.4.
- `GERT/Thm31.lean`, `GERT/Thm32.lean`, `GERT/Thm35.lean`, `GERT/Thm41.lean`,
  `GERT/Thm42.lean` — one module per theorem.
- `core/` — the no-Mathlib track from the second version: self-contained files
  with zero imports, each stating its algebraic axioms explicitly (satisfied by
  ℝ, with a ℤ instance certifying consistency) and carrying `#eval` smoke
  tests. They cover the entrywise identities of Theorems 3.1, 3.2 and 4.2 and
  are checkable with a bare `lean` binary.
- `retired/` — the second version's files for the theorems that v3 removed
  (the M-matrix properties of C and D) and for the proof routes that v3
  replaced (the single-row and maximum-principle proofs of Theorem 3.5, and
  the reciprocity-based proof of Theorem 4.1). They are kept as the record of
  how the gaps in v2 were found; they are not part of the build.

## Theorem ↔ module map

| Paper | Content | Module | Main declarations |
|---|---|---|---|
| Thm 3.1 | A + R = F; 𝟙 is an eigenvector with eigenvalue 1 | `GERT/Thm31.lean` | `thm31_step1`, `thm31_step2`, `thm31_eigvec` |
| Thm 3.2 | ρ(R) ≤ bmax; D invertible for bmax < 1 | `GERT/Thm32.lean` | `thm32_eigenvalue_bound`, `thm32_uniform_attained`, `thm32_D_isUnit` |
| Lemma 3.4 | Masked minimum principle | `GERT/MinPrinciple.lean` | `masked_min_principle` |
| Thm 3.5 | Non-singularity of mixed boundaries | `GERT/Thm35.lean` | `thm35_isUnit`, `Dmat_isUnit_of_irreducible` |
| Thm 4.1 | Non-negative radiation (M⁻¹ ≥ 0, j ≥ 0, r ≥ 0, e ≥ 0) | `GERT/Thm41.lean` | `thm41_inv_nonneg`, `thm41_solution_nonneg`, `thm41_reflected_nonneg`, `thm41_emissive_nonneg` |
| Thm 4.2 | Unconditional energy conservation | `GERT/Thm42.lean` | `thm42_energy_conservation` |

Every module docstring restates the theorem as the paper gives it, and every
declaration's docstring names the equation numbers it corresponds to.

Hypotheses consumed, in every result: F non-negative with rows summing to
one, and, where connectivity is needed (Lemma 3.4, Theorems 3.5 and 4.1),
strong connectivity of the positive-entry relation. Theorem 4.2 uses only the
row sums. None of the modules assumes reciprocity.

Not formalized: the spectral radius statement of Theorem 3.1 (Perron–Frobenius),
the conditioning bound of Corollary 3.3, the albedo invariance of Proposition
4.3, and the equivalence with Hottel's zonal balances of Section 5.

## History

The formalization of the second version of the paper found two problems,
both corrected in the third version and recorded in its Appendix B: the proof
of non-singularity of the mixed-boundary system had a gap, and the proof of
non-negativity inferred the M-matrix property from the sign pattern and
non-singularity alone, which does not follow in general. The repair in v2's
files used a maximum principle for non-singularity and a minimum principle in
reciprocity-scaled variables for non-negativity. The third version of the
paper derives both from one minimum principle without reciprocity, and the
modules in `GERT/` follow that proof.

## Building

Requires `elan`. From this directory:

```
lake exe cache get
lake build
```

The first command downloads the compiled Mathlib matching
`lake-manifest.json`; the second checks every module. `lean-toolchain` pins
the compiler version the modules were checked with. Mathlib evolves, so a
later snapshot may need occasional one-line fixes to lemma names.

The core files need no project:

```
lake env lean core/Thm31_basic.lean
```

## Provenance and scope of the verification claim

The Lean files, both the formal statements and the proofs, were drafted by
Anthropic's Claude (Fable 5 for the second version, Fable 5.1 for the third)
in collaboration with the author, who compiled and machine-checked every file
and iterated with the model until the project built without errors or
`sorry`.

**Scope of the claim.** The author is not proficient in Lean and does not
independently certify that the formal statements are faithful renderings of
the paper's theorems. The precise claim is therefore conditional: *if* the
formal statements correspond to the paper's theorem statements, *then* those
theorems are machine verified. This is the standard trusted base of any
formalization: a proof checker verifies proofs against statements, and
statement faithfulness is a human-audit question.

**Auditing faithfulness.** The audit surface is small: only the theorem
statements need review, not the proofs. Each statement is a few lines,
cross-referenced in its docstring to the paper's equation numbers. Readers
proficient in Lean are invited to audit the statements against the paper;
corrections and issue reports are welcome.