/-
  Theorem 3.1 (Steps 1–2), Bielefeld (2026), arXiv:2512.22157v2 — Mathlib version.

  ✅ STATUS: compiled clean on user's Mathlib (this revision only rescopes
  typeclass assumptions to silence the unused-instance linters; proofs unchanged).

  Instance hygiene, which the linter helpfully surfaces as mathematical content:
  * Step 1 (A + R = F) is pure entrywise algebra — needs NO finiteness at all;
    it holds for infinite index types.
  * Step 2 (row sums / eigenvector) needs [Fintype n] for Σ, nothing more.
  * [DecidableEq n] is needed nowhere in Thm 3.1 (it first appears in Thm 3.2c
    for the identity matrix / determinant).

  NB the recurring pitfall: the paper's "1" in (1 − B) is the ALL-ONES matrix.
  Mathlib's `(1 : Matrix n n ℝ)` is the IDENTITY. All definitions below are
  therefore entrywise, never via `(1 : Matrix _ _ _) - B`.
-/
import Mathlib

open Matrix Finset

variable {n : Type*}

variable (F : Matrix n n ℝ) (b : n → ℝ)

/-- Column-constant reflection-scattering matrix B (paper eq. 8): `B i j = b j`. -/
def Bcol : Matrix n n ℝ := Matrix.of fun _ j => b j

/-- Absorption matrix A = F ∘ (𝟙 − B) entrywise (paper eq. 9). -/
def Amat : Matrix n n ℝ := Matrix.of fun i j => F i j * (1 - b j)

/-- Reflection-scattering matrix R = F ∘ B entrywise (paper eq. 10). -/
def Rmat : Matrix n n ℝ := Matrix.of fun i j => F i j * b j

/-- Faithfulness check: `Amat` is literally the Hadamard product of the paper. -/
theorem Amat_eq_hadamard :
    Amat F b = F.hadamard (Matrix.of fun _ j => 1 - b j) := by
  ext i j
  simp [Amat, Matrix.hadamard_apply]

/-- **Theorem 3.1, Step 1**: A + R = F (paper eq. 20).
    Pure entrywise algebra — no finiteness assumption needed. -/
theorem thm31_step1 : Amat F b + Rmat F b = F := by
  ext i j
  simp only [Matrix.add_apply, Amat, Rmat, Matrix.of_apply]
  ring

section RowSums
variable [Fintype n]

/-- Row-stochasticity: every row sums to 1. -/
def RowStochastic (M : Matrix n n ℝ) : Prop := ∀ i, ∑ j, M i j = 1

/-- **Theorem 3.1, Step 2a**: row-stochasticity transfers to A + R (eq. 21). -/
theorem thm31_step2 (hF : RowStochastic F) :
    RowStochastic (Amat F b + Rmat F b) := by
  rw [thm31_step1]
  exact hF

/-- **Theorem 3.1, Step 2b**: 𝟙 is an eigenvector of A + R with eigenvalue 1,
    i.e. (A + R) *ᵥ 𝟙 = 𝟙, given F row stochastic. -/
theorem thm31_eigvec (hF : RowStochastic F) :
    (Amat F b + Rmat F b) *ᵥ (fun _ => (1:ℝ)) = fun _ => (1:ℝ) := by
  rw [thm31_step1]
  funext i
  simp only [Matrix.mulVec, dotProduct, mul_one]
  exact hF i

end RowSums

#check @thm31_step1
#check @thm31_step2
#check @thm31_eigvec
