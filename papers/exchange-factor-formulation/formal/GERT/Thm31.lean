import GERT.Defs

/-!
# Theorem 3.1: A + R = F, and 𝟙 is the Perron vector

Bielefeld (2026), arXiv:2512.22157v3, Section 3.1.

* `thm31_step1`  — A + R = F (eq. 19); pure entrywise algebra, no finiteness
* `thm31_step2`  — row-stochasticity transfers to A + R (eq. 20)
* `thm31_eigvec` — (A + R) 𝟙 = 𝟙, i.e. 𝟙 is an eigenvector with eigenvalue 1

The spectral part of the theorem, ρ(A + R) = 1 by Perron–Frobenius, is not
formalized; it is not used by any later result.
-/

open Matrix Finset

namespace GERT

variable {n : Type*} {F : Matrix n n ℝ} {b : n → ℝ}

/-- **Theorem 3.1, Step 1**: A + R = F (eq. 19). -/
theorem thm31_step1 : Amat F b + Rmat F b = F := by
  ext i j
  simp only [Matrix.add_apply, Amat, Rmat, Matrix.of_apply]
  ring

section RowSums

variable [Fintype n]

/-- **Theorem 3.1, Step 2a**: row-stochasticity transfers to A + R (eq. 20). -/
theorem thm31_step2 (hF1 : ∀ i, ∑ j, F i j = 1) :
    ∀ i, ∑ j, (Amat F b + Rmat F b) i j = 1 := by
  rw [thm31_step1]
  exact hF1

/-- **Theorem 3.1, Step 2b**: 𝟙 is an eigenvector of A + R with eigenvalue 1. -/
theorem thm31_eigvec (hF1 : ∀ i, ∑ j, F i j = 1) :
    (Amat F b + Rmat F b) *ᵥ (fun _ => (1:ℝ)) = fun _ => (1:ℝ) := by
  rw [thm31_step1]
  funext i
  simp only [Matrix.mulVec, dotProduct, mul_one]
  exact hF1 i

end RowSums

end GERT
