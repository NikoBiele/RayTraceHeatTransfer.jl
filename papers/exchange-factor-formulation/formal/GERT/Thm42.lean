import GERT.Defs

/-!
# Theorem 4.2: unconditional energy conservation

Bielefeld (2026), arXiv:2512.22157v3, Section 4.2.

The identity 𝟙ᵀ(Cj) = 0 holds for EVERY vector j, independently of the
mixed system M, the prescribed quantities h and the coefficients b, because
every column of C = I − Fᵀ sums to zero. Only row-stochasticity of F enters,
which is the formal expression of "unconditional".

* `col_sum_zero`              — Σₖ Cₖᵢ = 0 for every i
* `thm42_energy_conservation` — Σₖ (C *ᵥ j)ₖ = 0 for every j
-/

open Matrix Finset

namespace GERT

variable {n : Type*} [Fintype n] [DecidableEq n] {F : Matrix n n ℝ}

/-- Every column of C sums to zero: Σₖ Cₖᵢ = Σₖ δₖᵢ − Σₖ Fᵢₖ = 1 − 1 = 0. -/
theorem col_sum_zero (hF1 : ∀ i, ∑ k, F i k = 1) (i : n) :
    ∑ k, Cmat F k i = 0 := by
  have h : ∀ k, Cmat F k i = (1 : Matrix n n ℝ) k i - F i k := by
    intro k
    simp [Cmat, Matrix.sub_apply, Matrix.transpose_apply]
  calc ∑ k, Cmat F k i
      = ∑ k, ((1 : Matrix n n ℝ) k i - F i k) := Finset.sum_congr rfl fun k _ => h k
    _ = (∑ k, (1 : Matrix n n ℝ) k i) - ∑ k, F i k := by
        rw [Finset.sum_sub_distrib]
    _ = 1 - 1 := by
        rw [hF1 i]
        congr 1
        simp [Matrix.one_apply, Finset.sum_ite_eq', Finset.mem_univ]
    _ = 0 := sub_self 1

/-- **Theorem 4.2 (Unconditional Energy Conservation)**: 𝟙ᵀ(Cj) = 0 for
    EVERY total radiant power vector j (eq. 29):
    Σₖ (C *ᵥ j)ₖ = Σᵢ (Σₖ Cₖᵢ) jᵢ = 0. -/
theorem thm42_energy_conservation (hF1 : ∀ i, ∑ k, F i k = 1) (j : n → ℝ) :
    ∑ k, (Cmat F *ᵥ j) k = 0 := by
  have hexpand : ∑ k, (Cmat F *ᵥ j) k = ∑ k, ∑ i, Cmat F k i * j i := by
    refine Finset.sum_congr rfl fun k _ => ?_
    simp [Matrix.mulVec, dotProduct]
  rw [hexpand, Finset.sum_comm]
  have h : ∀ i, ∑ k, Cmat F k i * j i = 0 := by
    intro i
    rw [← Finset.sum_mul, col_sum_zero hF1 i, zero_mul]
  rw [Finset.sum_congr rfl fun i _ => h i]
  exact Finset.sum_const_zero

end GERT
