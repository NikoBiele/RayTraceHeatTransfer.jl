/-
  Theorem 4.2, Bielefeld (2026), arXiv:2512.22157v2 — Mathlib version.
  "Unconditional Energy Conservation"     𝟙ᵀCj = 0

  STATUS: rev 2 — sole rev-1 error fixed: `Finset.sum_sub_distrib` has
  explicit function arguments, so it works under `rw` but not as a bare
  calc-step term; now applied via `rw`.

  The formal statement captures the paper's "unconditional" precisely:
  the identity holds for EVERY j : n → ℝ, with no hypotheses on b, no
  mixed system, no boundary conditions — only row-stochasticity of F.
-/
import Mathlib

open Matrix Finset

variable {n : Type*} [Fintype n]

section Thm42

variable (F : Matrix n n ℝ)

def Cmat [DecidableEq n] : Matrix n n ℝ := 1 - Fᵀ

variable {F} [DecidableEq n]

/-- Columns of C sum to zero (Thm 3.3 part 1; restated standalone). -/
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
    EVERY total radiant power vector j (eq. 37):
    Σₖ (C *ᵥ j)ₖ = Σᵢ (Σₖ Cₖᵢ)·jᵢ = 0. -/
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

end Thm42

#check @col_sum_zero
#check @thm42_energy_conservation
