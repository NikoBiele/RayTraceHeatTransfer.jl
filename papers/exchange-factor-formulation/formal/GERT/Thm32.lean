import GERT.Defs

/-!
# Theorem 3.2: upper bound on the spectral radius of R

Bielefeld (2026), arXiv:2512.22157v3, Section 3.1.

* `eigenvalue_norm_le_of_row_norm_sum_le` — general lemma (any normed field):
  rows satisfy Σⱼ ‖Mᵢⱼ‖ ≤ c ⟹ every eigenvalue μ has ‖μ‖ ≤ c
  (the max-coordinate argument behind ρ(M) ≤ ‖M‖∞)
* `Rmat_row_norm_sum_le` — rows of R = F ∘ B satisfy Σⱼ |Rᵢⱼ| ≤ bmax (eq. 21)
* `thm32_eigenvalue_bound` — part (a): complex eigenvalues of R have ‖μ‖ ≤ bmax
* `thm32_uniform_attained` — part (b): b ≡ bmax ⟹ R 𝟙 = bmax 𝟙 (bound attained)
* `thm32_D_isUnit` — part (c): bmax < 1 ⟹ D = I − Rᵀ invertible

Theorem 3.5 (`Dmat_isUnit_of_irreducible`) sharpens part (c) for irreducible
F: a single bₖ₀ < 1 then suffices.
-/

open Matrix Finset

namespace GERT

variable {n : Type*} [Fintype n]

section GeneralLemma

variable [Nonempty n] {𝕜 : Type*} [NormedField 𝕜]

/-- If every row of `M` satisfies `Σⱼ ‖M i j‖ ≤ c`, then any eigenvalue `μ`
    of `M` satisfies `‖μ‖ ≤ c`. -/
theorem eigenvalue_norm_le_of_row_norm_sum_le
    {M : Matrix n n 𝕜} {c : ℝ} (hM : ∀ i, ∑ j, ‖M i j‖ ≤ c)
    {μ : 𝕜} {v : n → 𝕜} (hv : v ≠ 0) (heig : M.mulVec v = μ • v) :
    ‖μ‖ ≤ c := by
  obtain ⟨i, -, hi⟩ :=
    Finset.exists_max_image (Finset.univ : Finset n) (fun k => ‖v k‖)
      ⟨Classical.arbitrary n, Finset.mem_univ _⟩
  have hvi : 0 < ‖v i‖ := by
    obtain ⟨k, hk⟩ := Function.ne_iff.mp hv
    calc (0:ℝ) < ‖v k‖ := norm_pos_iff.mpr hk
      _ ≤ ‖v i‖ := hi k (Finset.mem_univ k)
  have hrow : ∑ j, M i j * v j = μ * v i := by
    have h := congrFun heig i
    simpa [Matrix.mulVec, dotProduct, Pi.smul_apply, smul_eq_mul] using h
  have hbound : ‖μ‖ * ‖v i‖ ≤ c * ‖v i‖ := by
    calc ‖μ‖ * ‖v i‖ = ‖μ * v i‖ := (norm_mul μ (v i)).symm
      _ = ‖∑ j, M i j * v j‖ := by rw [hrow]
      _ ≤ ∑ j, ‖M i j * v j‖ := norm_sum_le _ _
      _ = ∑ j, ‖M i j‖ * ‖v j‖ := by simp [norm_mul]
      _ ≤ ∑ j, ‖M i j‖ * ‖v i‖ :=
          Finset.sum_le_sum fun j _ =>
            mul_le_mul_of_nonneg_left (hi j (Finset.mem_univ j)) (norm_nonneg _)
      _ = (∑ j, ‖M i j‖) * ‖v i‖ := (Finset.sum_mul ..).symm
      _ ≤ c * ‖v i‖ := mul_le_mul_of_nonneg_right (hM i) (norm_nonneg _)
  exact le_of_mul_le_mul_right hbound hvi

end GeneralLemma

section Thm32

variable {F : Matrix n n ℝ} {b : n → ℝ} {bmax : ℝ}

/-- Row sums of |R| are bounded by `bmax` (eq. 21): Σⱼ bⱼ Fᵢⱼ ≤ bmax Σⱼ Fᵢⱼ = bmax. -/
theorem Rmat_row_norm_sum_le
    (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : ∀ i, ∑ j, F i j = 1)
    (hb0 : ∀ j, 0 ≤ b j) (hbm : ∀ j, b j ≤ bmax) (i : n) :
    ∑ j, ‖Rmat F b i j‖ ≤ bmax := by
  have habs : ∀ j, ‖Rmat F b i j‖ = F i j * b j := fun j => by
    simp only [Rmat, Matrix.of_apply, Real.norm_eq_abs]
    exact abs_of_nonneg (mul_nonneg (hF0 i j) (hb0 j))
  calc ∑ j, ‖Rmat F b i j‖ = ∑ j, F i j * b j := by simp only [habs]
    _ ≤ ∑ j, F i j * bmax :=
        Finset.sum_le_sum fun j _ => mul_le_mul_of_nonneg_left (hbm j) (hF0 i j)
    _ = (∑ j, F i j) * bmax := (Finset.sum_mul ..).symm
    _ = bmax := by rw [hF1 i, one_mul]

/-- **Theorem 3.2, part (a)**: every complex eigenvalue μ of (the
    complexification of) R satisfies ‖μ‖ ≤ bmax, i.e. ρ(R) ≤ bmax. -/
theorem thm32_eigenvalue_bound [Nonempty n]
    (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : ∀ i, ∑ j, F i j = 1)
    (hb0 : ∀ j, 0 ≤ b j) (hbm : ∀ j, b j ≤ bmax)
    {μ : ℂ} {v : n → ℂ} (hv : v ≠ 0)
    (heig : ((Rmat F b).map (↑· : ℝ → ℂ)).mulVec v = μ • v) :
    ‖μ‖ ≤ bmax := by
  refine eigenvalue_norm_le_of_row_norm_sum_le (fun i => ?_) hv heig
  have hmap : ∀ j, ‖((Rmat F b).map (↑· : ℝ → ℂ)) i j‖ = ‖Rmat F b i j‖ := fun j => by
    simp [Matrix.map_apply, Complex.norm_real]
  simpa [hmap] using Rmat_row_norm_sum_le hF0 hF1 hb0 hbm i

/-- **Theorem 3.2, part (b)**: with uniform coefficients `b ≡ bmax` the bound
    is attained: 𝟙 is an eigenvector of R with eigenvalue bmax. -/
theorem thm32_uniform_attained
    (hF1 : ∀ i, ∑ j, F i j = 1) (hb : ∀ j, b j = bmax) :
    (Rmat F b).mulVec (fun _ => 1) = bmax • (fun _ => (1:ℝ)) := by
  funext i
  simp only [Matrix.mulVec, dotProduct, Rmat, Matrix.of_apply, hb,
    mul_one, Pi.smul_apply, smul_eq_mul]
  rw [← Finset.sum_mul, hF1 i, one_mul]

/-- **Theorem 3.2, part (c)**: if bmax < 1 then D = I − Rᵀ is invertible.
    Proof: det(I − Rᵀ) = det(I − R); if it vanished, some v ≠ 0 would satisfy
    R v = v, so 1 would be an eigenvalue of R, contradicting ‖μ‖ ≤ bmax < 1.
    The transpose step matters: the ROW sums of Rᵀ are the COLUMN sums of R,
    which are not bounded by bmax; eigenvalue invariance under transposition
    carries the argument. -/
theorem thm32_D_isUnit [DecidableEq n] [Nonempty n]
    (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : ∀ i, ∑ j, F i j = 1)
    (hb0 : ∀ j, 0 ≤ b j) (hbm : ∀ j, b j ≤ bmax) (hlt : bmax < 1) :
    IsUnit (Dmat F b) := by
  unfold Dmat
  rw [Matrix.isUnit_iff_isUnit_det, isUnit_iff_ne_zero]
  intro hdet
  have htr : (1 - (Rmat F b)ᵀ) = (1 - Rmat F b)ᵀ := by
    rw [Matrix.transpose_sub, Matrix.transpose_one]
  rw [htr, Matrix.det_transpose] at hdet
  obtain ⟨v, hvne, hv0⟩ := (Matrix.exists_mulVec_eq_zero_iff).mpr hdet
  have heig : (Rmat F b).mulVec v = (1:ℝ) • v := by
    rw [Matrix.sub_mulVec, Matrix.one_mulVec, sub_eq_zero] at hv0
    simpa using hv0.symm
  have hle : ‖(1:ℝ)‖ ≤ bmax :=
    eigenvalue_norm_le_of_row_norm_sum_le
      (Rmat_row_norm_sum_le hF0 hF1 hb0 hbm) hvne heig
  rw [norm_one] at hle
  linarith

end Thm32

end GERT
