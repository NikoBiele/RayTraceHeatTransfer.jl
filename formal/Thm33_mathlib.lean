/-
  Theorem 3.3, Bielefeld (2026), arXiv:2512.22157v2 — Mathlib version, rev 2.
  "M-matrix Property of C"        C = I − Aᵀ − Rᵀ

  rev 2 changes:
  * `Matrix.smul_mulVec_assoc` does not exist on this Mathlib → replaced by a
    local helper `smul_one_mulVec : (c • 1) *ᵥ u = c • u`, proven entrywise.
  * fixed the `⊢ 1 ≤ 1` residue in the row-sum bound (calc proved `=`, goal
    was `≤` → now `le_of_eq`).
  * [Fintype n] rescoped: the entrywise results need only [DecidableEq n].

  Contents beyond the verified core file:
  * `thm33_diag_eq_erase`   — eq. 28 literally: Cᵢᵢ = Σ_{j≠i} Fᵢⱼ.
  * `thm33_det_zero`        — C is singular (0 is an eigenvalue; witness 𝟙ᵀ).
  * `thm33_eigen_re_nonneg` — every complex eigenvalue μ of C has 0 ≤ Re μ.
-/
import Mathlib

open Matrix Finset

variable {n : Type*}

/-! ### General eigenvalue bound (shared with Thm32; files are standalone) -/

section GeneralLemma
variable [Fintype n] [Nonempty n] {𝕜 : Type*} [NormedField 𝕜]

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

/-! ### Theorem 3.3 -/

section Thm33

variable (F : Matrix n n ℝ) (b : n → ℝ)

def Amat : Matrix n n ℝ := Matrix.of fun i j => F i j * (1 - b j)
def Rmat : Matrix n n ℝ := Matrix.of fun i j => F i j * b j

/-- C = I − Aᵀ − Rᵀ exactly as in the paper (eq. 11).
    Here Mathlib's `(1 : Matrix n n ℝ)` IS the correct object (the identity). -/
def Cmat [DecidableEq n] : Matrix n n ℝ := 1 - (Amat F b)ᵀ - (Rmat F b)ᵀ

variable {F b} [DecidableEq n]

/-- **C = I − Fᵀ** (via A + R = F, Thm 3.1). Entrywise — no finiteness. -/
theorem thm33_C_eq : Cmat F b = 1 - Fᵀ := by
  have hAR : Amat F b + Rmat F b = F := by
    ext i j
    simp only [Matrix.add_apply, Amat, Rmat, Matrix.of_apply]
    ring
  calc Cmat F b = 1 - ((Amat F b)ᵀ + (Rmat F b)ᵀ) := by
        rw [Cmat, sub_sub]
    _ = 1 - (Amat F b + Rmat F b)ᵀ := by rw [Matrix.transpose_add]
    _ = 1 - Fᵀ := by rw [hAR]

/-- **Theorem 3.3, part 3**: off-diagonal entries are non-positive. -/
theorem thm33_offdiag_nonpos (hF0 : ∀ i j, 0 ≤ F i j)
    {i j : n} (hij : i ≠ j) : Cmat F b i j ≤ 0 := by
  rw [thm33_C_eq]
  simp only [Matrix.sub_apply, Matrix.one_apply_ne hij, Matrix.transpose_apply,
    zero_sub]
  exact neg_nonpos.mpr (hF0 j i)

/-- Diagonal entries: Cᵢᵢ = 1 − Fᵢᵢ. -/
theorem thm33_diag_eq (i : n) : Cmat F b i i = 1 - F i i := by
  rw [thm33_C_eq]
  simp only [Matrix.sub_apply, Matrix.one_apply_eq, Matrix.transpose_apply]

section RowSums
variable [Fintype n]

/-- Row-stochasticity. -/
def RowStochastic (M : Matrix n n ℝ) : Prop := ∀ i, ∑ j, M i j = 1

/-- **Eq. 28**: the diagonal entry is the off-diagonal row sum,
    Cᵢᵢ = Σ_{j≠i} Fᵢⱼ. -/
theorem thm33_diag_eq_erase (hF1 : RowStochastic F) (i : n) :
    Cmat F b i i = ∑ j ∈ Finset.univ.erase i, F i j := by
  rw [thm33_diag_eq]
  have h := Finset.add_sum_erase Finset.univ (F i) (Finset.mem_univ i)
  have h1 := hF1 i
  linarith

/-- **Theorem 3.3, part 2**: diagonal entries are non-negative. -/
theorem thm33_diag_nonneg (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : RowStochastic F)
    (i : n) : 0 ≤ Cmat F b i i := by
  rw [thm33_diag_eq_erase hF1 i]
  exact Finset.sum_nonneg fun j _ => hF0 i j

/-- **Theorem 3.3, part 1** (eq. 25): 𝟙ᵀC = 0, stated as Cᵀ *ᵥ 𝟙 = 0. -/
theorem thm33_left_eigvec (hF1 : RowStochastic F) :
    (Cmat F b)ᵀ *ᵥ (fun _ => (1:ℝ)) = 0 := by
  have hCt : (Cmat F b)ᵀ = 1 - F := by
    rw [thm33_C_eq, Matrix.transpose_sub, Matrix.transpose_one,
      Matrix.transpose_transpose]
  rw [hCt, Matrix.sub_mulVec, Matrix.one_mulVec]
  have hF : F *ᵥ (fun _ => (1:ℝ)) = fun _ => (1:ℝ) := by
    funext i
    simp only [Matrix.mulVec, dotProduct, mul_one]
    exact hF1 i
  rw [hF, sub_self]

/-- **Singularity** (completing "singular M-matrix"): det C = 0, hence 0 is an
    eigenvalue of C — witnessed by the left eigenvector 𝟙ᵀ of part 1. -/
theorem thm33_det_zero [Nonempty n] (hF1 : RowStochastic F) :
    (Cmat F b).det = 0 := by
  rw [← Matrix.det_transpose]
  refine Matrix.exists_mulVec_eq_zero_iff.mp ⟨fun _ => (1:ℝ), ?_, ?_⟩
  · intro h
    have := congrFun h (Classical.arbitrary n)
    norm_num at this
  · exact thm33_left_eigvec hF1

/-- Helper (replaces the nonexistent `Matrix.smul_mulVec_assoc`):
    a scalar multiple of the identity acts on vectors as the scalar. -/
theorem smul_one_mulVec (c : ℂ) (u : n → ℂ) :
    ((c • (1 : Matrix n n ℂ)) *ᵥ u) = c • u := by
  funext i
  simp [Matrix.mulVec, dotProduct, Matrix.smul_apply, Matrix.one_apply,
    mul_ite, ite_mul, mul_one, mul_zero, zero_mul, smul_eq_mul,
    Finset.sum_ite_eq, Finset.mem_univ]

/-- **Theorem 3.3, part 1 (spectral half)**: every complex eigenvalue μ of
    (the complexification of) C satisfies 0 ≤ Re μ.  With `thm33_det_zero`
    this makes zero an eigenvalue of smallest real part.

    Proof: C·v = μ·v ⟹ det(C − μI) = 0 ⟹ det((1−μ)I − Fᶜ) = 0 (transpose)
    ⟹ Fᶜ has eigenvalue 1−μ ⟹ ‖1−μ‖ ≤ 1 (row-sum bound) ⟹ Re μ ≥ 0. -/
theorem thm33_eigen_re_nonneg [Nonempty n]
    (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : RowStochastic F)
    {μ : ℂ} {v : n → ℂ} (hv : v ≠ 0)
    (heig : ((Cmat F b).map (↑· : ℝ → ℂ)).mulVec v = μ • v) :
    0 ≤ μ.re := by
  set Fc : Matrix n n ℂ := F.map (↑· : ℝ → ℂ) with hFc
  -- complexified C is 1 − Fcᵀ
  have hCc : (Cmat F b).map (↑· : ℝ → ℂ) = 1 - Fcᵀ := by
    rw [thm33_C_eq]
    ext i j
    by_cases h : i = j
    · subst h
      simp [Matrix.map_apply, Matrix.sub_apply, Matrix.one_apply_eq,
        Matrix.transpose_apply, hFc]
    · simp [Matrix.map_apply, Matrix.sub_apply, Matrix.one_apply_ne h,
        Matrix.transpose_apply, hFc]
  -- v kills (Cᶜ − μ·1), so its determinant vanishes
  have hdet1 : ((Cmat F b).map (↑· : ℝ → ℂ) - μ • 1).det = 0 := by
    refine Matrix.exists_mulVec_eq_zero_iff.mp ⟨v, hv, ?_⟩
    rw [Matrix.sub_mulVec, smul_one_mulVec, heig, sub_self]
  -- rewrite as ((1−μ)·1 − Fᶜ)ᵀ and drop the transpose
  have hdet2 : ((1 - μ) • (1 : Matrix n n ℂ) - Fc).det = 0 := by
    have hshape : (Cmat F b).map (↑· : ℝ → ℂ) - μ • 1
        = ((1 - μ) • (1 : Matrix n n ℂ) - Fc)ᵀ := by
      rw [hCc, Matrix.transpose_sub, Matrix.transpose_smul,
        Matrix.transpose_one]
      ext i j
      simp only [Matrix.sub_apply, Matrix.smul_apply, Matrix.transpose_apply,
        smul_eq_mul]
      ring
    rw [hshape, Matrix.det_transpose] at hdet1
    exact hdet1
  -- hence Fᶜ has eigenvalue 1 − μ
  obtain ⟨w, hw, hw0⟩ := Matrix.exists_mulVec_eq_zero_iff.mpr hdet2
  have heigF : Fc.mulVec w = (1 - μ) • w := by
    rw [Matrix.sub_mulVec, smul_one_mulVec, sub_eq_zero] at hw0
    exact hw0.symm
  -- row-sum bound on Fᶜ gives ‖1 − μ‖ ≤ 1
  have hrow : ∀ i, ∑ j, ‖Fc i j‖ ≤ (1:ℝ) := by
    intro i
    have hnorm : ∀ j, ‖Fc i j‖ = F i j := fun j => by
      simp only [hFc, Matrix.map_apply, Complex.norm_real, Real.norm_eq_abs]
      exact abs_of_nonneg (hF0 i j)
    have hsum : ∑ j, ‖Fc i j‖ = 1 := by
      calc ∑ j, ‖Fc i j‖ = ∑ j, F i j := by simp only [hnorm]
        _ = 1 := hF1 i
    exact le_of_eq hsum
  have hle : ‖(1:ℂ) - μ‖ ≤ 1 :=
    eigenvalue_norm_le_of_row_norm_sum_le hrow hw heigF
  -- Re μ = 1 − Re(1−μ) ≥ 1 − ‖1−μ‖ ≥ 0
  have hre : ((1:ℂ) - μ).re ≤ 1 := le_trans (Complex.re_le_norm _) hle
  have : ((1:ℂ) - μ).re = 1 - μ.re := by simp [Complex.sub_re, Complex.one_re]
  linarith

end RowSums

end Thm33

#check @thm33_C_eq
#check @thm33_offdiag_nonpos
#check @thm33_diag_eq_erase
#check @thm33_diag_nonneg
#check @thm33_left_eigvec
#check @thm33_det_zero
#check @thm33_eigen_re_nonneg
