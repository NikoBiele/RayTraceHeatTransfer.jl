/-
  Theorem 3.4, Bielefeld (2026), arXiv:2512.22157v2 — Mathlib version.
  "M-matrix Property of D"        D = I − Rᵀ

  ✅ STATUS: rev 2 — compiled clean on user's Mathlib on first attempt; this
  revision only adds `omit [DecidableEq n]` on the row-bound lemma to silence
  the unused-instance linter. Companion to the verified core Thm34.lean.

  Contents:
  * `thm34_D_eq_C_add_At`   — D = C + Aᵀ (the paper's Step-2 identity, eq. 30)
  * `thm34_offdiag_nonpos`  — part 3: Dᵢⱼ ≤ 0 for i ≠ j
  * `thm34_diag_eq`         — Dᵢᵢ = 1 − Fᵢᵢbᵢ
  * `thm34_diag_nonneg`     — part 2: 0 ≤ Dᵢᵢ  (needs only b ≤ 1 and F ≥ 0)
  * `thm34_eigen_re_bound`  — part 1: every complex eigenvalue μ of D has
                              1 − bmax ≤ Re μ  (the paper's Re λ_D ≥ 1 − ρ(R)
                              through ρ(R) ≤ bmax of Thm 3.2)
  Non-singularity for bmax < 1 is `thm32_D_isUnit` (Thm32_mathlib.lean, verified).
-/
import Mathlib

open Matrix Finset

variable {n : Type*}

/-! ### General eigenvalue bound (shared machinery; files are standalone) -/

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

/-! ### Theorem 3.4 -/

section Thm34

variable (F : Matrix n n ℝ) (b : n → ℝ)

def Amat : Matrix n n ℝ := Matrix.of fun i j => F i j * (1 - b j)
def Rmat : Matrix n n ℝ := Matrix.of fun i j => F i j * b j

def Cmat [DecidableEq n] : Matrix n n ℝ := 1 - Fᵀ
def Dmat [DecidableEq n] : Matrix n n ℝ := 1 - (Rmat F b)ᵀ

variable {F b} [DecidableEq n]

/-- **D = C + Aᵀ** (the paper's Step-2 identity, eq. 30). Entrywise algebra
    with `(1 : Matrix) i j` as a shared atom — no case split needed. -/
theorem thm34_D_eq_C_add_At : Dmat F b = Cmat F + (Amat F b)ᵀ := by
  ext i j
  simp only [Dmat, Cmat, Amat, Rmat, Matrix.sub_apply, Matrix.add_apply,
    Matrix.transpose_apply, Matrix.of_apply]
  ring

/-- **Theorem 3.4, part 3**: off-diagonal entries are non-positive. -/
theorem thm34_offdiag_nonpos (hF0 : ∀ i j, 0 ≤ F i j) (hb0 : ∀ j, 0 ≤ b j)
    {i j : n} (hij : i ≠ j) : Dmat F b i j ≤ 0 := by
  simp only [Dmat, Matrix.sub_apply, Matrix.one_apply_ne hij,
    Matrix.transpose_apply, Rmat, Matrix.of_apply, zero_sub]
  exact neg_nonpos.mpr (mul_nonneg (hF0 j i) (hb0 i))

/-- Diagonal entries: Dᵢᵢ = 1 − Fᵢᵢ·bᵢ. -/
theorem thm34_diag_eq (i : n) : Dmat F b i i = 1 - F i i * b i := by
  simp only [Dmat, Matrix.sub_apply, Matrix.one_apply_eq,
    Matrix.transpose_apply, Rmat, Matrix.of_apply]

section RowSums
variable [Fintype n]

def RowStochastic (M : Matrix n n ℝ) : Prop := ∀ i, ∑ j, M i j = 1

/-- **Theorem 3.4, part 2**: diagonal entries are non-negative.
    NB: needs only b ≤ 1 (with F ≥ 0 and row sums 1), NOT 0 ≤ b —
    a hypothesis the paper's D = C + Aᵀ route also implicitly relies on. -/
theorem thm34_diag_nonneg (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : RowStochastic F)
    (hb1 : ∀ j, b j ≤ 1) (i : n) : 0 ≤ Dmat F b i i := by
  rw [thm34_diag_eq]
  have h1 : F i i * b i ≤ F i i := by
    calc F i i * b i ≤ F i i * 1 :=
          mul_le_mul_of_nonneg_left (hb1 i) (hF0 i i)
      _ = F i i := mul_one _
  have h2 : F i i ≤ 1 := by
    have := Finset.single_le_sum (f := F i)
      (fun j _ => hF0 i j) (Finset.mem_univ i)
    rwa [hF1 i] at this
  linarith

omit [DecidableEq n] in
/-- Row abs-sums of R bounded by bmax (from Thm 3.2; restated for standalone use). -/
theorem Rmat_row_norm_sum_le {bmax : ℝ}
    (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : RowStochastic F)
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

/-- Helper: a scalar multiple of the identity acts on vectors as the scalar. -/
theorem smul_one_mulVec (c : ℂ) (u : n → ℂ) :
    ((c • (1 : Matrix n n ℂ)) *ᵥ u) = c • u := by
  funext i
  simp [Matrix.mulVec, dotProduct, Matrix.smul_apply, Matrix.one_apply,
    mul_ite, ite_mul, mul_one, mul_zero, zero_mul, smul_eq_mul,
    Finset.sum_ite_eq, Finset.mem_univ]

/-- **Theorem 3.4, part 1**: every complex eigenvalue μ of (the
    complexification of) D satisfies 1 − bmax ≤ Re μ — the paper's
    Re λ_D ≥ 1 − ρ(R), through ρ(R) ≤ bmax (Thm 3.2).

    Proof: D·v = μ·v ⟹ det(D − μI) = 0 ⟹ det((1−μ)I − Rᶜ) = 0 (transpose)
    ⟹ Rᶜ has eigenvalue 1−μ ⟹ ‖1−μ‖ ≤ bmax ⟹ Re μ ≥ 1 − bmax. -/
theorem thm34_eigen_re_bound [Nonempty n] {bmax : ℝ}
    (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : RowStochastic F)
    (hb0 : ∀ j, 0 ≤ b j) (hbm : ∀ j, b j ≤ bmax)
    {μ : ℂ} {v : n → ℂ} (hv : v ≠ 0)
    (heig : ((Dmat F b).map (↑· : ℝ → ℂ)).mulVec v = μ • v) :
    1 - bmax ≤ μ.re := by
  set Rc : Matrix n n ℂ := (Rmat F b).map (↑· : ℝ → ℂ) with hRc
  -- complexified D is 1 − Rcᵀ
  have hDc : (Dmat F b).map (↑· : ℝ → ℂ) = 1 - Rcᵀ := by
    ext i j
    by_cases h : i = j
    · subst h
      simp [Dmat, Matrix.map_apply, Matrix.sub_apply, Matrix.one_apply_eq,
        Matrix.transpose_apply, hRc]
    · simp [Dmat, Matrix.map_apply, Matrix.sub_apply, Matrix.one_apply_ne h,
        Matrix.transpose_apply, hRc]
  -- v kills (Dᶜ − μ·1), so its determinant vanishes
  have hdet1 : ((Dmat F b).map (↑· : ℝ → ℂ) - μ • 1).det = 0 := by
    refine Matrix.exists_mulVec_eq_zero_iff.mp ⟨v, hv, ?_⟩
    rw [Matrix.sub_mulVec, smul_one_mulVec, heig, sub_self]
  -- rewrite as ((1−μ)·1 − Rᶜ)ᵀ and drop the transpose
  have hdet2 : ((1 - μ) • (1 : Matrix n n ℂ) - Rc).det = 0 := by
    have hshape : (Dmat F b).map (↑· : ℝ → ℂ) - μ • 1
        = ((1 - μ) • (1 : Matrix n n ℂ) - Rc)ᵀ := by
      rw [hDc, Matrix.transpose_sub, Matrix.transpose_smul,
        Matrix.transpose_one]
      ext i j
      simp only [Matrix.sub_apply, Matrix.smul_apply, Matrix.transpose_apply,
        smul_eq_mul]
      ring
    rw [hshape, Matrix.det_transpose] at hdet1
    exact hdet1
  -- hence Rᶜ has eigenvalue 1 − μ
  obtain ⟨w, hw, hw0⟩ := Matrix.exists_mulVec_eq_zero_iff.mpr hdet2
  have heigR : Rc.mulVec w = (1 - μ) • w := by
    rw [Matrix.sub_mulVec, smul_one_mulVec, sub_eq_zero] at hw0
    exact hw0.symm
  -- row-sum bound on Rᶜ gives ‖1 − μ‖ ≤ bmax
  have hrow : ∀ i, ∑ j, ‖Rc i j‖ ≤ bmax := by
    intro i
    have hnorm : ∀ j, ‖Rc i j‖ = ‖Rmat F b i j‖ := fun j => by
      simp only [hRc, Matrix.map_apply, Complex.norm_real]
    calc ∑ j, ‖Rc i j‖ = ∑ j, ‖Rmat F b i j‖ := by simp only [hnorm]
      _ ≤ bmax := Rmat_row_norm_sum_le hF0 hF1 hb0 hbm i
  have hle : ‖(1:ℂ) - μ‖ ≤ bmax :=
    eigenvalue_norm_le_of_row_norm_sum_le hrow hw heigR
  -- Re μ = 1 − Re(1−μ) ≥ 1 − ‖1−μ‖ ≥ 1 − bmax
  have hre : ((1:ℂ) - μ).re ≤ bmax := le_trans (Complex.re_le_norm _) hle
  have : ((1:ℂ) - μ).re = 1 - μ.re := by simp [Complex.sub_re, Complex.one_re]
  linarith

end RowSums

end Thm34

#check @thm34_D_eq_C_add_At
#check @thm34_offdiag_nonpos
#check @thm34_diag_nonneg
#check @thm34_eigen_re_bound
