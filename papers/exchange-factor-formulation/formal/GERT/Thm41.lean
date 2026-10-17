import GERT.Thm35

/-!
# Theorem 4.1: non-negative radiation

Bielefeld (2026), arXiv:2512.22157v3, Section 4.1.

Let F be non-negative, row stochastic and irreducible, let 0 ≤ bₖ ≤ 1 for
all k, and let M be assembled as in Theorem 3.5, with bₖ₀ < 1 for at least
one k₀ ∈ S. Then M⁻¹ ≥ 0 entrywise, and for any h ≥ 0 the unique solution
of Mj = h satisfies j ≥ 0, r = Rᵀj ≥ 0 and e = Dj = q + Aᵀj ≥ 0.

Part 1 is the paper's argument: for each i, the solution x of Mᵀx = eᵢ is
the i-th row of M⁻¹, and Mᵀx = eᵢ is eq. 25 with c = eᵢ ≥ 0, so Lemma 3.4
gives x ≥ 0.

* `thm41_transpose_nonneg` — Mᵀx = c with c ≥ 0 implies x ≥ 0
  (the general form; the paper applies it with c = eᵢ)
* `thm41_inv_nonneg`       — M⁻¹ ≥ 0 entrywise
* `thm41_solution_nonneg`  — Part 1: Mj = h, h ≥ 0 ⟹ j ≥ 0
* `thm41_reflected_nonneg` — Part 2: r = Rᵀj ≥ 0
* `thm41_emissive_nonneg`  — Part 3: e = q + Aᵀj ≥ 0

No reciprocity is assumed; the v2 formalization of this theorem needed it.
-/

open Matrix Finset

namespace GERT

variable {n : Type*} [Fintype n] [DecidableEq n] {F : Matrix n n ℝ} {b : n → ℝ}

/-- The unit vector eₖ. -/
def unitVec (k : n) : n → ℝ := fun l => if l = k then 1 else 0

omit [Fintype n] in
theorem unitVec_nonneg (k l : n) : 0 ≤ unitVec k l := by
  unfold unitVec
  split_ifs <;> norm_num

theorem dotProduct_unitVec (v : n → ℝ) (k : n) : v ⬝ᵥ unitVec k = v k := by
  simp [unitVec, dotProduct, mul_ite, mul_one, mul_zero,
    Finset.sum_ite_eq', Finset.mem_univ]

/-- **Theorem 4.1, the engine**: every solution of Mᵀx = c with c ≥ 0 is
    non-negative (Lemma 3.4 applied to eq. 25). -/
theorem thm41_transpose_nonneg {S : Finset n}
    (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : ∀ i, ∑ j, F i j = 1)
    (hirr : StronglyConnected F)
    (hb0 : ∀ k ∈ S, 0 ≤ b k) (hb1 : ∀ k ∈ S, b k ≤ 1)
    {k₀ : n} (hk₀S : k₀ ∈ S) (hk₀ : b k₀ < 1)
    {x c : n → ℝ} (hc : ∀ k, 0 ≤ c k)
    (hx : (MmatS F b S)ᵀ *ᵥ x = c) :
    ∀ k, 0 ≤ x k := by
  have hg0 : ∀ k, 0 ≤ bmask b S k := bmask_nonneg hb0
  have hg1 : ∀ k, bmask b S k ≤ 1 := bmask_le_one hb1
  have hgk₀ : bmask b S k₀ < 1 := by
    rw [bmask_of_mem hk₀S]
    exact hk₀
  exact masked_min_principle hF0 hF1 hirr hg0 hg1 hgk₀ hc
    (MmatS_fixed_of_transpose_mulVec hx)

/-- Mᵀ has a unit determinant, from Theorem 3.5. -/
theorem MmatS_transpose_det_isUnit {S : Finset n}
    (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : ∀ i, ∑ j, F i j = 1)
    (hirr : StronglyConnected F)
    (hb0 : ∀ k ∈ S, 0 ≤ b k) (hb1 : ∀ k ∈ S, b k ≤ 1)
    {k₀ : n} (hk₀S : k₀ ∈ S) (hk₀ : b k₀ < 1) :
    IsUnit (MmatS F b S)ᵀ.det := by
  rw [Matrix.det_transpose]
  exact (Matrix.isUnit_iff_isUnit_det _).mp
    (thm35_isUnit hF0 hF1 hirr hb0 hb1 hk₀S hk₀)

/-- **Theorem 4.1, M⁻¹ ≥ 0 entrywise**: the i-th row of M⁻¹ is the solution
    of Mᵀx = eᵢ, which is non-negative by `thm41_transpose_nonneg`. -/
theorem thm41_inv_nonneg {S : Finset n}
    (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : ∀ i, ∑ j, F i j = 1)
    (hirr : StronglyConnected F)
    (hb0 : ∀ k ∈ S, 0 ≤ b k) (hb1 : ∀ k ∈ S, b k ≤ 1)
    {k₀ : n} (hk₀S : k₀ ∈ S) (hk₀ : b k₀ < 1) (i k : n) :
    0 ≤ (MmatS F b S)⁻¹ i k := by
  have hunit := MmatS_transpose_det_isUnit hF0 hF1 hirr hb0 hb1 hk₀S hk₀
  -- x = Mᵀ⁻¹ eᵢ solves Mᵀx = eᵢ
  have hx : (MmatS F b S)ᵀ *ᵥ ((MmatS F b S)ᵀ⁻¹ *ᵥ unitVec i) = unitVec i := by
    rw [Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ hunit, Matrix.one_mulVec]
  have hxnn := thm41_transpose_nonneg hF0 hF1 hirr hb0 hb1 hk₀S hk₀
    (unitVec_nonneg i) hx
  -- and its k-th entry is (Mᵀ⁻¹)ₖᵢ = (M⁻¹)ᵢₖ
  have hxk : ((MmatS F b S)ᵀ⁻¹ *ᵥ unitVec i) k = (MmatS F b S)⁻¹ i k := by
    rw [← Matrix.transpose_nonsing_inv]
    simp [unitVec, Matrix.mulVec, dotProduct, Matrix.transpose_apply, mul_ite,
      mul_one, mul_zero, Finset.sum_ite_eq', Finset.mem_univ]
  rw [← hxk]
  exact hxnn k

/-- **Theorem 4.1, Part 1**: the solution of Mj = h is non-negative for h ≥ 0.
    Proof: with x the solution of Mᵀx = eₖ, jₖ = eₖ·j = (Mᵀx)·j = x·(Mj) = x·h ≥ 0. -/
theorem thm41_solution_nonneg {S : Finset n}
    (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : ∀ i, ∑ j, F i j = 1)
    (hirr : StronglyConnected F)
    (hb0 : ∀ k ∈ S, 0 ≤ b k) (hb1 : ∀ k ∈ S, b k ≤ 1)
    {k₀ : n} (hk₀S : k₀ ∈ S) (hk₀ : b k₀ < 1)
    {j h : n → ℝ} (hh : ∀ k, 0 ≤ h k)
    (hsol : MmatS F b S *ᵥ j = h) :
    ∀ k, 0 ≤ j k := by
  intro k
  have hunit := MmatS_transpose_det_isUnit hF0 hF1 hirr hb0 hb1 hk₀S hk₀
  obtain ⟨x, hx⟩ : ∃ x : n → ℝ, (MmatS F b S)ᵀ *ᵥ x = unitVec k :=
    ⟨(MmatS F b S)ᵀ⁻¹ *ᵥ unitVec k, by
      rw [Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ hunit, Matrix.one_mulVec]⟩
  have hxnn := thm41_transpose_nonneg hF0 hF1 hirr hb0 hb1 hk₀S hk₀
    (unitVec_nonneg k) hx
  have hpair : j ⬝ᵥ ((MmatS F b S)ᵀ *ᵥ x) = x ⬝ᵥ (MmatS F b S *ᵥ j) :=
    Matrix.dotProduct_transpose_mulVec (MmatS F b S) j x
  rw [hx, hsol, dotProduct_unitVec] at hpair
  rw [hpair]
  change 0 ≤ ∑ i, x i * h i
  exact Finset.sum_nonneg fun i _ => mul_nonneg (hxnn i) (hh i)

omit [DecidableEq n] in
/-- **Theorem 4.1, Part 2**: the reflected-scattered power r = Rᵀj is
    non-negative when j ≥ 0. -/
theorem thm41_reflected_nonneg
    (hF0 : ∀ i j, 0 ≤ F i j) (hb0 : ∀ k, 0 ≤ b k)
    {j : n → ℝ} (hj : ∀ i, 0 ≤ j i) (k : n) :
    0 ≤ ((Rmat F b)ᵀ *ᵥ j) k := by
  simp only [Matrix.mulVec, dotProduct, Matrix.transpose_apply, Rmat,
    Matrix.of_apply]
  exact Finset.sum_nonneg fun i _ =>
    mul_nonneg (mul_nonneg (hF0 i k) (hb0 k)) (hj i)

omit [DecidableEq n] in
/-- **Theorem 4.1, Part 3**: the emissive power e = q + Aᵀj is non-negative
    when q ≥ 0 and j ≥ 0 (needs b ≤ 1). -/
theorem thm41_emissive_nonneg
    (hF0 : ∀ i j, 0 ≤ F i j) (hb1 : ∀ k, b k ≤ 1)
    {q j : n → ℝ} (hq : ∀ k, 0 ≤ q k) (hj : ∀ i, 0 ≤ j i) (k : n) :
    0 ≤ q k + ((Amat F b)ᵀ *ᵥ j) k := by
  refine add_nonneg (hq k) ?_
  simp only [Matrix.mulVec, dotProduct, Matrix.transpose_apply, Amat,
    Matrix.of_apply]
  exact Finset.sum_nonneg fun i _ =>
    mul_nonneg (mul_nonneg (hF0 i k) (by linarith [hb1 k])) (hj i)

end GERT
