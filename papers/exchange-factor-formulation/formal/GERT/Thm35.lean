import GERT.MinPrinciple

/-!
# Theorem 3.5: non-singularity of mixed boundaries

Bielefeld (2026), arXiv:2512.22157v3, Section 3.2.

Let F be non-negative, row stochastic and irreducible, and let M be assembled
with the rows in S taken from D and the remaining rows from C. If 0 ≤ bₖ ≤ 1
for all k ∈ S, and bₖ₀ < 1 for at least one k₀ ∈ S, then M is non-singular.

The proof is the paper's: Mᵀx = 0 is eq. 25 with c = 0, so Lemma 3.4 gives
x ≥ 0; since −x satisfies the same equation, also −x ≥ 0, hence x = 0.

* `thm35_transpose_null_trivial` — Mᵀx = 0 ⟹ x = 0
* `thm35_isUnit`                 — M is invertible
* `Dmat_isUnit_of_irreducible`   — the remark after the theorem: with S = all
  elements, M = D, so a single element with bₖ₀ < 1 suffices for D to be
  invertible (sharpening Theorem 3.2 for irreducible F)

The hypothesis ∃ k₀ ∈ S with bₖ₀ < 1 cannot be dropped: with b ≡ 1 on S,
M = C, which is singular. No reciprocity and no spectral input is used.
-/

open Matrix Finset

namespace GERT

variable {n : Type*} [Fintype n] [DecidableEq n] {F : Matrix n n ℝ} {b : n → ℝ}

/-- **Theorem 3.5, core statement**: the transpose system Mᵀx = 0 has only
    the trivial solution. -/
theorem thm35_transpose_null_trivial {S : Finset n}
    (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : ∀ i, ∑ j, F i j = 1)
    (hirr : StronglyConnected F)
    (hb0 : ∀ k ∈ S, 0 ≤ b k) (hb1 : ∀ k ∈ S, b k ≤ 1)
    {k₀ : n} (hk₀S : k₀ ∈ S) (hk₀ : b k₀ < 1) :
    ∀ x : n → ℝ, (MmatS F b S)ᵀ *ᵥ x = 0 → x = 0 := by
  intro x hx
  -- the mask b̃ = bmask b S satisfies the hypotheses of Lemma 3.4
  have hg0 : ∀ k, 0 ≤ bmask b S k := bmask_nonneg hb0
  have hg1 : ∀ k, bmask b S k ≤ 1 := bmask_le_one hb1
  have hgk₀ : bmask b S k₀ < 1 := by
    rw [bmask_of_mem hk₀S]
    exact hk₀
  -- x satisfies eq. 25 with c = 0
  have hfix0 : ∀ k, x k = ∑ l, F k l * (bmask b S l * x l) := by
    intro k
    have h := MmatS_fixed_of_transpose_mulVec hx k
    rw [Pi.zero_apply, add_zero] at h
    exact h
  have hfixc : ∀ k, x k = ∑ l, F k l * (bmask b S l * x l) + (fun _ => (0:ℝ)) k := by
    intro k
    change x k = ∑ l, F k l * (bmask b S l * x l) + 0
    rw [add_zero]
    exact hfix0 k
  -- and so does −x
  have hfixnegc : ∀ k, (-x) k = ∑ l, F k l * (bmask b S l * (-x) l) + (fun _ => (0:ℝ)) k := by
    intro k
    change -x k = ∑ l, F k l * (bmask b S l * -x l) + 0
    rw [add_zero, hfix0 k, ← Finset.sum_neg_distrib]
    exact Finset.sum_congr rfl fun l _ => by ring
  -- Lemma 3.4 twice
  have hpos := masked_min_principle hF0 hF1 hirr hg0 hg1 hgk₀
    (c := fun _ => (0:ℝ)) (fun _ => le_refl 0) hfixc
  have hneg := masked_min_principle hF0 hF1 hirr hg0 hg1 hgk₀
    (c := fun _ => (0:ℝ)) (fun _ => le_refl 0) hfixnegc
  funext k
  change x k = 0
  have h1 := hpos k
  have h2 := hneg k
  rw [Pi.neg_apply] at h2
  exact le_antisymm (by linarith) h1

/-- **Theorem 3.5 (Non-singularity of Mixed Boundaries)**: M is invertible. -/
theorem thm35_isUnit {S : Finset n}
    (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : ∀ i, ∑ j, F i j = 1)
    (hirr : StronglyConnected F)
    (hb0 : ∀ k ∈ S, 0 ≤ b k) (hb1 : ∀ k ∈ S, b k ≤ 1)
    {k₀ : n} (hk₀S : k₀ ∈ S) (hk₀ : b k₀ < 1) :
    IsUnit (MmatS F b S) := by
  rw [Matrix.isUnit_iff_isUnit_det, isUnit_iff_ne_zero]
  intro hdet
  rw [← Matrix.det_transpose] at hdet
  obtain ⟨v, hv, hv0⟩ := Matrix.exists_mulVec_eq_zero_iff.mpr hdet
  exact hv (thm35_transpose_null_trivial hF0 hF1 hirr hb0 hb1 hk₀S hk₀ v hv0)

/-- **Remark after Theorem 3.5**: when all emissive powers are prescribed,
    M = D, so for irreducible F a single element with bₖ₀ < 1 suffices for
    D to be invertible. This sharpens Theorem 3.2, which needs bmax < 1. -/
theorem Dmat_isUnit_of_irreducible
    (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : ∀ i, ∑ j, F i j = 1)
    (hirr : StronglyConnected F)
    (hb0 : ∀ k, 0 ≤ b k) (hb1 : ∀ k, b k ≤ 1)
    {k₀ : n} (hk₀ : b k₀ < 1) :
    IsUnit (Dmat F b) := by
  rw [← MmatS_univ]
  exact thm35_isUnit hF0 hF1 hirr
    (fun k _ => hb0 k) (fun k _ => hb1 k) (Finset.mem_univ k₀) hk₀

end GERT
