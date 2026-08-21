/-
  Theorem 3.5, Bielefeld (2026), arXiv:2512.22157v2 — Mathlib version.
  "Non-singularity of Mixed Boundaries"  (single-row case, the paper's WLOG)

  STATUS: verified (compiled clean by the author). Superseded by
  Thm35multi_mathlib.lean for the theorem itself; retained as the record
  of the original proof's Step-3 gap and its reciprocity-based repair
  (hPF_pos_of_reciprocity).

  ── FORMALIZATION FINDING ────────────────────────────────────────────────
  The paper's Step 3 infers  "𝟙 ∉ null(Mᵀ) ⟹ nullity(Mᵀ) < nullity(Cᵀ)".
  That inference presupposes null(Mᵀ) ⊆ null(Cᵀ), which is not established.
  Formally, x ∈ null(Mᵀ) gives  Cᵀx = −xⱼ·A[:,j],  and the proof splits:
   * xⱼ = 0 :  Cᵀx = 0 → x ∈ span 𝟙 → x = 0.   (the paper's argument)
   * xⱼ ≠ 0 :  needs an EXTRA Perron–Frobenius ingredient: the right null
     vector w of C (Fᵀw = w) is strictly positive; pairing w against
     Mᵀx = 0 yields xⱼ·(wᵀA[:,j]) = 0 with wᵀA[:,j] > 0 — contradiction.
  Both PF ingredients appear below as named hypotheses (hPF_null, hPF_pos),
  mirroring the paper's Horn & Johnson citation.
  ─────────────────────────────────────────────────────────────────────────

  Hypothesis minimality (surfaced by the formalization):
  * only column j of F needs non-negativity and a positive entry (Step 2's
    irreducibility consequence, stated directly as `hcol` this round);
  * only b j < 1 is needed — no bound on other components of b, and no
    row-stochasticity hypothesis at all (it is subsumed by hPF_null/hPF_pos).
-/
import Mathlib

open Matrix Finset

variable {n : Type*} [Fintype n] [DecidableEq n]

section Thm35

variable (F : Matrix n n ℝ) (b : n → ℝ)

def Amat : Matrix n n ℝ := Matrix.of fun i j => F i j * (1 - b j)
def Rmat : Matrix n n ℝ := Matrix.of fun i j => F i j * b j
def Cmat : Matrix n n ℝ := 1 - Fᵀ
def Dmat : Matrix n n ℝ := 1 - (Rmat F b)ᵀ

/-- The mixed-boundary matrix with row j of C replaced by row j of D
    (the paper's WLOG single-row case of eq. 12 / eq. 31). -/
def Mmat (j : n) : Matrix n n ℝ := (Cmat F).updateRow j (Dmat F b j)

variable {F b}

omit [Fintype n] in
/-- Row j of D equals row j of C plus row j of Aᵀ  (D = C + Aᵀ, Thm 3.4). -/
theorem Dmat_row_eq (j k : n) :
    Dmat F b j k = Cmat F j k + Amat F b k j := by
  simp only [Dmat, Cmat, Amat, Rmat, Matrix.sub_apply, Matrix.transpose_apply,
    Matrix.of_apply]
  ring

/-- The key entrywise computation:  (Mᵀx)ₖ = (Cᵀx)ₖ + Aₖⱼ·xⱼ
    (the paper's eq. 31, M = C + eⱼ(Aᵀ)ⱼ, read through the transpose). -/
theorem Mmat_transpose_mulVec (j : n) (x : n → ℝ) (k : n) :
    ((Mmat F b j)ᵀ *ᵥ x) k = ((Cmat F)ᵀ *ᵥ x) k + Amat F b k j * x j := by
  simp only [Matrix.mulVec, dotProduct, Matrix.transpose_apply]
  rw [← Finset.add_sum_erase Finset.univ
        (fun i => Mmat F b j i k * x i) (Finset.mem_univ j),
      ← Finset.add_sum_erase Finset.univ
        (fun i => Cmat F i k * x i) (Finset.mem_univ j)]
  have h1 : ∑ i ∈ Finset.univ.erase j, Mmat F b j i k * x i
          = ∑ i ∈ Finset.univ.erase j, Cmat F i k * x i := by
    refine Finset.sum_congr rfl fun i hi => ?_
    -- [RISK] `Finset.ne_of_mem_erase : a ∈ s.erase b → a ≠ b`
    have hij : i ≠ j := Finset.ne_of_mem_erase hi
    rw [Mmat, Matrix.updateRow_ne hij]
  have h2 : Mmat F b j j k = Cmat F j k + Amat F b k j := by
    rw [Mmat, Matrix.updateRow_self]
    exact Dmat_row_eq j k
  rw [h1, h2]
  ring

/-- **Theorem 3.5 (single-row case): M has trivial left null space.**

    PF inputs as named hypotheses (both from the paper's Horn & Johnson
    citation, applied to the irreducible row-stochastic F):
    * `hPF_null` : null(Cᵀ) is contained in the constants (geometric
      multiplicity 1 of eigenvalue 1 of F, with eigenvector 𝟙);
    * `hPF_pos`  : C has a strictly positive right null vector
      (the Perron vector of Fᵀ) — the ingredient the paper's Step 3
      does not invoke but the xⱼ ≠ 0 case requires. -/
theorem thm35_left_null_trivial (j : n)
    (hF0 : ∀ k, 0 ≤ F k j) (hbj : b j < 1) (hcol : ∃ i, 0 < F i j)
    (hPF_null : ∀ y : n → ℝ, (Cmat F)ᵀ *ᵥ y = 0 → ∃ c, y = fun _ => c)
    (hPF_pos : ∃ w : n → ℝ, (∀ i, 0 < w i) ∧ Cmat F *ᵥ w = 0) :
    ∀ x : n → ℝ, (Mmat F b j)ᵀ *ᵥ x = 0 → x = 0 := by
  intro x hx
  have hA : ∀ k, Amat F b k j = F k j * (1 - b j) := fun k => rfl
  by_cases hxj : x j = 0
  · -- the paper's case:  Cᵀx = 0 → x constant → x = 0
    have hCx : (Cmat F)ᵀ *ᵥ x = 0 := by
      funext k
      have hk := congrFun hx k
      rw [Mmat_transpose_mulVec] at hk
      simpa [hxj] using hk
    obtain ⟨c, hc⟩ := hPF_null x hCx
    have hc0 : c = 0 := by rw [hc] at hxj; exact hxj
    funext k
    rw [hc, hc0]
    rfl
  · -- the missing case:  pair against the positive null vector w
    exfalso
    obtain ⟨w, hwpos, hwC⟩ := hPF_pos
    have hsum0 : ∑ k, w k * (((Mmat F b j)ᵀ *ᵥ x) k) = 0 := by
      simp [hx]
    have hexp : ∑ k, w k * (((Mmat F b j)ᵀ *ᵥ x) k)
        = (∑ k, w k * (((Cmat F)ᵀ *ᵥ x) k))
          + (∑ k, w k * Amat F b k j) * x j := by
      calc ∑ k, w k * (((Mmat F b j)ᵀ *ᵥ x) k)
          = ∑ k, (w k * (((Cmat F)ᵀ *ᵥ x) k) + w k * Amat F b k j * x j) := by
            refine Finset.sum_congr rfl fun k _ => ?_
            rw [Mmat_transpose_mulVec]
            ring
        -- [RISK] `Finset.sum_add_distrib : ∑ (f + g) = ∑ f + ∑ g`
        _ = (∑ k, w k * (((Cmat F)ᵀ *ᵥ x) k))
              + ∑ k, w k * Amat F b k j * x j := Finset.sum_add_distrib
        _ = (∑ k, w k * (((Cmat F)ᵀ *ᵥ x) k))
              + (∑ k, w k * Amat F b k j) * x j := by rw [Finset.sum_mul]
    have hpair : ∑ k, w k * (((Cmat F)ᵀ *ᵥ x) k) = 0 := by
      have h := Matrix.dotProduct_transpose_mulVec (Cmat F) w x
      rw [hwC] at h
      simpa [dotProduct] using h
    have hpos : 0 < ∑ k, w k * Amat F b k j := by
      have h1mb : 0 < 1 - b j := by linarith
      refine Finset.sum_pos' (fun k _ => ?_) ?_
      · rw [hA k]
        exact mul_nonneg (le_of_lt (hwpos k))
          (mul_nonneg (hF0 k) (le_of_lt h1mb))
      · obtain ⟨i, hi⟩ := hcol
        refine ⟨i, Finset.mem_univ i, ?_⟩
        rw [hA i]
        exact mul_pos (hwpos i) (mul_pos hi h1mb)
    rw [hexp, hpair, zero_add] at hsum0
    rcases mul_eq_zero.mp hsum0 with h | h
    · exact absurd h (ne_of_gt hpos)
    · exact hxj h

/-- **Theorem 3.5 (single-row case): M is non-singular.** -/
theorem thm35_isUnit (j : n)
    (hF0 : ∀ k, 0 ≤ F k j) (hbj : b j < 1) (hcol : ∃ i, 0 < F i j)
    (hPF_null : ∀ y : n → ℝ, (Cmat F)ᵀ *ᵥ y = 0 → ∃ c, y = fun _ => c)
    (hPF_pos : ∃ w : n → ℝ, (∀ i, 0 < w i) ∧ Cmat F *ᵥ w = 0) :
    IsUnit (Mmat F b j) := by
  rw [Matrix.isUnit_iff_isUnit_det, isUnit_iff_ne_zero]
  intro hdet
  rw [← Matrix.det_transpose] at hdet
  obtain ⟨v, hv, hv0⟩ := Matrix.exists_mulVec_eq_zero_iff.mpr hdet
  exact hv (thm35_left_null_trivial j hF0 hbj hcol hPF_null hPF_pos v hv0)

/-- **Reciprocity discharges `hPF_pos` with an explicit witness** (user's
    observation): if F satisfies reciprocity (paper 1, eq. 7; w-notation per
    paper 2) with strictly positive equivalent emission areas w, and F is row
    stochastic, then w itself is a strictly positive right null vector of C —
    no Perron–Frobenius needed.  (Fᵀw)ᵢ = Σⱼ Fⱼᵢwⱼ = Σⱼ wᵢFᵢⱼ = wᵢ.
    Physically: the radiative-equilibrium radiant power field is proportional
    to the emission capacities.  Paper 2's smoothing enforces exactly these
    two hypotheses, so every smoothed F̂ satisfies them by construction. -/
theorem hPF_pos_of_reciprocity (w : n → ℝ) (hw : ∀ i, 0 < w i)
    (hrec : ∀ i j, w i * F i j = F j i * w j)
    (hF1 : ∀ i, ∑ j, F i j = 1) :
    ∃ v : n → ℝ, (∀ i, 0 < v i) ∧ Cmat F *ᵥ v = 0 := by
  refine ⟨w, hw, ?_⟩
  funext i
  simp only [Cmat, Matrix.sub_mulVec, Matrix.one_mulVec, Pi.sub_apply,
    Pi.zero_apply]
  have hFt : ((Fᵀ) *ᵥ w) i = w i := by
    calc ((Fᵀ) *ᵥ w) i = ∑ j, F j i * w j := by
          simp only [Matrix.mulVec, dotProduct, Matrix.transpose_apply]
      _ = ∑ j, w i * F i j := by
          exact Finset.sum_congr rfl fun j _ => (hrec i j).symm
      _ = w i * ∑ j, F i j := by rw [Finset.mul_sum]
      _ = w i := by rw [hF1 i, mul_one]
  rw [hFt, sub_self]

end Thm35

#check @thm35_left_null_trivial
#check @thm35_isUnit
#check @hPF_pos_of_reciprocity
