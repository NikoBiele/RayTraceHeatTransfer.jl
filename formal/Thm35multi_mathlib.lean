/-
  Theorem 3.5, Bielefeld (2026), arXiv:2512.22157v2 — Mathlib version.
  MULTI-ROW CASE (the theorem as actually stated: "at least one row"),
  via the masked-albedo reduction and a maximum principle.

  STATUS: rev 3. Rev-3 fixes: the true source of the hgm mismatch was the
  `rw [one_mul, ← hxj]` in hgj rewriting BOTH occurrences of m (now a
  targeted calc), and `le_of_not_lt` is renamed on this Mathlib (now the
  stable `not_lt.mp`). All three name gambles of rev 1 confirmed to exist
  (`Relation.ReflTransGen.cases_tail`, `Finset.sum_sub_distrib`,
  `Finset.sum_ite_eq'`); rev 2 fixes two rewrite-direction bugs of mine
  (exists_pos_in_column, the hxj step), removes deprecated `push_neg`,
  and adds `omit [DecidableEq n]` where the linter requested.

  ── STRUCTURE ────────────────────────────────────────────────────────────
  * `StronglyConnected`  — irreducibility, formalized as: between any two
    elements there is a chain of strictly positive entries (verbatim the
    paper's §2.1 connectivity requirement), via Relation.ReflTransGen.
  * `MmatS`              — mixed matrix for an arbitrary replacement set
    S : Finset n (rows in S from D, rows outside S from C).
  * `gmask`              — the masked albedo: gⱼ = bⱼ on S, gⱼ = 1 off S.
  * `MmatS_null_fixed`   — the reduction: Mᵀx = 0 ⟹ x is a fixed vector
    of the masked reflection operator, xₖ = Σⱼ Fₖⱼ gⱼ xⱼ.
  * `masked_no_fixed_vector` — THE MAIN LEMMA (maximum principle): an
    irreducible row-stochastic F with mask 0 ≤ g ≤ 1, g(j₀) < 1, and an
    in-edge into j₀, admits no nonzero masked fixed vector.
  * `thm35_multirow_isUnit` — Theorem 3.5 for arbitrary nonempty S with
    the (necessary!) hypothesis ∃ j₀ ∈ S, b j₀ < 1.
  * `thm32c_of_multirow` — unification check: S = univ recovers Thm 3.2(c),
    now needing only ONE b j₀ < 1 (instead of bmax < 1) at the price of
    irreducibility.

  ── FINDINGS ─────────────────────────────────────────────────────────────
  1. NO Perron–Frobenius input remains: neither hPF_null nor hPF_pos.  The
     maximum principle re-derives the nullity structure elementarily, so the
     Horn & Johnson citation becomes background.  Theorem 3.5 is
     machine-checked end-to-end from {F ≥ 0, rows = 1, strongly connected,
     0 ≤ b ≤ 1 on S, ∃ j₀ ∈ S with b j₀ < 1}.
  2. The hypothesis ∃ j₀ ∈ S, b j₀ < 1 is NECESSARY: if b ≡ 1 on S, every
     replaced D-row equals the C-row it replaces, M = C, singular.  The
     paper's statement needs this condition added (it is the abstract's
     standing "max reflection-scattering coefficient < 1", localized to S).
  3. Physical reading: gmask says prescribed elements keep their true
     albedo while free elements act as perfect mirrors; the maximum
     principle is "the radiant power field has no interior maximum unless
     it is flat, and flatness is incompatible with absorption into a
     prescribed element."
-/
import Mathlib

open Matrix Finset

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- Irreducibility of F, formalized as strong connectivity of the
    positive-entry relation: from every i there is a chain
    i → k₁ → … → j of strictly positive entries (paper §2.1). -/
def StronglyConnected (F : Matrix n n ℝ) : Prop :=
  ∀ i j, Relation.ReflTransGen (fun a c => 0 < F a c) i j

section Thm35Multi

variable (F : Matrix n n ℝ) (b : n → ℝ)

def Rmat : Matrix n n ℝ := Matrix.of fun i j => F i j * b j
def Cmat : Matrix n n ℝ := 1 - Fᵀ
def Dmat : Matrix n n ℝ := 1 - (Rmat F b)ᵀ

/-- Mixed-boundary matrix for replacement set S: rows in S from D,
    rows outside S from C (paper eq. 12, general case). -/
def MmatS (S : Finset n) : Matrix n n ℝ :=
  Matrix.of fun i j => if i ∈ S then Dmat F b i j else Cmat F i j

/-- Masked albedo: prescribed elements keep b, free elements reflect fully. -/
def gmask (S : Finset n) : n → ℝ := fun j => if j ∈ S then b j else 1

variable {F b}

/-- **The reduction**: a left null vector of M is a fixed vector of the
    masked reflection operator:  xₖ = Σⱼ Fₖⱼ · g ⱼ · xⱼ. -/
theorem MmatS_null_fixed {S : Finset n} {x : n → ℝ}
    (hx : (MmatS F b S)ᵀ *ᵥ x = 0) (k : n) :
    x k = ∑ j, F k j * (gmask b S j * x j) := by
  have hk := congrFun hx k
  simp only [Matrix.mulVec, dotProduct, Matrix.transpose_apply,
    Pi.zero_apply] at hk
  -- hk : ∑ i, MmatS F b S i k * x i = 0
  have hM : ∀ i, MmatS F b S i k
      = (1 : Matrix n n ℝ) i k - F k i * gmask b S i := by
    intro i
    by_cases hiS : i ∈ S
    · simp [MmatS, hiS, Dmat, Rmat, gmask, Matrix.sub_apply,
        Matrix.transpose_apply]
    · simp [MmatS, hiS, Cmat, gmask, Matrix.sub_apply,
        Matrix.transpose_apply]
  rw [Finset.sum_congr rfl (fun i _ => by rw [hM i, sub_mul])] at hk
  -- [RISK] `Finset.sum_sub_distrib : ∑ (f - g) = ∑ f - ∑ g`
  rw [Finset.sum_sub_distrib] at hk
  have hid : ∑ i, (1 : Matrix n n ℝ) i k * x i = x k := by
    -- `Finset.sum_ite_eq'` (bound variable on the LEFT of the ite eq; confirmed)
    simp [Matrix.one_apply, ite_mul, one_mul, zero_mul,
      Finset.sum_ite_eq', Finset.mem_univ]
  rw [hid] at hk
  rw [sub_eq_zero.mp hk]
  exact Finset.sum_congr rfl fun i _ => by ring

omit [DecidableEq n] in
/-- Every column of a strongly connected row-stochastic matrix has a
    positive entry (an in-edge): the row of j₀ has a positive entry
    F j₀ j₁ (rows summing to 1 cannot be all non-positive — note that
    non-negativity of F is NOT needed here), and the return chain
    j₁ → j₀ supplies the in-edge (its last step; or the self-edge
    when j₁ = j₀). -/
theorem exists_pos_in_column
    (hF1 : ∀ i, ∑ j, F i j = 1) (hirr : StronglyConnected F) (j₀ : n) :
    ∃ k, 0 < F k j₀ := by
  have hout : ∃ j₁, 0 < F j₀ j₁ := by
    by_contra h
    have h' : ∀ j, F j₀ j ≤ 0 :=
      fun j => not_lt.mp (fun hj => h ⟨j, hj⟩)
    have hle : ∑ j, F j₀ j ≤ 0 :=
      Finset.sum_nonpos (fun j _ => h' j)
    rw [hF1 j₀] at hle
    linarith
  obtain ⟨j₁, hj₁⟩ := hout
  -- `Relation.ReflTransGen.cases_tail` (confirmed on this Mathlib):
  --          ReflTransGen r a b → b = a ∨ ∃ c, ReflTransGen r a c ∧ r c b
  rcases (hirr j₁ j₀).cases_tail with heq | ⟨c, -, hc⟩
  · rw [← heq] at hj₁
    exact ⟨j₀, hj₁⟩
  · exact ⟨c, hc⟩

omit [DecidableEq n] in
/-- **MAIN LEMMA (maximum principle)**: an irreducible row-stochastic
    non-negative F, masked by 0 ≤ g ≤ 1 with g j₀ < 1 and an in-edge into
    j₀, admits no nonzero fixed vector x = (F∘g)·x.

    Proof: at a coordinate of maximal modulus m > 0, the fixed-point
    identity forces every F-successor j to have g j = 1 and |x j| = m; the
    max-set floods along chains (strong connectivity) to the whole domain;
    the in-neighbor of j₀ then forces g j₀ = 1 — contradiction. -/
theorem masked_no_fixed_vector
    (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : ∀ i, ∑ j, F i j = 1)
    (hirr : StronglyConnected F)
    {g : n → ℝ} (hg0 : ∀ j, 0 ≤ g j) (hg1 : ∀ j, g j ≤ 1)
    {j₀ : n} (hj₀ : g j₀ < 1) (hcol : ∃ k, 0 < F k j₀)
    {x : n → ℝ} (hfix : ∀ k, x k = ∑ j, F k j * (g j * x j)) :
    x = 0 := by
  by_contra hx0
  obtain ⟨kmax, -, hmax⟩ :=
    Finset.exists_max_image (Finset.univ : Finset n) (fun k => |x k|)
      ⟨j₀, Finset.mem_univ j₀⟩
  set m := |x kmax| with hm_def
  have hm : 0 < m := by
    obtain ⟨k0, hk0⟩ := Function.ne_iff.mp hx0
    calc (0:ℝ) < |x k0| := abs_pos.mpr hk0
      _ ≤ m := hmax k0 (Finset.mem_univ k0)
  -- Step A: at a max coordinate, every F-successor is max and has g = 1
  have step : ∀ k, |x k| = m → ∀ j, 0 < F k j → g j = 1 ∧ |x j| = m := by
    intro k hk j hFj
    have h1 : m ≤ ∑ j, F k j * (g j * |x j|) := by
      calc m = |x k| := hk.symm
        _ = |∑ j, F k j * (g j * x j)| := by rw [hfix k]
        _ ≤ ∑ j, |F k j * (g j * x j)| := Finset.abs_sum_le_sum_abs _ _
        _ = ∑ j, F k j * (g j * |x j|) := by
            refine Finset.sum_congr rfl fun j _ => ?_
            rw [abs_mul, abs_mul, abs_of_nonneg (hF0 k j),
              abs_of_nonneg (hg0 j)]
    have hub : ∀ j, F k j * (g j * |x j|) ≤ F k j * m := by
      intro j
      refine mul_le_mul_of_nonneg_left ?_ (hF0 k j)
      exact le_trans (mul_le_of_le_one_left (abs_nonneg _) (hg1 j))
        (hmax j (Finset.mem_univ j))
    have h3 : ∑ j, F k j * m = m := by
      rw [← Finset.sum_mul, hF1 k, one_mul]
    have hTm : ∑ j, F k j * (g j * |x j|) = m :=
      le_antisymm (le_of_le_of_eq (Finset.sum_le_sum fun j _ => hub j) h3) h1
    have hsub0 : ∑ j, (F k j * m - F k j * (g j * |x j|)) = 0 := by
      rw [Finset.sum_sub_distrib, h3, hTm, sub_self]
    have hzero := (Finset.sum_eq_zero_iff_of_nonneg
      (fun j _ => sub_nonneg.mpr (hub j))).mp hsub0 j (Finset.mem_univ j)
    have hgm : g j * |x j| = m :=
      mul_left_cancel₀ (ne_of_gt hFj) (by linarith)
    have hxj : |x j| = m := by
      have hle1 : g j * |x j| ≤ |x j| :=
        mul_le_of_le_one_left (abs_nonneg _) (hg1 j)
      have h2 : m ≤ |x j| := by
        rw [← hgm]
        exact hle1
      exact le_antisymm (hmax j (Finset.mem_univ j)) h2
    have hgj : g j = 1 := by
      have hgmm : g j * m = 1 * m := by
        rw [one_mul]
        calc g j * m = g j * |x j| := by rw [hxj]
          _ = m := hgm
      exact mul_right_cancel₀ (ne_of_gt hm) hgmm
    exact ⟨hgj, hxj⟩
  -- Step B: the max-set floods along chains to the whole domain
  have flood : ∀ j, |x j| = m := by
    intro j
    have hpath := hirr kmax j
    induction hpath with
    | refl => rfl
    | tail hpath hedge ih => exact (step _ ih _ hedge).2
  -- Step C: the in-neighbor of j₀ forces g j₀ = 1 — contradiction
  obtain ⟨k, hk⟩ := hcol
  have : g j₀ = 1 := (step k (flood k) j₀ hk).1
  linarith

/-- **Theorem 3.5, multi-row case**: for any replacement set S containing
    some j₀ with b j₀ < 1 (necessary hypothesis, see header), the mixed
    matrix has trivial left null space. NO Perron–Frobenius input. -/
theorem thm35_multirow_left_null_trivial {S : Finset n}
    (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : ∀ i, ∑ j, F i j = 1)
    (hirr : StronglyConnected F)
    (hb0 : ∀ j ∈ S, 0 ≤ b j) (hb1 : ∀ j ∈ S, b j ≤ 1)
    {j₀ : n} (hj₀S : j₀ ∈ S) (hj₀ : b j₀ < 1) :
    ∀ x : n → ℝ, (MmatS F b S)ᵀ *ᵥ x = 0 → x = 0 := by
  intro x hx
  refine masked_no_fixed_vector hF0 hF1 hirr
    (g := gmask b S) ?_ ?_ (j₀ := j₀) ?_
    (exists_pos_in_column hF1 hirr j₀) (fun k => MmatS_null_fixed hx k)
  · intro j
    by_cases hjS : j ∈ S
    · simpa [gmask, hjS] using hb0 j hjS
    · simp [gmask, hjS]
  · intro j
    by_cases hjS : j ∈ S
    · simpa [gmask, hjS] using hb1 j hjS
    · simp [gmask, hjS]
  · simpa [gmask, hj₀S] using hj₀

/-- **Theorem 3.5, multi-row case: M is non-singular.** -/
theorem thm35_multirow_isUnit {S : Finset n}
    (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : ∀ i, ∑ j, F i j = 1)
    (hirr : StronglyConnected F)
    (hb0 : ∀ j ∈ S, 0 ≤ b j) (hb1 : ∀ j ∈ S, b j ≤ 1)
    {j₀ : n} (hj₀S : j₀ ∈ S) (hj₀ : b j₀ < 1) :
    IsUnit (MmatS F b S) := by
  rw [Matrix.isUnit_iff_isUnit_det, isUnit_iff_ne_zero]
  intro hdet
  rw [← Matrix.det_transpose] at hdet
  obtain ⟨v, hv, hv0⟩ := Matrix.exists_mulVec_eq_zero_iff.mpr hdet
  exact hv (thm35_multirow_left_null_trivial hF0 hF1 hirr hb0 hb1
    hj₀S hj₀ v hv0)

/-- **Unification**: S = univ recovers Theorem 3.2(c) — D = I − Rᵀ is
    invertible — now under {irreducible + ONE b j₀ < 1} in place of
    {bmax < 1}: neither hypothesis set implies the other, and together
    the two results cover complementary regimes. -/
theorem thm32c_of_multirow
    (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : ∀ i, ∑ j, F i j = 1)
    (hirr : StronglyConnected F)
    (hb0 : ∀ j, 0 ≤ b j) (hb1 : ∀ j, b j ≤ 1)
    {j₀ : n} (hj₀ : b j₀ < 1) :
    IsUnit (Dmat F b) := by
  have hDM : MmatS F b Finset.univ = Dmat F b := by
    ext i j
    simp [MmatS, Finset.mem_univ]
  rw [← hDM]
  exact thm35_multirow_isUnit hF0 hF1 hirr
    (fun j _ => hb0 j) (fun j _ => hb1 j) (Finset.mem_univ j₀) hj₀

end Thm35Multi

#check @masked_no_fixed_vector
#check @thm35_multirow_left_null_trivial
#check @thm35_multirow_isUnit
#check @thm32c_of_multirow
