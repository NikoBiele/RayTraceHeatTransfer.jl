/-
  Theorem 4.1, Bielefeld (2026), arXiv:2512.22157v2 — Mathlib version.
  "Non-negative Radiation"

  STATUS: rev 2. `masked_min_principle` and all rev-1 name gambles compiled
  first-try (`exists_min_image`, `mul_le_mul_of_nonpos_right`,
  `div_neg_of_neg_of_pos` all confirmed). Rev-2 fixes the two field_simp
  sites: hjw needed `w i ≠ 0` in context, and hyfix must unfold only the
  standalone `y k` — unfolding y inside Σ demands (w k)⁻¹·w k = 1, which
  `ring` correctly cannot do. Also omits unused [DecidableEq n] on parts 2–3.

  ── PROOF ROUTE (differs from the paper's citation) ─────────────────────
  The paper derives j = M⁻¹h ≥ 0 from the classical inverse-positivity of
  non-singular M-matrices. That classical proof (Neumann series + weak
  Perron) is not in Mathlib. Formalized here instead: a MINIMUM PRINCIPLE
  in reciprocity-scaled variables —
    Mj = h  reduces to  j_k = g_k·(Σᵢ Fᵢₖ jᵢ) + h_k   (column sums!),
    and y := j/w with reciprocity  wᵢFᵢₖ = Fₖᵢwₖ  turns this into
    y_k = g_k·(Σᵢ Fₖᵢ yᵢ) + h_k/wₖ   (row-stochastic form),
  where a negative minimum of y is impossible: at it, h must vanish,
  g must equal 1, and the min floods along F-chains (strong connectivity)
  to j₀ — contradicting g j₀ = b j₀ < 1.
  So uniqueness (Thm 3.5, max principle) and nonnegativity (Thm 4.1, min
  principle) are the SAME argument — the discrete elliptic maximum
  principle — and reciprocity (eq. 7 / Paper 2's enforced structure) is
  again the physically-supplied ingredient replacing a classical citation.
  Physical reading: y is radiant power per unit emission capacity.
  ─────────────────────────────────────────────────────────────────────────

  Contents:
  * `masked_min_principle`      — the main lemma (part 1's engine)
  * `thm41_solution_nonneg`     — part 1: Mj = h, h ≥ 0 ⟹ j ≥ 0
  * `thm41_reflected_nonneg`    — part 2: r = Rᵀj ≥ 0
  * `thm41_emissive_nonneg`     — part 3: e = q + Aᵀj ≥ 0
-/
import Mathlib

open Matrix Finset

variable {n : Type*} [Fintype n] [DecidableEq n]

def StronglyConnected (F : Matrix n n ℝ) : Prop :=
  ∀ i j, Relation.ReflTransGen (fun a c => 0 < F a c) i j

section Thm41

variable (F : Matrix n n ℝ) (b : n → ℝ)

def Amat : Matrix n n ℝ := Matrix.of fun i j => F i j * (1 - b j)
def Rmat : Matrix n n ℝ := Matrix.of fun i j => F i j * b j
def Cmat : Matrix n n ℝ := 1 - Fᵀ
def Dmat : Matrix n n ℝ := 1 - (Rmat F b)ᵀ

def MmatS (S : Finset n) : Matrix n n ℝ :=
  Matrix.of fun i j => if i ∈ S then Dmat F b i j else Cmat F i j

def gmask (S : Finset n) : n → ℝ := fun j => if j ∈ S then b j else 1

variable {F b}

omit [DecidableEq n] in
/-- **MAIN LEMMA (minimum principle)**: for irreducible row-stochastic
    non-negative F, mask 0 ≤ g ≤ 1 with g j₀ < 1, and source h ≥ 0, any
    solution of  y = g∘(F·y) + h  is non-negative. -/
theorem masked_min_principle
    (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : ∀ i, ∑ j, F i j = 1)
    (hirr : StronglyConnected F)
    {g : n → ℝ} (hg0 : ∀ k, 0 ≤ g k) (hg1 : ∀ k, g k ≤ 1)
    {j₀ : n} (hj₀ : g j₀ < 1)
    {h y : n → ℝ} (hh : ∀ k, 0 ≤ h k)
    (hfix : ∀ k, y k = g k * (∑ i, F k i * y i) + h k) :
    ∀ k, 0 ≤ y k := by
  by_contra hneg
  obtain ⟨k0, hk0'⟩ := not_forall.mp hneg
  have hk0 : y k0 < 0 := not_le.mp hk0'
  -- [RISK] `Finset.exists_min_image` (dual of the confirmed exists_max_image)
  obtain ⟨km, -, hmin⟩ :=
    Finset.exists_min_image (Finset.univ : Finset n) (fun k => y k)
      ⟨j₀, Finset.mem_univ j₀⟩
  set m := y km with hm_def
  have hm : m < 0 := lt_of_le_of_lt (hmin k0 (Finset.mem_univ k0)) hk0
  -- Step A: at a min coordinate, h = 0, g = 1, and every F-successor is min
  have step : ∀ k, y k = m → g k = 1 ∧ (∀ i, 0 < F k i → y i = m) := by
    intro k hk
    have hs_lb : m ≤ ∑ i, F k i * y i := by
      calc m = (∑ i, F k i) * m := by rw [hF1 k, one_mul]
        _ = ∑ i, F k i * m := Finset.sum_mul ..
        _ ≤ ∑ i, F k i * y i :=
            Finset.sum_le_sum fun i _ =>
              mul_le_mul_of_nonneg_left (hmin i (Finset.mem_univ i)) (hF0 k i)
    have hfx := hfix k
    rw [hk] at hfx
    -- hfx : m = g k * (∑ i, F k i * y i) + h k
    have h_gm : m ≤ g k * m := by
      -- [RISK] `mul_le_mul_of_nonpos_right : a ≤ b → c ≤ 0 → b * c ≤ a * c`
      have := mul_le_mul_of_nonpos_right (hg1 k) (le_of_lt hm)
      rwa [one_mul] at this
    have h_gs : g k * m ≤ g k * (∑ i, F k i * y i) :=
      mul_le_mul_of_nonneg_left hs_lb (hg0 k)
    have hh0 : h k = 0 := le_antisymm (by linarith) (hh k)
    have h_gsm : g k * (∑ i, F k i * y i) = m := by linarith
    have h_gmm : g k * m = m := by linarith
    have hgk : g k = 1 := by
      have hgm1 : g k * m = 1 * m := by rw [one_mul]; exact h_gmm
      exact mul_right_cancel₀ (ne_of_lt hm) hgm1
    have hs_eq : ∑ i, F k i * y i = m := by
      rw [hgk, one_mul] at h_gsm
      exact h_gsm
    refine ⟨hgk, ?_⟩
    intro i hFi
    have hsub0 : ∑ i, (F k i * y i - F k i * m) = 0 := by
      rw [Finset.sum_sub_distrib, hs_eq]
      have h1 : ∑ i, F k i * m = m := by
        rw [← Finset.sum_mul, hF1 k, one_mul]
      rw [h1, sub_self]
    have hz := (Finset.sum_eq_zero_iff_of_nonneg
      (fun i _ => sub_nonneg.mpr
        (mul_le_mul_of_nonneg_left (hmin i (Finset.mem_univ i)) (hF0 k i)))).mp
      hsub0 i (Finset.mem_univ i)
    have heq : F k i * y i = F k i * m := by linarith
    exact mul_left_cancel₀ (ne_of_gt hFi) heq
  -- Step B: the min-set floods along chains to the whole domain
  have flood : ∀ i, y i = m := by
    intro i
    have hpath := hirr km i
    induction hpath with
    | refl => rfl
    | tail hp he ih => exact (step _ ih).2 _ he
  -- Step C: j₀ itself is a min element, so g j₀ = 1 — contradiction
  have : g j₀ = 1 := (step j₀ (flood j₀)).1
  linarith

/-- **Theorem 4.1, part 1**: any solution of the mixed system Mj = h with
    h ≥ 0 is non-negative, given reciprocity with positive capacities w
    (paper eq. 7), row-stochasticity, strong connectivity, 0 ≤ b ≤ 1 on S,
    and some j₀ ∈ S with b j₀ < 1. -/
theorem thm41_solution_nonneg {S : Finset n}
    (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : ∀ i, ∑ j, F i j = 1)
    (hirr : StronglyConnected F)
    {w : n → ℝ} (hw : ∀ i, 0 < w i)
    (hrec : ∀ i j, w i * F i j = F j i * w j)
    (hb0 : ∀ j ∈ S, 0 ≤ b j) (hb1 : ∀ j ∈ S, b j ≤ 1)
    {j₀ : n} (hj₀S : j₀ ∈ S) (hj₀ : b j₀ < 1)
    {j h : n → ℝ} (hh : ∀ k, 0 ≤ h k)
    (hsol : MmatS F b S *ᵥ j = h) :
    ∀ k, 0 ≤ j k := by
  -- the mixed rows, uniformly:  M k i = 1ₖᵢ − Fᵢₖ · g k
  have hM : ∀ (k i : n), MmatS F b S k i
      = (1 : Matrix n n ℝ) k i - F i k * gmask b S k := by
    intro k i
    by_cases hkS : k ∈ S
    · simp [MmatS, hkS, Dmat, Rmat, gmask, Matrix.sub_apply,
        Matrix.transpose_apply]
    · simp [MmatS, hkS, Cmat, gmask, Matrix.sub_apply,
        Matrix.transpose_apply]
  -- reduction:  j k = g k · (Σᵢ Fᵢₖ jᵢ) + h k
  have hrowk : ∀ k, j k = gmask b S k * (∑ i, F i k * j i) + h k := by
    intro k
    have hk := congrFun hsol k
    simp only [Matrix.mulVec, dotProduct] at hk
    rw [Finset.sum_congr rfl (fun i _ => by rw [hM k i, sub_mul])] at hk
    rw [Finset.sum_sub_distrib] at hk
    have hid : ∑ i, (1 : Matrix n n ℝ) k i * j i = j k := by
      simp [Matrix.one_apply, ite_mul, one_mul, zero_mul,
        Finset.sum_ite_eq, Finset.mem_univ]
    rw [hid] at hk
    have hpull : ∑ i, F i k * gmask b S k * j i
        = gmask b S k * ∑ i, F i k * j i := by
      rw [Finset.mul_sum]
      exact Finset.sum_congr rfl fun i _ => by ring
    rw [hpull] at hk
    linarith
  -- reciprocity-scaled variables y = j / w
  set y : n → ℝ := fun k => j k / w k with hy
  have hjw : ∀ i, j i = y i * w i := by
    intro i
    have hwne : w i ≠ 0 := (hw i).ne'
    simp only [hy]
    field_simp
  have hyfix : ∀ k, y k = gmask b S k * (∑ i, F k i * y i) + h k / w k := by
    intro k
    have hwne : w k ≠ 0 := (hw k).ne'
    have hswap : ∑ i, F i k * j i = w k * ∑ i, F k i * y i := by
      calc ∑ i, F i k * j i
          = ∑ i, F i k * (y i * w i) :=
            Finset.sum_congr rfl fun i _ => by rw [← hjw i]
        _ = ∑ i, w k * (F k i * y i) := by
            refine Finset.sum_congr rfl fun i _ => ?_
            have hri := hrec i k
            calc F i k * (y i * w i) = (w i * F i k) * y i := by ring
              _ = (F k i * w k) * y i := by rw [hri]
              _ = w k * (F k i * y i) := by ring
        _ = w k * ∑ i, F k i * y i := (Finset.mul_sum ..).symm
    have hr := hrowk k
    rw [hswap] at hr
    -- hr : j k = g k · (w k · Σ) + h k ; goal in y with Σ kept atomic
    have hyk : y k = j k / w k := by simp only [hy]
    rw [hyk, hr]
    field_simp
  -- min principle in y, then scale back
  have hg0' : ∀ k, 0 ≤ gmask b S k := by
    intro k
    by_cases hkS : k ∈ S
    · simpa [gmask, hkS] using hb0 k hkS
    · simp [gmask, hkS]
  have hg1' : ∀ k, gmask b S k ≤ 1 := by
    intro k
    by_cases hkS : k ∈ S
    · simpa [gmask, hkS] using hb1 k hkS
    · simp [gmask, hkS]
  have hj₀' : gmask b S j₀ < 1 := by simpa [gmask, hj₀S] using hj₀
  have hh' : ∀ k, 0 ≤ h k / w k := fun k =>
    div_nonneg (hh k) (le_of_lt (hw k))
  have hy_nonneg := masked_min_principle hF0 hF1 hirr hg0' hg1' hj₀' hh' hyfix
  intro k
  by_contra hjk
  have hjk' : j k < 0 := not_le.mp hjk
  -- [RISK] `div_neg_of_neg_of_pos`
  have hy_neg : j k / w k < 0 := div_neg_of_neg_of_pos hjk' (hw k)
  have := hy_nonneg k
  simp only [hy] at this
  linarith

omit [DecidableEq n] in
/-- **Theorem 4.1, part 2**: reflected-scattered power r = Rᵀj is
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
/-- **Theorem 4.1, part 3**: emissive power e = q + Aᵀj is non-negative
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

end Thm41

#check @masked_min_principle
#check @thm41_solution_nonneg
#check @thm41_reflected_nonneg
#check @thm41_emissive_nonneg
