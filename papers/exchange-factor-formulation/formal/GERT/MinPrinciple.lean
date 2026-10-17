import GERT.Defs

/-!
# Lemma 3.4: the masked minimum principle

Bielefeld (2026), arXiv:2512.22157v3, Section 3.2.

Let F be non-negative, row stochastic and irreducible, let 0 ≤ b̃ ≤ 1 with
b̃ₖ₀ < 1 for at least one element k₀, and let c ≥ 0. If x satisfies

  xₖ = Σₗ Fₖₗ b̃ₗ xₗ + cₖ   for all k   (eq. 25)

then x ≥ 0.

The proof is the paper's: at a negative minimum μ = min x every term of the
weighted average is at least μ, so equality holds throughout; this forces
cₖ = 0 and, on every out-neighbour l, xₗ = μ and b̃ₗ = 1. The set of minimum
elements is therefore closed under direct exchange, floods the whole domain
by strong connectivity, and the in-neighbour of k₀ then forces b̃ₖ₀ = 1,
a contradiction.

The mask is an arbitrary vector `g` here; the mixed-boundary system supplies
`g = bmask b S`, the paper's b̃, in the modules for Theorems 3.5 and 4.1.
No reciprocity is assumed.
-/

open Matrix Finset

namespace GERT

variable {n : Type*} [Fintype n] {F : Matrix n n ℝ}

/-- **Lemma 3.4 (Masked Minimum Principle).** -/
theorem masked_min_principle
    (hF0 : ∀ i j, 0 ≤ F i j) (hF1 : ∀ i, ∑ j, F i j = 1)
    (hirr : StronglyConnected F)
    {g : n → ℝ} (hg0 : ∀ k, 0 ≤ g k) (hg1 : ∀ k, g k ≤ 1)
    {k₀ : n} (hk₀ : g k₀ < 1)
    {x c : n → ℝ} (hc : ∀ k, 0 ≤ c k)
    (hfix : ∀ k, x k = ∑ l, F k l * (g l * x l) + c k) :
    ∀ k, 0 ≤ x k := by
  by_contra hneg
  obtain ⟨k1, hk1'⟩ := not_forall.mp hneg
  have hk1 : x k1 < 0 := not_le.mp hk1'
  -- μ = min x < 0, attained at km
  obtain ⟨km, -, hmin⟩ :=
    Finset.exists_min_image (Finset.univ : Finset n) (fun k => x k)
      ⟨k₀, Finset.mem_univ k₀⟩
  set μ := x km with hμ_def
  have hμ : μ < 0 := lt_of_le_of_lt (hmin k1 (Finset.mem_univ k1)) hk1
  -- the pointwise bound: μ ≤ g l · x l for every l
  have hlb : ∀ l, μ ≤ g l * x l := by
    intro l
    by_cases hxl : 0 ≤ x l
    · have : 0 ≤ g l * x l := mul_nonneg (hg0 l) hxl
      linarith
    · have hxl' : x l < 0 := not_le.mp hxl
      have h1 : 1 * x l ≤ g l * x l :=
        mul_le_mul_of_nonpos_right (hg1 l) (le_of_lt hxl')
      rw [one_mul] at h1
      exact le_trans (hmin l (Finset.mem_univ l)) h1
  -- Step A: at a minimum element k, c k = 0, and every out-neighbour l
  -- is a minimum element with g l = 1
  have step : ∀ k, x k = μ → c k = 0 ∧ ∀ l, 0 < F k l → x l = μ ∧ g l = 1 := by
    intro k hk
    have hterm : ∀ l, F k l * μ ≤ F k l * (g l * x l) :=
      fun l => mul_le_mul_of_nonneg_left (hlb l) (hF0 k l)
    have hsumμ : ∑ l, F k l * μ = μ := by
      rw [← Finset.sum_mul, hF1 k, one_mul]
    have hge : μ ≤ ∑ l, F k l * (g l * x l) :=
      calc μ = ∑ l, F k l * μ := hsumμ.symm
        _ ≤ ∑ l, F k l * (g l * x l) := Finset.sum_le_sum fun l _ => hterm l
    have hfx := hfix k
    rw [hk] at hfx
    -- hfx : μ = ∑ l, F k l * (g l * x l) + c k, with the sum ≥ μ and c k ≥ 0
    have hck : c k = 0 := le_antisymm (by linarith) (hc k)
    have hS : ∑ l, F k l * (g l * x l) = μ := by linarith
    refine ⟨hck, ?_⟩
    intro l hFl
    -- equality throughout: every term F k l (g l x l) equals F k l μ
    have hsub0 : ∑ l, (F k l * (g l * x l) - F k l * μ) = 0 := by
      rw [Finset.sum_sub_distrib, hS, hsumμ, sub_self]
    have hz := (Finset.sum_eq_zero_iff_of_nonneg
      (fun l _ => sub_nonneg.mpr (hterm l))).mp hsub0 l (Finset.mem_univ l)
    have hgx : g l * x l = μ :=
      mul_left_cancel₀ (ne_of_gt hFl) (by linarith)
    -- hence x l < 0, and μ ≤ x l ≤ g l x l = μ
    have hxl : x l = μ := by
      by_cases hxl0 : 0 ≤ x l
      · have : 0 ≤ g l * x l := mul_nonneg (hg0 l) hxl0
        linarith
      · have hxl' : x l < 0 := not_le.mp hxl0
        have h1 : 1 * x l ≤ g l * x l :=
          mul_le_mul_of_nonpos_right (hg1 l) (le_of_lt hxl')
        rw [one_mul] at h1
        exact le_antisymm (by linarith) (hmin l (Finset.mem_univ l))
    -- and g l μ = μ with μ ≠ 0 gives g l = 1
    have hgl : g l = 1 := by
      have h1 : g l * μ = 1 * μ := by
        rw [one_mul]
        calc g l * μ = g l * x l := by rw [hxl]
          _ = μ := hgx
      exact mul_right_cancel₀ (ne_of_lt hμ) h1
    exact ⟨hxl, hgl⟩
  -- Step B: the minimum set is closed under direct exchange, so by strong
  -- connectivity it floods the whole domain
  have flood : ∀ l, x l = μ := by
    intro l
    have hpath := hirr km l
    induction hpath with
    | refl => rfl
    | tail hp he ih => exact ((step _ ih).2 _ he).1
  -- Step C: the in-neighbour of k₀ forces g k₀ = 1, contradicting g k₀ < 1
  obtain ⟨i, hi⟩ := exists_pos_in_column hF1 hirr k₀
  have : g k₀ = 1 := ((step i (flood i)).2 k₀ hi).2
  linarith

end GERT
