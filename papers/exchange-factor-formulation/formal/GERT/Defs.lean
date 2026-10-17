import Mathlib

/-!
# Shared definitions

Formal verification of Bielefeld (2026), arXiv:2512.22157v3,
"Non-Negativity, Non-Singularity and Energy Conservation of the Zonal Method
with Non-Reciprocal Exchange Factors".

Every module of the development imports this file. It fixes, once, the
objects of the paper's Section 2 and the masked form of the mixed-boundary
matrix of Section 3.2:

* `Amat`, `Rmat` — A = F ∘ (𝟙 − B), R = F ∘ B, entrywise (eqs. 9–10)
* `Cmat`, `Dmat` — C = I − Fᵀ, D = I − Rᵀ (eq. 11, with A + R = F)
* `MmatS` — the mixed-boundary matrix for a prescribed set S : Finset n
  (rows in S from D, others from C; eq. 12)
* `bmask` — the masked coefficients b̃ of eq. 23 (`bmask b S k` is the paper's b̃ₖ)
* `MmatS_apply` — eq. 24 entrywise: Mₖᵢ = δₖᵢ − b̃ₖ Fᵢₖ
* `MmatS_fixed_of_transpose_mulVec` — eq. 25: Mᵀx = c ⟹ xₖ = Σₗ Fₖₗ b̃ₗ xₗ + cₖ
* `StronglyConnected` — irreducibility as strong connectivity of the
  positive-entry relation (§2.1)
* `exists_pos_in_column` — every element has an in-neighbour
* `Dmat_eq_Cmat_add_transpose` — D = C + Aᵀ (eq. 28)
* `MmatS_univ` — S = all elements gives M = D

NB: the paper's "𝟙" in (𝟙 − B) is the ALL-ONES matrix, while Mathlib's
`(1 : Matrix n n ℝ)` is the IDENTITY. A and R are therefore defined
entrywise, never via `(1 : Matrix _ _ _) - B`.
-/

open Matrix Finset

namespace GERT

variable {n : Type*}

/-- Irreducibility of F (paper §2.1): from every element i there is a chain
    i → k₁ → … → j of strictly positive entries to every element j. -/
def StronglyConnected (F : Matrix n n ℝ) : Prop :=
  ∀ i j, Relation.ReflTransGen (fun a c => 0 < F a c) i j

section Definitions

variable (F : Matrix n n ℝ) (b : n → ℝ)

/-- Absorption matrix A = F ∘ (𝟙 − B), entrywise `A i j = F i j * (1 - b j)` (eq. 9). -/
def Amat : Matrix n n ℝ := Matrix.of fun i j => F i j * (1 - b j)

/-- Reflection-scattering matrix R = F ∘ B, entrywise `R i j = F i j * b j` (eq. 10). -/
def Rmat : Matrix n n ℝ := Matrix.of fun i j => F i j * b j

variable [DecidableEq n]

/-- C = I − Fᵀ (eq. 11, using A + R = F). -/
def Cmat : Matrix n n ℝ := 1 - Fᵀ

/-- D = I − Rᵀ (eq. 11). -/
def Dmat : Matrix n n ℝ := 1 - (Rmat F b)ᵀ

/-- Mixed-boundary matrix for the prescribed set S: rows in S from D,
    rows outside S from C (eq. 12). -/
def MmatS (S : Finset n) : Matrix n n ℝ :=
  Matrix.of fun i j => if i ∈ S then Dmat F b i j else Cmat F i j

/-- Masked reflection-scattering coefficients b̃ (eq. 23):
    `bmask b S k = b k` for k ∈ S and `= 1` for k ∉ S. -/
def bmask (S : Finset n) : n → ℝ := fun k => if k ∈ S then b k else 1

end Definitions

section Connectivity

variable [Fintype n] {F : Matrix n n ℝ}

/-- Every element of a strongly connected row-stochastic matrix has an
    in-neighbour: the row of j₀ has a positive entry F j₀ j₁ (rows summing
    to 1 cannot be all non-positive), and the return chain j₁ → j₀ supplies
    the in-edge (its last step; or the self-edge when j₁ = j₀). -/
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
  rcases (hirr j₁ j₀).cases_tail with heq | ⟨c, -, hc⟩
  · rw [← heq] at hj₁
    exact ⟨j₀, hj₁⟩
  · exact ⟨c, hc⟩

end Connectivity

section MixedSystem

variable [Fintype n] [DecidableEq n] {F : Matrix n n ℝ} {b : n → ℝ}

omit [Fintype n] in
/-- **D = C + Aᵀ** (eq. 28). -/
theorem Dmat_eq_Cmat_add_transpose : Dmat F b = Cmat F + (Amat F b)ᵀ := by
  ext i j
  simp only [Dmat, Cmat, Amat, Rmat, Matrix.sub_apply, Matrix.add_apply,
    Matrix.transpose_apply, Matrix.of_apply]
  ring

omit [Fintype n] in
/-- Prescribed elements keep their coefficient: b̃ₖ = bₖ for k ∈ S. -/
theorem bmask_of_mem {S : Finset n} {k : n} (hkS : k ∈ S) : bmask b S k = b k := by
  simp [bmask, hkS]

omit [Fintype n] in
/-- Free elements re-emit everything: b̃ₖ = 1 for k ∉ S. -/
theorem bmask_of_not_mem {S : Finset n} {k : n} (hkS : k ∉ S) : bmask b S k = 1 := by
  simp [bmask, hkS]

omit [Fintype n] in
/-- 0 ≤ b̃ everywhere, given 0 ≤ b on S. -/
theorem bmask_nonneg {S : Finset n} (hb0 : ∀ k ∈ S, 0 ≤ b k) (k : n) :
    0 ≤ bmask b S k := by
  by_cases hkS : k ∈ S
  · simpa [bmask, hkS] using hb0 k hkS
  · simp [bmask, hkS]

omit [Fintype n] in
/-- b̃ ≤ 1 everywhere, given b ≤ 1 on S. -/
theorem bmask_le_one {S : Finset n} (hb1 : ∀ k ∈ S, b k ≤ 1) (k : n) :
    bmask b S k ≤ 1 := by
  by_cases hkS : k ∈ S
  · simpa [bmask, hkS] using hb1 k hkS
  · simp [bmask, hkS]

/-- When every emissive power is prescribed, b̃ = b. -/
theorem bmask_univ : bmask b (Finset.univ : Finset n) = b := by
  funext k
  simp [bmask]

/-- When every emissive power is prescribed, M = D. -/
theorem MmatS_univ : MmatS F b Finset.univ = Dmat F b := by
  ext i j
  simp [MmatS, Finset.mem_univ]

omit [Fintype n] in
/-- **Eq. 24, entrywise**: every row of M has the same form,
    Mₖᵢ = δₖᵢ − b̃ₖ Fᵢₖ, i.e. M = I − diag(b̃) Fᵀ. -/
theorem MmatS_apply (S : Finset n) (k i : n) :
    MmatS F b S k i = (1 : Matrix n n ℝ) k i - F i k * bmask b S k := by
  by_cases hkS : k ∈ S
  · simp [MmatS, hkS, Dmat, Rmat, bmask, Matrix.sub_apply,
      Matrix.transpose_apply]
  · simp [MmatS, hkS, Cmat, bmask, Matrix.sub_apply,
      Matrix.transpose_apply]

/-- The transpose system written entrywise: (Mᵀx)ₖ = xₖ − Σₗ Fₖₗ b̃ₗ xₗ. -/
theorem MmatS_transpose_mulVec_apply (S : Finset n) (x : n → ℝ) (k : n) :
    ((MmatS F b S)ᵀ *ᵥ x) k = x k - ∑ l, F k l * (bmask b S l * x l) := by
  simp only [Matrix.mulVec, dotProduct, Matrix.transpose_apply]
  rw [Finset.sum_congr rfl (fun l _ => by rw [MmatS_apply, sub_mul])]
  rw [Finset.sum_sub_distrib]
  have hid : ∑ l, (1 : Matrix n n ℝ) l k * x l = x k := by
    simp [Matrix.one_apply, ite_mul, one_mul, zero_mul,
      Finset.sum_ite_eq', Finset.mem_univ]
  rw [hid]
  congr 1
  exact Finset.sum_congr rfl fun l _ => by ring

/-- **Eq. 25**: a solution of Mᵀx = c satisfies xₖ = Σₗ Fₖₗ b̃ₗ xₗ + cₖ for every k.
    This is the form in which the minimum principle (Lemma 3.4) is applied. -/
theorem MmatS_fixed_of_transpose_mulVec {S : Finset n} {x c : n → ℝ}
    (hx : (MmatS F b S)ᵀ *ᵥ x = c) (k : n) :
    x k = ∑ l, F k l * (bmask b S l * x l) + c k := by
  have hk := congrFun hx k
  rw [MmatS_transpose_mulVec_apply] at hk
  linarith

end MixedSystem

end GERT
