/-
  Theorem 3.4, Bielefeld (2026), arXiv:2512.22157v2
  "M-matrix Property of D"        D = I − Rᵀ

  Self-contained core-Lean formalization (no Mathlib).

  Formalized content:
  * `thm34_D_eq_C_add_At` — the identity D = C + Aᵀ used in the paper's Step 2
  * `thm34_offdiag_nonpos` — part 3: Dᵢⱼ ≤ 0 for i ≠ j
  * `thm34_diag_eq`        — Dᵢᵢ = 1 − Fᵢᵢ·bᵢ
  * `thm34_diag_nonneg`    — part 2: 0 ≤ Dᵢᵢ (needs 0 ≤ b ≤ 1: Fᵢᵢbᵢ ≤ Fᵢᵢ ≤ 1)
  Part 1 (Re λ_D ≥ 1 − ρ(R)) needs ℂ → Mathlib companion.
  Non-singularity for bmax < 1 is already verified: `thm32_D_isUnit`
  (Thm32_mathlib.lean) / `thm32_no_unit_eigen` (Thm32_basic.lean).
-/

/-- Ordered-ring axioms used by Theorem 3.4 (all hold in ℝ). -/
class OrdRingD (α : Type) extends Add α, Mul α, Sub α, LE α, OfNat α 1 where
  zero : α
  le_refl : ∀ a : α, a ≤ a
  le_trans : ∀ {a b c : α}, a ≤ b → b ≤ c → a ≤ c
  add_le_add : ∀ {a b c d : α}, a ≤ b → c ≤ d → a + c ≤ b + d
  zero_add : ∀ a : α, zero + a = a
  sub_nonneg_of_le : ∀ {a b : α}, a ≤ b → zero ≤ b - a
  sub_nonpos_of_nonneg : ∀ {a : α}, zero ≤ a → zero - a ≤ zero
  le_add_of_nonneg_right : ∀ {a b : α}, zero ≤ b → a ≤ a + b
  le_add_of_nonneg_left : ∀ {a b : α}, zero ≤ b → a ≤ b + a
  sub_add_sub_cancel : ∀ a b c : α, (a - b) + (b - c) = a - c
  mul_sub : ∀ a b c : α, a * (b - c) = a * b - a * c
  mul_one : ∀ a : α, a * 1 = a
  mul_nonneg : ∀ {a b : α}, zero ≤ a → zero ≤ b → zero ≤ a * b
  mul_le_mul_of_nonneg_left : ∀ {a b c : α}, b ≤ c → zero ≤ a → a * b ≤ a * c

namespace OrdRingD

variable {α : Type} [OrdRingD α] {n : Nat}

abbrev Mat (α : Type) (n : Nat) := Fin n → Fin n → α

def rsum : {n : Nat} → (Fin n → α) → α
  | 0, _ => zero
  | _+1, f => f 0 + rsum (fun i => f i.succ)

theorem rsum_nonneg : ∀ {n : Nat} {f : Fin n → α},
    (∀ j, zero ≤ f j) → zero ≤ rsum f
  | 0, _, _ => le_refl zero
  | _+1, f, h => by
    have := add_le_add (h 0) (rsum_nonneg (fun j => h j.succ))
    rwa [zero_add] at this

/-- Each entry of a non-negative vector is at most the total sum. -/
theorem rsum_single_le : ∀ {n : Nat} {f : Fin n → α},
    (∀ j, zero ≤ f j) → ∀ k, f k ≤ rsum f := by
  intro n
  induction n with
  | zero => intro f _ k; exact k.elim0
  | succ m ih =>
    intro f h k
    cases k using Fin.cases with
    | zero =>
      exact le_add_of_nonneg_right (rsum_nonneg (fun j => h j.succ))
    | succ k' =>
      exact le_trans (ih (fun j => h j.succ) k') (le_add_of_nonneg_left (h 0))

def delta (i j : Fin n) : α := if i = j then (1 : α) else zero

/-! ### The paper's matrices and Theorem 3.4 -/

/-- Absorption matrix A = F ∘ (1 − B) (eq. 9). -/
def Amat (F : Mat α n) (b : Fin n → α) : Mat α n := fun i j => F i j * (1 - b j)

/-- C = I − Fᵀ (equal to the paper's I − Aᵀ − Rᵀ by Thm 3.3 / `thm33_C_eq`). -/
def Cmat (F : Mat α n) : Mat α n := fun i j => delta i j - F j i

/-- D = I − Rᵀ entrywise: Dᵢⱼ = δᵢⱼ − Fⱼᵢ·bᵢ (eq. 11 with eq. 10). -/
def Dmat (F : Mat α n) (b : Fin n → α) : Mat α n :=
  fun i j => delta i j - F j i * b i

/-- **The identity D = C + Aᵀ** (used in the paper's Step 2, eq. 30). -/
theorem thm34_D_eq_C_add_At (F : Mat α n) (b : Fin n → α) (i j : Fin n) :
    Dmat F b i j = Cmat F i j + Amat F b j i := by
  show delta i j - F j i * b i = (delta i j - F j i) + F j i * (1 - b i)
  rw [mul_sub, mul_one, sub_add_sub_cancel]

variable {F : Mat α n} {b : Fin n → α}

/-- **Theorem 3.4, part 3**: off-diagonal entries are non-positive. -/
theorem thm34_offdiag_nonpos (hF0 : ∀ i j, zero ≤ F i j)
    (hb0 : ∀ j, zero ≤ b j) {i j : Fin n} (hij : i ≠ j) :
    Dmat F b i j ≤ zero := by
  show delta i j - F j i * b i ≤ zero
  rw [show (delta i j : α) = zero from ite_eq_right hij]
  exact sub_nonpos_of_nonneg (mul_nonneg (hF0 j i) (hb0 i))

/-- Diagonal entries: Dᵢᵢ = 1 − Fᵢᵢ·bᵢ. -/
theorem thm34_diag_eq (i : Fin n) : Dmat F b i i = 1 - F i i * b i := by
  show delta i i - F i i * b i = 1 - F i i * b i
  rw [show (delta i i : α) = 1 from ite_eq_left rfl]

/-- **Theorem 3.4, part 2**: diagonal entries are non-negative
    (Fᵢᵢbᵢ ≤ Fᵢᵢ·1 = Fᵢᵢ ≤ Σⱼ Fᵢⱼ = 1). NB: needs only b ≤ 1, not 0 ≤ b —
    non-negativity of b is consumed only by the off-diagonal part. -/
theorem thm34_diag_nonneg (hF0 : ∀ i j, zero ≤ F i j)
    (hF1 : ∀ i, rsum (F i) = 1)
    (hb1 : ∀ j, b j ≤ 1) (i : Fin n) :
    zero ≤ Dmat F b i i := by
  rw [thm34_diag_eq]
  have h1 : F i i * b i ≤ F i i :=
    le_trans (mul_le_mul_of_nonneg_left (hb1 i) (hF0 i i))
      (by rw [mul_one]; exact le_refl _)
  have h2 : F i i ≤ (1 : α) := by
    have := rsum_single_le (f := F i) (fun j => hF0 i j) i
    rwa [hF1 i] at this
  exact sub_nonneg_of_le (le_trans h1 h2)

end OrdRingD

/-! ### Non-vacuity: ℤ satisfies all the axioms -/

instance : OrdRingD Int where
  zero := 0
  le_refl := fun a => by omega
  le_trans := fun h1 h2 => by omega
  add_le_add := fun h1 h2 => by omega
  zero_add := fun a => by omega
  sub_nonneg_of_le := fun h => by omega
  sub_nonpos_of_nonneg := fun h => by omega
  le_add_of_nonneg_right := fun h => by omega
  le_add_of_nonneg_left := fun h => by omega
  sub_add_sub_cancel := fun a b c => by omega
  mul_sub := Int.mul_sub
  mul_one := Int.mul_one
  mul_nonneg := Int.mul_nonneg
  mul_le_mul_of_nonneg_left := fun h ha => Int.mul_le_mul_of_nonneg_left h ha

open OrdRingD in
#check @thm34_D_eq_C_add_At
open OrdRingD in
#check @thm34_offdiag_nonpos
open OrdRingD in
#check @thm34_diag_nonneg

/-! ### Numerical smoke test over ℤ (Fin 2): F = antidiagonal, b = [1, 0]. -/
section SmokeTest
open OrdRingD

def Ftest : Mat Int 2 := fun i j => if i = j then 0 else 1
def btest : Fin 2 → Int := fun j => if j = 0 then 1 else 0

#eval decide (∀ i j : Fin 2,
  Dmat Ftest btest i j = Cmat Ftest i j + Amat Ftest btest j i)
#eval decide (∀ i : Fin 2, 0 ≤ Dmat Ftest btest i i)
#eval decide (∀ i j : Fin 2, i ≠ j → Dmat Ftest btest i j ≤ 0)
end SmokeTest
