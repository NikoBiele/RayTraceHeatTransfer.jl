/-
  Theorem 3.3, Bielefeld (2026), arXiv:2512.22157v2
  "M-matrix Property of C"        C = I − Aᵀ − Rᵀ

  Self-contained core-Lean formalization (no Mathlib), companion to
  Thm31.lean / Thm32.lean.  All results over an abstract ordered ring
  (axioms satisfied by ℝ; ℤ instance below certifies consistency).

  Formalized content:
  * `thm33_C_eq`           — C = I − Fᵀ entrywise (uses A + R = F, Thm 3.1)
  * `thm33_offdiag_nonpos` — part 3: Cᵢⱼ ≤ 0 for i ≠ j
  * `thm33_diag_eq`        — Cᵢᵢ = 1 − Fᵢᵢ
  * `thm33_diag_nonneg`    — part 2: 0 ≤ Cᵢᵢ (eq. 28, via Fᵢᵢ ≤ Σⱼ Fᵢⱼ = 1)
  * `thm33_col_sum_zero`   — part 1: 𝟙ᵀC = 0, i.e. every column of C sums
                             to zero (eq. 25, the left-eigenvector identity)
  The spectral half of part 1 ("zero has smallest real part") needs complex
  eigenvalues → Mathlib companion file.
-/

import Init.PropLemmas

/-- Ordered-ring axioms used by Theorem 3.3 (all hold in ℝ). -/
class OrdRingC (α : Type) extends Add α, Mul α, Sub α, LE α, OfNat α 1 where
  zero : α
  le_refl : ∀ a : α, a ≤ a
  le_trans : ∀ {a b c : α}, a ≤ b → b ≤ c → a ≤ c
  add_le_add : ∀ {a b c d : α}, a ≤ b → c ≤ d → a + c ≤ b + d
  zero_add : ∀ a : α, zero + a = a
  add_zero : ∀ a : α, a + zero = a
  sub_self : ∀ a : α, a - a = zero
  sub_sub : ∀ a b c : α, a - b - c = a - (b + c)
  sub_add_sub : ∀ a b c d : α, (a - b) + (c - d) = (a + c) - (b + d)
  sub_zero : ∀ a : α, a - zero = a
  sub_nonneg_of_le : ∀ {a b : α}, a ≤ b → zero ≤ b - a
  sub_nonpos_of_nonneg : ∀ {a : α}, zero ≤ a → zero - a ≤ zero
  le_add_of_nonneg_right : ∀ {a b : α}, zero ≤ b → a ≤ a + b
  le_add_of_nonneg_left : ∀ {a b : α}, zero ≤ b → a ≤ b + a
  mul_sub : ∀ a b c : α, a * (b - c) = a * b - a * c
  mul_one : ∀ a : α, a * 1 = a
  sub_add_cancel : ∀ a b : α, a - b + b = a

namespace OrdRingC

variable {α : Type} [OrdRingC α] {n : Nat}

abbrev Mat (α : Type) (n : Nat) := Fin n → Fin n → α

/-- Row sum by structural recursion. -/
def rsum : {n : Nat} → (Fin n → α) → α
  | 0, _ => zero
  | _+1, f => f 0 + rsum (fun i => f i.succ)

theorem rsum_congr {f g : Fin n → α} (h : ∀ j, f j = g j) : rsum f = rsum g :=
  congrArg rsum (funext h)

theorem rsum_zero_fun : ∀ {n : Nat}, rsum (fun _ : Fin n => (zero : α)) = zero
  | 0 => rfl
  | n+1 => by
    show zero + rsum (fun _ : Fin n => (zero : α)) = zero
    rw [rsum_zero_fun, add_zero]

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

/-- Sums split across subtraction. -/
theorem rsum_sub : ∀ {n : Nat} (f g : Fin n → α),
    rsum (fun j => f j - g j) = rsum f - rsum g
  | 0, _, _ => (sub_zero zero).symm
  | _+1, f, g => by
    show (f 0 - g 0) + rsum (fun j => f j.succ - g j.succ)
       = (f 0 + rsum (fun j => f j.succ)) - (g 0 + rsum (fun j => g j.succ))
    rw [rsum_sub (fun j => f j.succ) (fun j => g j.succ), sub_add_sub]

/-- Kronecker delta. -/
def delta (i j : Fin n) : α := if i = j then (1 : α) else zero

/-- A delta row sums to 1. -/
theorem rsum_delta : ∀ {n : Nat} (j : Fin n),
    rsum (fun i => (delta i j : α)) = 1 := by
  intro n
  induction n with
  | zero => intro j; exact j.elim0
  | succ m ih =>
    intro j
    cases j using Fin.cases with
    | zero =>
      show delta 0 0 + rsum (fun i : Fin m => (delta i.succ 0 : α)) = 1
      have h1 : (delta (0 : Fin (m+1)) 0 : α) = 1 := ite_eq_left rfl
      have h2 : (fun i : Fin m => (delta i.succ 0 : α)) = fun _ => zero := by
        funext i
        exact ite_eq_right (Fin.succ_ne_zero i)
      rw [h1, h2, rsum_zero_fun, add_zero]
    | succ k =>
      show delta 0 k.succ + rsum (fun i : Fin m => (delta i.succ k.succ : α)) = 1
      have h1 : (delta (0 : Fin (m+1)) k.succ : α) = zero :=
        ite_eq_right (fun h => Fin.succ_ne_zero k h.symm)
      have h2 : (fun i : Fin m => (delta i.succ k.succ : α))
              = fun i => (delta i k : α) := by
        funext i
        show (if i.succ = k.succ then (1:α) else zero) = if i = k then (1:α) else zero
        by_cases h : i = k
        · rw [ite_eq_left (by rw [h]), ite_eq_left h]
        · rw [ite_eq_right (fun hc => h (Fin.succ_inj.mp hc)), ite_eq_right h]
      rw [h1, h2, ih k, zero_add]

/-! ### The paper's matrices (eqs. 9–11) and Theorem 3.3 -/

/-- Absorption matrix A = F ∘ (1 − B) (eq. 9). -/
def Amat (F : Mat α n) (b : Fin n → α) : Mat α n := fun i j => F i j * (1 - b j)

/-- Reflection-scattering matrix R = F ∘ B (eq. 10). -/
def Rmat (F : Mat α n) (b : Fin n → α) : Mat α n := fun i j => F i j * b j

/-- C = I − Aᵀ − Rᵀ, exactly as defined in the paper (eq. 11). -/
def CmatDef (F : Mat α n) (b : Fin n → α) : Mat α n :=
  fun i j => delta i j - Amat F b j i - Rmat F b j i

/-- The simplified form C = I − Fᵀ. -/
def Cmat (F : Mat α n) : Mat α n := fun i j => delta i j - F j i

/-- **C = I − Fᵀ**: the paper's C (eq. 11) collapses via A + R = F (Thm 3.1). -/
theorem thm33_C_eq (F : Mat α n) (b : Fin n → α) :
    CmatDef F b = Cmat F := by
  funext i j
  show delta i j - F j i * (1 - b i) - F j i * b i = delta i j - F j i
  rw [sub_sub, mul_sub, mul_one, sub_add_cancel]

variable {F : Mat α n}

/-- **Theorem 3.3, part 3**: off-diagonal entries are non-positive. -/
theorem thm33_offdiag_nonpos (hF0 : ∀ i j, zero ≤ F i j)
    {i j : Fin n} (hij : i ≠ j) : Cmat F i j ≤ zero := by
  show delta i j - F j i ≤ zero
  rw [show (delta i j : α) = zero from ite_eq_right hij]
  exact sub_nonpos_of_nonneg (hF0 j i)

/-- Diagonal entries: Cᵢᵢ = 1 − Fᵢᵢ (the value inside eq. 28). -/
theorem thm33_diag_eq (i : Fin n) : Cmat F i i = 1 - F i i := by
  show delta i i - F i i = 1 - F i i
  rw [show (delta i i : α) = 1 from ite_eq_left rfl]

/-- **Theorem 3.3, part 2**: diagonal entries are non-negative (eq. 28):
    Fᵢᵢ ≤ Σⱼ Fᵢⱼ = 1, hence 0 ≤ 1 − Fᵢᵢ = Cᵢᵢ. -/
theorem thm33_diag_nonneg (hF0 : ∀ i j, zero ≤ F i j)
    (hF1 : ∀ i, rsum (F i) = 1) (i : Fin n) : zero ≤ Cmat F i i := by
  rw [thm33_diag_eq]
  have h : F i i ≤ (1 : α) := by
    have := rsum_single_le (f := F i) (fun j => hF0 i j) i
    rwa [hF1 i] at this
  exact sub_nonneg_of_le h

/-- **Theorem 3.3, part 1** (eq. 25): 𝟙ᵀC = 0, i.e. every column of C sums to
    zero.  (Σᵢ Cᵢⱼ = Σᵢ δᵢⱼ − Σᵢ Fⱼᵢ = 1 − 1 = 0.)  Combined with parts 2–3
    this is the singular-M-matrix structure; the spectral statement
    "zero has smallest real part" lives in the Mathlib companion. -/
theorem thm33_col_sum_zero (hF1 : ∀ i, rsum (F i) = 1) (j : Fin n) :
    rsum (fun i => Cmat F i j) = zero := by
  show rsum (fun i => (delta i j : α) - F j i) = zero
  rw [rsum_sub (fun i => (delta i j : α)) (fun i => F j i),
      rsum_delta j, hF1 j, sub_self]

end OrdRingC

/-! ### Non-vacuity: ℤ satisfies all the axioms -/

instance : OrdRingC Int where
  zero := 0
  le_refl := fun a => by omega
  le_trans := fun h1 h2 => by omega
  add_le_add := fun h1 h2 => by omega
  zero_add := fun a => by omega
  add_zero := fun a => by omega
  sub_self := fun a => by omega
  sub_sub := fun a b c => by omega
  sub_add_sub := fun a b c d => by omega
  sub_zero := fun a => by omega
  sub_nonneg_of_le := fun h => by omega
  sub_nonpos_of_nonneg := fun h => by omega
  le_add_of_nonneg_right := fun h => by omega
  le_add_of_nonneg_left := fun h => by omega
  mul_sub := Int.mul_sub
  mul_one := Int.mul_one
  sub_add_cancel := fun a b => by omega

open OrdRingC in
#check @thm33_C_eq
open OrdRingC in
#check @thm33_offdiag_nonpos
open OrdRingC in
#check @thm33_diag_nonneg
open OrdRingC in
#check @thm33_col_sum_zero

/-! ### Numerical smoke test over ℤ (Fin 2): F = antidiagonal (row sums 1). -/
section SmokeTest
open OrdRingC

def Ftest : Mat Int 2 := fun i j => if i = j then 0 else 1

#eval decide (∀ j : Fin 2, rsum (fun i => Cmat Ftest i j) = 0)
#eval decide (∀ i : Fin 2, 0 ≤ Cmat Ftest i i)
#eval decide (∀ i j : Fin 2, i ≠ j → Cmat Ftest i j ≤ 0)
end SmokeTest
