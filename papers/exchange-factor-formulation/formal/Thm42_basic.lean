/-
  Theorem 4.2, Bielefeld (2026), arXiv:2512.22157v2
  "Unconditional Energy Conservation"     𝟙ᵀCj = 0

  Self-contained core-Lean formalization (no Mathlib).

  The theorem as formalized is the paper's actual content: the identity
  𝟙ᵀ(Cj) = 0 holds for EVERY vector j — independent of the mixed system M,
  the boundary conditions h, or the coefficients b — because every column
  of C = I − Fᵀ sums to zero (Thm 3.3 part 1). Only equalities are
  involved: no order structure is needed at all (the axiom class below has
  no ≤), which is the formal expression of "unconditional".
-/

/-- Ring axioms used by Theorem 4.2 (all hold in ℝ). No order needed. -/
class RingE (α : Type) extends Add α, Mul α, Sub α, OfNat α 1 where
  zero : α
  zero_add : ∀ a : α, zero + a = a
  add_zero : ∀ a : α, a + zero = a
  add_add_add_comm : ∀ a b c d : α, (a + b) + (c + d) = (a + c) + (b + d)
  sub_self : ∀ a : α, a - a = zero
  sub_zero : ∀ a : α, a - zero = a
  sub_add_sub : ∀ a b c d : α, (a - b) + (c - d) = (a + c) - (b + d)
  zero_mul : ∀ a : α, zero * a = zero
  add_mul : ∀ a b c : α, (a + b) * c = a * c + b * c

namespace RingE

variable {α : Type} [RingE α] {n : Nat}

abbrev Mat (α : Type) (n : Nat) := Fin n → Fin n → α

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

theorem rsum_add : ∀ {n : Nat} (f g : Fin n → α),
    rsum (fun j => f j + g j) = rsum f + rsum g
  | 0, _, _ => (add_zero zero).symm
  | _+1, f, g => by
    show (f 0 + g 0) + rsum (fun j => f j.succ + g j.succ)
       = (f 0 + rsum (fun j => f j.succ)) + (g 0 + rsum (fun j => g j.succ))
    rw [rsum_add (fun j => f j.succ) (fun j => g j.succ), add_add_add_comm]

theorem rsum_sub : ∀ {n : Nat} (f g : Fin n → α),
    rsum (fun j => f j - g j) = rsum f - rsum g
  | 0, _, _ => (sub_zero zero).symm
  | _+1, f, g => by
    show (f 0 - g 0) + rsum (fun j => f j.succ - g j.succ)
       = (f 0 + rsum (fun j => f j.succ)) - (g 0 + rsum (fun j => g j.succ))
    rw [rsum_sub (fun j => f j.succ) (fun j => g j.succ), sub_add_sub]

theorem rsum_mul (c : α) : ∀ {n : Nat} (f : Fin n → α),
    rsum (fun j => f j * c) = rsum f * c
  | 0, _ => (zero_mul c).symm
  | _+1, f => by
    show f 0 * c + rsum (fun j => f j.succ * c)
       = (f 0 + rsum (fun j => f j.succ)) * c
    rw [rsum_mul c (fun j => f j.succ), add_mul]

/-- Fubini for finite double sums. -/
theorem rsum_comm : ∀ {m n : Nat} (f : Fin m → Fin n → α),
    rsum (fun k => rsum (fun i => f k i)) = rsum (fun i => rsum (fun k => f k i)) := by
  intro m
  induction m with
  | zero =>
    intro n f
    show zero = rsum (fun _ : Fin n => (zero : α))
    rw [rsum_zero_fun]
  | succ m ih =>
    intro n f
    show rsum (fun i => f 0 i) + rsum (fun k : Fin m => rsum (fun i => f k.succ i))
       = rsum (fun i => f 0 i + rsum (fun k : Fin m => f k.succ i))
    rw [ih (fun k i => f k.succ i),
      rsum_add (fun i => f 0 i) (fun i => rsum (fun k : Fin m => f k.succ i))]

def delta (i j : Fin n) : α := if i = j then (1 : α) else zero

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

/-! ### Theorem 4.2 -/

/-- C = I − Fᵀ (Thm 3.3 / `thm33_C_eq`). -/
def Cmat (F : Mat α n) : Mat α n := fun i j => delta i j - F j i

/-- Columns of C sum to zero (verified in Thm33.lean; restated here for a
    standalone file). -/
theorem col_sum_zero {F : Mat α n} (hF1 : ∀ i, rsum (F i) = 1) (j : Fin n) :
    rsum (fun i => Cmat F i j) = zero := by
  show rsum (fun i => (delta i j : α) - F j i) = zero
  rw [rsum_sub (fun i => (delta i j : α)) (fun i => F j i),
      rsum_delta j, hF1 j, sub_self]

/-- **Theorem 4.2 (Unconditional Energy Conservation)**: 𝟙ᵀ(Cj) = 0 for
    EVERY total radiant power vector j — regardless of the mixed system,
    boundary conditions, or reflection-scattering coefficients (eq. 37):
    Σₖ (Cj)ₖ = Σᵢ (Σₖ Cₖᵢ)·jᵢ = Σᵢ 0·jᵢ = 0. -/
theorem thm42_energy_conservation {F : Mat α n}
    (hF1 : ∀ i, rsum (F i) = 1) (j : Fin n → α) :
    rsum (fun k => rsum (fun i => Cmat F k i * j i)) = zero := by
  rw [rsum_comm (fun k i => Cmat F k i * j i)]
  have h : ∀ i, rsum (fun k => Cmat F k i * j i) = zero := by
    intro i
    rw [rsum_mul (j i) (fun k => Cmat F k i), col_sum_zero hF1 i, zero_mul]
  rw [rsum_congr h, rsum_zero_fun]

end RingE

/-! ### Non-vacuity: ℤ satisfies all the axioms -/

instance : RingE Int where
  zero := 0
  zero_add := fun a => by omega
  add_zero := fun a => by omega
  add_add_add_comm := fun a b c d => by omega
  sub_self := fun a => by omega
  sub_zero := fun a => by omega
  sub_add_sub := fun a b c d => by omega
  zero_mul := Int.zero_mul
  add_mul := Int.add_mul

open RingE in
#check @thm42_energy_conservation

/-! ### Numerical smoke test over ℤ (Fin 2), arbitrary j-vector. -/
section SmokeTest
open RingE

def Ftest : Mat Int 2 := fun i j => if i = j then 0 else 1
def jtest : Fin 2 → Int := fun k => if k = 0 then 7 else -3

#eval decide (rsum (fun k => rsum (fun i => Cmat Ftest k i * jtest i)) = 0)
end SmokeTest
