/-
  Theorem 3.2, Bielefeld (2026), arXiv:2512.22157v2
  "Upper Bound on Spectral Radius of the Reflection-Scattering Matrix"

  Self-contained core-Lean formalization (no Mathlib), in the style of Thm31.lean:
  all results proven over an abstract ordered ring with absolute value, whose
  axioms are satisfied by ℝ (and, as the non-vacuity instance below shows, by ℤ).

  Formalized content:
  * `rsum` lemmas          — finite-sum toolkit (monotonicity, factoring, triangle ineq.)
  * `eigenvalue_abs_le`    — general bound: row abs-sums ≤ c  ⟹  any eigenvalue |μ| ≤ c
                             (the max-coordinate argument behind ρ(M) ≤ ‖M‖∞)
  * `Rmat_row_abs_le`      — rows of R = F∘B satisfy Σⱼ|Rᵢⱼ| ≤ bmax   (inside eq. 22)
  * `thm32_eigen_bound`    — part (a), real form: eigenvalues of R have |μ| ≤ bmax
  * `thm32_uniform`        — part (b): b ≡ bmax  ⟹  R·𝟙 = bmax·𝟙 (bound attained)
  * `thm32_no_unit_eigen`  — part (c), eigen-content: bmax < 1 ⟹ R has no fixed
                             vector; the step "hence I − Rᵀ invertible" is classical
                             linear algebra (det theory), left to Mathlib.
-/

/-- Axioms of a linearly ordered commutative ring with absolute value.
    ℝ satisfies every field below. -/
class OrdRing (α : Type) extends Add α, Mul α, LE α, LT α, OfNat α 1 where
  zero : α
  abs : α → α
  le_refl : ∀ a : α, a ≤ a
  le_trans : ∀ {a b c : α}, a ≤ b → b ≤ c → a ≤ c
  le_total : ∀ a b : α, a ≤ b ∨ b ≤ a
  lt_of_lt_of_le : ∀ {a b c : α}, a < b → b ≤ c → a < c
  not_le_of_lt : ∀ {a b : α}, a < b → ¬ b ≤ a
  add_le_add : ∀ {a b c d : α}, a ≤ b → c ≤ d → a + c ≤ b + d
  add_mul : ∀ a b c : α, (a + b) * c = a * c + b * c
  zero_mul : ∀ a : α, zero * a = zero
  one_mul : ∀ a : α, 1 * a = a
  mul_one : ∀ a : α, a * 1 = a
  mul_nonneg : ∀ {a b : α}, zero ≤ a → zero ≤ b → zero ≤ a * b
  mul_le_mul_of_nonneg_left : ∀ {a b c : α}, b ≤ c → zero ≤ a → a * b ≤ a * c
  mul_le_mul_of_nonneg_right : ∀ {a b c : α}, a ≤ b → zero ≤ c → a * c ≤ b * c
  le_of_mul_le_mul_right : ∀ {a b c : α}, a * c ≤ b * c → zero < c → a ≤ b
  zero_le_one : zero ≤ (1 : α)
  abs_nonneg : ∀ a : α, zero ≤ abs a
  abs_of_nonneg : ∀ {a : α}, zero ≤ a → abs a = a
  abs_mul : ∀ a b : α, abs (a * b) = abs a * abs b
  abs_add_le : ∀ a b : α, abs (a + b) ≤ abs a + abs b
  abs_pos_of_ne_zero : ∀ {a : α}, a ≠ zero → zero < abs a

namespace OrdRing

variable {α : Type} [OrdRing α] {n : Nat}

/-- Square matrices as functions. -/
abbrev Mat (α : Type) (n : Nat) := Fin n → Fin n → α

/-- Row sum, by structural recursion (no library dependencies). -/
def rsum : {n : Nat} → (Fin n → α) → α
  | 0, _ => zero
  | _+1, f => f 0 + rsum (fun i => f i.succ)

theorem rsum_congr {f g : Fin n → α} (h : ∀ j, f j = g j) : rsum f = rsum g :=
  congrArg rsum (funext h)

theorem rsum_le_rsum : ∀ {n : Nat} {f g : Fin n → α},
    (∀ j, f j ≤ g j) → rsum f ≤ rsum g
  | 0, _, _, _ => le_refl zero
  | _+1, _, _, h => add_le_add (h 0) (rsum_le_rsum (fun j => h j.succ))

theorem rsum_mul (c : α) : ∀ {n : Nat} (f : Fin n → α),
    rsum (fun j => f j * c) = rsum f * c
  | 0, _ => (zero_mul c).symm
  | _+1, f => by
    show f 0 * c + rsum (fun j => f j.succ * c)
       = (f 0 + rsum (fun j => f j.succ)) * c
    rw [rsum_mul c (fun j => f j.succ), add_mul]

theorem abs_rsum_le : ∀ {n : Nat} (f : Fin n → α),
    abs (rsum f) ≤ rsum (fun j => abs (f j))
  | 0, _ => by
    show abs zero ≤ zero
    rw [abs_of_nonneg (le_refl zero)]
    exact le_refl zero
  | _+1, f => by
    show abs (f 0 + rsum (fun j => f j.succ))
       ≤ abs (f 0) + rsum (fun j => abs (f j.succ))
    exact le_trans (abs_add_le _ _) (add_le_add (le_refl _) (abs_rsum_le _))

/-- A finite nonempty index set has an argmax. -/
theorem exists_max : ∀ {m : Nat} (f : Fin (m+1) → α), ∃ i, ∀ j, f j ≤ f i := by
  intro m
  induction m with
  | zero =>
    intro f
    refine ⟨0, fun j => ?_⟩
    have hj : j = 0 := Fin.ext (Nat.lt_one_iff.mp j.isLt)
    rw [hj]; exact le_refl _
  | succ m ih =>
    intro f
    obtain ⟨i', hi'⟩ := ih (fun k => f k.succ)
    rcases le_total (f 0) (f i'.succ) with h | h
    · refine ⟨i'.succ, fun j => ?_⟩
      cases j using Fin.cases with
      | zero => exact h
      | succ k => exact hi' k
    · refine ⟨0, fun j => ?_⟩
      cases j using Fin.cases with
      | zero => exact le_refl _
      | succ k => exact le_trans (hi' k) h

/-- **General eigenvalue bound** (the argument behind ρ(M) ≤ ‖M‖∞, Thm 3.2 proof):
    if every row satisfies Σⱼ|Mᵢⱼ| ≤ c and M·v = μ·v with v ≠ 0, then |μ| ≤ c. -/
theorem eigenvalue_abs_le {M : Mat α n} {c : α}
    (hM : ∀ i, rsum (fun j => abs (M i j)) ≤ c)
    {μ : α} {v : Fin n → α} (hv : ∃ k, v k ≠ zero)
    (heig : ∀ i, rsum (fun j => M i j * v j) = μ * v i) :
    abs μ ≤ c := by
  cases n with
  | zero => obtain ⟨k, _⟩ := hv; exact k.elim0
  | succ m =>
    obtain ⟨i, hi⟩ := exists_max (fun k => abs (v k))
    obtain ⟨k, hk⟩ := hv
    have hvi : zero < abs (v i) :=
      lt_of_lt_of_le (abs_pos_of_ne_zero hk) (hi k)
    -- |μ|·|v i| = |Σⱼ Mᵢⱼ vⱼ| ≤ Σⱼ|Mᵢⱼ||vⱼ| ≤ (Σⱼ|Mᵢⱼ|)·|v i| ≤ c·|v i|
    have h1 : abs μ * abs (v i) = abs (rsum (fun j => M i j * v j)) := by
      rw [heig i, abs_mul]
    have h2 : abs (rsum (fun j => M i j * v j))
            ≤ rsum (fun j => abs (M i j) * abs (v j)) := by
      have := abs_rsum_le (fun j => M i j * v j)
      rwa [rsum_congr (fun j => abs_mul (M i j) (v j))] at this
    have h3 : rsum (fun j => abs (M i j) * abs (v j))
            ≤ rsum (fun j => abs (M i j) * abs (v i)) :=
      rsum_le_rsum (fun j => mul_le_mul_of_nonneg_left (hi j) (abs_nonneg _))
    have h4 : rsum (fun j => abs (M i j) * abs (v i))
            = rsum (fun j => abs (M i j)) * abs (v i) :=
      rsum_mul (abs (v i)) (fun j => abs (M i j))
    have h5 : rsum (fun j => abs (M i j)) * abs (v i) ≤ c * abs (v i) :=
      mul_le_mul_of_nonneg_right (hM i) (abs_nonneg _)
    have key : abs μ * abs (v i) ≤ c * abs (v i) := by
      rw [h1]
      exact le_trans h2 (le_trans h3 (h4 ▸ h5))
    exact le_of_mul_le_mul_right key hvi

/-! ### The paper's setting: R = F ∘ B, B column-constant (eqs. 8, 10) -/

/-- Reflection-scattering matrix, entrywise Rᵢⱼ = Fᵢⱼ·bⱼ (eq. 10 with eq. 8). -/
def Rmat (F : Mat α n) (b : Fin n → α) : Mat α n := fun i j => F i j * b j

variable {F : Mat α n} {b : Fin n → α} {bmax : α}

/-- Row abs-sums of R are bounded by bmax (the computation inside eq. 22):
    Σⱼ bⱼFᵢⱼ ≤ bmax·ΣⱼFᵢⱼ = bmax. -/
theorem Rmat_row_abs_le
    (hF0 : ∀ i j, zero ≤ F i j) (hF1 : ∀ i, rsum (F i) = 1)
    (hb0 : ∀ j, zero ≤ b j) (hbm : ∀ j, b j ≤ bmax) (i : Fin n) :
    rsum (fun j => abs (Rmat F b i j)) ≤ bmax := by
  have h1 : rsum (fun j => abs (Rmat F b i j)) = rsum (fun j => F i j * b j) :=
    rsum_congr fun j => abs_of_nonneg (mul_nonneg (hF0 i j) (hb0 j))
  have h2 : rsum (fun j => F i j * b j) ≤ rsum (fun j => F i j * bmax) :=
    rsum_le_rsum fun j => mul_le_mul_of_nonneg_left (hbm j) (hF0 i j)
  have h3 : rsum (fun j => F i j * bmax) = rsum (F i) * bmax :=
    rsum_mul bmax (F i)
  rw [h1]
  rw [h3, hF1 i, one_mul] at h2
  exact h2

/-- **Theorem 3.2, part (a)** (real form): any real eigenvalue μ of R
    satisfies |μ| ≤ bmax.  (This is ρ(R) ≤ bmax restricted to real
    eigenvalues; the complex case needs ℂ, i.e. Mathlib.) -/
theorem thm32_eigen_bound
    (hF0 : ∀ i j, zero ≤ F i j) (hF1 : ∀ i, rsum (F i) = 1)
    (hb0 : ∀ j, zero ≤ b j) (hbm : ∀ j, b j ≤ bmax)
    {μ : α} {v : Fin n → α} (hv : ∃ k, v k ≠ zero)
    (heig : ∀ i, rsum (fun j => Rmat F b i j * v j) = μ * v i) :
    abs μ ≤ bmax :=
  eigenvalue_abs_le (fun i => Rmat_row_abs_le hF0 hF1 hb0 hbm i) hv heig

/-- **Theorem 3.2, part (b)**: with uniform coefficients b ≡ bmax the bound is
    attained — 𝟙 is an eigenvector of R with eigenvalue bmax. -/
theorem thm32_uniform
    (hF1 : ∀ i, rsum (F i) = 1) (hb : ∀ j, b j = bmax) (i : Fin n) :
    rsum (fun j => Rmat F b i j * 1) = bmax * 1 := by
  have h : (fun j => Rmat F b i j * 1) = fun j => F i j * bmax := by
    funext j
    show F i j * b j * 1 = F i j * bmax
    rw [mul_one, hb j]
  rw [h, rsum_mul bmax (F i), hF1 i, one_mul, mul_one]

/-- **Theorem 3.2, part (c)** (eigen-content): if bmax < 1 then R has no
    fixed vector, i.e. 1 is not an eigenvalue of R.  The classical
    linear-algebra step "hence D = I − Rᵀ is invertible" (via
    det(I−Rᵀ) = det(I−R)) is left to Mathlib. -/
theorem thm32_no_unit_eigen
    (hF0 : ∀ i j, zero ≤ F i j) (hF1 : ∀ i, rsum (F i) = 1)
    (hb0 : ∀ j, zero ≤ b j) (hbm : ∀ j, b j ≤ bmax) (hlt : bmax < 1) :
    ¬ ∃ v : Fin n → α, (∃ k, v k ≠ zero) ∧
        (∀ i, rsum (fun j => Rmat F b i j * v j) = v i) := by
  intro h
  obtain ⟨v, hv, hfix⟩ := h
  have heig : ∀ i, rsum (fun j => Rmat F b i j * v j) = 1 * v i := fun i => by
    rw [hfix i, one_mul]
  have hle : abs (1 : α) ≤ bmax :=
    eigenvalue_abs_le (fun i => Rmat_row_abs_le hF0 hF1 hb0 hbm i) hv heig
  rw [abs_of_nonneg zero_le_one] at hle
  exact not_le_of_lt hlt hle

end OrdRing

/-! ### Non-vacuity: ℤ satisfies all the axioms -/

instance : OrdRing Int where
  zero := 0
  abs := fun a => (a.natAbs : Int)
  le_refl := fun a => by omega
  le_trans := fun h1 h2 => by omega
  le_total := fun a b => by omega
  lt_of_lt_of_le := fun h1 h2 => by omega
  not_le_of_lt := fun h => by omega
  add_le_add := fun h1 h2 => by omega
  add_mul := Int.add_mul
  zero_mul := Int.zero_mul
  one_mul := Int.one_mul
  mul_one := Int.mul_one
  mul_nonneg := Int.mul_nonneg
  mul_le_mul_of_nonneg_left := fun h ha => Int.mul_le_mul_of_nonneg_left h ha
  mul_le_mul_of_nonneg_right := fun h hc => Int.mul_le_mul_of_nonneg_right h hc
  le_of_mul_le_mul_right := fun h hc => Int.le_of_mul_le_mul_right h hc
  zero_le_one := by decide
  abs_nonneg := fun a => by omega
  abs_of_nonneg := fun h => by omega
  abs_mul := fun a b => by
    show ((a * b).natAbs : Int) = (a.natAbs : Int) * (b.natAbs : Int)
    rw [Int.natAbs_mul]
    omega
  abs_add_le := fun a b => by
    have := Int.natAbs_add_le a b
    omega
  abs_pos_of_ne_zero := fun h => by omega

open OrdRing in
#check @eigenvalue_abs_le
open OrdRing in
#check @thm32_eigen_bound
open OrdRing in
#check @thm32_uniform
open OrdRing in
#check @thm32_no_unit_eigen

/-! ### Numerical smoke test over ℤ (Fin 2): F = antidiagonal, uniform b = 3.
    Part (b): R·𝟙 should equal 3·𝟙, and row abs-sums should equal 3. -/
section SmokeTest
open OrdRing

def Ftest : Mat Int 2 := fun i j => if i = j then 0 else 1
def btest : Fin 2 → Int := fun _ => 3

#eval decide (∀ i : Fin 2, rsum (fun j => Rmat Ftest btest i j * 1) = 3 * 1)
#eval decide (∀ i : Fin 2, rsum (fun j => OrdRing.abs (Rmat Ftest btest i j)) = 3)
end SmokeTest
