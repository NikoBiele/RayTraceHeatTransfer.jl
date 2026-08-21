/-
  Theorem 3.1 (Steps 1–2), Bielefeld (2026), arXiv:2512.22157v2
  "Spectral Radius of Combined Absorption and Reflection-Scattering"

  Formalized over an arbitrary structure satisfying three ring axioms,
  so the result holds in any commutative ring (ℝ, ℚ, ℤ, ...).

  Step 1:  A + R = F        where A = F ∘ (1−B), R = F ∘ B, B column-constant
  Step 2:  F row-stochastic → (A+R) row-stochastic, and (A+R)·𝟙 = 𝟙.

  Step 3 of the paper (ρ(A+R)=1 via Perron–Frobenius) is NOT formalized here;
  it requires spectral theory beyond core Lean.
-/

/-- Minimal algebraic axioms needed. Any commutative ring satisfies these. -/
class RingAxioms (α : Type) extends Add α, Sub α, Mul α, OfNat α 1 where
  zero           : α
  mul_sub        : ∀ a b c : α, a * (b - c) = a * b - a * c
  mul_one        : ∀ a : α, a * 1 = a
  sub_add_cancel : ∀ a b : α, a - b + b = a

export RingAxioms (mul_sub mul_one sub_add_cancel)

variable {α : Type} [RingAxioms α] {n : Nat}

/-- Square matrices as functions. -/
abbrev Mat (α : Type) (n : Nat) := Fin n → Fin n → α

/-- Hadamard (entrywise) product. -/
def hadamard (X Y : Mat α n) : Mat α n := fun i j => X i j * Y i j

/-- Entrywise sum. -/
def matAdd (X Y : Mat α n) : Mat α n := fun i j => X i j + Y i j

/-- Column-constant reflection-scattering matrix B from vector b (paper eq. 8). -/
def Bmat (b : Fin n → α) : Mat α n := fun _ j => b j

/-- The matrix (1 − B), entrywise (paper's 𝟙 is the all-ones matrix). -/
def oneMinusB (b : Fin n → α) : Mat α n := fun _ j => (1 : α) - b j

/-- Absorption matrix A = F ∘ (1 − B)  (paper eq. 9). -/
def Amat (F : Mat α n) (b : Fin n → α) : Mat α n := hadamard F (oneMinusB b)

/-- Reflection-scattering matrix R = F ∘ B  (paper eq. 10). -/
def Rmat (F : Mat α n) (b : Fin n → α) : Mat α n := hadamard F (Bmat b)

/-- The entrywise identity behind Step 1:  f·(1−b) + f·b = f. -/
theorem entry_identity (f b : α) : f * ((1:α) - b) + f * b = f := by
  rw [mul_sub, sub_add_cancel, mul_one]

/-- **Theorem 3.1, Step 1**:  A + R = F  (paper eq. 20). -/
theorem thm31_step1 (F : Mat α n) (b : Fin n → α) :
    matAdd (Amat F b) (Rmat F b) = F := by
  funext i j
  exact entry_identity (F i j) (b j)

/-- Row sum of a vector of length n. -/
def rsum (f : Fin n → α) : α := Fin.foldr n (fun i acc => f i + acc) (RingAxioms.zero)

/-- A matrix is row stochastic when every row sums to 1. -/
def RowStochastic (F : Mat α n) : Prop := ∀ i, rsum (F i) = (1:α)

/-- **Theorem 3.1, Step 2a**: row-stochasticity transfers to A + R
    (paper eq. 21, (A+R)𝟙 = F𝟙 = 𝟙). -/
theorem thm31_step2 (F : Mat α n) (b : Fin n → α) (hF : RowStochastic F) :
    RowStochastic (matAdd (Amat F b) (Rmat F b)) := by
  rw [thm31_step1]
  exact hF

/-- Matrix-vector product: (X v) i = Σ_j X i j · v j. -/
def mulVec (X : Mat α n) (v : Fin n → α) : Fin n → α :=
  fun i => rsum (fun j => X i j * v j)

/-- The all-ones vector 𝟙. -/
def onesVec : Fin n → α := fun _ => (1:α)

/-- **Theorem 3.1, Step 2b**: 𝟙 is an eigenvector of A + R with eigenvalue 1,
    i.e. (A+R)·𝟙 = 𝟙, given F row stochastic. -/
theorem thm31_eigvec (F : Mat α n) (b : Fin n → α) (hF : RowStochastic F) :
    mulVec (matAdd (Amat F b) (Rmat F b)) (onesVec (α := α) (n := n)) = onesVec := by
  funext i
  show rsum (fun j => matAdd (Amat F b) (Rmat F b) i j * (1:α)) = (1:α)
  have h1 : (fun j => matAdd (Amat F b) (Rmat F b) i j * (1:α)) = F i := by
    funext j
    rw [mul_one, thm31_step1]
  rw [h1]
  exact hF i

/-- Sanity check: the axioms are satisfiable — ℤ is an instance
    (so the development is not vacuous). -/
instance : RingAxioms Int where
  zero           := 0
  mul_sub        := Int.mul_sub
  mul_one        := Int.mul_one
  sub_add_cancel := Int.sub_add_cancel

#check @thm31_step1
#check @thm31_step2
#check @thm31_eigvec

/-- Concrete numerical smoke test over ℤ-scaled entries is not meaningful for
    stochastic F; instead check step 1 numerically on a 2×2 integer example. -/
def Ftest : Mat Int 2 := fun i j => if i = j then 3 else 7
def btest : Fin 2 → Int := fun j => if j = 0 then 2 else 5

#eval decide (∀ i : Fin 2, ∀ j : Fin 2,
  matAdd (Amat Ftest btest) (Rmat Ftest btest) i j = Ftest i j)
