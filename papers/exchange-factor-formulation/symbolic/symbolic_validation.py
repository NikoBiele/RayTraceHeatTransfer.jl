"""
Symbolic validation of the exchange factor formulation against the textbook
closed-form solutions for two infinite gray-diffuse enclosures:
infinite parallel plates and infinite concentric cylinders.

Reproduces Section 6.1 (6.1.1 and 6.1.2) and Appendix A of
"A Sparsity-Preserving Radiation Exchange Factor Formulation with Proven
Non-Negativity, Non-Singularity and Energy Conservation".

Requires only SymPy. Every check is an exact symbolic identity, not a
numerical comparison; the script stops with an AssertionError if one fails.
"""

from sympy import Matrix, symbols, simplify, eye


def build_matrices(F, B):
    """Build A, R, C, D from the exchange factor matrix F and the
    column-constant reflection-scattering matrix B (single-step formulation)."""
    n = F.shape[0]                                                    # number of elements
    A = Matrix(n, n, lambda i, j: (1 - B[i, j]) * F[i, j])            # absorption matrix, A = F ∘ (1 − B)
    R = Matrix(n, n, lambda i, j: B[i, j] * F[i, j])                  # reflection-scattering matrix, R = F ∘ B
    C = simplify(eye(n) - A.T - R.T)                                  # source-term matrix, C = I − Aᵀ − Rᵀ
    D = simplify(eye(n) - R.T)                                        # emissive-power matrix, D = I − Rᵀ
    return simplify(A), simplify(R), C, D


def parallel_plates():
    """Two infinite parallel plates of equal area A and equal reflectivity rho."""
    rho, sigma, T1, T2, A = symbols("rho sigma T_1 T_2 A", positive=True)   # reflectivity, Stefan-Boltzmann constant, temperatures, area
    eps = 1 - rho                                                     # emissivity of an opaque surface
    Eb1, Eb2 = sigma * T1**4, sigma * T2**4                           # blackbody emissive powers [W/m²]
    e1, e2 = eps * Eb1 * A, eps * Eb2 * A                             # emitted powers [W], the prescribed right-hand side

    F = Matrix([[0, 1], [1, 0]])                                      # each plate sees only the other
    B = Matrix([[rho, rho], [rho, rho]])                              # column-constant reflectivities
    A_mat, R_mat, C_mat, D_mat = build_matrices(F, B)                 # system matrices

    j = simplify(D_mat.solve(Matrix([e1, e2])))                       # both temperatures prescribed: M = D, solve D j = e
    g_a = simplify(A_mat.T * j)                                       # absorbed powers
    q1 = simplify(e1 - g_a[0])                                        # source term of plate 1
    q2 = simplify(e2 - g_a[1])                                        # source term of plate 2

    q_textbook = (Eb1 - Eb2) * A / (1 / eps + 1 / eps - 1)            # textbook net exchange between parallel plates
    assert simplify(q1 - q_textbook) == 0                             # exact symbolic match with the textbook formula
    assert simplify(q1 + q2) == 0                                     # energy conservation, q1 + q2 = 0
    print("Parallel plates:      q1 matches the textbook formula exactly; q1 + q2 = 0")


def concentric_cylinders():
    """Two infinite concentric cylinders (or spheres) with area ratio alpha = A1/A2 < 1."""
    rho, sigma, T1, T2, A1, alpha = symbols("rho sigma T_1 T_2 A_1 alpha", positive=True)   # as above, inner area A1 and area ratio
    eps = 1 - rho                                                     # emissivity of an opaque surface
    A2 = A1 / alpha                                                   # outer area
    Eb1, Eb2 = sigma * T1**4, sigma * T2**4                           # blackbody emissive powers [W/m²]
    e1, e2 = eps * Eb1 * A1, eps * Eb2 * A2                           # emitted powers [W], the prescribed right-hand side

    F = Matrix([[0, 1], [alpha, 1 - alpha]])                          # inner sees only outer; outer sees inner and itself
    B = Matrix([[rho, rho], [rho, rho]])                              # column-constant reflectivities
    A_mat, R_mat, C_mat, D_mat = build_matrices(F, B)                 # system matrices

    j = simplify(D_mat.solve(Matrix([e1, e2])))                       # both temperatures prescribed: M = D, solve D j = e
    g_a = simplify(A_mat.T * j)                                       # absorbed powers
    q1 = simplify(e1 - g_a[0])                                        # source term of the inner cylinder
    q2 = simplify(e2 - g_a[1])                                        # source term of the outer cylinder

    q_textbook = A1 * sigma * (T1**4 - T2**4) / (1 / eps + alpha * (1 / eps - 1))   # textbook net exchange between concentric cylinders
    assert simplify(q1 - q_textbook) == 0                             # exact symbolic match with the textbook formula
    assert simplify(q1 + q2) == 0                                     # energy conservation, q1 + q2 = 0
    print("Concentric cylinders: q1 matches the textbook formula exactly; q1 + q2 = 0")


if __name__ == "__main__":
    parallel_plates()                                                 # Section 6.1.1
    concentric_cylinders()                                            # Section 6.1.2