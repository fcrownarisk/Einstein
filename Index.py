import sympy as sp
from sympy import symbols, Matrix, sin, cos, simplify, IndexedBase, Idx

# --- 1. Define basic symbols ---
t, r, theta, phi = symbols('t r theta phi')
M = symbols('M')
coords = [t, r, theta, phi]

# Define index ranges and dummy indices
i, j, k, l, m, n = sp.symbols('i j k l m n', cls=Idx)
# We'll use these dummy indices for display; the actual arrays are still computed with integer loops.

# Create indexed objects that display with Greek letters
Gamma_up = IndexedBase('Γ')          # Γ^i_j_k  (upper index is first)
Riemann_tensor = IndexedBase('R')    # R^i_j_k_l
Ricci_tensor = IndexedBase('R_{μν}') # Ricci tensor (we'll just display as 'R' with two indices)
Ricci_scalar_sym = symbols('R')      # Ricci scalar R
G_tensor = IndexedBase('G')          # Einstein tensor G_{μν}

# --- 2. Schwarzschild metric ---
g00 = -(1 - 2*M/r)
g11 = 1/(1 - 2*M/r)
g22 = r**2
g33 = r**2 * sin(theta)**2

g = Matrix([
    [g00, 0, 0, 0],
    [0, g11, 0, 0],
    [0, 0, g22, 0],
    [0, 0, 0, g33]
])

# --- 3. Inverse metric ---
g_inv = g.inv()

# --- 4. Christoffel symbols (stored as array) ---
Gamma_arr = sp.MutableDenseNDimArray.zeros(4, 4, 4)
for rho in range(4):
    for mu in range(4):
        for nu in range(4):
            sum_expr = 0
            for lam in range(4):
                term = (sp.diff(g[lam, mu], coords[nu]) +
                        sp.diff(g[lam, nu], coords[mu]) -
                        sp.diff(g[mu, nu], coords[lam]))
                sum_expr += g_inv[rho, lam] * term
            Gamma_arr[rho, mu, nu] = simplify(sum_expr / 2)

# --- 5. Riemann tensor ---
Riemann_arr = sp.MutableDenseNDimArray.zeros(4, 4, 4, 4)
for rho in range(4):
    for sigma in range(4):
        for mu in range(4):
            for nu in range(4):
                term1 = sp.diff(Gamma_arr[rho, sigma, nu], coords[mu])
                term2 = sp.diff(Gamma_arr[rho, sigma, mu], coords[nu])
                sum_terms = 0
                for lam in range(4):
                    sum_terms += (Gamma_arr[rho, mu, lam] * Gamma_arr[lam, sigma, nu] -
                                  Gamma_arr[rho, nu, lam] * Gamma_arr[lam, sigma, mu])
                Riemann_arr[rho, sigma, mu, nu] = simplify(term1 - term2 + sum_terms)

# --- 6. Ricci tensor ---
Ricci_arr = sp.MutableDenseNDimArray.zeros(4, 4)
for mu in range(4):
    for nu in range(4):
        Ricci_arr[mu, nu] = simplify(sum(Riemann_arr[rho, mu, rho, nu] for rho in range(4)))

# --- 7. Ricci scalar ---
R_scalar = sum(g_inv[mu, nu] * Ricci_arr[mu, nu] for mu in range(4) for nu in range(4))
R_scalar = simplify(R_scalar)

# --- 8. Einstein tensor ---
Einstein_arr = sp.MutableDenseNDimArray.zeros(4, 4)
for mu in range(4):
    for nu in range(4):
        Einstein_arr[mu, nu] = simplify(Ricci_arr[mu, nu] - sp.Rational(1,2)*g[mu, nu]*R_scalar)

# --- 9. Display results using Greek letter symbols ---
# We'll print non-zero components with indices in LaTeX-like style.
print("Non-zero Christoffel symbols Γ^ρ_μν:")
for rho in range(4):
    for mu in range(4):
        for nu in range(mu, 4):  # symmetry
            val = Gamma_arr[rho, mu, nu]
            if val != 0:
                display(sp.Eq(Gamma_up[rho, mu, nu], val))

print("\nNon-zero Riemann tensor components R^ρ_σμν (upper first index):")
for rho in range(4):
    for sigma in range(4):
        for mu in range(4):
            for nu in range(mu+1, 4):  # antisymmetry
                val = Riemann_arr[rho, sigma, mu, nu]
                if val != 0:
                    display(sp.Eq(Riemann_tensor[rho, sigma, mu, nu], val))

print("\nRicci tensor R_μν (all zero):")
for mu in range(4):
    for nu in range(mu, 4):
        if Ricci_arr[mu, nu] != 0:
            display(sp.Eq(Ricci_tensor[mu, nu], Ricci_arr[mu, nu]))
        else:
            print(f"R_{mu}{nu} = 0")

print("\nRicci scalar R =", R_scalar)

print("\nEinstein tensor G_μν (all zero):")
for mu in range(4):
    for nu in range(mu, 4):
        if Einstein_arr[mu, nu] != 0:
            display(sp.Eq(G_tensor[mu, nu], Einstein_arr[mu, nu]))
        else:
            print(f"G_{mu}{nu} = 0")
