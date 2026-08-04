import sympy as sp
from sympy import symbols, Matrix, sin, cos, simplify

# --- 1. Define basic symbols ---
# Define coordinate variables (spherical coordinates)
t, r, theta, phi = symbols('t r theta phi')
# Define metric parameters (mass)
M = symbols('M')

# --- 2. Define metric tensor g_{mu,nu} ---
# Schwarzschild metric (in natural units c=G=1)
g00 = -(1 - 2*M/r)           # time-time component
g11 = 1/(1 - 2*M/r)         # radial component
g22 = r**2                  # theta-theta component
g33 = r**2 * sin(theta)**2  # phi-phi component

# Construct metric matrix
g = Matrix([
    [g00, 0, 0, 0],
    [0, g11, 0, 0],
    [0, 0, g22, 0],
    [0, 0, 0, g33]
])
coords = [t, r, theta, phi]  # Coordinate list

# --- 3. Calculate inverse metric g^{mu,nu} ---
g_inv = g.inv()

# --- 4. Calculate Christoffel symbols Gamma^{rho}_{mu,nu} ---
# Initialize 4x4x4 zero tensor
Gamma = sp.MutableDenseNDimArray.zeros(4, 4, 4)

for rho in range(4):
    for mu in range(4):
        for nu in range(4):
            # Sum over lambda: 1/2 * g^{rho,lambda} * (∂g_{lambda,mu}/∂x^{nu} 
            # + ∂g_{lambda,nu}/∂x^{mu} - ∂g_{mu,nu}/∂x^{lambda})
            sum_expr = 0
            for lam in range(4):
                term = (sp.diff(g[lam, mu], coords[nu]) +
                        sp.diff(g[lam, nu], coords[mu]) -
                        sp.diff(g[mu, nu], coords[lam]))
                sum_expr += g_inv[rho, lam] * term
            Gamma[rho, mu, nu] = simplify(sum_expr / 2)

# --- 5. Calculate Riemann curvature tensor R^{rho}_{sigma,mu,nu} ---
Riemann = sp.MutableDenseNDimArray.zeros(4, 4, 4, 4)
for rho in range(4):
    for sigma in range(4):
        for mu in range(4):
            for nu in range(4):
                # ∂Gamma^{rho}_{sigma,nu}/∂x^{mu} - ∂Gamma^{rho}_{sigma,mu}/∂x^{nu}
                # + Gamma^{rho}_{mu,lambda} Gamma^{lambda}_{sigma,nu} 
                # - Gamma^{rho}_{nu,lambda} Gamma^{lambda}_{sigma,mu}
                term1 = sp.diff(Gamma[rho, sigma, nu], coords[mu])
                term2 = sp.diff(Gamma[rho, sigma, mu], coords[nu])
                
                sum_terms = 0
                for lam in range(4):
                    sum_terms += (Gamma[rho, mu, lam] * Gamma[lam, sigma, nu] -
                                  Gamma[rho, nu, lam] * Gamma[lam, sigma, mu])
                Riemann[rho, sigma, mu, nu] = simplify(term1 - term2 + sum_terms)

# --- 6. Contract to obtain Ricci tensor R_{mu,nu} = R^{rho}_{mu,rho,nu} ---
Ricci = sp.MutableDenseNDimArray.zeros(4, 4)
for mu in range(4):
    for nu in range(4):
        Ricci[mu, nu] = simplify(sum(Riemann[rho, mu, rho, nu] for rho in range(4)))

# --- 7. Calculate Ricci scalar R = g^{mu,nu} R_{mu,nu} ---
Ricci_scalar = 0
for mu in range(4):
    for nu in range(4):
        Ricci_scalar += g_inv[mu, nu] * Ricci[mu, nu]
Ricci_scalar = simplify(Ricci_scalar)

# --- 8. Construct Einstein tensor G_{mu,nu} = R_{mu,nu} - 1/2 * g_{mu,nu} * R ---
Einstein = sp.MutableDenseNDimArray.zeros(4, 4)
for mu in range(4):
    for nu in range(4):
        Einstein[mu, nu] = simplify(Ricci[mu, nu] - sp.Rational(1, 2) * g[mu, nu] * Ricci_scalar)

# Print results to verify (should be zero for vacuum solution)
print("Einstein tensor G_00 =", Einstein[0, 0])
