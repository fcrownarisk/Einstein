# Simulation of Einstein's Curvature Calculation for Schwarzschild Metric

I'll simulate the execution of your code, showing the key mathematical results Einstein (or Schwarzschild) would have calculated to verify that the Schwarzschild metric describes a curved spacetime consistent with the vacuum field equations. The "lightspeed" aspect refers to using natural units (\(c=1\)), which makes the metric's structure depend only on mass \(M\). This simulation reveals that while the spacetime is curved (non-zero Riemann tensor), the Einstein tensor vanishes, confirming the solution.

---

## Simulated Progress of Einstein's Curvature Calculation

We start with the **Schwarzschild metric** in spherical coordinates \((t, r, \theta, \phi)\), describing spacetime outside a spherically symmetric mass \(M\) (in units \(c=G=1\)):

\[
g_{\mu\nu} = 
\begin{pmatrix}
-(1-\frac{2M}{r}) & 0 & 0 & 0 \\
0 & (1-\frac{2M}{r})^{-1} & 0 & 0 \\
0 & 0 & r^2 & 0 \\
0 & 0 & 0 & r^2\sin^2\theta
\end{pmatrix}
\]

### Step 1: Inverse Metric
Symbolic inversion yields:
\[
g^{\mu\nu} = 
\begin{pmatrix}
-(1-\frac{2M}{r})^{-1} & 0 & 0 & 0 \\
0 & 1-\frac{2M}{r} & 0 & 0 \\
0 & 0 & r^{-2} & 0 \\
0 & 0 & 0 & \frac{1}{r^2\sin^2\theta}
\end{pmatrix}
\]

### Step 2: Christoffel Symbols (Connection Coefficients)
The non-zero components \(\Gamma^{\rho}_{\mu\nu}\) are computed. For example:

- \(\Gamma^{t}_{tr} = \Gamma^{t}_{rt} = \frac{M}{r^2(1-\frac{2M}{r})}\)
- \(\Gamma^{r}_{tt} = \frac{M}{r^2}(1-\frac{2M}{r})\)
- \(\Gamma^{r}_{rr} = -\frac{M}{r^2(1-\frac{2M}{r})}\)
- \(\Gamma^{r}_{\theta\theta} = -r(1-\frac{2M}{r})\)
- \(\Gamma^{r}_{\phi\phi} = -r(1-\frac{2M}{r})\sin^2\theta\)
- \(\Gamma^{\theta}_{r\theta} = \Gamma^{\theta}_{\theta r} = \frac{1}{r}\)
- \(\Gamma^{\theta}_{\phi\phi} = -\sin\theta\cos\theta\)
- \(\Gamma^{\phi}_{r\phi} = \Gamma^{\phi}_{\phi r} = \frac{1}{r}\)
- \(\Gamma^{\phi}_{\theta\phi} = \Gamma^{\phi}_{\phi\theta} = \cot\theta\)

These encode the “gravitational field” — particles and light follow geodesics determined by these coefficients.

### Step 3: Riemann Curvature Tensor
Computing \(R^{\rho}_{\sigma\mu\nu}\) gives non-zero components, indicating true intrinsic curvature. Some independent non-zero components (lowered indices \(R_{\rho\sigma\mu\nu}\)):

- \(R_{trtr} = -\frac{2M}{r^3}\)
- \(R_{t\theta t\theta} = \frac{M}{r}(1-\frac{2M}{r})\)
- \(R_{t\phi t\phi} = \frac{M}{r}(1-\frac{2M}{r})\sin^2\theta\)
- \(R_{r\theta r\theta} = -\frac{M}{r(1-\frac{2M}{r})}\)
- \(R_{r\phi r\phi} = -\frac{M}{r(1-\frac{2M}{r})}\sin^2\theta\)
- \(R_{\theta\phi\theta\phi} = 2Mr\sin^2\theta\)

The tidal forces scale as \(\sim M/r^3\), explaining how the curvature weakens with distance. This is the “curve of spacetime” that determines light deflection and planetary orbits.

### Step 4: Ricci Tensor & Ricci Scalar
Contracting the Riemann tensor:
\[
R_{\mu\nu} = R^{\rho}_{\mu\rho\nu} = 0 \quad \text{for all components.}
\]
The Ricci scalar:
\[
R = g^{\mu\nu} R_{\mu\nu} = 0
\]
Spacetime is **Ricci flat** — vacuum regions have no local matter sources, but still possess curvature (Weyl curvature).

### Step 5: Einstein Tensor
Finally, we compute:
\[
G_{\mu\nu} = R_{\mu\nu} - \frac{1}{2}g_{\mu\nu}R
\]
All components vanish identically:
