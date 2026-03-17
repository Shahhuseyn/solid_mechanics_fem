# Solid Mechanics FEM — 1D Linear Finite Element Analysis

**Course:** Materials Simulation Practical | FAU Erlangen-Nürnberg  
**Tools:** Python · NumPy · Matplotlib

💻 [Notebook](https://github.com/Shahhuseyn/solid_mechanics_fem/blob/main/fem_1d_elastic_bar.ipynb)

---

## Overview

From-scratch implementation of the Finite Element Method for 1D linear elasticity in Python/NumPy. The project covers the full pipeline: constitutive theory and symmetry reduction of the stiffness tensor, weak-form derivation of the governing equations, linear shape functions and their derivatives, Newton–Cotes numerical integration, global stiffness matrix assembly, enforcement of Dirichlet and Neumann boundary conditions via a permutation-matrix reduction scheme, and element-wise stress recovery across four distinct loading cases.

---

## Part 1 — Theory: Constitutive Relations & Weak Form (Tasks 1–3)

Starting from the local balance of linear momentum in the static case, the strong form is reduced to:

$$\text{div}(\boldsymbol{\sigma}) = -\rho \boldsymbol{b}$$

Combined with Hooke's law σ = C:ε and the infinitesimal strain tensor, this yields the 1D governing equation. The full symmetry analysis of the stiffness tensor C_ijkl is carried out — minor symmetries (from stress/strain symmetry) reduce the 81 independent components to 36, and major symmetry (from the strain energy density) further reduces them to 21. The result is cast in Voigt notation. The 1D weak form is then derived by multiplying by a test function and integrating by parts:

$$\int_0^l EA \frac{du}{dx} \frac{d(\delta u)}{dx} \, dx = \int_0^l \rho b \, \delta u \, dx + \bar{t} A \, \delta u \Big|_{\Gamma_t}$$

---

## Part 2 — Shape Functions & Interpolation (Tasks 4.1–4.2)

Linear Lagrange shape functions N₁ and N₂ are defined over each element and plotted as overlapping hat functions across the 6-element mesh.

### Linear Shape Functions

![Shape functions](figures/shape_functions_linear.png)

The interpolation accuracy is tested on three functions — linear f(x)=2x+3, sinusoidal f(x)=sin(x), and quadratic f(x)=x² — evaluated on element 2 (nodes 2–3).

### Interpolation: Linear Function (element 2)

![Interpolation linear](figures/interpolation_linear.png)

The linear function is reproduced exactly — the interpolated curve overlaps the exact solution with zero error, confirming that linear shape functions exactly represent any function in the same polynomial space.

### Interpolation: Sinusoidal Function (element 2)

![Interpolation sinusoidal](figures/interpolation_sin.png)

For sin(x), the linear interpolation underestimates the curved function, producing a visible but bounded error that peaks near the midpoint of the element.

### Interpolation: Quadratic Function (element 2)

![Interpolation quadratic](figures/interpolation_quadratic.png)

Similarly, x² is only approximated — the interpolated line overshoots near the element endpoints and undershoots in the middle. The error profile is symmetric and bounded.

The same study is repeated on element 5 (nodes 5–6), demonstrating that the qualitative behavior is preserved regardless of element position while the magnitudes shift with local function curvature.

### Interpolation: Linear Function (element 5)

![Interpolation linear 2](figures/interpolation_linear_2.png)

### Interpolation: Sinusoidal Function (element 5)

![Interpolation sinusoidal 2](figures/interpolation_sin_2.png)

### Interpolation: Quadratic Function (element 5)

![Interpolation quadratic 2](figures/interpolation_quadratic_2.png)

The strain matrix **B** (derivative of shape functions) is constant over each linear element, meaning the stress σ = E·B·û is elementwise constant — stress discontinuities across element boundaries are therefore expected.

---

## Part 3 — Numerical Integration (Tasks 5.1–5.3)

Three Newton–Cotes quadrature rules are implemented and tested:

| Rule | Points | Exact for degree |
|---|---|---|
| Trapezoidal | 2 | ≤ 1 |
| Simpson's | 3 | ≤ 2 |
| 3/8 Rule | 4 | ≤ 3 |

The element stiffness matrix for a 1D bar (E = 210,000 N/mm², A = 25 mm², L = 50 mm) is computed via trapezoidal quadrature and verified analytically. The displacement at the free end matches the analytical solution u(L) = FL/(EA) with negligible error.

---

## Part 4 — Global Assembly & Boundary Conditions (Tasks 6–7)

An element assembly routine loops over all elements, maps local stiffness matrices into the global system via node connectivity, and produces the symmetric tridiagonal global stiffness matrix K. The normalized matrix matches the classic 1D spring chain structure.

Dirichlet boundary conditions are enforced using a permutation-matrix reduction: a matrix P of shape (n−m)×n maps the full displacement vector to the reduced free-DOF subspace, solves the reduced system K_r·d_r = f_r, and reconstructs the full solution via back-substitution.

---

## Part 5 — Stress Recovery & Loading Cases (Task 8)

Element-wise stress is recovered from σ = E·B·û. Four loading cases are solved on the 6-element, 7-node bar (E = 210,000 N/mm², A = 25 mm², L = 50 mm):

### Case 1 — Pure Neumann (tip load F = 5 N)

![Displacement Neumann](figures/displacement_neumann.png)
![Stress Neumann](figures/stress_neumann.png)

Linear displacement field and uniform stress σ = 0.2 N/mm² across all elements, consistent with the analytical solution for a bar under axial end load with no body forces.

### Case 2 — Double Dirichlet (u(0) = 0, u(node 3) = 0.01 mm)

![Displacement double Dirichlet](figures/displacement_dirichlet_double.png)
![Stress double Dirichlet](figures/stress_dirichlet_double.png)

Displacement ramps linearly up to the prescribed midpoint value and plateaus thereafter. The stress is elevated in the constrained region and drops sharply in the free segment — physically consistent with the imposed displacement gradient.

### Case 3 — Internal + End Load (F_node3 = 30 N, F_end = 40 N)

![Displacement internal load](figures/displacement_internal_load.png)
![Stress internal load](figures/stress_internal_load.png)

A concentrated internal load at node 3 produces a stress jump at that location. Elements to the left carry the combined load while elements to the right carry only the end load, resulting in a piecewise-constant stress profile with a clear discontinuity.

### Case 4 — Combined BCs (u(0) = 0, u(6) = 0.0015 mm, internal loads)

![Displacement combined](figures/displacement_combined_bc.png)
![Stress combined](figures/stress_combined_bc.png)

With both ends constrained and internal loads of opposite sign applied, the displacement profile is nonlinear and the stress field exhibits a pronounced jump at the internal load node, with the stress level reversing across the bar.

---

## Requirements

```
pip install numpy matplotlib
```

---

## Usage

Open `fem_1d_elastic_bar.ipynb` in Jupyter and run cells sequentially. All functions — shape functions, B-matrix, quadrature, assembly, BC enforcement, and stress recovery — are defined incrementally and reused across tasks.
