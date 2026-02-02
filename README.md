# General algorithm of the Embedded Runge—Kutta method
A high-precision, adaptive step-size numerical integrator for Ordinary Differential Equations (ODEs) in MATLAB. 
This solver implements various embedded Runge-Kutta pairs, ranging from order 5 to 9.

## Implemented Methods
The solver includes a comprehensive suite of methods, primarily based on the work of [**Jim Verner**](https://www.sfu.ca/~jverner/), optimized for different tolerance levels and problem stiffness.

### Order 5
* `"Tsit5(4)7*"`: Tsitouras 5(4) pair. A modern improvement over Dormand-Prince, featuring minimized truncation error coefficients for higher efficiency.
* `"RK5(4)7*M"`: Dormand-Prince 5(4) pair. The standard method used in MATLAB's `ode45`.

### Order 6 (Verner IIIXb Series)
* `"IIIXb+6(5)9*"`: **Robust** variant. 
* `"IIIXb-6(5)9*"`: **Efficient** variant.

### Order 7 (Verner IIa1 Series)
* `"IIa1+7(6)10"`: **Robust** variant.
* `"IIa1-7(6)10"`: **Efficient** variant.
* `"IIa1-7(6)10M"`: **An even more efficient** variant

### Order 8 (Verner IIa Series)
* `"IIa+8(7)13"`: **Robust** variant. 
* `"IIa-8(7)13"`: **Efficient** variant. 

### Order 9 (Verner IIa Series)
* `"IIa+9(8)16"`: **Robust** variant.
* `"IIa-9(8)16"`: **Efficient** variant.

### Naming Convention Legend
* **`*` (FSAL):** First Same As Last. The last stage of the current step is reused as the first stage of the next step, saving 1 function evaluation per step.
* **`+` (Robust):** Designed with "conservative" error estimators. Better for problems with sharp turns or slight stiffness.
* **`-` (Efficient):** Designed for maximum speed per step. Best for smooth, well-behaved problems.

## Table of Contents
- [Embedded Runge—Kutta methods. Butcher tableau](#embedded-rungekutta-methods-butcher-tableau)
- [Description of the implemented algorithm](#description-of-the-implemented-algorithm)
- [Example](#example)
- [Notes](#notes)
  - [Syntax](#syntax)
  - [Input Arguments](#input-arguments)
  - [Output Arguments](#output-arguments)
- [Planned Features](#planned-features)
- [References](#references)

## Embedded Runge—Kutta methods. Butcher tableau
Let an initial value problem be specified as follows:
```math
\dot{\mathbf{z}}\left(t\right)=\mathbf{f}\left(t,\mathbf{z}\right),\quad t \in \left[t_0,t_\text{end}\right],\quad \mathbf{z}\left(t_0\right) = \mathbf{z}_0, 
```
where $`\mathbf{z}\left(t\right): \mathbb{R}\mapsto\mathbb{R}^m, \mathbf{f}\left(t,\mathbf{z}\right):\mathbb{R}\times\mathbb{R}^m\mapsto\mathbb{R}^m.`$

The method provides two approximations for the next step: the high-order solution $`\mathbf{z}_{n+1}`$ (order $p$) and a lower-order embedded solution $`\hat{\mathbf{z}}_{n+1}`$ (order $`\hat{p}`$, typically $`p-1`$):
```math
\begin{gather}
\mathbf{z}_{n+1} = \mathbf{z}_n+\tau_n\sum\limits_{i=1}^{s}b_i\mathbf{k}_{i}^{(n)},
\\
\hat{\mathbf{z}}_{n+1} = \mathbf{z}_n+\tau_n\sum\limits_{i=1}^{s}\hat{b}_i\mathbf{k}_{i}^{(n)}.
\end{gather}
```

The difference between them yields an estimate of the Local Truncation Error (LTE), which is used for adaptive step size control. The error is normalized using a standard mixed absolute-relative criterion:
```math
\begin{gather}
\textbf{LTE}_{n+1} \approx \hat{\mathbf{z}}_{n+1}- \mathbf{z}_{n+1}= \mathbf{z}_n+\tau_n\mathbf{K}^{(n)}\hat{\mathbf{b}}-\left(\mathbf{z}_n+\tau_n\mathbf{K}^{(n)}\mathbf{b}\right)= \tau_n\mathbf{K}^{(n)}\mathbf{d},
\\
\mathbf{w}_{n+1} = \mathrm{ATol} \cdot \mathbf{1} + \mathrm{RTol}\max\left(\left|\mathbf{z}_{n+1}\right|, \left|\mathbf{z}_{n}\right|\right),
\\
\mathrm{err}_{n+1} = \|\textbf{LTE}_{n+1}\oslash\mathbf{w}_{n+1}\|_\infty.
\end{gather}
```
Here, $`\oslash`$ denotes element-wise division. The step size $`\tau_{n+1}`$ is adapted based on the scalar error metric $`\mathrm{err}_{n+1}`$ using an Integral Controller:
```math
\tau_{n+1} = S\tau_n\mathrm{err}_{n+1}^{-\alpha}
```
where $S$ is a safety factor and $\alpha \approx 0.7/p$.

The approximation $\mathbf{z}_{n+1}$ is used to continue the integration.

**Butcher tableau** for the $s$-stage Embedded Runge—Kutta methods represented as follows:

```math
\begin{array}{r|c}
    \mathbf{c} & \mathbf{A} \\ \hline
    & \mathbf{b}^{\top} \\
    & \hat{\mathbf{b}}^{\top} \\ \hline
	& \mathbf{d}^{\top}
\end{array}
```
where $`\mathbf{d} = \hat{\mathbf{b}} - \mathbf{b}`$.

## Description of the implemented algorithm
The procedure for filling the matrix is identical as in [my algorithm for Explicit Runge—Kutta methods](https://github.com/whydenyscry/General-algorithm-of-the-explicit-Runge-Kutta-method).

```math
\begin{gather}
\mathbf{K}^{(n)}_{m\times s}=\left[\mathbf{k}_1^{(n)},\mathbf{k}_2^{(n)},\ldots,\mathbf{k}_s^{(n)}\right]=\mathbf{0}_{m\times s},
\\
\mathbf{A}_{s\times s} = 
	\begin{bmatrix}
		\mathbf{a}^{(1)\top}
		\\
		\mathbf{a}^{(2)\top}
		\\
		\vdots 
		\\
		\mathbf{a}^{(s)\top}
	\end{bmatrix}.
\end{gather}
```

Then the formulas for filling the matrix $\mathbf{K}^{(n)}$ can be represented as follows:

```math
\begin{gather}
\begin{cases}
		\mathbf{k}_{1}^{(n)} = \mathbf{f}\left(t_n,\mathbf{z}_n\right),\\
		\vdots\\
		\mathbf{k}_{i}^{(n)} = \mathbf{f}\left(t_n + c_i \tau, \mathbf{z}_n + \tau\mathbf{K}^{(n)}_{m\times i-1}\mathbf{a}_{i-1\times 1}^{(i)}\right),
	\end{cases}
\\
\mathbf{z}_{n+1} = \mathbf{z}_n+\tau_n\mathbf{K}^{(n)}\mathbf{b},
\\
\hat{\mathbf{z}}_{n+1} = \mathbf{z}_n+\tau_n\mathbf{K}^{(n)}\hat{\mathbf{b}}.
\end{gather}
```


## Example
The _ExampleOfUse.mlx_ file shows the obtaining of Arenstorf Orbit using Tsitouras (5)4 method.
```math
\begin{gather}
\begin{cases}
x'' = x + 2y' - (1 - \mu)\frac{x + \mu}{((x + \mu)^2 + y^2)^{3/2}} - \mu\frac{x - (1-\mu)}{((x - (1-\mu))^2 + y^2)^{3/2}},\\
y'' = y - 2x' - (1-\mu)\frac{y}{((x + \mu)^2 + y^2)^{3/2}} - \mu\frac{y}{((x - (1-\mu))^2 + y^2)^{3/2}}
\end{cases}
\\
x(0) = 0.994, \quad y(0) = 0, \quad x'(0) = 0, \quad y'(0) = -2.00158510637908252240537862224.
\end{gather}
```
The solution is periodic with period $T \approx 17.0652165601579625588917206249$.

<p align="center">
  <img src="images/ArenstorfOrbit.svg"/>
</p>

## Notes

### Syntax
`[t, zsol, dzdt_eval, stats] = odeEmbeddedSolvers(odefun, tspan, incond, options)`

### Input Arguments
- `odefun`: function handle defining the right-hand sides of the differential equations $`\dot{\mathbf{z}}\left(t\right)=\mathbf{f}\left(t,\mathbf{z}\right)`$. It must accept arguments (`t`, `z`) and return a column vector of derivatives;
- `tspan`: interval of integration, specified as a two-element vector;
- `incond`: vector of initial conditions.

#### Options (Name-Value Pairs)

You can customize the solver by passing `Name, Value` arguments after the required inputs.

| Option | Default | Description |
| :--- | :--- | :--- |
| **`Method`** | `"IIIXb+6(5)9*"` | The integration method to use. See **Supported Methods** below. |
| **`ATol`** | `1e-10` | Absolute tolerance for error control. |
| **`RTol`** | `1e-10` | Relative tolerance for error control. |
| **`SafetyFactor`** | `0.8` | Safety factor for step size prediction (prevents rejected steps). |
| **`MinStepSize`** | `1e-16` | Minimum allowed step size. |
| **`MaxStepSize`** | `1.0` | Maximum allowed step size. |
| **`MaxGrowth`** | `5.0` | Maximum factor by which the step size can increase in a single step. |
| **`MinGrowth`** | `0.2` | Minimum factor by which the step size can decrease. |


### Output Arguments
- `t`: vector of evaluation points used to perform the integration;
- `zsol`: solution matrix in which each row corresponds to a solution at the value returned in the corresponding row of `t`;
- `dzdt_eval`: matrix of derivatives $`\dot{\mathbf{z}}\left(t\right)`$ evaluated at the times in `t`; each row contains the derivative of the solution corresponding to the matching row of `t`;
- `stats`: structure containing solver performance statistics.

```matlab
stats = 
  struct with fields:
            n_feval: [integer]  % Total number of function evaluations
    n_success_steps: [integer]  % Number of accepted steps
     n_failed_steps: [integer]  % Number of rejected steps
        tau_history: [vector]   % History of step sizes used
      error_history: [vector]   % History of error ratios
```

## Planned Features
- Dense Output
- More Methods

## References
1. Butcher, J. (2016). Numerical methods for ordinary differential equations. https://doi.org/10.1002/9781119121534
2. Hairer, E., Nørsett, S. P., & Wanner, G. (1993). Solving Ordinary Differential Equations I: Nonstiff Problems (2nd ed.). Springer. https://doi.org/10.1007/978-3-540-78862-1
3. Dormand, J. R., & Prince, P. J. (1980). A family of embedded Runge-Kutta formulae. Journal of Computational and Applied Mathematics, 6(1), 19–26. https://doi.org/10.1016/0771-050x(80)90013-3
4. Tsitouras, C. (2011). Runge–Kutta pairs of order 5(4) satisfying only the first column simplifying assumption. Computers & Mathematics With Applications, 62(2), 770–775. https://doi.org/10.1016/j.camwa.2011.06.002
5. Verner, J. H. (1978). Explicit Runge–Kutta Methods with Estimates of the Local Truncation Error. SIAM Journal on Numerical Analysis, 15(4), 772–790. https://doi.org/10.1137/0715051
6. Verner, J. H. (1993). Differentiable interpolants for High-Order Runge–Kutta Methods. SIAM Journal on Numerical Analysis, 30(5), 1446–1466. https://doi.org/10.1137/0730075
7. Verner, J. H. (2009). Numerically optimal Runge–Kutta pairs with interpolants. Numerical Algorithms, 53(2–3), 383–396. https://doi.org/10.1007/s11075-009-9290-3
8. Verner, J. H. (n.d.). Jim Verner’s Refuge for Runge–Kutta Pairs. Simon Fraser University. https://www.sfu.ca/~jverner/

