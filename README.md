# General algorithm of the Embedded Runge—Kutta method
Below I will describe an algorithm for solving IVP by any Embedded Runge—Kutta method of any order for any dimensionality of the system based on algorithm I have implemented in [General-algorithm-of-the-explicit-Runge-Kutta-method](https://github.com/whydenyscry/General-algorithm-of-the-explicit-Runge-Kutta-method).
I haven't implemented adaptive step size here, maybe I'll do it in a while. Here are two scripts — for a method with and without FSAL property.  The FSAL (First Same As Last) property means that if the same function evaluations (or calculations) are required at the start point of a step and at the end point of the next step, the function value calculated at the end point of the current step can be reused as the initial value for the next step without having to recalculate it. 
The output, compared to my algorithm for explicit methods, only adds an array with the estimates of local errors (ELE).

## Table of Contents

- [Embedded Runge—Kutta methods. Butcher tableau](#embedded-rungekutta-methods-butcher-tableau)
- [Description of the implemented algorithm](#description-of-the-implemented-algorithm)
- [Example](#example)
	- [Verner's method of order 6(5) (DVERK)](#verners-method-of-order-65-dverk)
	- [Dormand—Prince method of order 5(4) (RK5(4)7M)](#dormandprince-method-of-order-54-rk547m)
- [Notes](#notes)
  - [Input Arguments](#input-arguments)
  - [Output Arguments](#output-arguments)
  - [About Optimized Scripts](#about-optimized-scripts)
- [Planned Features](#planned-features)
- [References](#references)

## Embedded Runge—Kutta methods. Butcher tableau

Let an initial value problem (IVP) be specified as follows:

$$
\dot{\mathbf{x}}=\mathbf{f}\left(t,\mathbf{x}\right),\quad t \in \left[t_0,t_\text{end}\right],\quad \mathbf{x}\left(t_0\right) = \mathbf{x}_0 \in \mathbb{R}^m.
$$
 
where $\mathbf{x}=\left[x_1,\dots,x_m\right]^\mathbf{T},\quad
	\mathbf{f}\left(t,\mathbf{x}\right)=\left[f_1\left(t,x_1,\dots,x_n\right),\dots,f_m\left(t,x_1,\dots,x_n\right)\right]^\mathbf{T}.$
	
The idea behind the Embedded Runge—Kutta methods is to provide two approximations: $\mathbf{x}\_{n+1}$ and $\mathbf{\hat{x}}\_{n+1}$, such that 

$$
\mathbf{x}_{n+1} = \mathbf{x}_n+\tau_n\mathbf{K}^{(n)}\mathbf{b},
$$

is of order $p$, and

$$
\mathbf{\hat{x}}_{n+1} = \mathbf{x}_n+\tau_n\mathbf{K}^{(n)}\mathbf{\hat{b}},
$$

is of order $\hat{p}$. 

The difference between them gives an estimate of the local error for a less accurate result and can be used to control the step size:

$$
\text{ELE}\_{n+1} = \varepsilon\_{n+1} = \left\lVert\mathbf{\hat{x}}\_{n+1}- \mathbf{x}\_{n+1}\right\rVert\_\infty=\left\lVert \mathbf{x}\_n+\tau\_n\mathbf{K}^{(n)}\mathbf{\hat{b}}-\left(\mathbf{x}\_n+\tau\_n\mathbf{K}^{(n)}\mathbf{b}\right)\right\rVert\_\infty=\left\lVert \tau\_n\mathbf{K}^{(n)}\mathbf{d}\right\rVert\_\infty.
$$


The approximation $\mathbf{x}\_{n+1}$ is used to continue the integration.

**Butcher tableau** for the $s$-stage Embedded Runge—Kutta methods represented as follows:

$$
\begin{array}{r|c}
		\mathbf{c} & \mathbf{A} \\
		\hline
		& \mathbf{b}^{\mathbf{T}} \\
		& \mathbf{\hat{b}}^{\mathbf{T}} \\
		\hline
		& \mathbf{d}^{\mathbf{T}}
	\end{array} 	
	\quad \Rightarrow \quad 
	\begin{array}{r|ccccc}
		0     &         &         &         & \\
		c_2   & a\_{2,1}  &         &         & \\
		c_3   & a\_{3,1}  & a\_{3,2}  &         & \\
		\vdots& \vdots  & \vdots  & \ddots  & \\
		c_s   & a\_{s,1}  & a\_{s,2}  & \cdots  & a\_{s,s-1} \\
		\hline
		& b_1     & b_2     & \cdots  & b\_{s-1} & b_s \\
		& \hat{b}_1     & \hat{b}_2     & \cdots  & \hat{b}\_{s-1} & \hat{b}_s  \\
		\hline
		& d_1     & d_2     & \cdots  & d\_{s-1} & d_s
	\end{array},
$$

$$
\mathbf{d} = \mathbf{\hat{b}} - \mathbf{b},
$$

$$
\mathbf{c},\mathbf{b},\mathbf{\hat{b}}, \mathbf{d}\in \mathbb{R}^s,\quad \mathbf{A} \in \mathbb{R}^{s\times s}.
$$

## Description of the implemented algorithm
The procedure for filling the matrix is identical as in [my algorithm for Explicit Runge—Kutta methods](https://github.com/whydenyscry/General-algorithm-of-the-explicit-Runge-Kutta-method). At each iteration ([also take a look here](#AboutKmatrix)) we need to initialize the matrix $\mathbf{K}^{(n)}$ of the corresponding size as a zero matrix and this matrix is interpreted as follows:

$$
	\mathbf{K}^{(n)}\_{m\times s}=\left[\mathbf{k}_1^{(n)},\mathbf{k}_2^{(n)},\ldots,\mathbf{k}_s^{(n)}\right]=\mathbf{0}\_{m\times s},
$$

and the matrix $\mathbf{A}$ as

$$
\mathbf{A}_{s\times s} = 
	\begin{bmatrix}
		\mathbf{a}^{(1)\mathbf{T}}
		\\
		\mathbf{a}^{(2)\mathbf{T}}
		\\
		\vdots 
		\\
		\mathbf{a}^{(s)\mathbf{T}}
	\end{bmatrix}.
$$

Then the formulas for filling the matrix $\mathbf{K}^{(n)}$ can be represented as follows:

$$
\begin{cases}
		\mathbf{k}\_{1}^{(n)} = \mathbf{f}\left(t_n,\mathbf{x}_n\right),\\
		\vdots\\
		\mathbf{k}\_{i}^{(n)} = \mathbf{f}\left(t_n + c_i \tau, \mathbf{x}_n + \tau\mathbf{K}^{(n)}\_{m\times i-1}\mathbf{a}\_{i-1\times 1}^{(i)}\right),
	\end{cases}
$$

$$
	i=\overline{2,s}.
$$
## Example
The _ExampleOfUse.mlx_ file shows the obtaining of The Lotka—Volterra Attractor

$$ 
\begin{cases}
			\frac{\mathrm{d}x}{\mathrm{d}t}=x-xy+\varsigma x^2-\alpha z x^2, \\
			\frac{\mathrm{d}y}{\mathrm{d}t}=-y+xy, \\
			\frac{\mathrm{d}z}{\mathrm{d}t}=-\beta z +\alpha z x^2,
\end{cases}
$$
 
$$
\begin{bmatrix}
			\alpha\\
			\beta\\
			\varsigma
\end{bmatrix}=
		\begin{bmatrix}
			2.9851\\
			3\\
			2
\end{bmatrix},
$$

and The TSUCS2 Attractor

$$ 
\begin{cases}
			\frac{\mathrm{d}x}{\mathrm{d}t} = \alpha\left(y-x\right)+\delta xz, \\
			\frac{\mathrm{d}y}{\mathrm{d}t} = \varsigma x-xz+\xi y, \\
			\frac{\mathrm{d}z}{\mathrm{d}t} = \beta z+xy-\varepsilon x^2,
		\end{cases}
$$

$$ 
\begin{bmatrix}
			\alpha\\
			\beta\\
			\varsigma\\
			\delta\\
			\varepsilon\\
			\xi
		\end{bmatrix} = 
		\begin{bmatrix}
			40\\
			1.833\\
			55\\
			0.16\\
			20\\
			0.65
		\end{bmatrix}.
$$

### Verner's method of order 6(5) (DVERK)
The Lotka—Volterra Attractor has been obtained using DVERK method, which has no FSAL propety, so it's important to use the _odeEmbeddedGeneral_ (_odeEmbeddedGeneral\_optimized_) function.

![The Lotka—Volterra Attractor](https://github.com/whydenyscry/General-algorithm-of-the-embedded-Runge-Kutta-method/blob/main/images/The_Lotka_Volterra_Attractor.png)

### Dormand—Prince method of order 5(4) (RK5(4)7M)
The TSUCS2 Attractor has been obtained using RK5(4)7M method, which has FSAL propety, so it's recommended to use _odeFSALEmbeddedGeneral_ (_odeFSALEmbeddedGeneral\_optimized_), which takes this property into account, but it's also possible to use _odeEmbeddedGeneral_ (_odeEmbeddedGeneral\_optimized_) function, only in this case FSAL won't be taken into account.

![The TSUCS2 Attractor](https://github.com/whydenyscry/General-algorithm-of-the-embedded-Runge-Kutta-method/blob/main/images/The_TSUCS2_Attractor.png)

## Notes

### Input Arguments
- `c_vector`: vector of coefficients $\mathbf{c}$ of Butcher tableau for the selected method;
- `A_matrix`: matrix of coefficients $\mathbf{A}$ of Butcher tableau for the selected method;
- `b_vector`: vector of coefficients $\mathbf{b}$ of Butcher tableau for the selected method;
- `b_hat_vector`: vector of coefficients $\mathbf{\hat{b}}$ of Butcher tableau for the selected method;
- `odefun`: functions to solve, specified as a function handle that defines the functions to be integrated;
- `tspan`: interval of integration, specified as a two-element vector;
- `tau`: time discretization step;
- `incond`: vector of initial conditions.

### Output Arguments
- `t`: vector of evaluation points used to perform the integration;
- `xsol`: solution matrix in which each row corresponds to a solution at the value returned in the corresponding row of `t`.
- `ELE`: local errors vector in which each row corresponds to a estimated local error at the value returned in the corresponding row of `t`.

### About Optimized Scripts

The codes from the _odeEmbeddedGeneral.m_ & _odeFSALEmbeddedGeneral.m_ scripts shows a more illustrative integration procedure, for understanding from a theoretical point of view. The optimized versions of these scripts _odeEmbeddedGeneral_optimized.m_ & _odeFSALEmbeddedGeneral_optimized.m_ look as follows:
```MATLAB
function [t, xsol, ELE] = odeEmbeddedGeneral_optimized(c_vector, A_matrix, b_vector, b_hat_vector, odefun, tspan, tau, incond)

s_stages = length(c_vector);
m = length(incond);

c_vector = reshape(c_vector, [s_stages 1]);
b_vector = reshape(b_vector, [s_stages 1]);
b_hat_vector = reshape(b_hat_vector, [s_stages 1]);
d_vector = b_hat_vector - b_vector;
incond = reshape(incond, [m 1]);

t = (tspan(1) : tau : tspan(2))';
xsol = zeros(length(incond), length(t));
xsol(:, 1) = incond(:);
K_matrix = zeros(m, s_stages);
ELE = zeros(length(t), 1);

for n = 1:length(t)-1
    K_matrix(:, 1) = odefun(t(n), xsol(:, n));   
        for i = 2:s_stages
            K_matrix(:, i) = odefun(t(n) + tau * c_vector(i), xsol(:, n) + tau * K_matrix(:, 1:i-1) * A_matrix(i, 1:i-1)');
        end
    xsol(:, n+1) = xsol(:, n) + tau * K_matrix * b_vector;
    ELE(n+1) = norm(tau * K_matrix * d_vector, "inf");
end
xsol = xsol';
end
```

```MATLAB
function [t, xsol, ELE] = odeFSALEmbeddedGeneral_optimized(c_vector, A_matrix, b_vector, b_hat_vector, odefun, tspan, tau, incond)

s_stages = length(c_vector);
m = length(incond);

c_vector = reshape(c_vector, [s_stages 1]);
b_vector = reshape(b_vector, [s_stages 1]);
b_hat_vector = reshape(b_hat_vector, [s_stages 1]);
d_vector = b_hat_vector - b_vector;
incond = reshape(incond, [m 1]);

t = (tspan(1) : tau : tspan(2))';
xsol = zeros(length(incond), length(t));
xsol(:, 1) = incond(:);
K_matrix = zeros(m, s_stages);
ELE = zeros(length(t), 1);
K_matrix(:, s_stages) = odefun(t(1), xsol(:, 1));

for n = 1:length(t)-1
    K_matrix(:, 1) = K_matrix(:, s_stages);
        for i = 2:s_stages
            K_matrix(:, i) = odefun(t(n) + tau * c_vector(i), xsol(:, n) + tau * K_matrix(:, 1:i-1) * A_matrix(i, 1:i-1)');
        end
    xsol(:, n+1) = xsol(:, n) + tau * K_matrix * b_vector;
    ELE(n+1) = norm(tau * K_matrix * d_vector, "inf");
end
xsol = xsol';
end
```
With only 27 & 28 lines for such a powerful instrument, it looks awesome, doesn't it?

<span id="AboutKmatrix"></span> Here no unnecessary variables are created, and the `K_matrix` is initialized as zero matrix only once, because the algorithm allows not to fill it with zeros at each iteration, but just to overwrite the columns at this iteration without using the columns with coeficients from the previous one: 
```MATLAB
K_matrix(:, i) = odefun(t(n) + tau * c_vector(i), xsol(:, n) + tau * K_matrix(:, 1:i-1) * A_matrix(i, 1:i-1)')
```

## Planned Features
- based on this script, add specific integrators, as I did [here](https://github.com/whydenyscry/Dynamics-of-Nonlinear-Attractors/tree/main/scripts/odeExplicitSolvers).
- add adaptive step size.

## References
1. Butcher, J. (2016). Numerical methods for ordinary differential equations. https://doi.org/10.1002/9781119121534
2. Hairer, E., Nørsett, S. P., & Wanner, G. (1993). Solving Ordinary Differential Equations I: Nonstiff Problems (2nd ed.). Springer. https://doi.org/10.1007/978-3-540-78862-1
3. Dormand, J. R., & Prince, P. J. (1980). A family of embedded Runge-Kutta formulae. Journal of Computational and Applied Mathematics, 6(1), 19–26. https://doi.org/10.1016/0771-050x(80)90013-3
4. Samardzija, N., & Greller, L. D. (1988). Explosive route to chaos through a fractal torus in a generalized Lotka-Volterra model. Bulletin of Mathematical Biology, 50(5), 465–491. https://doi.org/10.1007/BF02458847
5. Li, D. (2008). A three-scroll chaotic attractor. Physics Letters A, 372(4), 387–393. https://doi.org/10.1016/j.physleta.2007.07.045
