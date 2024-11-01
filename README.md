# HyperbolicCrossApprox.jl

Introduction
============

This package implements a hyperbolic cross sparse-grid method for approximating multivariate continuous functions with Chebyshev polynomials as basis functions.  Flexibility is given over the resulting grid through parameters that govern the maximum number of points along each spatial dimension and the number of layers along each spatial dimension.  The package also allows for isotropic and anisotropic grids.

Hyperbolic cross approximation is a sparse-grid alternative to Smolyak's method.  Where existing treatments of the hyperbolic cross method are based on non-equispaced fast Fourier methods, the treatment in this code is based on Chebyshev polynomials with the coefficients in the approximating function constructed using Lagrange interpolation.

This package draws on and is related to others I have written, specifically: ChebyshevApprox.jl and SmolyakApprox.jl. 

To install this package you need to type in the REPL

```julia
using Pkg
Pkg.add("HyperbolicCrossApprox")
```

Then the package can be used by typing

```julia
using HyperbolicCrossApprox
```

Approximating grid
------------------

The nodes used to form the approximating grid can be computed using either the Chebyshev nodes (points of the first kind), Chebyshev extrema (points of the second kind), or the extended Chebyshev nodes, with the approximating grid and the multi-index computed by

```julia
grid, multi_ind = hyperbolic_cross_grid(chebyshev_nodes,d,k,domain)
```

where `d` is the dimension of the function, `k` is the number of layers along each dimension, and domain is a 2d-array (2xd) containing the upper and lower bound on each variable.  If domain is not provided, then it is assumed that the variables reside on the [-1,1] interval.  If `k` is an integer, then an isotropic grid is computed whereas if `k` is an array of integers with length `d`, then an anisotropic grid is computed.  For the function above, the number of grid points along each spacial dimension is given by `n = 2k+1`.  Alternatively, the values for `k` and `n` can be provided separately (and hence not linked by `n=2k+1`) by

```julia
grid, multi_ind = hyperbolic_cross_grid(chebyshev_nodes,d,k,n,domain)
```

where `n` is either an interger or a 1d array of integers.

In the functions above, `chebyshev_nodes` can be replaced with `chebyshev_extrema`, or `chebyshev_extended`.

Polynomial coefficients
-----------------------

With the grid and multi-index in hand, we can compute the weights, or coefficients in the approximation, according to

```julia
weights = hyperbolic_cross_weights(y,grid,multi_ind,domain)
```

where `y` is a 1d-array containing the evaluations at each grid point of the function being approximated.  The weights can be computed more efficiently by computing the inverse interpolation matrix (this generally needs to be done only once, outside any loops)

```julia
inv_interp_mat = hyperbolic_cross_inverse_interpolation_matrix(grid,multi_ind,domain)
```

with the weights then computed through

```julia
weights = hyperbolic_cross_weights(y,inv_interp_mat)
```

Interpolation
-------------

We can evaluate the hyperbolic cross approximation of the function at any point in the domain by

```julia
y_hat = hyperbolic_cross_evaluate(weights,point,multi_ind,domain)
```

where `point` (a 1d-array) is the point in the domain where the approximation is to be evaluated.

Derivatives, gradients, and hessians
------------------------------------

The package can also be used to compute derivatives of the approximating function, as per

```Julia
deriv = hyperbolic_cross_derivative(weights,point,multi_ind,domain,pos)
```

where `pos` is an integer reflecting the number of the variable being differentiated with respect to.  The package can compute gradients according to 

```julia
gradient = hyperbolic_cross_gradient(weights,point,multi_ind,domain)
```

and hessians according to

```julia
hess = hyperbolic_cross_hessian(weights,point,multi_ind,domain)
```

Integration
-----------

Lastly, the package implements a Clenshaw-Curtis integration and a Guass-Chebyshev integration scheme.  To integrate a function with repsect to all variables over the approximating domain use

```julia
integral = hyperbolic_cross_integrate(f,plan,:clenshaw_curtis)
integral = hyperbolic_cross_integrate(f,plan,:gauss_chebyshev_quad)
```

where `f` is the function to be integrated and `plan` is a `HCApproxPlan` structure, described below. The `:gauss_chebyshev_quad` approach constructs the quadrature weights using a variant on "designed quadrature",  (Keshavarzzadeh, et al, (2018, SIAM J. SCI. COMPUT)). The `:clenshaw_curtis` method approximates the function to be integrated with a hyperbolic cross and then numerically integrates the approximation.

Multi-threading
---------------

There are multi-threaded functions to compute the polynomial weights and the interpolation matrix.  These multi-threaded functions are accessed by adding `_threaded` to the end of the funtion, as per

```julia
weights = hyperbolic_cross_weights_threaded(y,inv_interp_mat)
```

Useful structures
-----------------

The key structure to be aware of is the HCApproxPlan, which contains the key information needed to approximate a function.

```julia
d = 3
k = 5
domain = [2.0 2.0 2.0; -2.0 -2.0 -2.0]
grid, mi = hyperbolic_cross_grid(chebyshev_nodes,d,k,domain)
plan = HCApproxPlan(grid,mi,domain)
```
or
```julia
plan = hyperbolic_cross_plan(chebyshev_nodes,d,k,n,domain)
plan = hyperbolic_cross_plan(chebyshev_nodes,d,k,domain)
```

Once the approximation plan has been constructed it can be used to create functions to interpolate and to compute gradients and hessians.

```julia
f = hyperbolic_cross_interp(y,plan)
g = hyperbolic_cross_gradient(y,plan)
h = hyperbolic_cross_hessian(y,plan)

point = [1.0, 1.0, 1.0]

f(point)
g(point)
h(point)
```

There are multi-threaded versions of `hyperbolic_cross_interp`, `hyperbolic_cross_gradient`, and `hyperbolic_cross_hessian`; just add `_threaded` to the end of the function name.

Related packages
----------------

- ChebyshevApprox.jl
- SmolyakApprox.jl
- PiecewiseLinearApprox.jl

References
----------

Some useful references are:

Dennis, R., (2024), "Using a Hyperbolic Cross to Solve Non-linear Macroeconomic Models", Journal of Economic Dynamics and Control, 163, 104860.

Döhler, M., Kunis, S., and D. Potts, (2010), "Nonequispaced Hyperbolic Cross Fast Fourier Transform", *SIAM Journal on Numerical Analysis*, 47, (6), pp. 4415-4428.

Dũng, D., Temlyakov, V., and T. Ullrich, (2018), "Hyperbolic Cross Approximation", *Advanced Courses in Mathematics*, CRM Barcelona, Centre de Recerca Matemàtica, Birkhäuser.

Shen, J., L-L. Wang, (2010), "Sparse Spectral Approximations of High-Dimensional Problems Based on Hyperbolic Cross" *SIAM Journal on Numerical Analysis*, 48, (3), pp. 1087-1109.
