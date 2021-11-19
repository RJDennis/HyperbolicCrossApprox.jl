# HyperbolicCrossApprox.jl
 
This package implements a hyperbolicbolic cross sparse-grid method for approximating multivariate continuous functions with Chebyshev polynomials as basis functions.  Flexibility is given over the resulting grid through parameters that govern the maximum number of points along each spacial dimension and the number of layers along each spacial dimension.  The package also allows for isotropic and ansiotropic grids.  The associated weights for the hyperbolic grid are also computed allowing the package to be used to perform cubature.  Hyperbolic cross approximation is a sparse-grid alternative to Smolyak's method.  Where existing treatments of the hyperbolic cross method are based on non-equispaced fast Fourier methods, the treatment in this code is based on Chebyshev polynomials with the coefficients in the approximating function constructed using Lagrange interpolation.

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

Chebyshev polynomials
---------------------

The nodes used to form the approximating grid can be computed using either the Chebyshev nodes (points of the first kind) or Chebyshev extrema (points of the second kind), with the approximating grid and the multi-index computed by

```julia
grid, multi_ind = hyperbolic_cross_grid(chebyshev_nodes,d,k,domain)
```

where `d` is the dimension of the function, `k` is the number of layers along each dimension, and domain is a 2d-array (2xd) containing the upper and lower bound on each variable.  If domain is not provided, then it is assumed that the variables reside on the [-1,1] interval.  If `k` is an integer, then an isotropic grid is computed whereas if `k` is an array of integers with length `d`, then an ansiotropic grid is computed.  For the function above, the number of grid points along each spacial dimension is given by `n = 2k+1`.  Alternatively, the values for `k` and `n` can be provided separately (and hence not linked by `n=2k+1`) by

```julia
grid, multi_ind = hyperbolic_cross_grid(chebyshev_nodes,d,k,n,domain)
```

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

Lastly, we can evaluate the hyperbolic cross approximation of the function at any point in the domain by

```julia
y_hat = hyperbolic_cross_evaluate(weights,point,multi_ind,domain)
```

where `point` (a 1d-array) is the point in the domain where the approximation is to be evaluated.

The package can also be used to compute derivatives and gradients of the approximating function, as per

```Julia
deriv = hyperbolic_cross_derivative(weights,point,multi_ind,domain,pos)
```

where `pos` is an integer reflecting the position of the variable being differentiated with respect to, and 

```julia
gradient = hyperbolic_cross_gradient(weights,point,multi_ind,domain)
```

Some useful references are:

Dennis, R., (2021), "Using a Hyperbolic Cross to Solve Non-linear Macroeconomic Models," CAMA working paper number 93/2021, Australian National University.

Döhler, M., Kunis, S., and D. Potts, (2010), "Nonequispaced Hyperbolic Cross Fast Fourier Transform," SIAM Journal on Numerical Analysis ,47 (6), pp. 4415-4428.

Dũng, D., Temlyakov, V., and T. Ullrich, (2018), "Hyperbolic Cross Approximation," Advanced Courses in Mathematics, CRM Barcelona, Centre de Recerca Matemàtica, Birkhäuser.

Shen, J., L-L. Wang, (2010), "Sparse Spectral Approximations of High-Dimensional Problems Based on Hyperbolic Cross," SIAM Journal on Numerical Analysis, 48 (3), pp. 1087-1109.