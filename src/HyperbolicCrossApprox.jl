module HyperbolicCrossApprox

import SmolyakApprox: scale_nodes,
                      combine_nodes

import ChebyshevApprox: normalize_node,
                        chebyshev_polynomial,
                        chebyshev_polynomial_deriv,
                        chebyshev_polynomial_sec_deriv,
                        chebyshev_nodes,
                        chebyshev_extrema,
                        chebyshev_extended,
                        vertesi_nodes,
                        legendre_nodes

using ThreadPools
using LinearAlgebra

include("hyperbolic_cross_functions.jl")

export HCApproxPlan

export chebyshev_nodes,
       chebyshev_extrema,
       chebyshev_extended,
       vertesi_nodes,
       legendre_nodes,
       hyperbolic_cross_grid,
       hyperbolic_cross_weights,
       hyperbolic_cross_weights_threaded,
       hyperbolic_cross_inverse_interpolation_matrix,
       hyperbolic_cross_inverse_interpolation_matrix_threaded,
       hyperbolic_cross_polynomial,
       hyperbolic_cross_evaluate,
       hyperbolic_cross_derivative,
       hyperbolic_cross_gradient,
       hyperbolic_cross_hessian,
       hyperbolic_cross_interp

end
