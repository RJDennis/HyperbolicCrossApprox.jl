module HyperbolicCrossApprox

import SmolyakApprox: scale_nodes!, # Not exposed by SmolyakApprox
                      combine_nodes # Not exposed by SmolyakApprox

import ChebyshevApprox: normalize_node, # Not exposed by ChebyshevApprox
                        chebyshev_polynomial,
                        chebyshev_polynomial_deriv,
                        chebyshev_polynomial_sec_deriv,
                        chebyshev_nodes,
                        chebyshev_extrema,
                        chebyshev_extended

using LinearAlgebra

include("hyperbolic_cross_functions.jl")

export HCApproxPlan

export hyperbolic_cross_plan,
       chebyshev_nodes,
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
       hyperbolic_cross_integrate,
       hyperbolic_cross_interp

end
