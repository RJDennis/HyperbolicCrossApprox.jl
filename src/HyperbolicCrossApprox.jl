module HyperbolicCrossApprox

import SmolyakApprox: normalize_node, scale_nodes, combine_nodes
import ChebyshevApprox: chebyshev_polynomial, derivative_of_chebyshev_polynomial, chebyshev_nodes, chebyshev_extrema, chebyshev_extended, vertesi_nodes

using ThreadPools
using FLoops

include("hyperbolic_cross_weights.jl")
include("hyperbolic_cross_evaluate.jl")
include("hyperbolic_cross_grid.jl")
include("generate_multi_index.jl")
include("hyperbolic_cross_derivative.jl")
include("hyperbolic_cross_gradient.jl")
include("determine_grid_size.jl")
include("hyperbolic_cross_inverse_interpolation_matrix.jl")

export chebyshev_nodes,
       chebyshev_extrema,
       chebyshev_extended,
       vertesi_nodes,
       hyperbolic_cross_grid,
       hyperbolic_cross_weights,
       hyperbolic_cross_weights_threaded,
       hyperbolic_cross_inverse_interpolation_matrix,
       hyperbolic_cross_inverse_interpolation_matrix_threaded,
       hyperbolic_cross_evaluate,
       hyperbolic_cross_derivative,
       hyperbolic_cross_gradient

end
