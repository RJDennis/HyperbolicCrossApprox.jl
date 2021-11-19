# This code presents an example to illustrate how HyperbolicCrossApprox can be used

using HyperbolicCrossApprox

function test_hyperbolic_cross_approx()

  d  = 5  # Set the number of dimensions
  k = 7  # Set the level of approximation

  g1, m1 = hyperbolic_cross_grid(chebyshev_nodes,d,k)
  g2, m2 = hyperbolic_cross_grid(chebyshev_extrema,d,k)  # Construct the hyperbolic cross grid and the multi index

  # An arbitrary test function

  function test(grid)

    y_value = (grid[:,1].+1).^0.1.*exp.(grid[:,2]).*log.(grid[:,3].+2).^0.2.*(grid[:,4].+2).^0.8.*(grid[:,5].+7).^0.1

    return y_value

  end

  y1 = test(g1)  # Evaluate the test function on the hyperbolic cross grid
  y2 = test(g2)

  point = [0.75, 0.45, 0.82, -0.15, -0.95]  # Choose a point to evaluate the approximated function

  # Chebyshev nodes

  # One way of computing the weights and evaluating the approximated function

  w1       = hyperbolic_cross_weights(y1,g1,m1)    # Compute the hyperbolic cross weights
  y1_hat   = hyperbolic_cross_evaluate(w1,point,m1)  # Evaluate the approximated function

  w2       = hyperbolic_cross_weights(y2,g2,m2)    # Compute the hyperbolic cross weights
  y2_hat   = hyperbolic_cross_evaluate(w2,point,m2)  # Evaluate the approximated function

  #= A second way of computing the weights and evaluating the approximated function that
  computes the interpolation matrix just once. =#

  interp_mat = hyperbolic_cross_inverse_interpolation_matrix(g1,m1)  # Compute the interpolation matrix
  w1b = hyperbolic_cross_weights(y1,interp_mat)             # Compute the hyperbolic cross weights
  y1b_hat = hyperbolic_cross_evaluate(w1b,point,m1)  # Evaluate the approximated function

  # Evaluate the exact function at point

  y_actual = test(point')

  # Now consider the ansiotropic case

  k = [7, 5, 5, 5, 7]
  g3, m3 = hyperbolic_cross_grid(chebyshev_nodes,d,k)  # Construct the hyperbolic cross grid and the multi index
  y3 = test(g3)
  w3 = hyperbolic_cross_weights(y3,g3,m3)                   # Compute the hyperbolic cross weights
  y3_hat_ansio = hyperbolic_cross_evaluate(w3,point,m3)       # Evaluate the approximated function
  w3_th = hyperbolic_cross_weights_threaded(y3,grid,m3)                   # Compute the hyperbolic cross weights

  g4, m4 =  hyperbolic_cross_grid(chebyshev_extrema,d,k)
  y4 = test(g4)
  w4 = hyperbolic_cross_weights(y4,g4,m4)
  w4_th = hyperbolic_cross_weights_threaded(y4,g4,m4)
  y4_hat_ansio = hyperbolic_cross_evaluate(w4,point,g4,m4)

  return y_actual, y1_hat, y2_hat, y1b_hat, y3_hat_ansio, y4_hat_ansio

end

test_hyperbolic_cross_approx()
