# This code presents an example to illustrate how HyperbolicCrossApprox can be used

using HyperbolicCrossApprox

function test_hyperbolic_cross_derivative()

  d  = 5  # Set the number of dimensions
  k = 5  # Set the level

  grid, multi_ind = hyperbolic_cross_grid(chebyshev_nodes,d,k)  # Construct the hyperbolic cross grid and the multi index

  # An arbitrary test function

  function test(grid)

    y_value = (grid[:,1].+1).^0.1.*exp.(grid[:,2]).*log.(grid[:,3].+2).^0.2.*(grid[:,4].+2).^0.8.*(grid[:,5].+7).^0.1

    return y_value

  end

  y = test(grid)  # Evaluate the test function on the Smolyak grid

  point = [0.75, 0.45, 0.82, -0.15, -0.95]  # Choose a point to evaluate the approximated function

  # One way of computing the weights and evaluating the approximated function

  weights = hyperbolic_cross_weights(y,grid,multi_ind)        # Compute the hyperbolic cross weights
  y_hat = hyperbolic_cross_evaluate(weights,point,multi_ind)  # Evaluate the approximated function

  #= A second way of computing the weights and evaluating the approximated function that
  computes the interpolation matrix just once. =#

  interp_mat = hyperbolic_cross_inverse_interpolation_matrix(grid,multi_ind)  # Compute the interpolation matrix
  w = hyperbolic_cross_weights(y,interp_mat)             # Compute the hyperbolic cross weights
  y_hatt = hyperbolic_cross_evaluate(w,point,multi_ind)  # Evaluate the approximated function

  # Evaluate the exact function at point

  y_actual = test(point')

  dom = [ones(1,5); -ones(1,5)]

  derivatives_2 = hyperbolic_cross_gradient(weights,point,multi_ind,dom)
  derivatives_3 = hyperbolic_cross_derivative(weights,point,multi_ind,dom,1)
  derivatives_4 = hyperbolic_cross_derivative(weights,point,multi_ind,dom,3)

  return derivatives_2, derivatives_3, derivatives_4

end

test_hyperbolic_cross_derivative()
