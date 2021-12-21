function hyperbolic_cross_derivative(weights::Array{T,1},node::Array{R,1},multi_index::Array{S,2},pos::S) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  n = 2*maximum(multi_index,dims=1).+1
  d = length(n) # d is the number of dimensions
    
  base_polynomials                   = Array{Array{R,2},1}(undef,d)
  base_polynomial_derivatives        = Array{Array{R,2},1}(undef,d)
  unique_base_polynomials            = Array{Array{R,1},2}(undef,size(multi_index))
  unique_base_polynomial_derivatives = Array{Array{R,1},2}(undef,size(multi_index))
     
  # Construct the base polynomials
      
  for i = 1:d
    base_polynomials[i] = chebyshev_polynomial(n[i]-1,node[i])
    base_polynomial_derivatives[i] = derivative_of_chebyshev_polynomial(n[i],node[i])
  end
      
  # Compute the unique polynomial terms from the base polynomials
      
  for i = 1:size(multi_index,1)
   for j = 1:d
      if multi_index[i,j] == 0
        unique_base_polynomials[i,j]            = [base_polynomials[j][1]]
        unique_base_polynomial_derivatives[i,j] = [base_polynomial_derivatives[j][1]]
      else
        unique_base_polynomials[i,j]            = [base_polynomials[j][2*multi_index[i,j]], base_polynomials[j][2*multi_index[i,j]+1]]
        unique_base_polynomial_derivatives[i,j] = [base_polynomial_derivatives[j][2*multi_index[i,j]], base_polynomial_derivatives[j][2*multi_index[i,j]+1]]
      end
    end
  end
        
  # Construct the first row of the interplation matrix
      
  polynomials = Array{R,1}(undef,length(weights))

  l = 1
  @inbounds for j = 1:size(multi_index,1)
    if pos == 1
      new_polynomials = unique_base_polynomial_derivatives[j,1]
    else
      new_polynomials = unique_base_polynomials[j,1]
    end
    for i = 2:size(multi_index,2)
      if pos == i
        new_polynomials = kron(new_polynomials,unique_base_polynomial_derivatives[j,i])
      else
        new_polynomials = kron(new_polynomials,unique_base_polynomials[j,i])
      end
    end
    m = length(new_polynomials)
    polynomials[l:l+m-1] = new_polynomials
    l += m
  end
  
  evaluated_derivative = zero(T)
  
  for i = 1:length(polynomials)
    evaluated_derivative += polynomials[i]*weights[i]
  end
  
  return evaluated_derivative
  
end
  
function hyperbolic_cross_derivative(weights::Array{T,1},node::Array{R,1},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}},pos::S) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  node = copy(node)

  if size(domain,2) != length(node)
    error("domain is inconsistent with the number of dimensions")
  end
  
  d = length(node)
  for i = 1:d
    node[i] = normalize_node(node[i],domain[:,i])
  end
  
  evaluated_derivative = hyperbolic_cross_derivative(weights,node,multi_index,pos)
  
  return evaluated_derivative
  
end
