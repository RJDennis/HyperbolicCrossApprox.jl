function hyperbolic_cross_polynomial(node::AbstractArray{R,1},multi_index::Array{S,2}) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  n = 2*maximum(multi_index,dims=1).+1
  d = length(n) # d is the number of dimensions
  
  base_polynomials        = Array{Array{R,2},1}(undef,d)
  unique_base_polynomials = Array{Array{R,1},2}(undef,size(multi_index))
   
  # Construct the base polynomials
    
  @inbounds for i = 1:d
    base_polynomials[i] = chebyshev_polynomial(n[i]-1,node[i])
  end
    
  # Compute the unique polynomial terms from the base polynomials
    
  @inbounds for i = 1:size(multi_index,1)
    for j = 1:d
      if multi_index[i,j] == 0
        unique_base_polynomials[i,j] = [base_polynomials[j][1]]
      else
        unique_base_polynomials[i,j] = [base_polynomials[j][2*multi_index[i,j]], base_polynomials[j][2*multi_index[i,j]+1]]
      end
    end
  end
      
  # Construct the first row of the interplation matrix

  n = determine_grid_size(multi_index)
  polynomial = Array{R,1}(undef,n[1])
  l = 1
  @inbounds for j = 1:size(multi_index,1)
    new_polynomials = unique_base_polynomials[j,1]
    for i = 2:d
      new_polynomials = kron(new_polynomials,unique_base_polynomials[j,i])
    end
    m = length(new_polynomials)
    polynomial[l:l+m-1] = new_polynomials
    l += m
  end
    
  return polynomial
  
end
  
function hyperbolic_cross_polynomial(node::AbstractArray{R,1},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  if size(domain,2) != length(node)
    error("domain is inconsistent with the number of dimensions")
  end
  
  node = copy(node)

  d = length(node)
  for i = 1:d
    node[i] = normalize_node(node[i],domain[:,i])
  end
  
  polynomial = hyperbolic_cross_polynomial(node,multi_index)

  return polynomial
  
end

function hyperbolic_cross_evaluate(weights::Array{T,1},node::AbstractArray{R,1},multi_index::Array{S,2}) where {T<:AbstractFloat,R<:Number,S<:Integer}

  poly = hyperbolic_cross_polynomial(node,multi_index)

  estimate = weights'poly

  return estimate

end

function hyperbolic_cross_evaluate(weights::Array{T,1},node::AbstractArray{R,1},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,R<:Number,S<:Integer}

if size(domain,2) != length(node)
  error("domain is inconsistent with the number of dimensions")
end

node = copy(node)

d = length(node)
for i = 1:d
  node[i] = normalize_node(node[i],domain[:,i])
end

estimate = hyperbolic_cross_evaluate(weights,node,multi_index)

return estimate

end

function hyperbolic_cross_evaluate(weights::Array{T,1},polynomial::Array{R,1}) where {T<:AbstractFloat,R<:Number}

  estimate = weights'polynomial

  return estimate

end

function hyperbolic_cross_evaluate(weights::Array{T,1},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

function goo(x::Array{R,1}) where {R<:Number}

  return hyperbolic_cross_evaluate(weights,x,multi_index)

end

return goo

end

function hyperbolic_cross_evaluate(weights::Array{T,1},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

function goo(x::Array{R,1}) where {R<:Number}

  return hyperbolic_cross_evaluate(weights,x,multi_index,domain)

end

return goo

end