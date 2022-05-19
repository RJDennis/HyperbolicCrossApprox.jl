function hyperbolic_cross_gradient(weights::Array{T,1},node::Array{R,1},multi_index::Array{S,2}) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  d = length(node)
  gradient = Array{R,2}(undef,1,d)
  
  for i = 1:d
    gradient[i] = hyperbolic_cross_derivative(weights,node,multi_index,i)
  end
  
  return gradient
  
end
  
function hyperbolic_cross_gradient(weights::Array{T,1},node::Array{R,1},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  d = length(node)
  gradient = Array{R,2}(undef,1,d)
  
  for i = 1:d
    gradient[i] = hyperbolic_cross_derivative(weights,node,multi_index,domain,i)
  end
  
  return gradient

end