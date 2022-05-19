function hyperbolic_cross_weights(y::AbstractArray{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2}) where {T <: AbstractFloat, S <: Integer}
  
  interpolation_matrix = zeros(size(nodes,1),size(nodes,1))
  
  n = 2*maximum(multi_index)+1
  
  unique_multi_index = sort(unique(multi_index))

  base_polynomials        = Array{Array{T,2},1}(undef,length(unique_multi_index))
  unique_base_polynomials = Array{Array{T,2},1}(undef,length(unique_multi_index))
 
  @inbounds for k in axes(nodes,1)
  
    # Construct the base polynomials
  
    for i in eachindex(unique_multi_index)
      base_polynomials[i] = chebyshev_polynomial(n-1,nodes[k,:])
    end
    
    # Compute the unique polynomial terms from the base polynomials
  
    for i in eachindex(unique_multi_index)
      if unique_multi_index[i] == 0
        unique_base_polynomials[i] = base_polynomials[i][:,1:1]
      else
        unique_base_polynomials[i] = base_polynomials[i][:,2*unique_multi_index[i]:2*unique_multi_index[i]+1]
      end
    end
  
    # Construct a row of the interplation matrix
  
    l = 1
    @inbounds @views for j in axes(multi_index,1)
      new_polynomials = unique_base_polynomials[multi_index[j,1]+1][1,:]
      for i = 2:size(multi_index,2)
        new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]+1][i,:])
      end
      m = length(new_polynomials)
      interpolation_matrix[k,l:l+m-1] = new_polynomials
      l += m
    end
  
  end
  
  weights = interpolation_matrix\y
  
  return weights
  
end
  
function hyperbolic_cross_weights(y::AbstractArray{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}
  
  # Normalize nodes to the [-1.0 1.0] interval
  
  d = size(multi_index,2)
  grid = similar(nodes)
  for i = 1:d
    grid[:,i] = normalize_node(nodes[:,i],domain[:,i])
  end
  
  weights = hyperbolic_cross_weights(y,grid,multi_index)
  
  return weights
  
end
  
function hyperbolic_cross_weights_threaded(y::AbstractArray{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}
  
  interpolation_matrix = zeros(size(nodes,1),size(nodes,1))
  
  n = 2*maximum(multi_index)+1
  
  unique_multi_index = sort(unique(multi_index))

  @inbounds @sync @qthreads for k in axes(nodes,1)
  
    base_polynomials        = Array{Array{T,2},1}(undef,length(unique_multi_index))
    unique_base_polynomials = Array{Array{T,2},1}(undef,length(unique_multi_index))
  
    # Construct the base polynomials
  
    for i in eachindex(unique_multi_index)
      base_polynomials[i] = chebyshev_polynomial(n-1,nodes[k,:])
    end

    # Compute the unique polynomial terms from the base polynomials
  
    for i in eachindex(unique_multi_index)
      if unique_multi_index[i] == 0
        unique_base_polynomials[i] = base_polynomials[i][:,1:1]
      else
        unique_base_polynomials[i] = base_polynomials[i][:,2*unique_multi_index[i]:2*unique_multi_index[i]+1]
      end
    end
  
    # Construct a row of the interplation matrix
  
    l = 1
    @inbounds @views for j in axes(multi_index,1)
      new_polynomials = unique_base_polynomials[multi_index[j,1]+1][1,:]
      for i = 2:size(multi_index,2)
        new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]+1][i,:])
      end
      m = length(new_polynomials)
      interpolation_matrix[k,l:l+m-1] = new_polynomials
      l += m
    end
  
  end

  weights = interpolation_matrix\y
  
  return weights

end
  
function hyperbolic_cross_weights_threaded(y::AbstractArray{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}
  
  # Normalize nodes to the [-1.0 1.0] interval
  
  d = size(multi_index,2)
  grid = similar(nodes)
  for i = 1:d
    grid[:,i] = normalize_node(nodes[:,i],domain[:,i])
  end
  
  weights = hyperbolic_cross_weights_threaded(y,grid,multi_index)
  
  return weights
  
end

function hyperbolic_cross_weights(y::AbstractArray{T,1},inverse_interpolation_matrix::Array{T,2}) where {T<:AbstractFloat}
  
  weights = inverse_interpolation_matrix*y
  
  return weights
  
end
  