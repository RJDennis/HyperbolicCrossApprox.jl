function generate_multi_index(d::S,k::S) where {S<:Integer} # Recursive function

  if d < 1
    error("d must be positive")
  end
  
  if k < 0
    error("k must be non-negative")
  end
  
  kprime = k+1
  
  if d == 1
    mi = zeros(S,kprime)
    for i = 1:k+1
      if i <= kprime
        mi[i] = i-1
      end
    end
    return mi
  else
    mi_base = generate_multi_index(d-1,k)
    N = size(mi_base,1)
    mi = zeros(S,N*kprime,d)
    pos = 0
    @inbounds @views for j = 1:N
      for i = 1:k+1
        if prod(mi_base[j,:].+1)*i <= kprime
          pos += 1
          mi[pos,1:d-1] .= mi_base[j,:]
          mi[pos,d] = i-1
        end
      end
    end
    return mi[1:pos,:]
  end
end
       
function generate_multi_index(d::S,k::S,n::S) where {S<:Integer} # Recursive function

  if d < 1
    error("d must be positive")
  end
      
  if k < 0
    error("k must be non-negative")
  end
  
  if n < 1
    error("n must be positive")
  end
  
  if !isodd(n)
    error("n must be odd")
  end
  
  if n >= 2*k+1
    return generate_multi_index(d,k)
  end

  kprime = k+1
  #ki = Int((n-1)/2) + 1
  ki = div(n-1,2) + 1
  
  if d == 1
    mi = zeros(S,ki)
    for i = 1:ki
      if i <= kprime
        mi[i] = i-1
      end
    end
    return mi
  else
    mi_base = generate_multi_index(d-1,k,n)
    N = size(mi_base,1)
    mi = zeros(S,N*ki,d)
    pos = 0
    @inbounds @views for j = 1:N
      for i = 1:ki
        if prod(mi_base[j,:].+1)*i <= kprime
          pos += 1
          mi[pos,1:d-1] .= mi_base[j,:]
          mi[pos,d] = i-1
        end
      end
    end
    return mi[1:pos,:]
  end
    
end
  
function generate_multi_index(d::S,k::S,n::Array{S,1}) where {S<:Integer} # Recursive function

  if d < 1
    error("d must be positive")
  end
  
  if k < 0
    error("k must be non-negative")
  end
    
  if length(n) != d
    error("node number for each dimension must be given")
  end
  
  for i = 1:d
    if n[i] < 1
      error("n must be positive in all dimensions")
    end
    if !isodd(n[i]) 
      error("each element in n must be odd")
    end
  end
  
  kprime = k + 1
  #ki = Int((n[end]-1)/2) + 1
  ki = div(n[end]-1,2) + 1
 
  if d == 1
    mi = zeros(S,ki)
    for i = 1:ki
      if i <= kprime
        mi[i] = i-1
      end
    end
    return mi
  else
    mi_base = generate_multi_index(d-1,k,n[1:d-1])
    nn = size(mi_base,1)
    mi = zeros(S,nn*ki,d)
    pos = 0
    @inbounds @views for j = 1:nn
      for i = 1:ki
        if prod(mi_base[j,:].+1)*i <= kprime
          pos += 1
          mi[pos,1:d-1] .= mi_base[j,:]
          mi[pos,d] = i-1
        end
      end
    end
    return mi[1:pos,:]
  end
end
  
function determine_grid_size(p::Array{S,2}) where {S<:Integer}
    
  q = similar(p)

  for i in eachindex(p)
    if p[i] > 0
      q[i] = 2
    else
      q[i] = 1
    end
  end

  s = 0
  for i in axes(q,1)
    t = 1
    for j in axes(q,2)
      t *= q[i,j]
    end
    s += t
  end
    
  return (s, size(q,2))
   
end

function hyperbolic_cross_grid(node_type::Function,d::S,ind::Array{S,2}) where {S<:Integer}
  
  multi_index = generate_multi_index(d,ind)
    
  # Create base nodes to be used in the sparse grid
    
  base_nodes = node_type(n)
  T = eltype(base_nodes)
    
  # Determine the unique nodes introduced at each higher level

  unique_base_nodes = Array{Array{T,1},2}(undef,size(multi_index))
  @inbounds for i in eachindex(multi_index)
    if multi_index[i] == 0
      unique_base_nodes[i] = [base_nodes[Int((n-1)/2)+1]]
    else
      unique_base_nodes[i] = [base_nodes[Int((n-1)/2)+1-multi_index[i]], base_nodes[Int((n-1)/2)+1+multi_index[i]]]
    end
  end

  # Construct the hyperbolic cross from the unique nodes
    
  nodes = Array{T,2}(undef,determine_grid_size(multi_index))
  l = 1
  @inbounds for j in axes(multi_index,1)
    new_nodes = unique_base_nodes[j,1]  # Here new_nodes is a 1d array
    for i = 2:d
      new_nodes = combine_nodes(new_nodes,unique_base_nodes[j,i])  # Here new_nodes becomes a 2d array
    end
    m = size(new_nodes,1)
    nodes[l:l+m-1,:] = new_nodes
    l += m
  end
    
  return nodes, multi_index
    
end

function hyperbolic_cross_grid(node_type::Function,d::S,k::S) where {S<:Integer}
  
  multi_index = generate_multi_index(d,k)
  
  n = 2*k+1
  
  # Create base nodes to be used in the sparse grid
    
  base_nodes = node_type(n)
  T = eltype(base_nodes)
    
  # Determine the unique nodes introduced at each higher level

  unique_base_nodes = Array{Array{T,1},2}(undef,size(multi_index))
  @inbounds for i in eachindex(multi_index)
    if multi_index[i] == 0
      unique_base_nodes[i] = [base_nodes[Int((n-1)/2)+1]]
    else
      unique_base_nodes[i] = [base_nodes[Int((n-1)/2)+1-multi_index[i]], base_nodes[Int((n-1)/2)+1+multi_index[i]]]
    end
  end

  # Construct the hyperbolic cross from the unique nodes
    
  nodes = Array{T,2}(undef,determine_grid_size(multi_index))
  l = 1
  @inbounds for j in axes(multi_index,1)
    new_nodes = unique_base_nodes[j,1]  # Here new_nodes is a 1d array
    for i = 2:d
      new_nodes = combine_nodes(new_nodes,unique_base_nodes[j,i])  # Here new_nodes becomes a 2d array
    end
    m = size(new_nodes,1)
    nodes[l:l+m-1,:] = new_nodes
    l += m
  end
    
  return nodes, multi_index
    
end

function hyperbolic_cross_grid(node_type::Function,d::S,k::S,n::S) where {S<:Integer}
  
  multi_index = generate_multi_index(d,k,n)
    
  # Create base nodes to be used in the sparse grid
    
  base_nodes = node_type(n)
  T = eltype(base_nodes)
    
  # Determine the unique nodes introduced at each higher level

  unique_base_nodes = Array{Array{T,1},2}(undef,size(multi_index))
  @inbounds for i in eachindex(multi_index)
    if multi_index[i] == 0
      unique_base_nodes[i] = [base_nodes[Int((n-1)/2)+1]]
    else
      unique_base_nodes[i] = [base_nodes[Int((n-1)/2)+1-multi_index[i]], base_nodes[Int((n-1)/2)+1+multi_index[i]]]
    end
  end

  # Construct the hyperbolic cross from the unique nodes
    
  nodes = Array{T,2}(undef,determine_grid_size(multi_index))
  l = 1
  @inbounds for j in axes(multi_index,1)
    new_nodes = unique_base_nodes[j,1]  # Here new_nodes is a 1d array
    for i = 2:d
      new_nodes = combine_nodes(new_nodes,unique_base_nodes[j,i])  # Here new_nodes becomes a 2d array
    end
    m = size(new_nodes,1)
    nodes[l:l+m-1,:] = new_nodes
    l += m
  end
    
  return nodes, multi_index
    
end

function hyperbolic_cross_grid(node_type::Function,d::S,k::S,n::Array{S,1}) where {S<:Integer}
  
  multi_index = generate_multi_index(d,k,n)
      
  # Create base nodes to be used in the sparse grid
      
  T = Float64
  base_nodes = Array{Array{T,1},1}(undef,d)

  for i = 1:d
    base_nodes[i] = node_type(n[i])
  end
      
  # Determine the unique nodes introduced at each higher level
  
  unique_base_nodes = Array{Array{T,1},2}(undef,size(multi_index))
  for i in axes(multi_index,1)
    for j = 1:d
      if multi_index[i,j] == 0
        unique_base_nodes[i,j] = [base_nodes[j][Int((n[j]-1)/2)+1]]
      else
        unique_base_nodes[i,j] = [base_nodes[j][Int((n[j]-1)/2)+1-multi_index[i,j]], base_nodes[j][Int((n[j]-1)/2)+1+multi_index[i,j]]]
      end
    end
  end
  
  # Construct the hyperbolic cross from the unique nodes
      
  nodes = Array{T,2}(undef,determine_grid_size(multi_index))
  l = 1
  @inbounds for j in axes(multi_index,1)
    new_nodes = unique_base_nodes[j,1]  # Here new_nodes is a 1d array
    for i = 2:d
      new_nodes = combine_nodes(new_nodes,unique_base_nodes[j,i])  # Here new_nodes becomes a 2d array
    end
    m = size(new_nodes,1)
    nodes[l:l+m-1,:] = new_nodes
    l += m
  end
      
  return nodes, multi_index
      
end

function hyperbolic_cross_grid(node_type::Function,d::S,ind::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {S<:Integer,T<:AbstractFloat}
  
  if size(domain,2) != d
    error("domain is inconsistent with the number of dimensions")
  end
    
  (nodes, multi_index) = hyperbolic_cross_grid(node_type,d,ind)
    
  # Now scale the nodes to the desired domain
    
  nodes .= scale_nodes(nodes,domain)
    
  return nodes, multi_index
    
end

function hyperbolic_cross_grid(node_type::Function,d::S,k::S,domain::Union{Array{T,1},Array{T,2}}) where {S<:Integer,T<:AbstractFloat}
  
  if size(domain,2) != d
    error("domain is inconsistent with the number of dimensions")
  end
    
  (nodes, multi_index) = hyperbolic_cross_grid(node_type,d,k)
    
  # Now scale the nodes to the desired domain
    
  nodes .= scale_nodes(nodes,domain)
    
  return nodes, multi_index
    
end

function hyperbolic_cross_grid(node_type::Function,d::S,k::S,n::Union{S,Array{S,1}},domain::Union{Array{T,1},Array{T,2}}) where {S<:Integer,T<:AbstractFloat}
  
  if size(domain,2) != d
    error("domain is inconsistent with the number of dimensions")
  end
    
  (nodes, multi_index) = hyperbolic_cross_grid(node_type,d,k,n)
    
  # Now scale the nodes to the desired domain
    
  nodes .= scale_nodes(nodes,domain)
    
  return nodes, multi_index
    
end

function hyperbolic_cross_weights(y::AbstractArray{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}
  
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
  
function hyperbolic_cross_inverse_interpolation_matrix(nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}
  
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
  
  inverse_interpolation_matrix = inv(interpolation_matrix)
  
  return inverse_interpolation_matrix
  
end
  
function hyperbolic_cross_inverse_interpolation_matrix(nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}
  
  # Normalize nodes to the [-1.0 1.0] interval
  
  d = size(multi_index,2)
  grid = similar(nodes)
  for i = 1:d
    grid[:,i] = normalize_node(nodes[:,i],domain[:,i])
  end
  
  inverse_interpolation_matrix = hyperbolic_cross_inverse_interpolation_matrix(grid,multi_index)
  
  return inverse_interpolation_matrix
  
end
  
function hyperbolic_cross_inverse_interpolation_matrix_threaded(nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}
  
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
  
  inverse_interpolation_matrix = inv(interpolation_matrix)
  
  return inverse_interpolation_matrix
  
end
  
function hyperbolic_cross_inverse_interpolation_matrix_threaded(nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}
  
  # Normalize nodes to the [-1.0 1.0] interval
  
  d = size(multi_index,2)
  grid = similar(nodes)
  for i = 1:d
    grid[:,i] = normalize_node(nodes[:,i],domain[:,i])
  end
  
  inverse_interpolation_matrix = hyperbolic_cross_inverse_interpolation_matrix_threaded(grid,multi_index)
  
  return inverse_interpolation_matrix
  
end
  
function hyperbolic_cross_weights(y::AbstractArray{T,1},inverse_interpolation_matrix::Array{T,2}) where {T<:AbstractFloat}
  
  weights = inverse_interpolation_matrix*y
  
  return weights
  
end

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
      
    @inbounds for i in axes(multi_index,1)
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
    @inbounds for j in axes(multi_index,1)
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
      
  for i in axes(multi_index,1)
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
  @inbounds for j in axes(multi_index,1)
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
  
  for i in eachindex(polynomials)
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