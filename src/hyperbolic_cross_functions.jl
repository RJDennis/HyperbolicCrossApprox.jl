abstract type HCApproximationPlan end

struct HCApproxPlan{S<:Integer,T<:AbstractFloat} <: HCApproximationPlan

  grid::Union{Array{T,1},Array{T,2}}
  multi_index::Union{Array{S,1},Array{S,2}}
  domain::Union{Array{T,1},Array{T,2}}

end

function generate_multi_index(d::S,k::S) where {S<:Integer} # Recursive function

  if d < 1
    error("d must be positive")
  end
  
  if k < 0
    error("k must be non-negative")
  end
  
  kprime = k+1
  
  if d == 1
    mi = [0:1:k;]
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
  ki = div(n-1,2) + 1
  
  if d == 1
    mi = [0:1:ki-1;]
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
  
function generate_multi_index(d::S,k::S,n::Union{NTuple{N,S},Array{S,1}}) where {S<:Integer,N} # Recursive function

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
  ki = div(n[end]-1,2) + 1
 
  if d == 1
    mi = [0:1:ki-1;]
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
  
function determine_grid_size(p::Union{Array{S,1},Array{S,2}}) where {S<:Integer}
    
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
    
  return (s,size(q,2))
   
end

function hyperbolic_cross_grid(node_type::Function,d::S,k::S,domain=[ones(1,d);-ones(1,d)]) where {S<:Integer}
  
  multi_index = generate_multi_index(d,k)
  
  n = 2*k+1
  
  # Create base nodes to be used in the sparse grid
    
  base_nodes = node_type(n)
  T = eltype(base_nodes)
    
  # Determine the unique nodes introduced at each higher level

  if d == 1
    unique_base_nodes = Array{Array{T,1},2}(undef,(length(multi_index),1))
  else
    unique_base_nodes = Array{Array{T,1},2}(undef,size(multi_index))
  end
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
  
  # Now scale the nodes to the desired domain
    
  nodes .= scale_nodes(nodes,domain)

  return nodes, multi_index
    
end

function hyperbolic_cross_grid(node_type::Function,d::S,k::S,n::S,domain=[ones(1,d);-ones(1,d)]) where {S<:Integer}
  
  multi_index = generate_multi_index(d,k,n)
    
  # Create base nodes to be used in the sparse grid
    
  base_nodes = node_type(n)
  T = eltype(base_nodes)
    
  # Determine the unique nodes introduced at each higher level

  if d == 1
    unique_base_nodes = Array{Array{T,1},2}(undef,(length(multi_index),1))
  else
    unique_base_nodes = Array{Array{T,1},2}(undef,size(multi_index))
  end
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
    
  # Now scale the nodes to the desired domain

  nodes .= scale_nodes(nodes,domain)

  return nodes, multi_index
    
end

function hyperbolic_cross_grid(node_type::Function,d::S,k::S,n::Union{NTuple{N,S},Array{S,1}},domain=[ones(1,d);-ones(1,d)]) where {S<:Integer,N}
  
  multi_index = generate_multi_index(d,k,n)

  # Create base nodes to be used in the sparse grid
      
  T = Float64
  base_nodes = Array{Array{T,1},1}(undef,d)

  for i = 1:d
    base_nodes[i] = node_type(n[i])
  end
      
  # Determine the unique nodes introduced at each higher level
  
  if d == 1
    unique_base_nodes = Array{Array{T,1},2}(undef,(length(multi_index),1))
  else
    unique_base_nodes = Array{Array{T,1},2}(undef,size(multi_index))
  end
  @inbounds for i in axes(multi_index,1)
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
  
  # Now scale the nodes to the desired domain

  nodes .= scale_nodes(nodes,domain)

  return nodes, multi_index

end

function hyperbolic_cross_plan(node_type::Function,d::S,k::S,n::Union{NTuple{N,S},Array{S,1}},domain=[ones(1,d);-ones(1,d)]) where {S<:Integer,N}

    g, mi = hyperbolic_cross_grid(node_type,d,k,n,domain)
    plan = HCApproxPlan(g,mi,domain)

    return plan

end

function hyperbolic_cross_plan(node_type::Function,d::S,k::S,domain=[ones(1,d);-ones(1,d)]) where {S<:Integer}

  n = [2k+1 for _ in 1:d]

  g, mi = hyperbolic_cross_grid(node_type,d,k,n,domain)
  plan = HCApproxPlan(g,mi,domain)

  return plan

end


function hyperbolic_cross_weights(y::AbstractArray{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}}) where {T<:AbstractFloat,S<:Integer}
  
  interpolation_matrix = zeros(size(nodes,1),size(nodes,1))
  
  n = 2*maximum(multi_index)+1
  
  unique_multi_index = sort(unique(multi_index))

  base_polynomials        = Array{Array{T,2},1}(undef,length(unique_multi_index))
  unique_base_polynomials = Array{Array{T,2},1}(undef,length(unique_multi_index))
 
  @inbounds for k in axes(nodes,1)
  
    # Construct the base polynomials
  
    @inbounds for i in eachindex(unique_multi_index)
      base_polynomials[i] = chebyshev_polynomial(n-1,nodes[k,:])
    end
    
    # Compute the unique polynomial terms from the base polynomials
  
    @inbounds for i in eachindex(unique_multi_index)
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
  
function hyperbolic_cross_weights(y::AbstractArray{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}
  
  # Normalize nodes to the [-1.0 1.0] interval
  
  d = size(multi_index,2)
  grid = similar(nodes)
  for i = 1:d
    grid[:,i] = normalize_node(nodes[:,i],domain[:,i])
  end
  
  weights = hyperbolic_cross_weights(y,grid,multi_index)
  
  return weights
  
end
  
function hyperbolic_cross_weights_threaded(y::AbstractArray{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}}) where {T<:AbstractFloat,S<:Integer}
  
  interpolation_matrix = zeros(size(nodes,1),size(nodes,1))
  
  n = 2*maximum(multi_index)+1
  
  unique_multi_index = sort(unique(multi_index))

  @inbounds @sync @qthreads for k in axes(nodes,1)
  
    base_polynomials        = Array{Array{T,2},1}(undef,length(unique_multi_index))
    unique_base_polynomials = Array{Array{T,2},1}(undef,length(unique_multi_index))
  
    # Construct the base polynomials
  
    @inbounds for i in eachindex(unique_multi_index)
      base_polynomials[i] = chebyshev_polynomial(n-1,nodes[k,:])
    end

    # Compute the unique polynomial terms from the base polynomials
  
    @inbounds for i in eachindex(unique_multi_index)
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
  
function hyperbolic_cross_weights_threaded(y::AbstractArray{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}
  
  # Normalize nodes to the [-1.0 1.0] interval
  
  d = size(multi_index,2)
  grid = similar(nodes)
  for i = 1:d
    grid[:,i] = normalize_node(nodes[:,i],domain[:,i])
  end
  
  weights = hyperbolic_cross_weights_threaded(y,grid,multi_index)
  
  return weights
  
end
  
function hyperbolic_cross_inverse_interpolation_matrix(nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}}) where {T<:AbstractFloat,S<:Integer}
  
  interpolation_matrix = zeros(size(nodes,1),size(nodes,1))
  
  n = 2*maximum(multi_index)+1
  
  unique_multi_index = sort(unique(multi_index))

  base_polynomials        = Array{Array{T,2},1}(undef,length(unique_multi_index))
  unique_base_polynomials = Array{Array{T,2},1}(undef,length(unique_multi_index))
 
  @inbounds for k in axes(nodes,1)
  
    # Construct the base polynomials
  
    @inbounds for i in eachindex(unique_multi_index)
      base_polynomials[i] = chebyshev_polynomial(n-1,nodes[k,:])
    end
    
    # Compute the unique polynomial terms from the base polynomials
  
    @inbounds for i in eachindex(unique_multi_index)
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
  
function hyperbolic_cross_inverse_interpolation_matrix(nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}
  
  # Normalize nodes to the [-1.0 1.0] interval
  
  d = size(multi_index,2)
  grid = similar(nodes)
  for i = 1:d
    grid[:,i] = normalize_node(nodes[:,i],domain[:,i])
  end
  
  inverse_interpolation_matrix = hyperbolic_cross_inverse_interpolation_matrix(grid,multi_index)
  
  return inverse_interpolation_matrix
  
end
  
function hyperbolic_cross_inverse_interpolation_matrix_threaded(nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}}) where {T<:AbstractFloat,S<:Integer}
  
  interpolation_matrix = zeros(size(nodes,1),size(nodes,1))
  
  n = 2*maximum(multi_index)+1
  
  unique_multi_index = sort(unique(multi_index))
 
  @inbounds @sync @qthreads for k in axes(nodes,1)

    base_polynomials        = Array{Array{T,2},1}(undef,length(unique_multi_index))
    unique_base_polynomials = Array{Array{T,2},1}(undef,length(unique_multi_index))
  
    # Construct the base polynomials
  
    @inbounds for i in eachindex(unique_multi_index)
      base_polynomials[i] = chebyshev_polynomial(n-1,nodes[k,:])
    end
    
    # Compute the unique polynomial terms from the base polynomials
  
    @inbounds for i in eachindex(unique_multi_index)
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
  
function hyperbolic_cross_inverse_interpolation_matrix_threaded(nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}
  
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

function hyperbolic_cross_polynomial(node::AbstractArray{R,1},multi_index::Union{Array{S,1},Array{S,2}}) where {R<:Number,S<:Integer}
  
    n = 2*maximum(multi_index,dims=1).+1
    d = length(n) # d is the number of dimensions
    
    base_polynomials        = Array{Array{R,2},1}(undef,d)
    if d == 1
      unique_base_polynomials = Array{Array{R,1},1}(undef,length(multi_index))
    else
      unique_base_polynomials = Array{Array{R,1},2}(undef,size(multi_index))
    end
     
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
    
function hyperbolic_cross_polynomial(node::AbstractArray{R,1},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,R<:Number,S<:Integer}
    
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

function hyperbolic_cross_evaluate(weights::Array{T,1},node::AbstractArray{R,1},multi_index::Union{Array{S,1},Array{S,2}}) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
    poly = hyperbolic_cross_polynomial(node,multi_index)

    estimate = weights'poly

    return estimate

end
  
function hyperbolic_cross_evaluate(weights::Array{T,1},node::AbstractArray{R,1},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
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

function hyperbolic_cross_interp(y::AbstractArray{T,1},plan::P) where {T<:AbstractFloat,P<:HCApproximationPlan}

  weights = hyperbolic_cross_weights(y,plan.grid,plan.multi_index,plan.domain)
  
  function hcross_interp(x::Array{R,1}) where {R<:Number}
  
    return hyperbolic_cross_evaluate(weights,x,plan.multi_index,plan.domain)
  
  end
  
  return hcross_interp
  
end
  
function hyperbolic_cross_interp_threaded(y::AbstractArray{T,1},plan::P) where {T<:AbstractFloat,P<:HCApproximationPlan}
  
  weights = hyperbolic_cross_weights_threaded(y,plan.grid,plan.multi_index,plan.domain)
  
  function hcross_interp(x::Array{R,1}) where {R<:Number}
  
    return hyperbolic_cross_evaluate(weights,x,plan.multi_index,plan.domain)
  
  end
  
  return hcross_interp
  
end

function hyperbolic_cross_derivative(weights::Array{T,1},point::Array{R,1},multi_index::Union{Array{S,1},Array{S,2}},pos::S) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  n = 2*maximum(multi_index,dims=1).+1
  d = length(n) # d is the number of dimensions
    
  base_polynomials                   = Array{Array{R,2},1}(undef,d)
  base_polynomial_derivatives        = Array{Array{R,2},1}(undef,d)
  if d == 1
    unique_base_polynomials            = Array{Array{R,1},1}(undef,length(multi_index))
    unique_base_polynomial_derivatives = Array{Array{R,1},1}(undef,length(multi_index))
  else
    unique_base_polynomials            = Array{Array{R,1},2}(undef,size(multi_index))
    unique_base_polynomial_derivatives = Array{Array{R,1},2}(undef,size(multi_index))
  end

  # Construct the base polynomials
      
  for i = 1:d
    base_polynomials[i] = chebyshev_polynomial(n[i]-1,point[i])
    base_polynomial_derivatives[i] = chebyshev_polynomial_deriv(n[i]-1,point[i])
  end
      
  # Compute the unique polynomial terms from the base polynomials
      
  @inbounds for i in axes(multi_index,1)
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
  
function hyperbolic_cross_derivative(weights::Array{T,1},point::Array{R,1},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}},pos::S) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  point = copy(point)

  if size(domain,2) != length(point)
    error("domain is inconsistent with the number of dimensions")
  end
  
  d = length(point)
  for i = 1:d
    point[i] = normalize_node(point[i],domain[:,i])
  end
  
  evaluated_derivative = hyperbolic_cross_derivative(weights,point,multi_index,pos)
  
  return evaluated_derivative*(2.0/(domain[1,pos]-domain[2,pos]))
  
end
  
function hyperbolic_cross_gradient(weights::Array{T,1},point::Array{R,1},multi_index::Union{Array{S,1},Array{S,2}}) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  d = length(point)
  gradient = Array{R,2}(undef,1,d)
  
  for i = 1:d
    gradient[i] = hyperbolic_cross_derivative(weights,point,multi_index,i)
  end
  
  return gradient
  
end

function hyperbolic_cross_gradient(weights::Array{T,1},point::Array{R,1},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  d = length(point)
  gradient = Array{R,2}(undef,1,d)
  
  for i = 1:d
    gradient[i] = hyperbolic_cross_derivative(weights,point,multi_index,domain,i)
  end
  
  return gradient

end

function hyperbolic_cross_gradient(y::AbstractArray{T,1},plan::P) where {T<:AbstractFloat,P<:HCApproximationPlan}

  weights = hyperbolic_cross_weights(y,plan.grid,plan.multi_index,plan.domain)
  
  function hcross_grad(x::Array{R,1}) where {R<:Number}
  
    return hyperbolic_cross_gradient(weights,x,plan.multi_index,plan.domain)
  
  end
  
  return hcross_grad
  
end
  
function hyperbolic_cross_gradient_threaded(y::AbstractArray{T,1},plan::P) where {T<:AbstractFloat,P<:HCApproximationPlan}
  
  weights = hyperbolic_cross_weights_threaded(y,plan.grid,plan.multi_index,plan.domain)
  
  function hcross_grad(x::Array{R,1}) where {R<:Number}
  
    return hyperbolic_cross_gradient(weights,x,plan.multi_index,plan.domain)
  
  end
  
  return hcross_grad
  
end

function hyperbolic_cross_hessian(weights::Array{T,1},point::Array{R,1},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  point = copy(point)

  if size(domain,2) != length(point)
    error("domain is inconsistent with the number of dimensions")
  end
  
  d = length(point)
  for i = 1:d
    point[i] = normalize_node(point[i],domain[:,i])
  end
  
  n = 2*maximum(multi_index,dims=1).+1
  d = length(n) # d is the number of dimensions
  
  hess = Array{T,2}(undef,d,d)

  base_polynomials                   = Array{Array{R,2},1}(undef,d)
  base_polynomial_derivatives        = Array{Array{R,2},1}(undef,d)
  base_polynomial_sec_derivatives    = Array{Array{R,2},1}(undef,d)
  if d == 1
    unique_base_polynomials                = Array{Array{R,1},1}(undef,length(multi_index))
    unique_base_polynomial_derivatives     = Array{Array{R,1},1}(undef,length(multi_index))
    unique_base_polynomial_sec_derivatives = Array{Array{R,1},1}(undef,length(multi_index))
  else
    unique_base_polynomials                = Array{Array{R,1},2}(undef,size(multi_index))
    unique_base_polynomial_derivatives     = Array{Array{R,1},2}(undef,size(multi_index))
    unique_base_polynomial_sec_derivatives = Array{Array{R,1},2}(undef,size(multi_index))
  end

  # Construct the base polynomials
      
  for i = 1:d
    base_polynomials[i] = chebyshev_polynomial(n[i]-1,point[i])
    base_polynomial_derivatives[i] = chebyshev_polynomial_deriv(n[i]-1,point[i])
    base_polynomial_sec_derivatives[i] = chebyshev_polynomial_sec_deriv(n[i]-1,point[i])
  end
      
  # Compute the unique polynomial terms from the base polynomials
      
  @inbounds for i in axes(multi_index,1)
   for j = 1:d
      if multi_index[i,j] == 0
        unique_base_polynomials[i,j]                = [base_polynomials[j][1]]
        unique_base_polynomial_derivatives[i,j]     = [base_polynomial_derivatives[j][1]]
        unique_base_polynomial_sec_derivatives[i,j] = [base_polynomial_sec_derivatives[j][1]]
      else
        unique_base_polynomials[i,j]                = [base_polynomials[j][2*multi_index[i,j]], base_polynomials[j][2*multi_index[i,j]+1]]
        unique_base_polynomial_derivatives[i,j]     = [base_polynomial_derivatives[j][2*multi_index[i,j]], base_polynomial_derivatives[j][2*multi_index[i,j]+1]]
        unique_base_polynomial_sec_derivatives[i,j] = [base_polynomial_sec_derivatives[j][2*multi_index[i,j]], base_polynomial_sec_derivatives[j][2*multi_index[i,j]+1]]
      end
    end
  end
        
  # Construct the first row of the interplation matrix
      
  polynomials = Array{R,1}(undef,length(weights))

  @inbounds for c in CartesianIndices(hess)
    l = 1
    for j in axes(multi_index,1)
      if 1 == c[1] == c[2]
        new_polynomials = unique_base_polynomial_sec_derivatives[j,1]
      elseif 1 == c[1] || 1 == c[2]
        new_polynomials = unique_base_polynomial_derivatives[j,1]
      else
        new_polynomials = unique_base_polynomials[j,1]
      end
      for i = 2:size(multi_index,2)
        if i == c[1] == c[2]
          new_polynomials = kron(new_polynomials,unique_base_polynomial_sec_derivatives[j,i])
        elseif i == c[1] || i == c[2]
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
  
    hess[c] = evaluated_derivative*(2.0/(domain[1,c[1]]-domain[2,c[1]]))*(2.0/(domain[1,c[2]]-domain[2,c[2]]))

  end
  
  return hess
  
end

function hyperbolic_cross_hessian(y::AbstractArray{T,1},plan::P) where {T<:AbstractFloat,P<:HCApproximationPlan}

  weights = hyperbolic_cross_weights(y,plan.grid,plan.multi_index,plan.domain)
  
  function hcross_hess(x::Array{R,1}) where {R<:Number}
  
    return hyperbolic_cross_hessian(weights,x,plan.multi_index,plan.domain)
  
  end
  
  return hcross_hess
  
end
  
function hyperbolic_cross_hessian_threaded(y::AbstractArray{T,1},plan::P) where {T<:AbstractFloat,P<:HCApproximationPlan}
  
  weights = hyperbolic_cross_weights_threaded(y,plan.grid,plan.multi_index,plan.domain)
  
  function hcross_hess(x::Array{R,1}) where {R<:Number}
  
    return hyperbolic_cross_hessian(weights,x,plan.multi_index,plan.domain)
  
  end
  
  return hcross_hess
  
end

function integrate_cheb_polys(order::S) where {S <: Integer}

  # Integrates Chebyshev polynomials over the domain [-1,1]

  p = zeros(order+1)
    
  for i in 1:order+1
    if i == 2
      p[i] = 0.0
    else
      p[i] = ((-1)^(i-1)+1)/(1-(i-1)^2)
    end
  end

  return p

end

function hyperbolic_cross_integrate(f::Function,plan::HCApproxPlan,method::Symbol)

    if method == :clenshaw_curtis
        integral = hyperbolic_cross_clenshaw_curtis(f,plan)
    elseif method == :gauss_chebyshev_quad
        integral = hyperbolic_cross_gauss_chebyshev_quad(f,plan)
    else
        error("Integration not implemented for that method")
    end

    return integral

end

#function hyperbolic_cross_integrate(f::Function,plan::HCApproxPlan,method::Symbol,pos::S) where {S <: Integer}

#    if method == :clenshaw_curtis
#        integral = hyperbolic_cross_clenshaw_curtis(f,plan,pos)
#    elseif method == :gauss_chebyshev_quad
#        integral = hyperbolic_cross_gauss_chebyshev_quad(f,plan,pos)
#    else
#        error("Integration not implemented for that method")
#   end

#    return integral

#end

function hyperbolic_cross_clenshaw_curtis(f::Function,plan::HCApproxPlan)
  
  # Uses Clenshaw-Curtis to integrate over all dimensions

  grid        = plan.grid
  multi_index = plan.multi_index
  domain      = plan.domain

  y = zeros(size(grid,1))
  for i in eachindex(y)
    y[i] = f(grid[i,:])
  end

  weights = hyperbolic_cross_weights(y,grid,multi_index,domain)

  n = 2*maximum(multi_index,dims=1).+1
  d = length(n) # d is the number of dimensions
    
  T = eltype(grid)

  base_polynomial_integrals        = Array{Array{T,1},1}(undef,d)
  if d == 1
    unique_base_polynomial_integrals = Array{Array{T,1},1}(undef,length(multi_index))
  else
    unique_base_polynomial_integrals = Array{Array{T,1},2}(undef,size(multi_index))
  end

  # Construct the base polynomials
      
  for i = 1:d
    base_polynomial_integrals[i] = integrate_cheb_polys(n[i]-1)
  end

  # Compute the unique polynomial terms from the base polynomials
      
  @inbounds for i in axes(multi_index,1)
   for j = 1:d
      if multi_index[i,j] == 0
        unique_base_polynomial_integrals[i,j] = [base_polynomial_integrals[j][1]]
      else
        unique_base_polynomial_integrals[i,j] = [base_polynomial_integrals[j][2*multi_index[i,j]], base_polynomial_integrals[j][2*multi_index[i,j]+1]]
      end
    end
  end
        
  # Construct the first row of the interplation matrix
      
  polynomials = Array{T,1}(undef,length(weights))

  l = 1
  @inbounds for j in axes(multi_index,1)
    new_polynomials = unique_base_polynomial_integrals[j,1]
    for i = 2:size(multi_index,2)
        new_polynomials = kron(new_polynomials,unique_base_polynomial_integrals[j,i])
    end
    m = length(new_polynomials)
    polynomials[l:l+m-1] = new_polynomials
    l += m
  end
  
  evaluated_integral = zero(T)
  
  for i in eachindex(polynomials)
    evaluated_integral += polynomials[i]*weights[i]
  end

  scale_factor = (domain[1,1]-domain[2,1])/2
  for i in 2:size(multi_index,2)
    scale_factor = scale_factor*(domain[1,i]-domain[2,i])/2
  end
  
  return evaluated_integral*scale_factor
  
end

function hyperbolic_cross_clenshaw_curtis(f::Function,plan::HCApproxPlan,pos::S) where {S <: Integer}
  
  # Uses Clenshaw-Curtis to integrate over all dimensions except for pos

  grid        = plan.grid
  multi_index = plan.multi_index
  domain      = plan.domain

  y = zeros(size(grid,1))
  for i in eachindex(y)
    y[i] = f(grid[i,:])
  end

  weights = hyperbolic_cross_weights(y,grid,multi_index,domain)

  function hyperbolic_cross_int(point::R) where {R <: Number}

    point = normalize_node(point,domain[:,pos])
  
    n = 2*maximum(multi_index,dims=1).+1
    d = length(n) # d is the number of dimensions
    
    T = eltype(grid)

    base_polynomials                 = Array{Array{T,1},1}(undef,d)
    base_polynomial_integrals        = Array{Array{T,1},1}(undef,d)
    if d == 1
      unique_base_polynomials          = Array{Array{T,1},1}(undef,length(multi_index))
      unique_base_polynomial_integrals = Array{Array{T,1},1}(undef,length(multi_index))    
    else
      unique_base_polynomials          = Array{Array{T,1},2}(undef,size(multi_index))
      unique_base_polynomial_integrals = Array{Array{T,1},2}(undef,size(multi_index))
    end

    # Construct the base polynomials
      
    for i = 1:d
      base_polynomials[i]          = chebyshev_polynomial(n[i]-1,point)[:]
      base_polynomial_integrals[i] = integrate_cheb_polys(n[i]-1)
    end
      
    # Compute the unique polynomial terms from the base polynomials
      
    @inbounds for i in axes(multi_index,1)
      for j = 1:d
        if multi_index[i,j] == 0
          unique_base_polynomials[i,j]          = [base_polynomials[j][1]]
          unique_base_polynomial_integrals[i,j] = [base_polynomial_integrals[j][1]]
        else
          unique_base_polynomials[i,j]          = [base_polynomials[j][2*multi_index[i,j]], base_polynomials[j][2*multi_index[i,j]+1]]
          unique_base_polynomial_integrals[i,j] = [base_polynomial_integrals[j][2*multi_index[i,j]], base_polynomial_integrals[j][2*multi_index[i,j]+1]]
        end
      end
    end
        
    # Construct the first row of the interplation matrix
      
    polynomials = Array{T,1}(undef,length(weights))

    l = 1
    @inbounds for j in axes(multi_index,1)
      if pos == 1
        new_polynomials = unique_base_polynomials[j,1]
      else
        new_polynomials = unique_base_polynomial_integrals[j,1]
      end
      for i = 2:size(multi_index,2)
        if pos == i
          new_polynomials = kron(new_polynomials,unique_base_polynomials[j,i])
        else
          new_polynomials = kron(new_polynomials,unique_base_polynomial_integrals[j,i])
        end
      end
      m = length(new_polynomials)
      polynomials[l:l+m-1] = new_polynomials
      l += m
    end
  
    evaluated_integral = zero(T)
  
    for i in eachindex(polynomials)
      evaluated_integral += polynomials[i]*weights[i]
    end
  
    scale_factor = 1.0
    for i in 1:size(multi_index,2)
      if pos != i
        scale_factor = scale_factor*(domain[1,i]-domain[2,i])/2
      end
    end
    
    return evaluated_integral*scale_factor
    
  end

  return hyperbolic_cross_int

end

function hyperbolic_cross_gauss_chebyshev_quad(f::Function,plan::HCApproxPlan)
  
  # Uses Gauss-Chebyshev quadrature to integrate over all dimensions
  
  grid        = plan.grid
  multi_index = plan.multi_index
  domain      = plan.domain

  iim = hyperbolic_cross_inverse_interpolation_matrix(grid,multi_index,domain) 

  d = size(grid,2)

  e = zeros(1,size(grid,1))
  e[1] = π^d
  w = e*iim

  y = zeros(size(grid,1))
  for i in eachindex(y)
    integrating_weights = sqrt(1.0-normalize_node(grid[i,1],domain[:,1])^2)
    for j = 2:d
      integrating_weights *= sqrt(1.0-normalize_node(grid[i,j],domain[:,j])^2)
    end
    y[i] = f(grid[i,:])*integrating_weights
  end

  scale_factor = (domain[1,1]-domain[2,1])/2
  for i in 2:d
    scale_factor = scale_factor*(domain[1,i]-domain[2,i])/2
  end

  return (w*y)[1]*scale_factor

end