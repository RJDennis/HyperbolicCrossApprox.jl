abstract type HCApproximationPlan end

"""
HApproxPlan is an immutable struct that contains the information used to approximate a multi-dimensional function.
HApproxPlan has three fields: the approximation grid, the multi-index associated with the polynomial, and the 
approximation domain.
""" 
struct HCApproxPlan{S<:Integer,T<:AbstractFloat} <: HCApproximationPlan

  grid::Array{T,2}
  multi_index::Array{S,2}
  domain::Array{T,2}

end

"""
Creates the multi-index for a ```d```-variable ansiotropic grid with order given by ```k```.  Returns
a matrix of integers.

Signature
=========

m_index = generate_multi_index(d,k)

Example
=======
```
julia> m_index = generate_multi_index(2,3)
[0  0
 0  1
 0  2
 0  3
 1  0
 1  1
 2  0
 3  0]
 ```
"""
function generate_multi_index(d::S,k::S) where {S<:Integer} # Recursive function

  if d < 1
    error("d must be positive")
  end
  
  if k < 0
    error("k must be non-negative")
  end
  
  kprime = k+1
  
  if d == 1
    mi = [0:1:k;;] # mi is always a matrix
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

"""
Creates the multi-index for a ```d```-variable ansiotropic grid with order given by ```k``` and the maximum number
nodes along each dimension given by ```n```.  Returns a matrix of integers.

Signature
=========

m_index = generate_multi_index(d,k,n)

Example
=======
```
julia> m_index = generate_multi_index(2,3,5)
[0  0
 0  1
 0  2
 1  0
 1  1
 2  0]
```
"""
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
    mi = [0:1:ki-1;;] # mi is always a matrix
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
  
"""
Creates the multi-index for a ```d```-variable ansiotropic grid with order given by ```k``` and the maximum number
nodes along each dimension given by the vector ```n```.  Returns a matrix of integers.

Signature
=========

m_index = generate_multi_index(d,k,n)

Example
=======
```
julia> m_index = generate_multi_index(2,3,[5,3])
[0  0
 0  1
 1  0
 1  1
 2  0]
```
"""
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
  
  for i in eachindex(n)
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
    mi = [0:1:ki-1;;] # mi is always a matrix
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
  
function determine_grid_size(p::Array{S,2}) where {S<:Integer} # p is a multi-index, internal function, not exported
    
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

"""
Ensures the approximation domain has the expected form.

Signatures
==========

domain = check_domain(dom)

Example
=======

```
julia> dom = [-1.0 2.0; 2.0 -1.0]
julia> domain = check_domain(dom)
```
"""
function check_domain(d::S,domain::Union{Array{T,1},Array{T,2}}) where {S<:Integer, T<:AbstractFloat}

  if ndims(domain) == 1 # domain is a vector, so convert to a matrix
    dom = reshape(domain,2,1)
  else
    dom = copy(domain)
  end

  n = size(dom)
  if d != n[2]
    error("The size of 'domain' is inconsistent with the number of variables entering the function.")
  end

  for i in 1:n[2]
    if dom[1,i] < dom[2,i]
      dom[1,i] = maximum(domain[:,i])
      dom[2,i] = minimum(domain[:,i])
    end
  end

  return dom

end

"""
Uses the ```node_type``` function to construct the ```d```-dimensional hyperbolic cross grid with approximation layer
```k``` and ```domain``` (defaults to [-1.0,1.0]^d).  Returns the approximating grid and the associated multi index.

Signatures
==========

grid, multi_index = hyperbolic_cross_grid(node_type,d,k)
grid, multi_index = hyperbolic_cross_grid(chebyshev_extrema,d,k,domain)

Examples
========
```
julia> grid, m_index = hyperbolic_cross_grid(chebyshev_extrema,2,2)
julia> grid, m_index = hyperbolic_cross_grid(chebyshev_extrema,2,2,[3.0 1.5; 2.0 0.5])
```
"""
function hyperbolic_cross_grid(node_type::Function,d::S,k::S,domain=[ones(1,d);-ones(1,d)]) where {S<:Integer}
  
  dom = check_domain(d,domain)

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
      unique_base_nodes[i] = [base_nodes[k+1]]
    else
      unique_base_nodes[i] = [base_nodes[k+1-multi_index[i]], base_nodes[k+1+multi_index[i]]]
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
  
  scale_nodes!(nodes,dom)

  return nodes, multi_index
    
end

"""
Uses the ```node_type``` function to construt the ```d```-dimensional hyperbolic cross grid with approximation layer
```k```, and spacial node-maximum, ```n```, and ```domain``` (defaults to [-1.0,1.0]^d).  Returns the approximating 
grid and the associated multi index.

Signatures
==========

grid, multi_index = hyperbolic_cross_grid(node_type,d,k,n)
grid, multi_index = hyperbolic_cross_grid(chebyshev_extrema,d,k,n,domain)

Examples
========
```
julia> grid, m_index = hyperbolic_cross_grid(chebyshev_extrema,2,3,5)
julia> grid, m_index = hyperbolic_cross_grid(chebyshev_extrema,2,3,5,[3.0 1.5; 2.0 0.5])
```
"""
function hyperbolic_cross_grid(node_type::Function,d::S,k::S,n::S,domain=[ones(1,d);-ones(1,d)]) where {S<:Integer}
  
  dom = check_domain(d,domain)

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
      unique_base_nodes[i] = [base_nodes[div(n-1,2)+1]]
    else
      unique_base_nodes[i] = [base_nodes[div(n-1,2)+1-multi_index[i]], base_nodes[div(n-1,2)+1+multi_index[i]]]
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

  scale_nodes!(nodes,dom)

  return nodes, multi_index
    
end

"""
Uses the ```node_type``` function to construt the ```d```-dimensional hyperbolic cross grid with approximation layer
```k```, and spacial node-maximum, ```n```, and ```domain``` (defaults to [-1.0,1.0]^d).  Returns the approximating 
grid and the associated multi index.

Signatures
==========

grid, multi_index = hyperbolic_cross_grid(node_type,d,k,n)
grid, multi_index = hyperbolic_cross_grid(chebyshev_extrema,d,k,n,domain)

Examples
========
```
julia> grid, m_index = hyperbolic_cross_grid(chebyshev_extrema,2,3,[5,3])
julia> grid, m_index = hyperbolic_cross_grid(chebyshev_extrema,2,3,[7,5],[3.0 1.5; 2.0 0.5])
```
"""
function hyperbolic_cross_grid(node_type::Function,d::S,k::S,n::Union{NTuple{N,S},Array{S,1}},domain=[ones(1,d);-ones(1,d)]) where {S<:Integer,N}
  
  dom = check_domain(d,domain)

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
        unique_base_nodes[i,j] = [base_nodes[j][div(n[j]-1,2)+1]]
      else
        unique_base_nodes[i,j] = [base_nodes[j][div(n[j]-1,2)+1-multi_index[i,j]], base_nodes[j][div(n[j]-1,2)+1+multi_index[i,j]]]
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

  scale_nodes!(nodes,dom)

  return nodes, multi_index

end

"""
Constructs a hyperbolic cross approximation plan, given the ```node_type``` function, the number of dimensions,
```d```, the approximation layers, ```k```, the maximum number of spacial nodes, ```n``` (defaults to 2k+1), and the approximation
```domain``` (defaults to [-1.0,1.0]^d).  Returns an HApproxPlan struct.

Signatures
==========

hplan = hyperbolic_cross_plan(node_type,d,k)
hplan = hyperbolic_cross_plan(node_type,d,k,n)
hplan = hyperbolic_cross_plan(node_type,d,k,domain)
hplan = hyperbolic_cross_plan(node_type,d,k,n,domain)

Examples
========
```
julia> hplan = hyperbolic_cross_plan(chebyshev_extrema,2,2)
julia> hplan = hyperbolic_cross_plan(chebyshev_nodes,2,3,[5,3])
julia> hplan = hyperbolic_cross_plan(chebyshev_extrema,2,3,[3.0 1.5; 2.0 0.5])
julia> hplan = hyperbolic_cross_plan(chebyshev_nodes,2,3,[7,5],[3.0 1.5; 2.0 0.5])
```
"""
function hyperbolic_cross_plan(node_type::Function,d::S,k::S,n::Union{NTuple{N,S},Array{S,1}},domain=[ones(1,d);-ones(1,d)]) where {S<:Integer,N}

  dom = check_domain(d,domain)

  g, mi = hyperbolic_cross_grid(node_type,d,k,n,dom)
  plan = HCApproxPlan(g,mi,dom)

  return plan

end

function hyperbolic_cross_plan(node_type::Function,d::S,k::S,domain=[ones(1,d);-ones(1,d)]) where {S<:Integer}

  dom = check_domain(d,domain)

  n = [2k+1 for _ in 1:d]

  g, mi = hyperbolic_cross_grid(node_type,d,k,n,dom)
  plan = HCApproxPlan(g,mi,dom)

  return plan

end

"""
Uses Chebyshev polynomials as basis functions to compute the weights in a hyperbolic cross polynomial approximation
given the approximation sample, ```y```, the approximation ```grid```, the ```multi_index```, and the approximation 
```domain``` (defaults to [-1.0,1.0]^d).  Returns a vector containing the weights in the hyperbolic cross polynomial.

Signatures
==========

w = hyperbolic_cross_weights(y,grid,multi_index)
w = hyperbolic_cross_weights(y,grid,multi_index,domain)

Example
=======
```
julia> g,mi = hyperbolic_cross_grid(chebyshev_extrema,2,3,[7,5],[1.0 1.0; 0.0 0.0])
julia> f(x) = sum(x.^2)
julia> y = [f(g[i,:]) for i in axes(g,1)]
julia> w = hyperbolic_cross_weights(y,g,mi,[1.0 1.0; 0.0 0.0])
```
"""
function hyperbolic_cross_weights(y::AbstractArray{T,1},grid::Array{T,2},multi_index::Array{S,2},domain=[ones(1,size(grid,2));-ones(1,size(grid,2))]) where {T<:AbstractFloat,S<:Integer}
  
  dom = check_domain(size(multi_index,2),domain)

  # Normalize nodes to the [-1.0 1.0] interval
  
  grid = copy(grid)
  for i in axes(grid,2)
    grid[:,i] = normalize_node(grid[:,i],dom[:,i])
  end
   
  interpolation_matrix = zeros(size(grid,1),size(grid,1))
  
  n = 2*maximum(multi_index)+1
  
  unique_multi_index = sort(unique(multi_index))

  base_polynomials        = Array{Array{T,2},1}(undef,length(unique_multi_index))
  unique_base_polynomials = Array{Array{T,2},1}(undef,length(unique_multi_index))
 
  @inbounds for k in axes(grid,1)
  
    # Construct the base polynomials
  
    @inbounds for i in eachindex(unique_multi_index)
      base_polynomials[i] = chebyshev_polynomial(n-1,grid[k,:])
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

"""
Uses Chebyshev polynomials as basis functions to compute using multi-threading the weights in a hyperbolic cross 
polynomial approximation given the approximation sample, ```y```, the approximation ```grid```, the 
```multi_index```, and the approximation ```domain``` (defaults to [-1.0,1.0]^d).  Returns a vector containing
the weights in the hyperbolic cross polynomial.

Signatures
==========

w = hyperbolic_cross_weights_threaded(y,grid,multi_index)
w = hyperbolic_cross_weights_threaded(y,grid,multi_index,domain)

Example
=======
```
julia> g,mi = hyperbolic_cross_grid(chebyshev_extrema,2,3,[7,5],[1.0 1.0; 0.0 0.0])
julia> f(x) = sum(x.^2)
julia> y = [f(g[i,:]) for i in axes(g,1)]
julia> w = hyperbolic_cross_weights_threaded(y,g,mi,[1.0 1.0; 0.0 0.0])
```
"""
function hyperbolic_cross_weights_threaded(y::AbstractArray{T,1},grid::Array{T,2},multi_index::Array{S,2},domain=[ones(1,size(grid,2));-ones(1,size(grid,2))]) where {T<:AbstractFloat,S<:Integer}
  
  dom = check_domain(size(multi_index,2),domain)

  # Normalize nodes to the [-1.0 1.0] interval
  
  grid = copy(grid)
  for i in axes(grid,2)
    grid[:,i] = normalize_node(grid[:,i],dom[:,i])
  end

  interpolation_matrix = zeros(size(grid,1),size(grid,1))
  
  n = 2*maximum(multi_index)+1
  
  unique_multi_index = sort(unique(multi_index))

  @inbounds @sync Threads.@threads for k in axes(grid,1)
  
    base_polynomials        = Array{Array{T,2},1}(undef,length(unique_multi_index))
    unique_base_polynomials = Array{Array{T,2},1}(undef,length(unique_multi_index))
  
    # Construct the base polynomials
  
    @inbounds for i in eachindex(unique_multi_index)
      base_polynomials[i] = chebyshev_polynomial(n-1,grid[k,:])
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

"""
Uses Chebyshev polynomials as basis functions to compute the inverse interpolation matrix for a hyperbolic cross
polynomial approximation given the approximation grid, ```grid```, the ```multi_index```, and the approximation 
```domain``` (defaults to [-1.0,1.0]^d).  Returns a matrix containing the inverse interoplation matrix.

Signatures
==========

iim = hyperbolic_cross_inverse_interpolation_matrix(grid,multi_index)
iim = hyperbolic_cross_inverse_interpolation_matrix(grid,multi_index,domain)

Example
=======
```
julia> g,mi = hyperbolic_cross_grid(chebyshev_extrema,2,3,5,[1.0 1.0; 0.0 0.0])
julia> iim = hyperbolic_cross_inverse_interpolation_matrix(g,mi)
```
"""
function hyperbolic_cross_inverse_interpolation_matrix(grid::Array{T,2},multi_index::Array{S,2},domain=[ones(1,size(grid,2));-ones(1,size(grid,2))]) where {T<:AbstractFloat,S<:Integer}
  
  dom = check_domain(size(multi_index,2),domain)

  # Normalize nodes to the [-1.0 1.0] interval
  
  grid = copy(grid)
  for i in axes(grid,2)
    grid[:,i] = normalize_node(grid[:,i],dom[:,i])
  end

  interpolation_matrix = zeros(size(grid,1),size(grid,1))
  
  n = 2*maximum(multi_index)+1
  
  unique_multi_index = sort(unique(multi_index))

  base_polynomials        = Array{Array{T,2},1}(undef,length(unique_multi_index))
  unique_base_polynomials = Array{Array{T,2},1}(undef,length(unique_multi_index))
 
  @inbounds for k in axes(grid,1)
  
    # Construct the base polynomials
  
    @inbounds for i in eachindex(unique_multi_index)
      base_polynomials[i] = chebyshev_polynomial(n-1,grid[k,:])
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

"""
Uses Chebyshev polynomials as basis functions to compute using multi-threading the inverse interpolation matrix 
for a hyperbolic cross polynomial approximation given the approximation grid, ```grid```, the ```multi_index```,
and the approximation ```domain``` (defaults to [-1.0,1.0]^d).  Returns a matrix containing the inverse 
interoplation matrix.

Signatures
==========

iim = hyperbolic_cross_inverse_interpolation_matrix_threaded(grid,multi_index)
iim = hyperbolic_cross_inverse_interpolation_matrix_threaded(grid,multi_index,domain)

Example
=======
```
julia> g,mi = hyperbolic_cross_grid(chebyshev_extrema,2,3,5,[1.0 1.0; 0.0 0.0])
julia> iim = hyperbolic_cross_inverse_interpolation_matrix_threaded(g,mi)
```
"""
function hyperbolic_cross_inverse_interpolation_matrix_threaded(grid::Array{T,2},multi_index::Array{S,2},domain=[ones(1,size(grid,2));-ones(1,size(grid,2))]) where {T<:AbstractFloat,S<:Integer}
  
  dom = check_domain(size(multi_index,2),domain)

  # Normalize grid to the [-1.0 1.0] interval
  
  grid = copy(grid)
  for i in axes(grid,2)
    grid[:,i] = normalize_node(grid[:,i],dom[:,i])
  end

  interpolation_matrix = zeros(size(grid,1),size(grid,1))
  
  n = 2*maximum(multi_index)+1
  
  unique_multi_index = sort(unique(multi_index))
 
  @inbounds @sync Threads.@threads for k in axes(grid,1)

    base_polynomials        = Array{Array{T,2},1}(undef,length(unique_multi_index))
    unique_base_polynomials = Array{Array{T,2},1}(undef,length(unique_multi_index))
  
    # Construct the base polynomials
  
    @inbounds for i in eachindex(unique_multi_index)
      base_polynomials[i] = chebyshev_polynomial(n-1,grid[k,:])
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

"""
Uses Chebyshev polynomials as basis functions to compute the weights in a hyperbolic cross polynomial
approximation given the approximation sample, ```y```, the ```inverse interpolation matrix```.
Returns a vector containing the weights in the hyperbolic cross polynomial.

Signature
=========

w = hyperbolic_cross_weights(y,inverse_interpolation_matrix)

Example
=======
```
julia> g,mi = hyperbolic_cross_grid(chebyshev_extrema,2,3,5,[1.0 1.0; 0.0 0.0])
julia> f(x) = sum(x.^2)
julia> y = [f(g[i,:]) for i in axes(g,1)]
julia> iim = hyperbolic_cross_inverse_interpolation_matrix(g,mi)
julia> w = hyperbolic_cross_weights(y,iim)
```
"""
function hyperbolic_cross_weights(y::AbstractArray{T,1},inverse_interpolation_matrix::Array{T,2}) where {T<:AbstractFloat}
  
  weights = inverse_interpolation_matrix*y
  
  return weights
  
end

"""
Computes a hyperbolic cross polynomial given the ```point```, the ```multi-index```, and the approximation ```domain```
(defaults to [-1.0,1.0]^d).  Returns a vector containing the basis functions in the polynomial evaluated at ```point```.

Signatures
==========

hpoly = hyperbolic_cross_polynomial(point,multi_index)
hpoly = hyperbolic_cross_polynomial(point,multi_index,domain)

Example
=======
```
julia> g,mi = hyperbolic_cross_grid(chebyshev_extrema,2,3,5,[1.0 1.0; 0.0 0.0])
julia> point = g[5,:]
julia> hpoly = hyperbolic_cross_polynomial(point,mi,[1.0 1.0; 0.0 0.0])
```
"""
function hyperbolic_cross_polynomial(point::AbstractArray{R,1},multi_index::Array{S,2},domain=[ones(1,length(point));-ones(1,length(point))]) where {R<:Number,S<:Integer}
  
  if length(point) != size(domain,2)
    error("Inconsistency between the length of 'point' and the size of 'domain'.")
  end

  dom = check_domain(size(multi_index,2),domain)

  # Normalize grid to the [-1.0 1.0] interval
  
    point = copy(point)
    for i in eachindex(point)
      point[i] = normalize_node(point[i],dom[:,i])
    end  

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
      base_polynomials[i] = chebyshev_polynomial(n[i]-1,point[i])
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
        
    # Construct the hyperbolic cross polynomial

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

"""
Evaluates a hyperbolic cross polynomial formed using Chebyshev basis functions, given the ```weights```, the 
```point``` at which to evaluate the polynomial, the ```multi_index```, and the approximation ```domain``` 
(defaults to [-1.0,1.0]^d).  Returns a scalar.

Signatures
==========

yhat = hyperbolic_cross_evaluate(weights,point,multi_index)
yhat = hyperbolic_cross_evaluate(weights,point,multi_index,domain)

Example
=======
```
julia> g,mi = hyperbolic_cross_grid(chebyshev_extrema,2,3,5,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> w = hyperbolic_cross_weights(y,g,mi,[1.0 1.0; 0.0 0.0])
julia> yhat = hyperbolic_cross_evaluate(w,[0.37,0.71],mi,[1.0 1.0; 0.0 0.0])
0.6100559317314179
```
"""
function hyperbolic_cross_evaluate(weights::Array{T,1},point::AbstractArray{R,1},multi_index::Array{S,2},domain=[ones(1,length(point));-ones(1,length(point))]) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  if length(point) != size(domain,2)
    error("Inconsistency between the length of 'point' and the size of 'domain'.")
  end

  dom = check_domain(size(multi_index,2),domain)

  point = copy(point)

  for i in eachindex(point)
    point[i] = normalize_node(point[i],dom[:,i])
  end

  poly = hyperbolic_cross_polynomial(point,multi_index)

  estimate = weights'poly

  return estimate

end

"""
Evaluates a hyperbolic cross polynomial formed using Chebyshev basis functions, given the ```weights``` and the
hyperbolic cross polynomial.  Returns a scalar.

Signature
=========

yhat = hyperbolic_cross_evaluate(weights,poly)

Example
=======
```
julia> g,mi = hyperbolic_cross_grid(chebyshev_extrema,2,3,5,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> w = hyperbolic_cross_weights(y,g,mi,[1.0 1.0; 0.0 0.0])
julia> hpoly = hyperbolic_cross_polynomial([0.37,0.71],mi,[1.0 1.0; 0.0 0.0])
julia> yhat = hyperbolic_cross_evaluate(w,hpoly)
0.6100559317314179
```
"""
function hyperbolic_cross_evaluate(weights::Array{T,1},polynomial::Array{R,1}) where {T<:AbstractFloat,R<:Number}

    estimate = weights'polynomial

    return estimate

end

"""
Creates an interpolating function for a hyperbolic cross approximation given the sampling points, ```y```, and the
approximation ```plan```.  Returns an interpolating function.

Signature
=========

f = hyperbolic_cross_interp(y,plan)

Example
=======
```
julia> g,mi = hyperbolic_cross_grid(chebyshev_extrema,2,3,7,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> hplan = hyperbolic_cross_plan(chebyshev_extrema,2,3,[7,7],[1.0 1.0; 0.0 0.0])
julia> f = hyperbolic_cross_interp(y,hplan)
julia> f([0.37,0.71])
0.5638586549437409
```
"""
function hyperbolic_cross_interp(y::AbstractArray{T,1},plan::P) where {T<:AbstractFloat,P<:HCApproximationPlan}

  weights = hyperbolic_cross_weights(y,plan.grid,plan.multi_index,plan.domain)
  
  function hcross_interp(x::Array{R,1}) where {R<:Number}
  
    return hyperbolic_cross_evaluate(weights,x,plan.multi_index,plan.domain)
  
  end
  
  return hcross_interp
  
end
  
"""
Uses multi-threading to create an interpolating function for a hyperbolic cross approximation given the sampling
points, ```y```, and the approximation ```plan```.  Returns an interpolating function.

Signature
=========

f = hyperbolic_cross_interp_threaded(y,plan)

Example
=======
```
julia> g,mi = hyperbolic_cross_grid(chebyshev_extrema,2,3,7,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> hplan = hyperbolic_cross_plan(chebyshev_extrema,2,3,[7,7],[1.0 1.0; 0.0 0.0])
julia> f = hyperbolic_cross_interp_threaded(y,hplan)
julia> f([0.37,0.71])
0.5638586549437409
```
"""
function hyperbolic_cross_interp_threaded(y::AbstractArray{T,1},plan::P) where {T<:AbstractFloat,P<:HCApproximationPlan}
  
  weights = hyperbolic_cross_weights_threaded(y,plan.grid,plan.multi_index,plan.domain)
  
  function hcross_interp(x::Array{R,1}) where {R<:Number}
  
    return hyperbolic_cross_evaluate(weights,x,plan.multi_index,plan.domain)
  
  end
  
  return hcross_interp
  
end

function _hyperbolic_cross_derivative(weights::Array{T,1},point::Array{R,1},multi_index::Array{S,2},pos::S) where {T<:AbstractFloat,R<:Number,S<:Integer} # Internal function, not exported
  
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

"""
Computes the partial derivative of a hyperbolic cross polynomial formed using Chebyshev basis functions, given the 
```weights```, the ```point``` at which to evaluate the polynomial, the ```multi_index```, the approximation 
```domain```, and the index of the variable to differentiate with respect to, ```pos```.  Returns a scalar.

Signature
=========

deriv = hyperbolic_cross_derivative(weights,point,multi_index,domain,pos)

Example
=======
```
julia> g,mi = hyperbolic_cross_grid(chebyshev_extrema,2,3,7,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> w = hyperbolic_cross_weights(y,g,mi,[1.0 1.0; 0.0 0.0])
julia> deriv1 = hyperbolic_cross_derivative(w,[0.37,0.71],mi,[1.0 1.0; 0.0 0.0],1)
0.5157820877272242
julia> deriv2 = hyperbolic_cross_derivative(w,[0.37,0.71],mi,[1.0 1.0; 0.0 0.0],2)
0.5514971149820203
```
"""
function hyperbolic_cross_derivative(weights::Array{T,1},point::Array{R,1},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}},pos::S) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  if length(point) != size(domain,2)
    error("Inconsistency between the length of 'point' and the size of 'domain'.")
  end

  dom = check_domain(size(multi_index,2),domain)

  point = copy(point)
 
  d = length(point)
  for i = 1:d
    point[i] = normalize_node(point[i],dom[:,i])
  end
  
  evaluated_derivative = _hyperbolic_cross_derivative(weights,point,multi_index,pos)
  
  return evaluated_derivative*(2.0/(dom[1,pos]-dom[2,pos]))
  
end

"""
Computes the gradient a hyperbolic cross polynomial formed using Chebyshev basis functions, given the
```weights```, the ```point``` at which to evaluate the polynomial, the ```multi_index```, and the
approximation ```domain``` (defaults to [-1.0,1.0]^d).  Returns a one-row matrix.

Signatures
==========

grad = hyperbolic_cross_gradient(weights,point,multi_index)
grad = hyperbolic_cross_gradient(weights,point,multi_index,domain)

Example
=======
```
julia> g,mi = hyperbolic_cross_grid(chebyshev_extrema,2,3,7,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> w = hyperbolic_cross_weights(y,g,mi,[1.0 1.0; 0.0 0.0])
julia> grad = hyperbolic_cross_gradient(w,[0.37,0.71],mi,[1.0 1.0; 0.0 0.0])
[0.515782  0.551497]
```
"""
function hyperbolic_cross_gradient(weights::Array{T,1},point::Array{R,1},multi_index::Array{S,2},domain=[ones(1,length(point));-ones(1,length(point))]) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  if length(point) != size(domain,2)
    error("Inconsistency between the length of 'point' and the size of 'domain'.")
  end

  dom = check_domain(size(multi_index,2),domain)

  d = length(point)
  gradient = Array{R,2}(undef,1,d)
  
  for i = 1:d
    gradient[i] = hyperbolic_cross_derivative(weights,point,multi_index,dom,i)
  end
  
  return gradient

end

"""
Creates an interpolating function to compute the gradient of a hyperbolic cross polynomial given the sampling
points, ```y```, and the approximation ```plan```.  Returns an interpolating function.

Signature
=========

grad = hyperbolic_cross_gradient(y,plan)

Example
=======
```
julia> g,mi = hyperbolic_cross_grid(chebyshev_extrema,2,3,7,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> hplan = hyperbolic_cross_plan(chebyshev_extrema,2,3,[7,7],[1.0 1.0; 0.0 0.0])
julia> grad = hyperbolic_cross_gradient(y,hplan)
julia> grad([0.37,0.71])
[0.515782  0.551497]
```
"""
function hyperbolic_cross_gradient(y::AbstractArray{T,1},plan::P) where {T<:AbstractFloat,P<:HCApproximationPlan}

  weights = hyperbolic_cross_weights(y,plan.grid,plan.multi_index,plan.domain)
  
  function hcross_grad(x::Array{R,1}) where {R<:Number}
  
    return hyperbolic_cross_gradient(weights,x,plan.multi_index,plan.domain)
  
  end
  
  return hcross_grad
  
end

"""
Creates an interpolating function that uses multi-threading to compute the gradient of a hyperbolic cross polynomial
  given the sampling points, ```y```, and the approximation ```plan```.  Returns an interpolating function.

Signature
=========

grad = hyperbolic_cross_gradient_threaded(y,plan)

Example
=======
```
julia> g,mi = hyperbolic_cross_grid(chebyshev_extrema,2,3,7,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> hplan = hyperbolic_cross_plan(chebyshev_extrema,2,3,[7,7],[1.0 1.0; 0.0 0.0])
julia> grad = hyperbolic_cross_gradient_threaded(y,hplan)
julia> grad([0.37,0.71])
[0.515782  0.551497]
```
"""
function hyperbolic_cross_gradient_threaded(y::AbstractArray{T,1},plan::P) where {T<:AbstractFloat,P<:HCApproximationPlan}
  
  weights = hyperbolic_cross_weights_threaded(y,plan.grid,plan.multi_index,plan.domain)
  
  function hcross_grad(x::Array{R,1}) where {R<:Number}
  
    return hyperbolic_cross_gradient(weights,x,plan.multi_index,plan.domain)
  
  end
  
  return hcross_grad
  
end

"""
Computes the hessian of a hyperbolic cross polynomial formed using Chebyshev basis functions, given the 
```weights```, the ```point``` at which to evaluate the polynomial, the ```multi_index```, and the
approximation ```domain``` (defaults to [-1.0,1.0]^d).  Returns a matrix.

Signatures
==========

hess = hyperbolic_cross_hessian(weights,point,multi_index)
hess = hyperbolic_cross_hessian(weights,point,multi_index,domain)

Example
=======
```
julia> g,mi = hyperbolic_cross_grid(chebyshev_extrema,2,11,23,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> w = hyperbolic_cross_weights(y,g,mi,[1.0 1.0; 0.0 0.0])
julia> hess = hyperbolic_cross_hessian(w,[0.37,0.71],mi,[1.0 1.0; 0.0 0.0])
[-4.0189     0.466767
  0.466767  -0.20577]
```
"""
function hyperbolic_cross_hessian(weights::Array{T,1},point::Array{R,1},multi_index::Array{S,2},domain=[ones(1,length(point));-ones(1,length(point))]) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  if length(point) != size(domain,2)
    error("Inconsistency between the length of 'point' and the size of 'domain'.")
  end

  dom = check_domain(size(multi_index,2),domain)

  point = copy(point)

  if size(domain,2) != length(point)
    error("domain is inconsistent with the number of dimensions")
  end
  
  for i in eachindex(point)
    point[i] = normalize_node(point[i],dom[:,i])
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
  
    hess[c] = evaluated_derivative*(2.0/(dom[1,c[1]]-dom[2,c[1]]))*(2.0/(dom[1,c[2]]-dom[2,c[2]]))

  end
  
  return hess
  
end

"""
Creates an interpolating function to compute the hessian of a hyperbolic cross polynomial given the sampling
points, ```y```, and the approximation ```plan```.  Returns an interpolating function.

Signature
=========

grad = hyperbolic_cross_hessian(y,plan)

Example
=======
```
julia> g,mi = hyperbolic_cross_grid(chebyshev_extrema,2,11,23,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> hplan = hyperbolic_cross_plan(chebyshev_extrema,2,11,[23,23],[1.0 1.0; 0.0 0.0])
julia> hess = hyperbolic_cross_hessian(y,hplan)
julia> hess([0.37,0.71])
[-4.0189     0.466767
  0.466767  -0.20577]
```
"""
function hyperbolic_cross_hessian(y::AbstractArray{T,1},plan::P) where {T<:AbstractFloat,P<:HCApproximationPlan}

  weights = hyperbolic_cross_weights(y,plan.grid,plan.multi_index,plan.domain)
  
  function hcross_hess(x::Array{R,1}) where {R<:Number}
  
    return hyperbolic_cross_hessian(weights,x,plan.multi_index,plan.domain)
  
  end
  
  return hcross_hess
  
end

"""
Creates an interpolating function that uses multi-threading to compute the hessian of a hyperbolic cross polynomial
  given the sampling points, ```y```, and the approximation ```plan```.  Returns an interpolating function.

Signature
=========

grad = hyperbolic_cross_hessian_threaded(y,plan)

Example
=======
```
julia> g,mi = hyperbolic_cross_grid(chebyshev_extrema,2,11,23,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> hplan = hyperbolic_cross_plan(chebyshev_extrema,2,11,[23,23],[1.0 1.0; 0.0 0.0])
julia> hess = hyperbolic_cross_hessian_threaded(y,hplan)
julia> hess([0.37,0.71])
[-4.0189     0.466767
  0.466767  -0.20577]
```
"""
function hyperbolic_cross_hessian_threaded(y::AbstractArray{T,1},plan::P) where {T<:AbstractFloat,P<:HCApproximationPlan}
  
  weights = hyperbolic_cross_weights_threaded(y,plan.grid,plan.multi_index,plan.domain)
  
  function hcross_hess(x::Array{R,1}) where {R<:Number}
  
    return hyperbolic_cross_hessian(weights,x,plan.multi_index,plan.domain)
  
  end
  
  return hcross_hess
  
end

function integrate_cheb_polys(order::S) where {S <: Integer} # Internal function, not exported

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

"""
Numerically integrates a function, ```f```, over all dimensions by approximating the function with a hyperbolic
cross polynomial according to the approximation ```plan```, using either :clenshaw_curtis or gauss_Chebyshev_quad 
as ```method```.  Returns a scalar.

Signature
=========

integral = hyperbolic_cross_integrate(f,plan,method)

Example
=======
```
julia> f(x) = sum(x.^2)
julia> hplan = hyperbolic_cross_plan(chebyshev_extrema,4,8,[17,17,17,17],[ones(1,4); zeros(1,4)])
julia> integral = hyperbolic_cross_integrate(f,hplan,:clenshaw_curtis)
1.3333333333336967
julia> integral = hyperbolic_cross_integrate(f,hplan,:gauss_chebyshev_quad)
0.6897405247069276
```
"""
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

"""
Uses the Clenshaw-Curtis method to numerically integrates a function, ```f```, over all dimensions by approximating
the function with a hyperbolic cross polynomial according to the approximation ```plan```.  Returns a scalar.

Signature
=========

integral = hyperbolic_cross_clenshaw_curtis(f,plan)

Example
=======
```
julia> f(x) = sum(x.^2)
julia> hplan = hyperbolic_cross_plan(chebyshev_extrema,4,8,[17,17,17,17],[ones(1,4); zeros(1,4)])
julia> integral = hyperbolic_cross_clenshaw_curtis(f,hplan)
1.3333333333336967
```
"""
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

"""
Uses the Clenshaw-Curtis method to numerically integrates a function, ```f```, over all dimensions except ```pos```
by approximating the function with a hyperbolic cross polynomial according to the approximation ```plan```.
Returns a function.

Signature
=========

integral = hyperbolic_cross_clenshaw_curtis(f,plan,pos)

Example
=======
```
julia> f(x) = sum(x.^2)
julia> hplan = hyperbolic_cross_plan(chebyshev_extrema,4,8,[17,17,17,17],[ones(1,4); zeros(1,4)])
julia> integral = hyperbolic_cross_clenshaw_curtis(f,hplan,4)
julia> integral(0.5)
1.2500000000001648
```
"""
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

"""
Uses Gauss-Chebyshev quadrature to numerically integrates a function, ```f```, over all dimensions by approximating
the function with a hyperbolic cross polynomial according to the approximation ```plan```.  Returns a scalar.

Signature
=========

integral = hyperbolic_cross_gauss_chebyshev_quad(f,plan)

Example
=======
```
julia> f(x) = sum(x.^2)
julia> hplan = hyperbolic_cross_plan(chebyshev_extrema,4,8,[17,17,17,17],[ones(1,4); zeros(1,4)])
julia> integral = hyperbolic_cross_gauss_chebyshev_quad(f,hplan)
0.6897405247069276
```
"""
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