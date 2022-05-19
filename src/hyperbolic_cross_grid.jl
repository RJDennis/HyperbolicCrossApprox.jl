function hyperbolic_cross_grid(node_type::Function,d::S,ind::Array{S,2}) where {S <: Integer}
  
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

function hyperbolic_cross_grid(node_type::Function,d::S,k::S) where {S <: Integer}
  
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

function hyperbolic_cross_grid(node_type::Function,d::S,k::S,n::S) where {S <: Integer}
  
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

function hyperbolic_cross_grid(node_type::Function,d::S,k::S,n::Array{S,1}) where {S <: Integer}
  
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

function hyperbolic_cross_grid(node_type::Function,ind::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {S <: Integer, T <: AbstractFloat}
  
  if size(domain,2) != d
    error("domain is inconsistent with the number of dimensions")
  end
    
  (nodes, multi_index) = hyperbolic_cross_grid(node_type,ind)
    
  # Now scale the nodes to the desired domain
    
  nodes .= scale_nodes(nodes,domain)
    
  return nodes, multi_index
    
end

function hyperbolic_cross_grid(node_type::Function,d::S,k::S,domain::Union{Array{T,1},Array{T,2}}) where {S <: Integer, T <: AbstractFloat}
  
  if size(domain,2) != d
    error("domain is inconsistent with the number of dimensions")
  end
    
  (nodes, multi_index) = hyperbolic_cross_grid(node_type,d,k)
    
  # Now scale the nodes to the desired domain
    
  nodes .= scale_nodes(nodes,domain)
    
  return nodes, multi_index
    
end

function hyperbolic_cross_grid(node_type::Function,d::S,k::S,n::Union{S,Array{S,1}},domain::Union{Array{T,1},Array{T,2}}) where {S <: Integer, T <: AbstractFloat}
  
  if size(domain,2) != d
    error("domain is inconsistent with the number of dimensions")
  end
    
  (nodes, multi_index) = hyperbolic_cross_grid(node_type,d,k,n)
    
  # Now scale the nodes to the desired domain
    
  nodes .= scale_nodes(nodes,domain)
    
  return nodes, multi_index
    
end
