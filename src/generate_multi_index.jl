function generate_multi_index(d::S,k::S) where {S <: Integer} # Recursive function

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
       
function generate_multi_index(d::S,k::S,n::S) where {S <: Integer} # Recursive function

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
  
function generate_multi_index(d::S,k::S,n::Array{S,1}) where {S <: Integer} # Recursive function

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
