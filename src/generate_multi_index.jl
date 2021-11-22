function generate_multi_index(d::S,k::S) where {S <: Integer}

  if d < 1
    error("d must be positive")
  end

  if k < 0
    error("k must be non-negative")
  end

  kprime = k+1
  n = 2*k+1

  # First determine the number of dimensions that will contain a non-zero element
  dd = 1
  for i = 2:d
    if Int(floor(kprime/(2^(i-1)))) > 1
      dd += 1
    else
      break
    end
  end

  N = kprime^dd
  
  nt = (kprime+1)*dd # This is intended to be an over-estimate
  ind = zeros(Int64,nt,d)
  ddd = [kprime for _ = 1:dd]
  ci = CartesianIndices(Tuple(ddd))
  pos = 1
  @floop for i = 2:N
    flag = true
    ii = collect(Tuple(ci[i]))
    if prod(ii) <= kprime
      for j = 1:pos
        if sort(ii.-1) == sort(ind[j,1:dd])
          flag = false
          break
        end
      end
      if flag == true
        pos += 1
        try
          ind[pos,1:dd] = ii .- 1
        catch
          ind = [ind; collect(ii.-1)']
        end
      end
    end
  end

  mi = generate_multi_index(d,ind[1:pos,:])

 return mi
    
end

function generate_multi_index(d::S,ind::Array{S,2}) where {S <: Integer}

  if d < 1
    error("d must be positive")
  end

  (n1,n2) = size(ind)

  if n2 != d
    error("'ind' must have 'd' columns")
  end

  nt = d*(maximum(ind)^2+1) # This is intended to be an over-estimate
  mi = zeros(Int64,nt,d)
  pos = 0
  @floop @inbounds for i = 1:n1
    p = unique(permutations(ind[i,:],d))
    for l = 1:length(p)
      pos += 1
      try
        mi[pos:pos,:] .= p[l]'
      catch
        mi = [mi; p[l]']
      end
    end
  end
  
  return mi[1:pos,:]
    
end

function generate_multi_index(d::S,k::S,n::S) where {S <: Integer}

  if n == 2*k+1
    return generate_multi_index(d,k)
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

  kprime = k + 1

  ki = Int((n-1)/2) + 1
  N  = ki^d

  nt = d*(ki^2+1) # This is intended to be an over-estimate
  mi = zeros(Int64,nt,d)
  pos = 1
  dd = [ki for _ = 1:d]
  ci = CartesianIndices(Tuple(dd))
  @floop @inbounds for i in 2:N
    ii = Tuple(ci[i])
    if prod(ii) <= kprime
      pos += 1
      try
        mi[pos,:] .= ii.-1
      catch
        mi = [mi; collect(ii.-1)']
      end
    end
  end
  
  return mi[1:pos,:]
    
end

function generate_multi_index(d::S,k::S,n::Array{S,1}) where {S <: Integer}

  if k < 0
    error("k must be non-negative")
  end

  if length(n) != d
    error("node number for each dimension must be given")
  end

  for i = 1:d
    if !isodd(n[i]) 
      error("each element in n must be odd")
    end
  end

  kprime = k + 1

  ki = Int.((n.-1)/2) .+ 1
  N  = prod(ki)

  nt = d*(maximum(ki)^2+1) # This is an approximation that is intended to be an over-estimate
  mi = zeros(Int64,nt,d)
  pos = 1
  ci = CartesianIndices(Tuple(ki))
  @floop @inbounds for i in 2:N
    ii = Tuple(ci[i])
    if prod(ii) <= kprime
      pos += 1
      try
        mi[pos,:] .= ii.-1
      catch
        mi = [mi; collect(ii.-1)']
      end
    end
  end

  return mi[1:pos,:]

end
