function determine_grid_size(p::Array{S,2}) where {S <: Integer}
    
    q = similar(p)
  
    for i in eachindex(p)
      if p[i] > 0
        q[i] = 2
      else
        q[i] = 1
      end
    end
  
    s = 0
    for i = 1:size(q,1)
      t = 1
      for j = 1:size(q,2)
        t *= q[i,j]
      end
      s += t
    end
      
    return (s, size(q,2))
     
  end
  