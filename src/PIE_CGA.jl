using LinearAlgebra
# C = matriz formada pelos centros dispostos em linha
# r = vetor dos raios
function PIE_CGA(C::Array{Float64,2}, r::Vector{Float64})
    n = length(r)
    d = copy(C[n, :])
    
    flag = :twosol

    for i = 1:n-1
        C[i, :] .= C[i, :] - d
    end
    
    b = zeros(n-1)
    xnpo = 0.5 * (r[n]^2)

    
    for k = 1:n-1
        b[k] = 0.5 * (dot(C[k, :], C[k, :]) - r[k]^2) + xnpo
    end
    
    A = copy(C[1:n-1,:])
    
    F = lu(A,check = false)
    
    if issuccess(F)
        freeidx = n 
    else
        freeidx = findfirst(i -> abs(F.U[i, i]) < 1e-12, 1:n)     
    end
    a = C[1:n-1, freeidx]
    cols = setdiff(1:n, freeidx)

    
    U = F.U[:,cols]
    w = U \ (F.L \ b[F.p])
    v = U \ (-F.U[:,freeidx])

    coefc = dot(w,w) - 2.0 * xnpo
    coefb = 2.0 * dot(w, v)
    coefa = dot(v, v) + 1.0

    Δ = coefb^2 - 4.0 * coefa * coefc
    if Δ<0.0
        flag = :nosol
        return flag, Nothing, Nothing
    end

    (Δ < sqrt(eps(Float64))) && (flag = :onesol)

    sq = sqrt(Δ)
    # --------------------------------------------------
    # Primeira solução
    # --------------------------------------------------

    xfree = (-coefb - sq) / (2.0 * coefa)
    x1 = zeros(n)
    x1[cols] = w + xfree * v
    x1[freeidx] = xfree

    x1 .+= d

    # --------------------------------------------------
    # Segunda solução
    # --------------------------------------------------

    xfree = (-coefb + sq) / (2.0 * coefa)

    x2 = zeros(n)

    x2[cols] = w + xfree * v
    x2[freeidx] = xfree

    x2 .+= d

    # Se as duas soluções forem praticamente iguais,
    # retorna apenas uma
    return flag, x1, x2

end
