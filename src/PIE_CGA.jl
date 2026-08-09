using LinearAlgebra
# C = matriz formada pelos centros dispostos em linha
# r = vetor dos raios
function PIE_CGA(C::Array{Float64,2}, r::Vector{Float64}, n::Int)

    d = copy(C[n, :])

    for i = 1:n-1
        C[i, :] .= C[i, :] - d
    end

    b = zeros(n-1)
    xnpo = 0.5 * (r[n]^2)

    for k = 1:n-1
        b[k] = 0.5 * (dot(C[k, :], C[k, :]) - r[k]^2) + xnpo
    end

    A = C[1:n-1, 1:n-1]
    a = C[1:n-1, n]

    F = lu(A)
    w = F \ b
    v = F \ -a

    coefc = dot(w,w) - 2.0 * xnpo
    coefb = 2.0 * dot(w, v)
    coefa = dot(v, v) + 1.0
    Δ = coefb^2 - 4.0 * coefa * coefc
    if Δ<0.0
        error("There is not an intersection")
    end
    # primeira solução
    sq =  sqrt(Δ)
    xn = (-coefb - sq) / (2.0 * coefa)
    x1 = [w + xn*v; xn] .+ d
    # segunda solução
    xn = (-coefb + sq) / (2.0 * coefa)
    x2 =  [w + xn*v; xn] .+ d

    if norm(x1 - x2) < 1.0e-10
        return x1
    else 
        return x1,x2
    end
end
