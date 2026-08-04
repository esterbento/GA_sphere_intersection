using LinearAlgebra
using DelimitedFiles

# C = matriz formada pelos centros dispostos em linha
# r = vetor dos raios
function PIE_CGA(filename::String)

    lines = readlines(filename)

    # n: dimensão do espaço
    n = parse(Int, strip(lines[2]))

    # C: matriz dos centros das esferas
    points_block = join(lines[4:3+n], "\n")
    C = Matrix{Float64}(readdlm(IOBuffer(points_block), ' '))

    # r: vetor dos raios
    r_line = replace(replace(lines[5+n], '[' => ""), ']' => "")
    r = parse.(Float64, split(r_line, ","))

    # sol: expected solution vector
    sol_line = replace(replace(lines[7+n], '[' => ""), ']' => "")
    sol = parse.(Float64, split(sol_line, ","))

    d = copy(C[n, :])
    x = zeros(n + 1)
    w = zeros(n-1)
    
    for i = 1:n-1
        C[i, :] .= C[i, :] - d
    end
    C[n, :] .= 0.0

    b = zeros(n)
    rn2 = r[n]^2
    b[n] = 0.5 * - rn2
    
    for k = 1:n-1
        b[k] = 0.5 * (dot(C[k, :], C[k, :]) - r[k]^2) - b[n]
    end
    x[n+1] = -b[n]

    A = C[1:n-1, 1:n-1]
    a = C[1:n-1, n]

    F = lu(A)
    x[1:n-1] .= F \ b[1:n-1]
    w .=x[1:n-1]
    v = F \ -a

    coefc = dot(x[1:n-1], x[1:n-1]) - rn2
    coefb = 2.0 * dot(x[1:n-1], v)
    coefa = dot(v, v) + 1.0
    
    Δ = coefb^2 - 4.0 * coefa * coefc
    if Δ<0.0
        error("There is not an intersection")
    end
    sqrtΔ=sqrt(Δ)
    den= 2.0 * coefa 

    # primeira solução
    x[n] = (-coefb - sqrtΔ) / den
    x[1:n-1] .= w + x[n] * v
    sol1 = x[1:n] + d

    # segunda solução 
    x[n] = (-coefb + sqrtΔ) / den
    x[1:n-1] .= w + x[n] * v
    sol2 = x[1:n] + d
    return sol1, sol2, sol
end

function PIE_CGA1(C,r)
    n=length(r)
    d = copy(C[n, :])
    x = zeros(n + 1)
    w = zeros(n-1)
    b = zeros(Float64, n)

    for i = 1:n-1
        C[i, :] .= C[i, :] - d
    end
    C[n, :] .= 0.0

    
    rn2 = r[n]^2
    b[n] = -0.5 * rn2

    for k = 1:n-1
        b[k] = 0.5 * (dot(C[k, :], C[k, :]) - r[k]^2) - b[n]
    end

    A = C[1:n-1, 1:n-1]
    a = C[1:n-1, n]

    F = lu(A)
    x[1:n-1] .= F \ b[1:n-1]
    w .= x[1:n-1]
    v = F \ -a
    coefc = dot(w,w) - rn2
    coefb = 2.0 * dot(w, v)
    coefa = dot(v, v) + 1.0

    Δ = coefb^2 - 4.0 * coefa * coefc
    if Δ < 0.0
        error("There is not an intersection")
    end

    sqrtΔ=sqrt(Δ)
    den= 2.0 * coefa 

    # primeira solução
    x[n] = (-coefb - sqrtΔ) / den
    x[1:n-1] .= w + x[n] * v
    sol1 = x[1:n] + d
    

    # segunda solução 
    x[n] = (-coefb + sqrtΔ) / den
    x[1:n-1] .= w + x[n] * v
    sol2 = x[1:n] + d
    return sol1, sol2
end
