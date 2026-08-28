using LinearAlgebra

function PIE_QR(P::Array{Float64,2}, r::Vector{Float64}, n::Int)
    # Constrói a matriz A = [a_1 - a_n, ..., a_(n-1) - a_n]
    A = zeros(n, n-1)

    for i = 1:n-1
        A[:, i] = P'[:, i] - P'[:, n]
    end

    # Calcula a decomposição QR da matriz A
    F = qr(A)

    # Determina a matriz Rp, conforme a equação (2.12)
    Rp = F.R[1:n-1, :]

    # Calcula o vetor c, conforme a equação (2.11)
    c = zeros(n-1)

    rn2=r[n]^2
    for i = 1:n-1
        c[i] = -0.5 * (r[i]^2 - rn2 - norm(A[:, i], 2)^2)
    end

    # Resolve o sistema linear da equação (2.13)
    y = (Rp') \ c

    delta = rn2 - norm(y, 2)^2
    # Verifica as condições de interseção
    if delta < 0.0
        println("A interseção é vazia.")
        return nothing

    elseif abs(delta) < 1e-10
        z = 0

    else
        # Se o quadrado do raio for maior que o quadrado da norma de y,
        # calcula o valor de z
        z = sqrt(delta)
    end

    # Calcula os pontos de interseção
    x1 = F.Q * [y; z] + P'[:, n]   # Primeiro ponto de interseção
    x2 = F.Q * [y; -z] + P'[:, n]  # Segundo ponto de interseção

    # Retorna uma ou duas soluções, dependendo de elas coincidirem ou não
    d = x1 - x2

    if dot(d, d) < 1.0e-20
        return x1
    else
        return x1, x2
    end
end

# Observação:
# Este algoritmo foi adaptado do artigo de Maioli.
# A numeração das equações, como (2.11), (2.12) e (2.13),
# segue a notação utilizada no artigo.

function PIE_LU(P::Array{Float64,2}, r::Vector{Float64}, n::Int)
   
    # Calcula o vetor b
    b = zeros(n)

    for k = 1:n
        b[k] = dot(P[k, :], P[k, :]) - r[k]^2 
    end

    e = ones(n)

    F = lu(P)
    u = F\e
    v = F\b

    # Produtos internos utilizados na equação quadrática.
    p = 2-dot(u, v)
    q = dot(u, u)
    s = dot(v, v)

    # Discriminante.
    Δ = (p)^2 - q * s

    if Δ < 0.0
        error("There is not an intersection")
    end

    sq = sqrt(Δ)

    # Valores de rho.
    rho1 = (p- sq) / q
    rho2 = (p - p + sq) / q

    # Pontos de interseção.
    x1 = 0.5 * (rho1 * u + v)
    x2 = 0.5 * (rho2 * u + v)

    # Verifica se as duas soluções coincidem.
    d = x1 - x2

    if dot(d, d) < 1.0e-20
        return x1
    else
        return x1, x2
    end
end

