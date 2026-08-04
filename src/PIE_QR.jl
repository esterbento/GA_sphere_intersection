using LinearAlgebra
using DelimitedFiles

function PIE_QR(nome_arquivo::String)
    linhas = readlines(nome_arquivo)

    # n: dimensão do espaço
    n = parse(Int, strip(linhas[2]))

    # P: matriz de pontos (cada coluna representa um ponto)
    bloco_pontos = join(linhas[4:3+n], "\n")
    P = readdlm(IOBuffer(bloco_pontos), ' ')

    # r: vetor dos raios
    linha_raios = replace(replace(linhas[5+n], '[' => ""), ']' => "")
    r = parse.(Float64, split(linha_raios, ","))

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

    for i = 1:n-1
        c[i] = -0.5 * (r[i]^2 - r[n]^2 - norm(A[:, i], 2)^2)
    end

    # Resolve o sistema linear da equação (2.13)
    y = (Rp') \ c

    # Verifica as condições de interseção
    if r[n]^2 - norm(y, 2)^2 < 0.0
        println("A interseção é vazia.")
        return nothing

    elseif abs(r[n]^2 - norm(y, 2)^2) < 1e-10
        z = 0

    else
        # Se o quadrado do raio for maior que o quadrado da norma de y,
        # calcula o valor de z
        z = sqrt(r[n]^2 - norm(y, 2)^2)
    end

    # Calcula os pontos de interseção
    x1 = F.Q * [y; z] + P'[:, n]   # Primeiro ponto de interseção
    x2 = F.Q * [y; -z] + P'[:, n]  # Segundo ponto de interseção

    # Retorna uma ou duas soluções, dependendo de elas coincidirem ou não
    if norm(x1 - x2) < 1.0e-10
        return x1
    else
        return x1, x2
    end
end

# Observação:
# Este algoritmo foi adaptado do artigo de Maioli.
# A numeração das equações, como (2.11), (2.12) e (2.13),
# segue a notação utilizada no artigo.