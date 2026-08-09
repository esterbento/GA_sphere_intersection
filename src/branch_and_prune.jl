using LinearAlgebra 

function viavel_BP(candidato, i, x, D)
    n = length(candidato)

    for j in 1:i-1

        # Calcula a distância sem criar o vetor
        # candidato - x[j].
        distancia_quadrada = 0.0

        for k in 1:n
            diferenca = candidato[k] - x[j][k]
            distancia_quadrada += diferenca * diferenca
        end

        distancia = sqrt(distancia_quadrada)

        if !isapprox(distancia, D[j, i];atol = 1.0e-3, rtol = 1.0e-8)
            return false
        end
    end

    return true
end


function BP!(i, n, N, D, x, solucoes, C, r)
    # Todos os pontos foram posicionados.
    if i > N
        push!(solucoes, deepcopy(x))
        return
    end

    # Preenche C e r com os n predecessores imediatos.
    primeiro_predecessor = i - n

    for k in 1:n
        predecessor = primeiro_predecessor + k - 1
        C[k, :] .= x[predecessor]
        r[k] = D[predecessor, i]
    end

    # Calcula os candidatos para o ponto i.
    sol1, sol2 = PIE_CGA(C, r)

    for candidato in (sol1, sol2)
        if viavel_BP(candidato, i, x, D)
            x[i] .= candidato
            BP!(i + 1, n, N, D, x, solucoes, C, r)
        end
    end

    return
end



function read_problem_BP(filename::String)

    lines = strip.(readlines(filename))

    dimension_line =
        findfirst(==("Dimension ="), lines)

    number_points_line =
        findfirst(==("Number of points ="), lines)

    initial_points_line =
        findfirst(==("Initial points ="), lines)

    matrix_line =
        findfirst(==("Distance matrix ="), lines)

    solution_line =
        findfirst(==("Reference solution ="), lines)

    n = parse(Int, lines[dimension_line + 1])
    N = parse(Int, lines[number_points_line + 1])

    # Primeiros n pontos usados para iniciar o BP.
    pontos_iniciais =
        Vector{Vector{Float64}}(undef, n)

    for i in 1:n
        valores = parse.( Float64, split(lines[initial_points_line + i]))

        if length(valores) != n
            error( "O ponto inicial $i possui $(length(valores)) " * "coordenadas, mas deveria possuir $n.")
        end

        pontos_iniciais[i] = valores
    end

    # Matriz de distâncias.
    D = Matrix{Float64}(undef, N, N)

    for i in 1:N
        valores = parse.( Float64, split(lines[matrix_line + i]))

        if length(valores) != N
            error("A linha $i da matriz possui $(length(valores)) " * "elementos, mas deveria possuir $N.")
        end

        D[i, :] .= valores
    end

    # Realização usada na criação do problema.
    Psol = Vector{Vector{Float64}}(undef, N)

    for i in 1:N
        valores = parse.(Float64,split(lines[solution_line + i]))

        if length(valores) != n
            error("O ponto $i da solução possui $(length(valores)) " * "coordenadas, mas deveria possuir $n.")
        end

        Psol[i] = valores
    end

    return n, N, pontos_iniciais, D, Psol
end


function solve_problem_BP(n, N, pontos_iniciais, D)

    x = [zeros(Float64, n) for _ in 1:N]

    # Fixamos os primeiros n pontos.
    for i in 1:n
        x[i] .= pontos_iniciais[i]
    end

    solucoes = Vector{Vector{Vector{Float64}}}()

    # C e r são criados apenas uma vez e reutilizados
    # durante toda a recursão.
    C = Matrix{Float64}(undef, n, n)
    r = Vector{Float64}(undef, n)

    BP!(n + 1, n, N, D, x, solucoes, C, r)

    return solucoes
end