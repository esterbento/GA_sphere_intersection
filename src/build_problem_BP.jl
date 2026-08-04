using LinearAlgebra
using Random

function build_problem_BP(n::Int)
    # n pontos iniciais e m pontos que serão encontrados pelo BP.
    m = 3*n
    number_of_points = n + m

    # Cada linha representa um ponto de Rⁿ.
    P = randn(number_of_points, n)

    # Matriz completa de distâncias reais.
    # Sua função intersecoes1 recebe raios e depois calcula r[i]^2.
    D = zeros(Float64, number_of_points, number_of_points)

    for i in 1:number_of_points-1
        for j in i+1:number_of_points
            D[i, j] = norm(P[i, :] - P[j, :])
            D[j, i] = D[i, j]
        end
    end

    filename = "BP_problem_$(n)_$(round(P[end, 1], digits=6))_$(round(P[end, 2], digits=6)).txt"


    open(filename, "w") do file
        write(file, "Dimension =\n")
        write(file, "$(n)\n")

        write(file, "Branching levels =\n")
        write(file, "$(m)\n")

        write(file, "Number of points =\n")
        write(file, "$(number_of_points)\n")

        # Os primeiros n pontos são necessários para iniciar o BP.
        write(file, "Initial points =\n")

        for i in 1:n
            write(file, join(P[i, :], " "))
            write(file, "\n")
        end

        # Matriz completa de distâncias.
        write(file, "Distance matrix =\n")

        for i in 1:number_of_points
            write(file, join(D[i, :], " "))
            write(file, "\n")
        end

        # Realização completa usada para gerar as distâncias.
        write(file, "Reference solution =\n")

        for i in 1:number_of_points
            write(file, join(P[i, :], " "))
            write(file, "\n")
        end
    end
    return P, D, filename
end