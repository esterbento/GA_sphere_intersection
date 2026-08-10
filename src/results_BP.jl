using LinearAlgebra
using BenchmarkTools
using Statistics
using Printf

# Lista apenas os arquivos das instâncias do BP
files = filter(
    file ->
        endswith(file, ".txt") &&
        startswith(file, "BP_problem_"),
    readdir(".")
)
println("\n=== Testes de erro e tempo do Branch-and-Prune ===\n")

# Menor erro obtido em cada instância
min_errors_BP = Float64[]

# Tempos medianos de cada instância, em segundos
median_times_BP = Float64[]

# Número de execuções realizadas pelo BenchmarkTools
num_runs_BP = Int[]

# Quantidade de soluções encontradas em cada instância
number_of_solutions = Int[]

for (i, file) in enumerate(files)

    caminho = joinpath(".", file)

    println("Arquivo $(i): $file")

    try
        # --------------------------------------------------
        # Cálculo do erro
        # --------------------------------------------------
        n, N, pontos_iniciais, D, Psol = read_problem_BP(caminho)

        solucoes  = solve_problem_BP(n, N, pontos_iniciais, D)

        quantidade = length(solucoes)

        push!(number_of_solutions, quantidade)

        println("Quantidade de soluções encontradas: ", quantidade)

        if isempty(solucoes)
            println("Nenhuma solução foi encontrada.")
        else
            errors = [
                norm(Psol - solution)^2
                for solution in solucoes
            ]

            min_error = minimum(errors)

            push!(min_errors_BP, min_error)

            @printf(
                "Menor erro da instância: %.2e\n",
                min_error
            )
        end
    catch error
        println("Erro ao resolver o arquivo:")
        println(error)
    end

    println("--------------------------------------------------\n")
end

for (i, file) in enumerate(files)

    caminho = joinpath(".", file)

    println("Arquivo $(i): $file")

    try
        # --------------------------------------------------
        # Benchmark
        # --------------------------------------------------

        n, N, pontos_iniciais, D, Psol = read_problem_BP(caminho)

        bench = @benchmark solve_problem_BP($n, $N, $pontos_iniciais, $D) seconds=180

        println("\nBenchmark:")
        show(stdout, MIME"text/plain"(), bench)
        println("\n")

        # bench.times está em nanossegundos
        median_time = median(bench.times) / 1e9
        executions = length(bench.times)

        push!(median_times_BP, median_time)
        push!(num_runs_BP, executions)

        @printf(
            "Tempo mediano: %.6e segundos\n",
            median_time
        )

        println("Número de execuções: ", executions)

    catch error
        println("Erro ao resolver o arquivo:")
        println(error)
    end

    println("--------------------------------------------------\n")
end

# ==================================================
# Estatísticas finais dos erros
# ==================================================

if !isempty(min_errors_BP)

    println("\n=== Estatísticas dos erros ===\n")

    @printf(
        "Menor erro obtido pelo BP: %.2e\n",
        minimum(min_errors_BP)
    )

    @printf(
        "Maior erro obtido pelo BP: %.2e\n",
        maximum(min_errors_BP)
    )

    @printf(
        "Erro médio do BP: %.2e\n",
        mean(min_errors_BP)
    )

else
    println("\nNenhum erro foi calculado.")
end

# ==================================================
# Estatísticas finais dos tempos
# ==================================================

if !isempty(median_times_BP)

    println("\n=== Estatísticas dos tempos medianos ===\n")

    @printf(
        "Menor tempo mediano: %.2e segundos\n",
        minimum(median_times_BP)
    )

    @printf(
        "Maior tempo mediano: %.2e segundos\n",
        maximum(median_times_BP)
    )

    @printf(
        "Tempo mediano médio: %.2e segundos\n",
        mean(median_times_BP)
    )

    @printf(
        "Menor número de execuções: %d\n",
        minimum(num_runs_BP)
    )
end

# ==================================================
# Estatísticas das soluções
# ==================================================

if !isempty(number_of_solutions)

    println("\n=== Estatísticas das soluções ===\n")

    @printf(
        "Menor número de soluções encontradas: %d\n",
        minimum(number_of_solutions)
    )

    @printf(
        "Maior número de soluções encontradas: %d\n",
        maximum(number_of_solutions)
    )

    @printf(
        "Número médio de soluções encontradas: %.2f\n",
        mean(number_of_solutions)
    )
end


