using LinearAlgebra
using BenchmarkTools
using Statistics
using Printf

# Lista apenas os arquivos das instâncias do PIE.
files = filter(
    file -> endswith(file, ".txt"), readdir("."))

println("\n=== Testes de erro e tempo da interseção de esferas ===\n")

# Menores erros obtidos em cada instância.
min_errors_AG = Float64[]
min_errors_QR = Float64[]

# Tempos medianos de cada instância, em segundos.
median_times_AG = Float64[]
median_times_QR = Float64[]

# Números de execuções realizadas pelo BenchmarkTools.
num_runs_AG = Int[]
num_runs_QR = Int[]

# ==================================================
# Cálculo dos erros
# ==================================================

for (i, file) in enumerate(files)
    caminho = joinpath(".", file)

    println("Arquivo $(i): $file")

    try
        x1, x2, sol = PIE_CGA(caminho)
        z1, z2 = PIE_QR(caminho)

        println("Solução exata:")
        println(sol)

        println("Solução 1 obtida pelo método AG:")
        println(x1)

        println("Solução 2 obtida pelo método AG:")
        println(x2)

        println("Solução 1 obtida pelo método QR:")
        println(z1)

        println("Solução 2 obtida pelo método QR:")
        println(z2)

        error_AG_1 = norm(sol - x1[1:length(sol)])^2
        error_AG_2 = norm(sol - x2[1:length(sol)])^2
        error_QR_1 = norm(sol - z1)^2
        error_QR_2 = norm(sol - z2)^2

        min_error_AG = min(error_AG_1, error_AG_2)
        min_error_QR = min(error_QR_1, error_QR_2)

        push!(min_errors_AG, min_error_AG)
        push!(min_errors_QR, min_error_QR)
    println("--------------------------------------------------\n")
        @printf("Menor erro do método AG: %.2e\n", min_error_AG)
        @printf("Menor erro do método QR: %.2e\n", min_error_QR)

    
    catch error
        println("Erro ao resolver o arquivo:")
        println(error)
    end

    println("--------------------------------------------------\n")
end

# ==================================================
# Benchmarks
# ==================================================

for (i, file) in enumerate(files)
    caminho = joinpath(".", file)

    println("Arquivo $(i): $file")

    try
        bench_AG = @benchmark PIE_CGA($caminho) seconds=10
        bench_QR = @benchmark PIE_QR($caminho) seconds=10

        println("\nBenchmark do método AG:")
        show(stdout, MIME"text/plain"(), bench_AG)
        println("\n")

        println("Benchmark do método QR:")
        show(stdout, MIME"text/plain"(), bench_QR)
        println("\n")

        # bench.times está em nanossegundos.
        median_time_AG = median(bench_AG.times) / 1e9
        median_time_QR = median(bench_QR.times) / 1e9

        executions_AG = length(bench_AG.times)
        executions_QR = length(bench_QR.times)

        push!(median_times_AG, median_time_AG)
        push!(median_times_QR, median_time_QR)
        push!(num_runs_AG, executions_AG)
        push!(num_runs_QR, executions_QR)

        @printf("Tempo mediano do método AG: %.6e segundos\n", median_time_AG)
        @printf("Tempo mediano do método QR: %.6e segundos\n", median_time_QR)

        println("Número de execuções do método AG: ", executions_AG)
        println("Número de execuções do método QR: ", executions_QR)

    catch error
        println("Erro ao executar o benchmark:")
        println(error)
    end

    println("--------------------------------------------------\n")
end

# ==================================================
# Estatísticas finais dos erros
# ==================================================

if !isempty(min_errors_AG)
    println("\n=== Estatísticas dos erros do método AG ===\n")

    @printf("Menor erro obtido: %.2e\n", minimum(min_errors_AG))
    @printf("Maior erro obtido: %.2e\n", maximum(min_errors_AG))
    @printf("Erro médio: %.2e\n", mean(min_errors_AG))
else
    println("\nNenhum erro foi calculado para o método AG.")
end

if !isempty(min_errors_QR)
    println("\n=== Estatísticas dos erros do método QR ===\n")

    @printf("Menor erro obtido: %.2e\n", minimum(min_errors_QR))
    @printf("Maior erro obtido: %.2e\n", maximum(min_errors_QR))
    @printf("Erro médio: %.2e\n", mean(min_errors_QR))
else
    println("\nNenhum erro foi calculado para o método QR.")
end

# ==================================================
# Estatísticas finais dos tempos
# ==================================================

if !isempty(median_times_AG)
    println("\n=== Estatísticas dos tempos medianos do método AG ===\n")

    @printf("Menor tempo mediano: %.2e segundos\n", minimum(median_times_AG))
    @printf("Maior tempo mediano: %.2e segundos\n", maximum(median_times_AG))
    @printf("Tempo mediano médio: %.2e segundos\n", mean(median_times_AG))
    @printf("Menor número de execuções: %d\n", minimum(num_runs_AG))
end

if !isempty(median_times_QR)
    println("\n=== Estatísticas dos tempos medianos do método QR ===\n")

    @printf("Menor tempo mediano: %.2e segundos\n", minimum(median_times_QR))
    @printf("Maior tempo mediano: %.2e segundos\n", maximum(median_times_QR))
    @printf("Tempo mediano médio: %.2e segundos\n", mean(median_times_QR))
    @printf("Menor número de execuções: %d\n", minimum(num_runs_QR))
end
