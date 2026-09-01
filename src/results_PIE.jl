using LinearAlgebra
using BenchmarkTools
using Statistics
using Printf
using DelimitedFiles

using Random

function read_PIE(filename::String)
    lines = readlines(filename)
    # n: dimension of the space
    n = parse(Int, strip(lines[2]))

    # A: matrix of sphere centers (each row is a center)
    points_block = join(lines[4:3+n], "\n")
    A = readdlm(IOBuffer(points_block), ' ')

    # r: vector of sphere radii
    r_line = replace(replace(lines[5+n], '[' => ""), ']' => "")
    r = parse.(Float64, split(r_line, ","))

    # sol: expected solution vector
    sol_line = replace(replace(lines[7+n], '[' => ""), ']' => "")
    sol = parse.(Float64, split(sol_line, ","))

    return n, A, r, sol 
end

function benchmark_PIE(A, r, n, solver)

    bench = @benchmark $solver(C, rc) setup=(C = copy($A); rc = copy($r)) evals=1 seconds=1

    return bench
end

function quality_test_PIE(n, N, solver; seed = 4114453573, rank::Bool = false)

    Random.seed!(seed)

    A = Matrix{Float64}(undef, n, n)
    r = Vector{Float64}(undef, n)
    s = Vector{Float64}(undef, n)

    min_errors = Vector{Float64}(undef, N)

    for i = 1:N

        build_problem_PIE!(A, r, s, rank)

        flag, x1, x2 = solver(A, r)

        (flag == :nosol) && println("Problema sem interseção.")

        min_errors[i] = min(norm(s - x1)^2, norm(s - x2)^2)

    end

    return min_errors

end

function time_test_PIE(n, N, solver; seed = 4114453573, rank::Bool = false)

    Random.seed!(seed)

    A = Matrix{Float64}(undef, n, n)
    r = Vector{Float64}(undef, n)
    s = Vector{Float64}(undef, n)

    benchmarks = Vector{Any}(undef, N)

    for i = 1:N

        build_problem_PIE!(A, r, s, rank)

        benchmarks[i] = benchmark_PIE(A, r, n, solver)

    end

    return benchmarks

end

# Lista apenas os arquivos das instâncias do PIE.
files = filter(file -> endswith(file, ".txt"), readdir("."))

println("\n=== Testes de erro e tempo da interseção de esferas ===\n")

# Menores erros obtidos em cada instância.
min_errors_AG = Float64[]
#min_errors_QR = Float64[]

# Tempos medianos de cada instância, em segundos.
median_times_AG = Float64[]
#median_times_QR = Float64[]

# Números de execuções realizadas pelo BenchmarkTools.
num_runs_AG = Int[]
#num_runs_QR = Int[]

# ==================================================
# Cálculo dos erros
# ==================================================

for (i, file) in enumerate(files)
    caminho = joinpath(".", file)

    println("Arquivo $(i): $file")

    try
        n, A, r, sol = read_PIE(caminho)

        x1, x2 = PIE_CGA(copy(A), copy(r), n)
        #z1, z2 = PIE_QR(copy(A), copy(r), n)

        println("Solução exata:")
        println(sol)

        println("Solução 1 obtida pelo método AG:")
        println(x1)

        println("Solução 2 obtida pelo método AG:")
        println(x2)

#        println("Solução 1 obtida pelo método QR:")
#        println(z1)

#        println("Solução 2 obtida pelo método QR:")
#        println(z2)

        error_AG_1 = norm(sol - x1[1:length(sol)])^2
        error_AG_2 = norm(sol - x2[1:length(sol)])^2
#        error_QR_1 = norm(sol - z1)^2
#        error_QR_2 = norm(sol - z2)^2

        min_error_AG = min(error_AG_1, error_AG_2)
#        min_error_QR = min(error_QR_1, error_QR_2)

        push!(min_errors_AG, min_error_AG)
#        push!(min_errors_QR, min_error_QR)
    println("--------------------------------------------------\n")
        @printf("Menor erro do método AG: %.2e\n", min_error_AG)
#        @printf("Menor erro do método QR: %.2e\n", min_error_QR)

    
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
        n, A, r, sol = read_PIE(caminho)

        bench_AG, bench_QR = benchmark_PIE(A, r, n)

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

