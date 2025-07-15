using BenchmarkTools
using Statistics
using Printf

# List all files in the current directory ending with ".txt"
files = filter(f -> endswith(f, ".txt"), readdir("."))

println("\n=== Benchmark for multiple files ===\n")

# Vectors to store median times (in seconds)
median_times_AG = Float64[]
median_times_QR = Float64[]

# Vectors to store the number of executions performed
num_runs_AG = Float64[]
num_runs_QR = Float64[]

for (i, file) in enumerate(files)
    println(" File $(i): $file")
    
    bench1 = @benchmark sphere_intersection($file)
    bench2 = @benchmark qr_factorization($file)
    # Print detailed benchmark results
    println("\nBenchmark sphere_intersection:")
    show(stdout, MIME"text/plain"(), bench1)
    println("\n")

    println("\nBenchmark qr_factorization:")
    show(stdout, MIME"text/plain"(), bench2)
    println("\n")

    # Calculate medians in seconds (bench.times is in nanoseconds)
    time_AG = median(bench1.times) / 1e9
    time_QR = median(bench2.times) / 1e9
    
    @printf("Median time sphere_intersection: %.6f seconds\n", time_AG)
    @printf("Median time qr_factorization: %.6f seconds\n", time_QR)
    
    # Store medians and number of executions for final statistics
    push!(median_times_AG, time_AG)
    push!(median_times_QR, time_QR)
    
    push!(num_runs_AG, length(bench1.times))
    push!(num_runs_QR, length(bench2.times))
    
    println("--------------------------------------------------\n")
end

# Final statistics
println(" Median times statistics:")
@printf("Minimum median time using AG: %.2e seconds\n", minimum(median_times_AG))
@printf("Maximum median time using AG: %.2e seconds\n", maximum(median_times_AG))
@printf("Average median time using AG: %.2e seconds\n", mean(median_times_AG))
@printf("Minimum number of runs for AG: %d\n",minimum(num_runs_AG))

@printf("Minimum median time using QR: %.2e seconds\n", minimum(median_times_QR))
@printf("Maximum median time using QR: %.2e seconds\n", maximum(median_times_QR))
@printf("Average median time using QR: %.2e seconds\n", mean(median_times_QR))
@printf("Minimum number of runs for QR: %d\n", minimum(num_runs_QR))
