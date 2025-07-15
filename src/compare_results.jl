using Statistics    # For statistical functions like mean()
using Printf        # For formatted output with @printf

# List all files with ".txt" extension in the current directory
files = filter(f -> endswith(f, ".txt"), readdir("."))

println("\n=== Running test comparisons ===\n")

# Vectors to store the smallest distances between solutions
min_distances_AG = Float64[]   # For results from the AG method
min_distances_QR = Float64[]   # For results from the QR method

# Loop over each test file
for (i, file) in enumerate(files)
    println("File $(i): $file")
    
    # Call sphere_intersection function which returns two solutions (x1, x2) and the exact solution sol
    x1, x2, sol = sphere_intersection(file)
    
    # Call intersecao_esferas function which returns two solutions (z1, z2) for the same problem
    z1, z2 = qr_factorization(file)
    
    # Display obtained solutions and the exact solution
    println("Exact solution:")
    println(sol)
    println("Obtained result 1 (AG):")
    println(x1)
    println("Obtained result 2 (AG):")
    println(x2)
    println("Obtained result 1 (QR):")
    println(z1)
    println("Obtained result 2 (QR):")
    println(z2)
    # Calculate squared distances between the exact solution and each obtained solution
    dist1 = norm(sol - x1[1:length(sol)])^2
    dist2 = norm(sol - x2[1:length(sol)])^2
    dist3 = norm(sol - z1)^2
    dist4 = norm(sol - z2)^2

    # Store the smallest distance for each method
    push!(min_distances_AG, min(dist1, dist2))
    push!(min_distances_QR, min(dist3, dist4))
    
    println("--------------------------------------------------\n")
end

# Print summary statistics of the minimum distances for each method
@printf("Minimum distance using AG: %.2e\n", minimum(min_distances_AG))
@printf("Maximum distance using AG: %.2e\n", maximum(min_distances_AG))
@printf("Mean of minimum distances in AG: %.2e\n", mean(min_distances_AG))

@printf("Minimum distance using QR: %.2e\n", minimum(min_distances_QR))
@printf("Maximum distance using QR: %.2e\n", maximum(min_distances_QR))
@printf("Mean of minimum distances in QR: %.2e\n", mean(min_distances_QR))
