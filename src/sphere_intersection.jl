using LinearAlgebra
using DelimitedFiles

function sphere_intersection(filename::String)
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
    
    
    B = copy(A) # Copy matrix A to build the augmented system
    u = -ones(n)  # auxiliary vector of -1's

    # Replace each radius by the constant ( (||a_i||^2 - r_i^2) / 2 )
    for i=1:n
        r[i]=0.5*(norm(B[i,:],2)^2 - r[i]^2)
    end
    # Construct the augmented matrix: [A | -1 | r]
    B = hcat(B, u, r)
    m = n+2
    piv_col = 1                 # pivot column index
    coluna_sem_pivo = 0.0       # tracks if a column without a useful pivot is found

    # Perform row-reduced echelon form reduction of matrix B
    for i in 1:min(n, m)
        # Search for a column with a non-zero pivot
        while piv_col <= m && all(abs.(B[i:end, piv_col]) .< eps(Float32))
            coluna_sem_pivo = piv_col
            piv_col += 1
        end
        if piv_col > m
            break
        end

        # Find the row with the largest absolute value for pivot in current column
        pivot_row = argmax(abs.(B[i:end, piv_col])) + (i - 1)

        # Swap the current row with the pivot row
        B[i, :], B[pivot_row, :] = B[pivot_row, :], B[i, :]

        # Normalize the pivot row
        B[i, :] /= B[i, piv_col]

        # Eliminate elements above and below the pivot in the current column
        for j in 1:n
            if j != i
                B[j, :] -= B[j, piv_col] * B[i, :]
            end
        end
        piv_col += 1  # move to the next column
    end
    
    # Initialize variables to store the two solutions
    Δ = 0.0
    x1 = zeros(n+1)
    x2 = zeros(n+1)

    if coluna_sem_pivo == 0.0 
        # Case of unique solution (full rank)

        # Solve the quadratic equation for the scalar parameter
        c = norm(B[1:n, m], 2)^2
        a = norm(B[1:n, m-1], 2)^2
        b = -2 * (sum(B[1:n, m-1] .* B[1:n, m]) + 1)

        sqrt_Δ = sqrt(b^2 - 4*a*c)

        x1[n+1] = (-b + sqrt_Δ) / (2*a)
        x2[n+1] = (-b - sqrt_Δ) / (2*a)

        # Compute the remaining coordinates of x1 and x2
        for i=1:n
            x1[i] = -B[i,m-1]*x1[n+1] + B[i,m]
            x2[i] = -B[i,m-1]*x2[n+1] + B[i,m]
        end
    else
        # Case with a column without a pivot (rank(A) = n - 1)

        # Set the known value in the last coordinate
        x1[n+1] = B[n, m]
        x2[n+1] = x1[n+1]

        # Solve the quadratic equation for the free coordinate
        c = norm(B[1:n-1, m], 2)^2 - 2.0*x1[n+1]
        a = norm(B[coluna_sem_pivo-1:-1:1, coluna_sem_pivo], 2)^2 + 1.0
        b = -2.0 * sum(B[1:coluna_sem_pivo - 1, coluna_sem_pivo] .* B[1:coluna_sem_pivo - 1, m])

        sqrt_Δ = sqrt(b^2 - 4.0*a*c)

        x1[coluna_sem_pivo] = (-b + sqrt_Δ) / (2.0*a)
        x2[coluna_sem_pivo] = (-b - sqrt_Δ) / (2.0*a)

        # Complete vectors x1 and x2 using the free coordinate
        d = 1
        for j in 1:n+1
            if j == coluna_sem_pivo
                continue
            end
            x1[j] = -B[d, coluna_sem_pivo]*x1[coluna_sem_pivo] + B[d, m]
            x2[j] = -B[d, coluna_sem_pivo]*x2[coluna_sem_pivo] + B[d, m]
            d += 1
        end
    end

    return x1, x2, sol # return the two solutions and the original solution for comparison
end
