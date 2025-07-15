using LinearAlgebra
using DelimitedFiles

function qr_factorization(filename::String)
    lines = readlines(filename)

    # n: dimension of the space
    n = parse(Int, strip(lines[2]))

    # P: matrix of points (each column is a point)
    points_block = join(lines[4:3+n], "\n")
    P = readdlm(IOBuffer(points_block), ' ')

    # r: vector of radii
    r_line = replace(replace(lines[5+n], '[' => ""), ']' => "")
    r = parse.(Float64, split(r_line, ","))

     

    # Construct matrix A = [a_1 - a_n, ..., a_(n-1) - a_n]
    A = zeros(n,n-1)
    for i=1:n-1
        A[:,i] = P'[:,i]-P'[:,n]
    end

    # Compute the QR decomposition of matrix A
    F=qr(A)

    # Determine the matrix Rp as in equation (2.12)
    Rp = F.R[1:n-1, :]

    # Compute vector c as in equation (2.11)
    c = zeros(n-1)
    for i=1:n-1
        c[i] = -0.5*(r[i]^2-r[n]^2-norm(A[:,i],2)^2)
    end

    # Solve the linear system (2.13)
    y = (Rp')\c

    # Check the intersection conditions
    if r[n]^2- norm(y,2)^2 < 0.0
        println("The intersection is empty.")
        return nothing
    elseif abs(r[n]^2 - norm(y,2)^2) < 1e-10
        z = 0
    else # If radius squared is larger than norm(y)^2, compute z
        z = sqrt(r[n]^2 - norm(y,2)^2)
    end
    
    # Calculate the solution points
    x1 = F.Q * [y; z] + P'[:, n]  # Compute first intersection point
    x2 = F.Q * [y; -z] + P'[:, n] # Compute second intersection point

    # Return one or two solutions depending on whether they coincide
    if norm(x1-x2)<1.0e-10
        return x1
    else
        return x1,x2
    end
end


# Note:
# This algorithm is adapted from Maioli's paper.
# The equation numbers (e.g., (2.11), (2.12), (2.13)) follow his notation.