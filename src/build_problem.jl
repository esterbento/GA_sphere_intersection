using LinearAlgebra

function build_problem(n::Int,rank::Bool=false)
    P = Vector{Vector{Float64}}(undef,n)
    r = Vector{Float64}(undef,n)
    Psol = rand(n)
    a = rand(1)

    file = open("problem_$(n)_first_comp_sol_$(Psol[1])_$(Psol[2]).txt","w")
    write(file,"Dimension = \n$(n) \n")
    write(file,"Points = \n")

    if rank == false
        for i in 1:n
            P[i] = rand(n)
            write(file, join(P[i], " "))
            write(file, i < n ? "\n" : "\n")  
            r[i] = norm(P[i] - Psol)
        end
    else
        for i in 1:n
            P[i] = rand(n)
        end
        P[3] = a[1]*P[1]
        for i in 1:n
            write(file, join(P[i], " "))
            write(file, i < n ? "\n" : "\n")  
            r[i] = norm(P[i] - Psol)  
        end      
    end

    write(file, "Radious = \n")
    write(file, "$(r) \n")
    write(file, "Solution = \n")
    write(file, "$(Psol) \n")
    if rank == true
        write(file, "λ = $(a)\n")
    end
    close(file)
end