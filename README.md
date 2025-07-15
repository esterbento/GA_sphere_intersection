# Sphere Intersection

This project implements an algorithm to solve the problem of the intersection of $n$ spheres in $\mathbb{R}^n$, where $n$ is a natural number with $n \geq 2$. The problem is modeled using **Conformal Geometric Algebra (CGA)**, and the implementation is written in the **Julia** programming language.
## Installation via Git

To use this project locally via Git:

1. **Clone the repository.**

Open your terminal (Git Bash on Windows or Terminal on macOS/Linux) and run:

```bash
git clone https://github.com/esterbento/GA_sphere_intersection.git
```

2. **Navigate to the project folder:**
```bash
cd GA_sphere_intersection
```

3. **Start the Julia REPL and load the files:**

```julia
include("src/build_problem.jl")
include("src/sphere_intersection.jl")
include("src/qr_factorization.jl")
```

You're all set! You can now generate examples and solve sphere intersection problems.

---


## Examples

The `build_problem` function generates test cases for sphere intersection problems in arbitrary dimensions.

### Syntax

```julia
build_problem(dimension::Int, rank::Bool)
```

- **`dimension`**: the problem's dimension (positive integer).
- **`rank`** (optional):  
  - `true` → generates an example with a **rank deficient** center matrix (rank $n-1$).  
  - `false` (default) → generates an example with a **full rank** center matrix.

Example files are saved in the current working directory.

### Generate an Example

To generate a problem with dimension 3 and a full rank matrix:

```julia
build_problem(3)
```

To generate a problem with dimension 5 and a rank deficient matrix:

```julia
build_problem(5, true)
```

---

## Algorithm Execution

If you generated a file using `build_problem`, it will be saved in the current directory with a name like:

```
problem_3_first_comp_sol_0.959977521807434_0.041049203988499317.txt
```

After that, you can solve the problem using either of the available methods.

To solve using the main approach:

```julia
sphere_intersection("problem_3_first_comp_sol_0.959977521807434_0.041049203988499317.txt")
```
Or, alternatively, using the implemented QR factorization:

```julia
qr_factorization("problem_3_first_comp_sol_0.959977521807434_0.041049203988499317.txt")
```

**Note:** Make sure to replace the filename with the actual name shown in your console or file explorer.

Both methods print the results to the console.

## Running the Test Scripts with the Repository Examples

After loading the main functions (`sphere_intersection` and `qr_factorization`) in the Julia REPL as previously instructed, you can run the scripts that test the existing example files.

**Important:** For the scripts to work correctly, navigate in your terminal to the folder containing the `.txt` files for the dimension and category you want to test. For example:

```julia
cd("examples/fullrank/dim_3")
```

Then, in the Julia REPL, run:

```julia
include("../../../src/compare_results.jl")
include("../../../src/benchmark_tests.jl")
```