# Sphere Intersection

This project implements an algorithm to solve the problem of the intersection of $n$ spheres in $\mathbb{R}^n$, with $n \geq 2$. The problem is modeled using **Conformal Geometric Algebra (CGA)**, and the implementation is written in the **Julia** programming language.

The repository also includes a method based on QR decomposition for comparison and an application of the proposed method to the Branch-and-Prune algorithm for the $K$-DMDGP.

The examples are organized by dimension in the `fullrank`, `rankdeficient`, and `tests_BP` folders, which contain subfolders in the format `dim_n`.

## Installation via Git

To use the project locally:

1. **Clone the repository:**

```bash
git clone https://github.com/esterbento/GA_sphere_intersection.git
```

2. **Navigate to the project folder:**

```bash
cd GA_sphere_intersection
```

3. **Install BenchmarkTools in the Julia REPL:**

```julia
using Pkg
Pkg.add("BenchmarkTools")
```

4. **Load the project files:**

```julia
include("src/build_problem_PIE.jl")
include("src/PIE_CGA.jl")
include("src/PIE_QR.jl")
include("src/build_problem_BP.jl")
include("src/branch_and_prune.jl")
```

You can now generate examples and run the sphere intersection and Branch-and-Prune algorithms.

---

## Sphere Intersection Examples

The `build_problem` function generates sphere intersection problem instances in arbitrary dimensions.

### Syntax

```julia
build_problem(dimension::Int, rank::Bool)
```

- **`dimension`**: dimension of the problem;
- **`rank`** (optional):
  - `true` generates a rank-deficient instance;
  - `false` (default) generates a full-rank instance.

To generate a full-rank problem in dimension 3:

```julia
build_problem(3)
```

To generate a rank-deficient problem in dimension 5:

```julia
build_problem(5, true)
```

The generated file is saved in the current working directory.

## Sphere Intersection Execution

To solve an instance using the CGA-based method:

```julia
PIE_CGA("file_name.txt")
```

To solve the same instance using the QR decomposition-based method:

```julia
PIE_QR("file_name.txt")
```

Replace `file_name.txt` with the name of the file you want to solve.

## Running the Sphere Intersection Tests

After loading `PIE_CGA.jl` and `PIE_QR.jl`, navigate to the folder containing the `.txt` files for the desired dimension and category. For example:

```julia
cd("examples/fullrank/dim_3")
```

Then run:

```julia
include("../../../src/results_PIE.jl")
```

To test the rank-deficient instances, navigate to the corresponding folder. For example:

```julia
cd("examples/rankdeficient/dim_3")
```

The `results_PIE.jl` file computes the errors and execution times of the CGA and QR methods.

---

## Branch-and-Prune

The Branch-and-Prune implementation uses the sphere intersection method to solve instances of the $K$-DMDGP.

To generate a Branch-and-Prune instance:

```julia
build_problem_BP(n)
```

To solve an existing instance:

```julia
solve_problem_BP("BP_file_name.txt")
```

## Running the Branch-and-Prune Tests

After loading `PIE_CGA.jl` and `branch_and_prune.jl`, navigate to the folder containing the instances for the desired dimension. For example:

```julia
cd("examples/tests_BP/dim_3")
```

Then run:

```julia
include("../../../src/results_BP.jl")
```

The `results_BP.jl` file computes the errors, execution times, and number of solutions found for all instances in the selected folder.
