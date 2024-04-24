using MPI
using BenchmarkTools
using BenchmarkPlots
using StatsPlots

PROCS = [1, 2, 4, 8]
FILE = "mandelbrot.jl"
timings = [;]
GRID_RES = 1_000 # Higher value improves resolution
MAX_ITER = 1_000_000 # Higher value improves accuracy

for (i, v) in enumerate(PROCS)
    println("procs: $v")
    t = @elapsed mpiexec(cmd->run(`$cmd -np $v julia --project=. $FILE $GRID_RES $MAX_ITER`));
    push!(timings, (v, t))
end
display(timings)
