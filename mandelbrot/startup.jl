using MPI
PROCS = 2
FILE = "mandelbrot.jl"
mpiexec(cmd->run(`$cmd -np $PROCS julia --project=. $FILE`));

# using BenchmarkTools
# @btime mpiexec(cmd->run(`$cmd -np $PROCS julia --project=. $FILE`));