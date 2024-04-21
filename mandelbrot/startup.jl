using MPI
using BenchmarkTools
using BenchmarkPlots
using StatsPlots

PROCS = [1, 2, 4]
FILE = "mandelbrot.jl"
# @time mpiexec(cmd->run(`$cmd -np $PROCS julia --project=. $FILE`));
# b = @benchmark mpiexec(cmd->run(`$cmd -np $(PROCS[3]) julia --project=. $FILE`));

suite = BenchmarkGroup()
suite["proc-1"] = @benchmarkable mpiexec(cmd->run(`$cmd -np $(PROCS[1]) julia --project=. $FILE`));
suite["proc-2"] = @benchmarkable mpiexec(cmd->run(`$cmd -np $(PROCS[2]) julia --project=. $FILE`));
suite["proc-4"] = @benchmarkable mpiexec(cmd->run(`$cmd -np $(PROCS[3]) julia --project=. $FILE`));
results = run(suite, verbose = true)
plot(results)
savefig("benchmarks.png")
BenchmarkTools.save("benchmarks.json", results)