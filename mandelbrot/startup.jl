using MPI
using BenchmarkTools
using BenchmarkPlots
using StatsPlots

PROCS = [1, 2, 4, 8]
FILE = "mandelbrot.jl"

suite = BenchmarkGroup()
suite["gres-1K"] = BenchmarkGroup()
suite["gres-1K"]["proc-1"] = @benchmarkable mpiexec(cmd->run(`$cmd -np $(PROCS[1]) julia --project=. $FILE 1000 100`));
suite["gres-1K"]["proc-2"] = @benchmarkable mpiexec(cmd->run(`$cmd -np $(PROCS[2]) julia --project=. $FILE 1000 100`));
suite["gres-1K"]["proc-4"] = @benchmarkable mpiexec(cmd->run(`$cmd -np $(PROCS[3]) julia --project=. $FILE 1000 100`));
suite["gres-1K"]["proc-8"] = @benchmarkable mpiexec(cmd->run(`$cmd -np $(PROCS[4]) julia --project=. $FILE 1000 100`));

# suite["gres-10K"] = BenchmarkGroup()
# suite["gres-10K"]["proc-1"] = @benchmarkable mpiexec(cmd->run(`$cmd -np $(PROCS[1]) julia --project=. $FILE 10000 100`));
# suite["gres-10K"]["proc-2"] = @benchmarkable mpiexec(cmd->run(`$cmd -np $(PROCS[2]) julia --project=. $FILE 10000 100`));
# suite["gres-10K"]["proc-4"] = @benchmarkable mpiexec(cmd->run(`$cmd -np $(PROCS[3]) julia --project=. $FILE 10000 100`));
# suite["gres-10K"]["proc-8"] = @benchmarkable mpiexec(cmd->run(`$cmd -np $(PROCS[4]) julia --project=. $FILE 10000 100`));

suite["gres-10K"] = BenchmarkGroup()
suite["gres-10K"]["proc-1"] = @benchmarkable mpiexec(cmd->run(`$cmd -np $(PROCS[1]) julia --project=. $FILE 10000 1000`));
suite["gres-10K"]["proc-2"] = @benchmarkable mpiexec(cmd->run(`$cmd -np $(PROCS[2]) julia --project=. $FILE 10000 1000`));
suite["gres-10K"]["proc-4"] = @benchmarkable mpiexec(cmd->run(`$cmd -np $(PROCS[3]) julia --project=. $FILE 10000 1000`));
suite["gres-10K"]["proc-8"] = @benchmarkable mpiexec(cmd->run(`$cmd -np $(PROCS[4]) julia --project=. $FILE 10000 1000`));

results = run(suite, verbose = true)
m = median(results)
println(m)

# println("suite.data ", suite.data)
# println("suite.tags ", suite.tags)

# plot(results)
# savefig("benchmarks.png")
BenchmarkTools.save("benchmarks.json", results)