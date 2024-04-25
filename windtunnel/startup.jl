using MPI

PROCS = [1, 2, 4]
FILE = "wind.jl"
timings = [;]

for (i, v) in enumerate(PROCS)
    println("procs: $v")
    t = @elapsed mpiexec(cmd->run(`$cmd -np $v julia --project=. $FILE`));
    push!(timings, (v, t))
end
display(timings)
