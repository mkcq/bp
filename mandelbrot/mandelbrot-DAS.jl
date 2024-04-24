using MPI

function mandel(x, y, max_iters)
    z = Complex(x, y)
    c = z
    threshold = 2
    for n in 1:max_iters
        if abs(z) > threshold
            return n - 1
        end
        z = z^2 + c
    end
    max_iters
end

function main(GRID_RES, MAX_ITER)
    MPI.Init()
    comm = MPI.Comm_dup(MPI.COMM_WORLD)
    nranks = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)
    
    ROOT = 0
    xmin = -1.7; xmax = 0.7; ymin = -1.2; ymax = 1.2

    # Divide rows to calculate by number of workers
    rows_w = GRID_RES ÷ nranks + (GRID_RES % nranks > rank)
    lb = 1 + (rank - 1) * rows_w + rows_w
    ub = rank * rows_w + rows_w

    rows = LinRange(ymin, ymax, GRID_RES)
    cols = LinRange(xmin, xmax, GRID_RES)

    t_comp = @elapsed snd = [mandel(c, r, MAX_ITER) for c in cols, r in rows[lb:ub]]
    t_comm = @elapsed rcv = MPI.Gather(snd, comm;root=ROOT)
    MPI.Barrier(comm)
    println("\n$MAX_ITER $rank $size $t_comp $t_comm")
end

#=
job.sh file

#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=16

. /etc/profile.d/lmod.sh
module load openmpi/gcc/64
. /etc/profile.d/lmod.sh
module load julia/1.7.3

pwd

# Parallel

mpiexec -np 8 julia --project=. -e '
    include("mandelbrot.jl")
    GRID_RES = 1000
    MAX_ITER = [100, 1000, 10000, 100000, 1000000]
    for (i, v) in enumerate(MAX_ITER)
    	main(GRID_RES, v)
    end
'
=#