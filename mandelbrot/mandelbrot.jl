#=
Exercise 3 from: https://www.francescverdugo.com/XM_40017/dev/julia_basics/

Function mandel estimates if a given point (x,y) in the complex plane belongs to the Mandelbrot set.

If the value of mandel is less than max_iters, the point is provably outside the Mandelbrot set.
If mandel is equal to max_iters, then the point is provably inside the set.
The larger max_iters, the better the quality of the estimate (the nicer will be your plot).

Plot the value of function mandel for each pixel in a 2D grid of the box.

(−1.7, 0.7) × (−1.2, 1.2)

Use a grid resolution of at least 1000 points in each direction and max_iters at least 10.
You can increase these values to get nicer plots.
To plot the values use function heatmap from the Julia package GLMakie.
Use LinRange to divide the horizontal and vertical axes into pixels.
See the documentation of these functions for help.
GLMakie is a GPU-accelerated plotting back-end for Julia.
It is a large package and it can take some time to install and to generate the first plot.
Be patient.
=#

using LinearAlgebra
using MPIClusterManagers
using Distributed

NUM_RANKS = 4

# Create workers that can call MPI functions
if procs() == workers()
    manager = MPIWorkerManager(NUM_RANKS)
    addprocs(manager)
end

@everywhere workers() begin
    using MPI
    using GLMakie

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

    MPI.Init()
    comm = MPI.Comm_dup(MPI.COMM_WORLD)
    nranks = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    # println("Hello, I am process $rank of $nranks processes!")
    root = 0

    GRID_RESOLUTION = 1000
    MAX_ITERS = 100
    xmin = -1.7
    xmax = 0.7
    ymin = -1.2
    ymax = 1.2
    
    # Divide rows to calculate by number of workers
    nRowsR = GRID_RESOLUTION ÷ nranks
    lb = 1 + (myid() - 2) * nRowsR
    ub = (myid() - 1) * nRowsR
    # println("Rank $rank will compute rows $lb until $ub")

    cols = LinRange(ymin, ymax, GRID_RESOLUTION)
    xAxis = LinRange(xmin, xmax, GRID_RESOLUTION)
    rows = LinRange(xAxis[lb], xAxis[ub], nRowsR)
    snd = [mandel(r, c, MAX_ITERS) for r in rows, c in cols]
    
    # @show typeof(snd), size(snd)
    
    # rcv = MPI.Gather(snd, comm;root)

    # MPI.Barrier(comm)

    # type of rcv should be a matrix Matrix{Int64} of (1000, 1000)

    if rank == root
        global res = snd
        # @show nranks
        rcv = Ref(0)
        @show typeof(res), size(res)
        
        MPI.Recv!(rcv, comm;source=1, tag=0)
        res = vcat(res, rcv)
        @show typeof(res), size(res)

        # for i in 1:(nranks-1)
        #     res = vcat(res, MPI.Recv!(rcv, comm;source=i, tag=0)) 
        # end
        @show typeof(res), size(res)

        # @show typeof(rcv)
        # @show size(rcv)
        # p = heatmap(rcv)
        # save("mandelbrot2.png", p)
    else
        MPI.Send(snd, comm;dest=root, tag=0)
    end
end
