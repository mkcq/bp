# using MPIClusterManagers
# using Distributed
using GLMakie
using LinearAlgebra
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

function plotMandel(r, s)
    m = Matrix{Int64}(undef, GRID_RES, GRID_RES)
    for (i, v) in enumerate(r)
        m[i] = v
    end
    p = heatmap(m, dpi=GRID_RES)
    filename = "mandelbrot-mpi-gres-$(ARGS[1])-maxi-$(ARGS[2])-pr-$s.png"
    save(filename, p)
end

function main()
    MPI.Init()
    comm = MPI.Comm_dup(MPI.COMM_WORLD)
    nranks = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)
    
    ROOT = 0
    GRID_RES = parse(Int, ARGS[1])
    MAX_ITER = parse(Int, ARGS[2])
    xmin = -1.7; xmax = 0.7; ymin = -1.2; ymax = 1.2

   # Divide rows to calculate by number of workers
    rows_w = GRID_RES ÷ nranks + (GRID_RES % nranks > rank)
    lb = 1 + (rank - 1) * rows_w + rows_w
    ub = rank * rows_w + rows_w

    rows = LinRange(ymin, ymax, GRID_RES)
    cols = LinRange(xmin, xmax, GRID_RES)

    snd = [mandel(c, r, MAX_ITER) for c in cols, r in rows[lb:ub]]
    rcv = MPI.Gather(snd, comm;root=ROOT)
    if rank == ROOT
        # temp = rcv
        # plotMandel(temp, nranks)
    end 
end

#=
Citations/References:

Simon Byrne, Lucas C. Wilcox, and Valentin Churavy (2021) "MPI.jl: Julia bindings for the Message Passing Interface". JuliaCon Proceedings, 1(1), 68, doi: 10.21105/jcon.00068

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
