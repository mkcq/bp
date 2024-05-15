using Test
using MPI
using MPIClusterManagers
using Distributed
# include("wind.jl")
include("wind2.jl")

@testset "tunnelStartRow" begin
    tunnel_rows = 10
    # Cartesian = [c, r]
    cart_dims = [2, 1]; cart_coord = [[0, 0], [1, 0]]
    @test tunnelStartRow(tunnel_rows, cart_dims, cart_coord[1]) == 1
    @test tunnelStartRow(tunnel_rows, cart_dims, cart_coord[2]) == 1

    cart_dims = [2, 2]; cart_coord = [[0, 0], [0, 1], [1, 0], [1, 1]]
    @test tunnelStartRow(tunnel_rows, cart_dims, cart_coord[1]) == 6
    @test tunnelStartRow(tunnel_rows, cart_dims, cart_coord[2]) == 1
    @test tunnelStartRow(tunnel_rows, cart_dims, cart_coord[3]) == 6
    @test tunnelStartRow(tunnel_rows, cart_dims, cart_coord[4]) == 1
end

@testset "tunnelStartColumn" begin
    tunnel_cols = 10

    cart_dims = [2, 1]; cart_coord = [[0, 0], [1, 0]]
    @test tunnelStartColumn(tunnel_cols, cart_dims, cart_coord[1]) == 1
    @test tunnelStartColumn(tunnel_cols, cart_dims, cart_coord[2]) == 6

    cart_dims = [2, 2]; cart_coord = [[0, 0], [0, 1], [1, 0], [1, 1]]
    @test tunnelStartColumn(tunnel_cols, cart_dims, cart_coord[1]) == 1
    @test tunnelStartColumn(tunnel_cols, cart_dims, cart_coord[2]) == 1
    @test tunnelStartColumn(tunnel_cols, cart_dims, cart_coord[3]) == 6
    @test tunnelStartColumn(tunnel_cols, cart_dims, cart_coord[4]) == 6
end

@testset "inCartGrid" begin
    tunnel_rows = 10; tunnel_cols = 10
    cart_dims = [2, 2]; cart_coord = [[0, 0], [0, 1], [1, 0], [1, 1]]
    own_rows = tunnel_rows ÷ cart_dims[2]; own_cols = tunnel_cols ÷ cart_dims[1]
    
    tsr = tunnelStartRow(tunnel_rows, cart_dims, cart_coord[1])
    tsc = tunnelStartColumn(tunnel_cols, cart_dims, cart_coord[1])
    @test inCartGrid(tsr, tsc, 6, 5, own_rows, own_cols) == true
    @test inCartGrid(tsr, tsc, 5, 5, own_rows, own_cols) == false
    @test inCartGrid(tsr, tsc, 6, 6, own_rows, own_cols) == false

    tsr = tunnelStartRow(tunnel_rows, cart_dims, cart_coord[2])
    tsc = tunnelStartColumn(tunnel_cols, cart_dims, cart_coord[2])
    @test inCartGrid(tsr, tsc, 5, 5, own_rows, own_cols) == true
    @test inCartGrid(tsr, tsc, 6, 5, own_rows, own_cols) == false
    @test inCartGrid(tsr, tsc, 5, 6, own_rows, own_cols) == false

    tsr = tunnelStartRow(tunnel_rows, cart_dims, cart_coord[3])
    tsc = tunnelStartColumn(tunnel_cols, cart_dims, cart_coord[3])
    @test inCartGrid(tsr, tsc, 6, 6, own_rows, own_cols) == true
    @test inCartGrid(tsr, tsc, 5, 6, own_rows, own_cols) == false
    @test inCartGrid(tsr, tsc, 6, 5, own_rows, own_cols) == false

    tsr = tunnelStartRow(tunnel_rows, cart_dims, cart_coord[4])
    tsc = tunnelStartColumn(tunnel_cols, cart_dims, cart_coord[4])
    @test inCartGrid(tsr, tsc, 5, 6, own_rows, own_cols) == true
    @test inCartGrid(tsr, tsc, 5, 5, own_rows, own_cols) == false
    @test inCartGrid(tsr, tsc, 6, 6, own_rows, own_cols) == false
end

@testset "updateFan 2c x 2r" begin
    maxiter = 10; fp = 2; fs = 7
    cart_dims = [2, 2]; cart_coord = [[0, 0], [0, 1], [1, 0], [1, 1]]
    tunnel_cols = 10; tunnel_rows = 6
    own_cols = tunnel_cols ÷ cart_dims[1]; own_rows = tunnel_rows ÷ cart_dims[2]
    
    flow0 = zeros(Int, own_rows + 2, own_cols + 2)
    flow1 = zeros(Int, own_rows + 2, own_cols + 2)
    flow2 = zeros(Int, own_rows + 2, own_cols + 2)
    flow3 = zeros(Int, own_rows + 2, own_cols + 2)

    for iter in 1:maxiter
        updateFan(iter, fp + 1, fs, flow0, cart_dims, cart_coord[1], tunnel_cols, own_cols)
        updateFan(iter, fp + 1, fs, flow1, cart_dims, cart_coord[2], tunnel_cols, own_cols)
        updateFan(iter, fp + 1, fs, flow2, cart_dims, cart_coord[3], tunnel_cols, own_cols)
        updateFan(iter, fp + 1, fs, flow3, cart_dims, cart_coord[4], tunnel_cols, own_cols)
    end

    lb = view(flow0, 2:own_rows + 1, 2:own_cols + 1)
    rb = view(flow2, 2:own_rows + 1, 2:own_cols + 1)
    lt = view(flow1, 2:own_rows + 1, 2:own_cols + 1)
    rt = view(flow3, 2:own_rows + 1, 2:own_cols + 1)
    b = hcat(lb, rb)
    t = hcat(lt, rt)
    res = vcat(t, b)
    # display(res)
    
    for j in 1:own_cols
        (fp <= j && j < fp + fs) ? (@test res[1, j] != 0) : @test res[1, j] == 0
    end
end

@testset "updateFan 2c x 1r" begin
    maxiter = 1; fp = 3; fs = 5
    cart_dims = [2, 1]; cart_coord = [[0, 0], [1, 0]]
    tunnel_cols = 10; tunnel_rows = 6
    own_cols = tunnel_cols ÷ cart_dims[1]; own_rows = tunnel_rows ÷ cart_dims[2]

    flow0 = zeros(Int, own_rows + 2, own_cols + 2)
    flow1 = zeros(Int, own_rows + 2, own_cols + 2)

    for iter in 1:maxiter
        updateFan(iter, fp + 1, fs, flow0, cart_dims, cart_coord[1], tunnel_cols, own_cols)
        updateFan(iter, fp + 1, fs, flow1, cart_dims, cart_coord[2], tunnel_cols, own_cols)
    end

    l = view(flow0, 2:own_rows + 1, 2:own_cols + 1)
    r = view(flow1, 2:own_rows + 1, 2:own_cols + 1)
    res = hcat(l, r)
    # display(res)

    for j in 1:own_cols
        (fp <= j && j < fp + fs) ? (@test res[1, j] != 0) : @test res[1, j] == 0
    end
end

append!(ARGS, ["tr", "tc", "mi", "th", "fp", "fs", "fbp", "fbs", "fd", "mbp", "mbs", "md", "rs1", "rs2", "rs3", "5", "5", "0.5", "5", "6", "0.5", "6", "5", "0.5", "6", "6","0.5"])

@testset "readFixedParticles" begin
    tunnel_rows = 10; tunnel_cols = 10
    cart_dims = [2, 2]; cart_coord = [[0, 0], [0, 1], [1, 0], [1, 1]]
    own_cols = tunnel_cols ÷ cart_dims[1]; own_rows = tunnel_rows ÷ cart_dims[2]

    lb = 1; ub = (length(ARGS) - 15) ÷ 3
    particles0::Vector{Particle} = []
    particles1::Vector{Particle} = []
    particles2::Vector{Particle} = []
    particles3::Vector{Particle} = []
    
    for i in 1:4    
        tsr = tunnelStartRow(tunnel_rows, cart_dims, cart_coord[i])
        tsc = tunnelStartColumn(tunnel_cols, cart_dims, cart_coord[i])
        if i == 1
            readFixedParticles(lb, ub, tsr, tsc, own_rows, own_cols, particles0)
        elseif i == 2
            readFixedParticles(lb, ub, tsr, tsc, own_rows, own_cols, particles1)
        elseif i == 3
            readFixedParticles(lb, ub, tsr, tsc, own_rows, own_cols, particles2)
        else
            readFixedParticles(lb, ub, tsr, tsc, own_rows, own_cols, particles3)
        end
    end

    @test particles0[1].position == [6, 5] * PRECISION
    @test particles1[1].position == [5, 5] * PRECISION
    @test particles2[1].position == [6, 6] * PRECISION
    @test particles3[1].position == [5, 6] * PRECISION
end

# # TODO: Make this test cleaner.
@testset "getNeighborRanks 2c x 2r" begin
    mpiexec(cmd->run(`$cmd -np 4 julia --project=. -e '
    using Test
    include("wind2.jl")
    tunnel_rows = 16; tunnel_cols = 16

    MPI.Init()
    comm = MPI.Comm_dup(MPI.COMM_WORLD)
    n_ranks = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)

    cart_dims = getCartDims(n_ranks)
    periodic = map(_ -> false, cart_dims)
    reorder = false
    cart_comm = MPI.Cart_create(comm, cart_dims; periodic, reorder)
    cart_coord = MPI.Cart_coords(cart_comm)

    n = getNeighborRanks(cart_dims, cart_coord, cart_comm)
    
    if rank == 0
        @test n["N"] == 1
        @test n["S"] == -1
        @test n["E"] == 2
        @test n["W"] == -1
        @test n["NE"] == 3
        @test n["NW"] == -1
        @test n["SE"] == -1
        @test n["SW"] == -1
    end; MPI.Barrier(cart_comm)
    if rank == 1
        @test n["N"] == -1
        @test n["S"] == 0
        @test n["E"] == 3
        @test n["W"] == -1
        @test n["NE"] == -1
        @test n["NW"] == -1
        @test n["SE"] == 2
        @test n["SW"] == -1
    end; MPI.Barrier(cart_comm)
    if rank == 2
        @test n["N"] == 3
        @test n["S"] == -1
        @test n["E"] == -1
        @test n["W"] == 0
        @test n["NE"] == -1
        @test n["NW"] == 1
        @test n["SE"] == -1
        @test n["SW"] == -1
    end; MPI.Barrier(cart_comm)
    if rank == 3
        @test n["N"] == -1
        @test n["S"] == 2
        @test n["E"] == -1
        @test n["W"] == 1
        @test n["NE"] == -1
        @test n["NW"] == -1
        @test n["SE"] == -1
        @test n["SW"] == 0
    end; MPI.Barrier(cart_comm)
    '`))
end

# @testset "getNeighborRanks 2c x 2r 2" begin
#     if procs() == workers()
#         n_ranks = 4
#         manager = MPIWorkerManager(n_ranks)
#         addprocs(manager)
#     end

#     @everywhere workers() begin
#         using MPI
#         include("wind2.jl")

#         tunnel_rows = 16; tunnel_cols = 16

#         MPI.Init()
#         comm = MPI.Comm_dup(MPI.COMM_WORLD)
#         n_ranks = MPI.Comm_size(comm)
#         rank = MPI.Comm_rank(comm)

#         cart_dims = getCartDims(n_ranks)
#         periodic = map(_ -> false, cart_dims)
#         reorder = false
#         cart_comm = MPI.Cart_create(comm, cart_dims; periodic, reorder)
#         cart_coord = MPI.Cart_coords(cart_comm)

#         n = getNeighborRanks(cart_dims, cart_coord, cart_comm)

#         if rank == 0
#             @show n["N"] == 1
#             @show n["S"] == -1
#             @show n["E"] == 2
#             @show n["W"] == -1
#         end; MPI.Barrier(cart_comm)
#         if rank == 1
#             @show n["N"] == -1
#             @show n["S"] == 0
#             @show n["E"] == 3
#             @show n["W"] == -1
#         end; MPI.Barrier(cart_comm)
#         if rank == 2
#             @show n["N"] == 3
#             @show n["S"] == -1
#             @show n["E"] == -1
#             @show n["W"] == 0
#         end; MPI.Barrier(cart_comm)
#         if rank == 3
#             @show n["N"] == -1
#             @show n["S"] == 2
#             @show n["E"] == -1
#             @show n["W"] == 1
#         end; MPI.Barrier(cart_comm)
#     end
# end

@testset "moveParticle" begin
    # TODO: With MPI.
    #=
    8 x 8
    -1 -1 -1 -1 -1 -1  -1 -1 -1 -1 -1 -1 
    -1  0  0  0  0 -1  -1  0  0  0  0 -1
    -1  0  0  0  0 -1  -1  0  0  0  0 -1
    -1  0  0  0  0 -1  -1  0  0  0  0 -1
    -1  0  0  0  0 -1  -1  0  0  0  0 -1
    -1 -1 -1 -1 -1 -1  -1 -1 -1 -1 -1 -1 

    -1 -1 -1 -1 -1 -1  -1 -1 -1 -1 -1 -1 
    -1  0  0  0  0 -1  -1  0  0  0  0 -1
    -1  0  0  0  0 -1  -1  0  0  0  0 -1
    -1  0  0  0  0 -1  -1  0  0  0  0 -1
    -1  0  0  0  0 -1  -1  0  0  0  0 -1
    -1 -1 -1 -1 -1 -1  -1 -1 -1 -1 -1 -1 
    =#
end