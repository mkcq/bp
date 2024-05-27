using Test
using MPI
using MPIClusterManagers
using Distributed
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
    maxiter = 8; fp = 2; fs = 7
    cart_dims = [2, 2]; cart_coord = [[0, 0], [0, 1], [1, 0], [1, 1]]
    tunnel_cols = 10; tunnel_rows = 4
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
    maxiter = 9; fp = 3; fs = 5
    cart_dims = [2, 1]; cart_coord = [[0, 0], [1, 0]]
    tunnel_cols = 10; tunnel_rows = 6
    own_cols = tunnel_cols ÷ cart_dims[1]; own_rows = tunnel_rows ÷ cart_dims[2]

    flow0 = zeros(Int, own_rows + 2, own_cols + 2)
    flow1 = zeros(Int, own_rows + 2, own_cols + 2)

    for iter in 1:maxiter
        updateFan(iter, fp + 1, fs, flow0, cart_dims, cart_coord[1], tunnel_cols, own_cols)
        updateFan(iter, fp + 1, fs, flow1, cart_dims, cart_coord[2], tunnel_cols, own_cols)
        # ll = view(flow0, 2:own_rows + 1, 2:own_cols + 1)
        # rr = view(flow1, 2:own_rows + 1, 2:own_cols + 1)
        # ress = hcat(ll, rr)
        # display(ress)
    end

    l = view(flow0, 2:own_rows + 1, 2:own_cols + 1)
    r = view(flow1, 2:own_rows + 1, 2:own_cols + 1)
    res = hcat(l, r)
    # println("\n###\nAfter $maxiter iterations.\n###\n")
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

@testset "cleanParticleLocations" begin 
    mpiexec(cmd->run(`$cmd -np 4 julia --project=. -e '
    using Test
    include("wind2.jl")

    tunnel_rows = 10
    tunnel_cols = 10
    
    MPI.Init()
    comm = MPI.Comm_dup(MPI.COMM_WORLD)
    n_ranks = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    cart_dims = getCartDims(n_ranks)
    periodic = map(_ -> false, cart_dims)
    reorder = false
    cart_comm = MPI.Cart_create(comm, cart_dims; periodic, reorder)
    cart_coord = MPI.Cart_coords(cart_comm)
    cart_neighbors = getNeighborRanks(cart_dims, cart_coord, cart_comm)
    
    own_cols = tunnel_cols ÷ cart_dims[1]
    own_rows = tunnel_rows ÷ cart_dims[2]

    p_locs = zeros(Int, own_rows + 2, own_cols + 2); p_locs[1, :] = p_locs[end, :] .= -1; p_locs[:, 1] = p_locs[:, end] .= -1

    p_locs[2:own_rows + 1, 2:own_cols + 1] .= 5

    for iter in 1:80
        cleanParticleLocations(p_locs, own_rows, own_cols, iter, cart_dims, cart_coord)
        ub = min(iter - own_rows * (cart_dims[2] - 1 - cart_coord[2]), own_rows) + 1
        if rank == 1
            for i in 2:own_rows + 1, j in 2:own_cols + 1
                i <= ub ? (@test p_locs[i, j] == 0) : (@test p_locs[i, j] != 0)
            end
        elseif rank == 0
            for i in 2:own_rows + 1, j in 2:own_cols + 1
                i <= ub ? (@test p_locs[i, j] == 0) : (@test p_locs[i, j] != 0)
            end
        end
        MPI.Barrier(cart_comm)
    end
    '`))
end

@testset "sendFlowToNeighbors and recvFlowFromNeighbors" begin
    mpiexec(cmd->run(`$cmd -np 4 julia --project=. -e '
    using Test
    include("wind2.jl")

    tunnel_rows = 10
    tunnel_cols = 10
    
    MPI.Init()
    comm = MPI.Comm_dup(MPI.COMM_WORLD)
    n_ranks = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    cart_dims = getCartDims(n_ranks)
    periodic = map(_ -> false, cart_dims)
    reorder = false
    cart_comm = MPI.Cart_create(comm, cart_dims; periodic, reorder)
    cart_coord = MPI.Cart_coords(cart_comm)
    cart_neighbors = getNeighborRanks(cart_dims, cart_coord, cart_comm)
    
    own_cols = tunnel_cols ÷ cart_dims[1]
    own_rows = tunnel_rows ÷ cart_dims[2]

    flow = zeros(Int, own_rows + 2, own_cols + 2)
    flow[1, :] = flow[end, :] .= -1
    flow[:, 1] = flow[:, end] .= -1
    flow[2:own_rows + 1, 2:own_cols + 1] .= rank

    sendFlowToNeighbors(flow, cart_neighbors, cart_comm)
    recvFlowFromNeighbors(flow, cart_neighbors, cart_comm)
    
    if rank == 0
        for c in 2:own_cols + 1
            @test flow[1, c] == 1
        end

        @test flow[1, own_cols + 2] == 3
        
        for r in 2:own_rows + 1
            @test flow[r, own_cols + 2] == 2
        end
    elseif rank == 1
        for r in 2:own_rows + 1
            @test flow[r, own_cols + 2] == 3
        end
    elseif rank == 2
        for c in 2:own_cols + 1
            @test flow[1, c] == 3
        end

        @test flow[1, 1] == 1

        for r in 2:own_rows + 1
            @test flow[r, 1] == 0
        end
    elseif rank == 3
        for r in 2:own_rows + 1
            @test flow[r, 1] == 1
        end
    end

    # if rank == 0; display(flow) end; MPI.Barrier(cart_comm)
    # if rank == 1; display(flow) end; MPI.Barrier(cart_comm)
    # if rank == 2; display(flow) end; MPI.Barrier(cart_comm)
    # if rank == 3; display(flow) end; MPI.Barrier(cart_comm)
    '`))
end

@testset "AnnotateParticleLocation" begin
    mpiexec(cmd->run(`$cmd -np 4 julia --project=. -e '
    using Test
    include("wind2.jl")

    tunnel_rows = 10
    tunnel_cols = 10
    
    MPI.Init()
    comm = MPI.Comm_dup(MPI.COMM_WORLD)
    n_ranks = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    cart_dims = getCartDims(n_ranks)
    periodic = map(_ -> false, cart_dims)
    reorder = false
    cart_comm = MPI.Cart_create(comm, cart_dims; periodic, reorder)
    cart_coord = MPI.Cart_coords(cart_comm)
    cart_neighbors = getNeighborRanks(cart_dims, cart_coord, cart_comm)
    
    own_cols = tunnel_cols ÷ cart_dims[1]
    own_rows = tunnel_rows ÷ cart_dims[2]
    tsc = tunnelStartColumn(tunnel_cols, cart_dims, cart_coord)
    tsr = tunnelStartRow(tunnel_rows, cart_dims, cart_coord)

    p_locs = zeros(Int, own_rows + 2, own_cols + 2)
    p_locs[1, :] = p_locs[end, :] .= -1
    p_locs[:, 1] = p_locs[:, end] .= -1

    ps::Vector{Particle} = []

    append!(ARGS, ["tr", "tc", "mi", "th", "fp", "fs", "fbp", "fbs", "fd", "mbp", "mbs", "md", "rs1", "rs2", "rs3"
    , "1", "1", "0.3", "5", "5", "0.3", "1", "5", "0.3", "5", "1", "0.3", "4", "4", "0.3"
    , "6", "1", "0.3", "10", "1", "0.3", "6", "5", "0.3", "10", "5", "0.3", "7", "4", "0.3"
    , "1", "6", "0.3", "5", "6", "0.3", "1", "10", "0.3", "5", "10", "0.3", "3", "8", "0.3"
    , "6", "6", "0.3", "10", "6", "0.3", "10", "10", "0.3", "6", "10", "0.3", "7", "8", "0.3", "7", "8", "0.6"])
    
    lb = 1
    ub = (length(ARGS) - 15) ÷ 3
    readFixedParticles(lb, ub, tsr, tsc, own_rows, own_cols, ps)
    MPI.Barrier(cart_comm)

    AnnotateParticleLocation(ps, own_rows, own_cols, p_locs)

    # if rank == 0; println(rank); display(p_locs) end; MPI.Barrier(cart_comm)
    # if rank == 1; println(rank); display(p_locs) end; MPI.Barrier(cart_comm)
    # if rank == 2; println(rank); display(p_locs) end; MPI.Barrier(cart_comm)
    # if rank == 3; println(rank); display(p_locs) end; MPI.Barrier(cart_comm)

    @test p_locs[2, 2] == 1
    @test p_locs[6, 6] == 1
    @test p_locs[2, 6] == 1
    @test p_locs[6, 2] == 1
    if rank == 0
        @test p_locs[3, 5] == 1
    elseif rank == 1
        @test p_locs[5, 5] == 1
    elseif rank == 2
        @test p_locs[3, 4] == 2
    elseif rank == 3
        @test p_locs[4, 4] == 1
    end

    '`))
end

@testset "sendParticlesToNeighbors S SE E" begin
    mpiexec(cmd->run(`$cmd -np 4 julia --project=. -e '
    using Test
    include("wind2.jl")

    tunnel_rows = 10
    tunnel_cols = 10
    
    MPI.Init()
    comm = MPI.Comm_dup(MPI.COMM_WORLD)
    n_ranks = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    cart_dims = getCartDims(n_ranks)
    periodic = map(_ -> false, cart_dims)
    reorder = false
    cart_comm = MPI.Cart_create(comm, cart_dims; periodic, reorder)
    cart_coord = MPI.Cart_coords(cart_comm)
    cart_neighbors = getNeighborRanks(cart_dims, cart_coord, cart_comm)
    
    own_cols = tunnel_cols ÷ cart_dims[1]
    own_rows = tunnel_rows ÷ cart_dims[2]
    tsr = tunnelStartRow(tunnel_rows, cart_dims, cart_coord)
    tsc = tunnelStartColumn(tunnel_cols, cart_dims, cart_coord)
    
    n_ps = (length(ARGS) - 15) ÷ 3
    ips::Vector{Particle} = []
    ops::Vector{Particle} = []

    if rank == 1
        push!(ops, Particle([6, 5], 1, 2, [3, 4], 7)) # Send to 0
        push!(ops, Particle([6, 6], 1, 2, [3, 4], 7)) # Send to 2
        push!(ops, Particle([5, 6], 1, 2, [3, 4], 7)) # Send to 3
        
        sendParticlesToNeighbors(ops, cart_comm, cart_neighbors, tsr, tsc, own_rows, own_cols)
    end

    # println(" $rank : $ops. Sending")
    MPI.Barrier(cart_comm)

    if rank == 1
        @test length(ops) == 0
    else
        rb = Array{Int}(undef, 7)
        rs = MPI.Irecv!(rb, cart_comm;source=1, tag=0)
        MPI.Wait(rs)
        push!(ips, Particle([rb[1], rb[2]], rb[3], rb[4], [rb[5], rb[6]], rb[7]))
    end
    # println(" $rank : $ips. Receiving. ")
    '`))
end

@testset "sendParticlesToNeighbors SW W" begin
    mpiexec(cmd->run(`$cmd -np 4 julia --project=. -e '
    using Test
    include("wind2.jl")

    tunnel_rows = 10
    tunnel_cols = 10
    
    MPI.Init()
    comm = MPI.Comm_dup(MPI.COMM_WORLD)
    n_ranks = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    cart_dims = getCartDims(n_ranks)
    periodic = map(_ -> false, cart_dims)
    reorder = false
    cart_comm = MPI.Cart_create(comm, cart_dims; periodic, reorder)
    cart_coord = MPI.Cart_coords(cart_comm)
    cart_neighbors = getNeighborRanks(cart_dims, cart_coord, cart_comm)
    
    own_cols = tunnel_cols ÷ cart_dims[1]
    own_rows = tunnel_rows ÷ cart_dims[2]
    tsr = tunnelStartRow(tunnel_rows, cart_dims, cart_coord)
    tsc = tunnelStartColumn(tunnel_cols, cart_dims, cart_coord)
    
    n_ps = (length(ARGS) - 15) ÷ 3
    ips::Vector{Particle} = []
    ops::Vector{Particle} = []

    src = 3

    if rank == src
        push!(ops, Particle([6, 5], 1, 2, [3, 4], 7)) # Send to 0, SW
        push!(ops, Particle([5, 5], 1, 2, [3, 4], 7)) # Send to 1, W
        
        sendParticlesToNeighbors(ops, cart_comm, cart_neighbors, tsr, tsc, own_rows, own_cols)
    end

    MPI.Barrier(cart_comm)

    if rank == src
        @test length(ops) == 0
    elseif rank == 1
        rb = Array{Int}(undef, 7)
        rs = MPI.Irecv!(rb, cart_comm;source=src, tag=0)
        MPI.Wait(rs)
        
        push!(ips, Particle([rb[1], rb[2]], rb[3], rb[4], [rb[5], rb[6]], rb[7]))
        @test ips[1].position == [5, 5]
    elseif rank == 0
        rb = Array{Int}(undef, 7)
        rs = MPI.Irecv!(rb, cart_comm;source=src, tag=0)
        MPI.Wait(rs)
        
        push!(ips, Particle([rb[1], rb[2]], rb[3], rb[4], [rb[5], rb[6]], rb[7]))
        @test ips[1].position == [6, 5]
    else
        @test length(ips) == 0
    end

    '`))
end

@testset "recvParticlesFromNeighbors" begin
    mpiexec(cmd->run(`$cmd -np 4 julia --project=. -e '
    using Test
    include("wind2.jl")

    tunnel_rows = 10
    tunnel_cols = 10
    
    MPI.Init()
    comm = MPI.Comm_dup(MPI.COMM_WORLD)
    n_ranks = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    cart_dims = getCartDims(n_ranks)
    periodic = map(_ -> false, cart_dims)
    reorder = false
    cart_comm = MPI.Cart_create(comm, cart_dims; periodic, reorder)
    cart_coord = MPI.Cart_coords(cart_comm)
    cart_neighbors = getNeighborRanks(cart_dims, cart_coord, cart_comm)
    
    own_cols = tunnel_cols ÷ cart_dims[1]
    own_rows = tunnel_rows ÷ cart_dims[2]
    tsr = tunnelStartRow(tunnel_rows, cart_dims, cart_coord)
    tsc = tunnelStartColumn(tunnel_cols, cart_dims, cart_coord)
    
    n_ps = (length(ARGS) - 15) ÷ 3
    ips::Vector{Particle} = []
    ops::Vector{Particle} = []

    if rank == 0
        push!(ops, Particle([8, 6], 1, 2, [3, 4], 7)) # Send to 2
    elseif rank == 1
        push!(ops, Particle([6, 5], 1, 2, [3, 4], 7)) # Send to 0
        push!(ops, Particle([7, 5], 1, 2, [3, 4], 7)) # Send to 0
        push!(ops, Particle([2, 6], 1, 2, [3, 4], 7)) # Send to 3
    elseif rank == 2
        push!(ops, Particle([9, 5], 1, 2, [3, 4], 7)) # Send to 0
    elseif rank == 3
        push!(ops, Particle([3, 5], 1, 2, [3, 4], 7)) # Send to 1
    end

    sendParticlesToNeighbors(ops, cart_comm, cart_neighbors, tsr, tsc, own_rows, own_cols)

    recvParticlesFromNeighbors(ips, cart_neighbors, cart_comm)
    
    @test length(ops) == 0

    if rank == 0; @test length(ips) == 3 end
    if rank == 1; @test length(ips) == 1 end
    if rank == 2; @test length(ips) == 1 end
    if rank == 3; @test length(ips) == 1 end

    # if rank == 0; println(" $rank : $(length(ips)) elements : $ips. ") end; MPI.Barrier(cart_comm)
    # if rank == 1; println(" $rank : $(length(ips)) elements : $ips. ") end; MPI.Barrier(cart_comm)
    # if rank == 2; println(" $rank : $(length(ips)) elements : $ips. ") end; MPI.Barrier(cart_comm)
    # if rank == 3; println(" $rank : $(length(ips)) elements : $ips. ") end; MPI.Barrier(cart_comm)
    '`))
end

@testset "particleMovements" begin
    mpiexec(cmd->run(`$cmd -np 4 julia --project=. -e '
    using Test
    include("wind2.jl")

    append!(ARGS, ["tr", "tc", "mi", "th", "fp", "fs", "fbp", "fbs", "fd", "mbp", "mbs", "md", "rs1", "rs2", "rs3", 
    "2", "3", "0.5", "2", "3", "0.5", "2", "3", "0.5",
    "3", "3", "0.5",
    "2", "2", "0.5", "2", "2", "0.5", 
    "2", "7", "0.5", "2", "7", "0.5", "2", "7", "0.5", "2", "7", "0.5",
    "2", "8", "0.5", "2", "8", "0.5", "2", "8", "0.5", 
    "2", "9", "0.5", "2", "9", "0.5", "2", "9", "0.5", "2", "9", "0.5", "2", "9", "0.5",])

    tunnel_rows = 10
    tunnel_cols = 10
    fp = 1
    fs = 10
    
    MPI.Init()
    comm = MPI.Comm_dup(MPI.COMM_WORLD)
    n_ranks = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    cart_dims = getCartDims(n_ranks)
    periodic = map(_ -> false, cart_dims)
    reorder = false
    cart_comm = MPI.Cart_create(comm, cart_dims; periodic, reorder)
    cart_coord = MPI.Cart_coords(cart_comm)
    cart_neighbors = getNeighborRanks(cart_dims, cart_coord, cart_comm)
    
    rrr = rank

    own_cols = tunnel_cols ÷ cart_dims[1]
    own_rows = tunnel_rows ÷ cart_dims[2]
    tsr = tunnelStartRow(tunnel_rows, cart_dims, cart_coord)
    tsc = tunnelStartColumn(tunnel_cols, cart_dims, cart_coord)
    
    n_ps = (length(ARGS) - 15) ÷ 3
    ips::Vector{Particle} = []
    ops::Vector{Particle} = []
    readFixedParticles(1, n_ps, tsr, tsc, own_rows, own_cols, ips)
    
    flow = zeros(Int, own_rows + 2, own_cols + 2)
    flow_copy = zeros(Int, own_rows + 2, own_cols + 2)
    p_locs = zeros(Int, own_rows + 2, own_cols + 2)
    
    flow[1, :] = flow[end, :] .= -1; flow[:, 1] = flow[:, end] .= -1
    flow_copy[1, :] = flow_copy[end, :] .= -1; flow_copy[:, 1] = flow_copy[:, end] .= -1
    p_locs[1, :] = p_locs[end, :] .= -1; p_locs[:, 1] = p_locs[:, end] .= -1


    for iter in 1:1
    
        updateFan(iter, fp + 1, fs, flow, cart_dims, cart_coord, tunnel_cols, own_cols)
    
        particleMovements(iter, tsr, tsc, tunnel_rows, tunnel_cols, p_locs, flow, ips, ops, cart_comm, cart_dims, cart_coord, cart_neighbors)        
        
        # particleEffects(iter, ips, flow, flow_copy, p_locs, tunnel_cols)
    end

    if rank == 1; println("$rank : flow"); display(flow) end; MPI.Barrier(cart_comm)
    if rank == 0; println("$rank : flow"); display(flow) end; MPI.Barrier(cart_comm)
    if rank == 3; println("$rank : flow"); display(flow) end; MPI.Barrier(cart_comm)
    if rank == 2; println("$rank : flow"); display(flow) end; MPI.Barrier(cart_comm)
    
    if rank == 1; println("$rank : p_locs"); display(p_locs) end; MPI.Barrier(cart_comm)
    if rank == 0; println("$rank : p_locs"); display(p_locs) end; MPI.Barrier(cart_comm)
    if rank == 3; println("$rank : p_locs"); display(p_locs) end; MPI.Barrier(cart_comm)
    if rank == 2; println("$rank : p_locs"); display(p_locs) end; MPI.Barrier(cart_comm)
    '`))
end

# TODO: updateFlow
@testset "updateFlow" begin 
end

# TODO: particleEffects
@testset "particleEffects" begin
end