using Test
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

@testset "updateFan 2cx2r" begin
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

@testset "updateFan 2cx1r" begin
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



# @testset "particleMovements" begin
#     iter = 1; rows = 10; cols = 20; dims = [2, 1]
#     coorl = [0, 0]; coordr = [1, 0]
#     coll = cols ÷ dims[1]; rowl = rows ÷ dims[2]
#     colr = cols ÷ dims[1]; rowr = rows ÷ dims[2]

#     flowl = zeros(Int, rowl + 2, coll + 2); flowr = zeros(Int, rowr + 2, colr + 2)
#     pll = zeros(Int, rowl + 2, coll + 2); plr = zeros(Int, rowr + 2, colr + 2)

#     particles = []

#     # particleMovements(iter, rows, cols, particle_locs, particles, flow)
# end
