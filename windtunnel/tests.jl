using Test
include("wind.jl")
#=
|0        |1        |
|:--:     |:--:     |
|0 (0, 0) | 1 (1, 0)|
dims[col, row]

cols = 30   own_col = cols ÷ dims[1] = 15
rows = 2    own_row = rows ÷ dims[2] = 2


=#

@testset "updateFan" begin
    maxiter = 10; fp = 4; fs = 17; dims = [2, 1]
    coordl = [0, 0]
    coordr = [1, 0]
    target_row = 2
    cols = 30
    coll = cols ÷ dims[1]
    colr = cols ÷ dims[1]
    flowl = zeros(Int, 3+2, coll+2)
    flowr = zeros(Int, 3+2, colr+2)
    for iter in 1:maxiter
        updateFan(iter, fp + 1, fs, flowl, dims, coordl, target_row, cols, coll)
        updateFan(iter, fp + 1, fs, flowr, dims, coordr, target_row, cols, colr)
    end
    res = vcat(view(flowl, target_row, 2:coll+1), view(flowr, target_row, 2:colr+1))
    @show res
    for (i,v) in enumerate(res)
        if fp <= i && i <= fp+fs
            @test v != 0
        end
    end
end

@testset "moveParticle" begin
    
    # moveParticle(flow, particles, p, rows, cols)
end

@testset "particleMovements" begin
    iter = 1; rows = 10; cols = 20; dims = [2, 1]
    coorl = [0, 0]; coordr = [1, 0]
    coll = cols ÷ dims[1]; rowl = rows ÷ dims[2]
    colr = cols ÷ dims[1]; rowr = rows ÷ dims[2]

    flowl = zeros(Int, rowl + 2, coll + 2); flowr = zeros(Int, rowr + 2, colr + 2)
    pll = zeros(Int, rowl + 2, coll + 2); plr = zeros(Int, rowr + 2, colr + 2)

    particles = []

    # particleMovements(iter, rows, cols, particle_locs, particles, flow)
end
