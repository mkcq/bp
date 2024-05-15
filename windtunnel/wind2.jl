using MPI
using Random

PRECISION = 10000
STEPS = 8

"""Position and speed as [row, column]."""
mutable struct Particle
    position::Vector
    mass::Int
    resistance::Int
    speed::Vector
    old_flow::Int
end


"""Given a size the function returns cartesian dimensions as [columns, rows].
|0        |1        |
|:--:     |:--:     |
|0 (0, 0) | 1 (1, 0)|

|0        |1        |
|:--:     |:--:     |
|1 (0, 1) | 3 (1, 1)|
|0 (0, 0) | 2 (1, 0)|

|0        |1        |
|:--:     |:--:     |
|3 (0, 3) | 7 (1, 3)|
|2 (0, 2) | 6 (1, 2)|
|1 (0, 1) | 5 (1, 1)|
|0 (0, 0) | 4 (1, 0)|"""
function getCartDims(n_ranks)
    if n_ranks == 1
        return [1, 1]
    elseif n_ranks == 2
        return [2, 1]
    elseif n_ranks == 4
        return [2, 2]
    elseif n_ranks == 8
        return [2, 4]
    elseif n_ranks == 16
        return [4, 4]
    end
end


"""tunnel_rows ÷ cart_dims[rows] * cart_coord[row]"""
tunnelStartRow(tunnel_rows, cart_dims, cart_coord) = tunnel_rows ÷ cart_dims[2] * (cart_dims[2] - 1 - cart_coord[2]) + 1


"""tunnel_cols ÷ cart_dims[columns] * cart_coord[column]"""
tunnelStartColumn(tunnel_cols, cart_dims, cart_coord) = tunnel_cols ÷ cart_dims[1] * cart_coord[1] + 1


"""Returns true if [r, c] is in own cartesian grid cell. Otherwise it returns false."""
inCartGrid(tsr, tsc, r, c, own_rows, own_cols) = (tsr <= r && r < tsr + own_rows && tsc <= c && c < tsc + own_cols) ? true : false


"""From lb until ub times, the function adds a particle to its own list if the particle is in own cartesian grid cell."""
function readFixedParticles(lb, ub, tsr, tsc, own_rows, own_cols, particles)
    for i in lb:ub
        row = parse(Int, ARGS[13 + i * 3])
        col = parse(Int, ARGS[14 + i * 3])
        if !inCartGrid(tsr, tsc, row, col, own_rows, own_cols)
            continue
        end
        mass = 0
        resistance = trunc(Int, parse(Float16, ARGS[15 + i * 3]) * PRECISION)
        speed = [0, 0]
        push!(particles, Particle([row * PRECISION, col * PRECISION], mass, resistance, speed, 0))
    end
end


"""Update only the first (NOT the ghostcells) row of the ranks at the top of the topology."""
function updateFan(iter, fan_pos, fan_size, flow, cart_dims, cart_coord, tunnel_cols, own_cols)
    if !(cart_dims[2] == cart_coord[2] + 1) || !(iter % STEPS == 1)
        # println("$(cart_dims[2]) == $(cart_coord[2] + 1)")
        return
    end

    tsc = tunnelStartColumn(tunnel_cols, cart_dims, cart_coord)

    for j in 2:own_cols + 1
        if (fan_pos < tsc + j) && (tsc + j <= fan_pos + fan_size)
            phase = iter ÷ STEPS * (π ÷ 4)
            phase_step = π ÷ 2 ÷ fan_size
            pressure_level = 9 + 2 * sin(phase + (j - fan_pos) * phase_step)
            noise = 0.5 - rand(1)[1]
            # Store level in the first row (NOT in the ghostcell) of the ancillary structure.
            flow[2, j] = trunc(Int, PRECISION * (pressure_level + noise))
        end
    end
end


"""Returns the neighbor ranks based on the cartesian topology of getCartDims."""
function getNeighborRanks(cart_dims, cart_coord, cart_comm)
    north = cart_coord[2] + 1 < cart_dims[2] ? cart_coord[2] + 1 : -1
    south = cart_coord[2] - 1 >= cart_dims[2] - 2 ? cart_coord[2] - 1 : -1
    east = cart_coord[1] + 1 < cart_dims[1] ? cart_coord[1] + 1 : -1
    west = cart_coord[1] - 1 >= cart_dims[1] - 2 ? cart_coord[1] - 1 : -1
    
    rank_north = north == -1 ? -1 : MPI.Cart_rank(cart_comm, [cart_coord[1], north])
    rank_south = south == -1 ? -1 : MPI.Cart_rank(cart_comm, [cart_coord[1], south])
    rank_east = east == -1 ? -1 : MPI.Cart_rank(cart_comm, [east, cart_coord[2]])
    rank_west = west == -1 ? -1 : MPI.Cart_rank(cart_comm, [west, cart_coord[2]])
    
    rank_ne = (north != -1 && east != -1) ? MPI.Cart_rank(cart_comm, [east, north]) : -1 
    rank_nw = (north != -1 && west != -1) ? MPI.Cart_rank(cart_comm, [west, north]) : -1 
    rank_se = (south != -1 && east != -1) ? MPI.Cart_rank(cart_comm, [east, south]) : -1 
    rank_sw = (south != -1 && west != -1) ? MPI.Cart_rank(cart_comm, [west, south]) : -1 

    Dict([("N", rank_north), ("S", rank_south), ("E", rank_east), ("W", rank_west), ("NE", rank_ne), ("NW", rank_nw), ("SE", rank_se), ("SW", rank_sw)])
end


""""""
function moveParticle(flow, particles, p, tunnel_rows, tunnel_cols, cart_dims, cart_coord, cart_comm)
    for i in 1:STEPS
        # TODO: I could combine particles and p into one.
    
        r = particles[p].position[1] ÷ PRECISION
        c = particles[p].position[2] ÷ PRECISION
        
        # TODO: Left here.
        neighbors = getNeighborRanks(cart_dims, cart_coord, cart_comm)

        # TODO: Data dependency here.
        pressure = flow[r - 1, c]
        
        # TODO: Data dependency here.
        left = 0; right = 0
        c == 2 ? left = 0 : left = pressure - flow[r - 1, c - 1]
        c == own_cols + 1 ? right = 0 : right = pressure - flow[r - 1, c + 1]

        flow_row = trunc(Int, pressure + particles[p].mass * PRECISION)
        flow_col = trunc(Int, (right - left) ÷ particles[p].mass * PRECISION)

        # Speed change.
        particles[p].speed[1] = (particles[p].speed[1] + flow_row) ÷ 2
        particles[p].speed[2] = (particles[p].speed[2] + flow_col) ÷ 2
        
        # Movement.
        particles[p].position[1] = particles[p].position[1] + particles[p].speed[1] ÷ STEPS ÷ 2
        particles[p].position[2] = particles[p].position[2] + particles[p].speed[2] ÷ STEPS ÷ 2
    
        # TODO: Control limits.
    end
end


""""""
function particleMovements(iter, tunnel_rows, tunnel_cols, p_locs, particles, flow)
    if !(iter % STEPS == 1); return end

    # TODO: Clean particle positions.
    
    # Move particles.
    ub = length(particles)
    for p in 1:ub
        (particles[p].mass == 0) ? continue : moveParticle(flow, particles, p, tunnel_rows, tunnel_cols)
    end
    
    # TODO: Annotate position.
end


"""Pushes a particle from the ARGS to particles only if row and column are in the cartesian grid cell."""
function readFixedParticles(lb, ub, tsr, tsc, own_rows, own_cols, particles)
    for i in lb:ub
        r = parse(Int, ARGS[13 + i * 3])
        c = parse(Int, ARGS[14 + i * 3])
        if !inCartGrid(tsr, tsc, r, c, own_rows, own_cols)
            continue
        end
        mass = 0
        resistance = trunc(Int, parse(Float16, ARGS[15 + i * 3]) * PRECISION)
        speed = [0, 0]
        push!(particles, Particle([r, c] * PRECISION, mass, resistance, speed, 0))
    end
end


"""For DAS the arguments are passed as parameters to main().
For local execution the arguments are in the command line.

Order of arguments: rows, cols, max_iter, threshold, fan_pos, fan_size, fixed_band_pos, fixed_band_size, fixed_density, moving_band_pos, moving_band_size, moving_density, rand_seq_1, rand_seq_2, rand_seq_3, particle_row, particle_col, particle_density, ...
"""
function main()
    # Simulation data.
    tunnel_rows = parse(Int, ARGS[1])
    tunnel_cols = parse(Int, ARGS[2])
    max_iter = parse(Int, ARGS[3])              # Max number of simulation steps.
    threshold = parse(Float16, ARGS[4])         # Threshold of variability to continue the simulation.
    fan_pos = parse(Int, ARGS[5])               # Fan position in the tunnel.
    fan_size = parse(Int, ARGS[6])              # Fan size in the tunnel.

    fixed_band_pos = parse(Int, ARGS[7])        # First position of the band where FIXED particles Start.
    fixed_band_size = parse(Int, ARGS[8])       # Size of the band where FIXED particles start.
    fixed_density = parse(Float16, ARGS[9])     # Density of starting FIXED particles.
    
    moving_band_pos = parse(Int, ARGS[10])      # First position of the band where MOVING particles start.
    moving_band_size = parse(Int, ARGS[11])     # Size of the band where MOVING particles start.
    moving_density = parse(Float16, ARGS[12])   # Density of starting MOVING particles.

    # Initialize random sequences.
    random_seq = [parse(Int, ARGS[i]) for i in 13:15]
    seed = rand(random_seq)
    Random.seed!(seed)

    # Initialize MPI.
    MPI.Init()
    comm = MPI.Comm_dup(MPI.COMM_WORLD)
    n_ranks = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)

    # Initialize MPI cartesian grid.
    cart_dims = getCartDims(n_ranks)
    periodic = map(_ -> false, cart_dims)
    reorder = false
    cart_comm = MPI.Cart_create(comm, cart_dims; periodic, reorder)
    cart_coord = MPI.Cart_coords(cart_comm)

    # Initialize own rows and cols.
    own_cols = tunnel_cols ÷ cart_dims[1]
    own_rows = tunnel_rows ÷ cart_dims[2]

    # Read particles form the arugments and only add them if the are within own cartesian grid cell.
    n_particles = (length(ARGS) - 15) ÷ 3
    particles::Vector{Particle} = []

    tsr = tunnelStartRow(tunnel_rows, cart_dims, cart_coord)
    tsc = tunnelStartColumn(tunnel_cols, cart_dims, cart_coord)
    readFixedParticles(lb, ub, tsr, tsc, own_rows, own_cols, particles)

    # Initialization for parallelization.
    flow = zeros(Int, own_rows + 2, own_cols + 2)          # With ghostcell. Tunnel air flow.
    p_locs = zeros(Int, own_rows + 2, own_cols + 2)        # With ghostcell. Quickly locate particle positions.
    
    # Initialize ghostcells with value -1.
    flow[1, :] = flow[end, :] .= -1; flow[:, 1] = flow[:, end] .= -1
    p_locs[1, :] = p_locs[end, :] .= -1; p_locs[:, 1] = p_locs[:, end] .= -1
    
    # Simulation
    for iter in 1:max_iter
        # Change the fan values each STEP iterations.
        updateFan(iter, fan_pos + 1, fan_size, flow, cart_dims, cart_coord, tunnel_cols, own_cols)    
        
        # Particles movement each STEPS iterations.
        particleMovements(iter, tunnel_rows, tunnel_cols, p_locs, particles, flow)
        
        # Effects due to particles each STEPS iterations.

        # Copy data in ancillary structure.

        # Propagation stage.
    end

    # TODO: Stop global timer.
    MPI.Barrier(cart_comm)

    # TODO: Print results.
end

# Comment out call to main() whenever running tests.jl.
# main()