using MPI
using Random

PRECISION = 10000
STEPS = 8

global rrr = 0

"""Position as [tunnel_row, tunnel_column] and speed as [row, column]."""
mutable struct Particle
    position::Vector
    mass::Int
    resistance::Int
    speed::Vector
    old_flow::Int
end


"""Given n_ranks ranks with the segments in matrices, the function returns a combined matrix.
The values of n_ranks can be 1, 2, 4, 8, 16.
The way in which the returned matrix is ordered depends on the value of getDims(n_ranks)."""
function combineMatrices(n_ranks, matrices, range_rows, range_cols)
    if n_ranks == 1
        return [view(matrices[1], range_rows, range_cols)]
    elseif n_ranks == 2
        return [view(matrices[1], range_rows, range_cols) view(matrices[2], range_rows, range_cols)]
    elseif n_ranks == 4
        return [view(matrices[2], range_rows, range_cols) view(matrices[4], range_rows, range_cols)
            ;view(matrices[1], range_rows, range_cols) view(matrices[3], range_rows, range_cols)]
    elseif n_ranks == 8
        return [view(matrices[4], range_rows, range_cols) view(matrices[8], range_rows, range_cols)
        ;view(matrices[3], range_rows, range_cols) view(matrices[7], range_rows, range_cols)
        ;view(matrices[2], range_rows, range_cols) view(matrices[6], range_rows, range_cols)
        ;view(matrices[1], range_rows, range_cols) view(matrices[5], range_rows, range_cols)]
    elseif n_ranks == 16
        return [view(matrices[4], range_rows, range_cols) view(matrices[8], range_rows, range_cols) view(matrices[12], range_rows, range_cols) view(matrices[16], range_rows, range_cols)
        ;view(matrices[3], range_rows, range_cols) view(matrices[7], range_rows, range_cols) view(matrices[11], range_rows, range_cols) view(matrices[15], range_rows, range_cols)
        ;view(matrices[2], range_rows, range_cols) view(matrices[6], range_rows, range_cols) view(matrices[10], range_rows, range_cols) view(matrices[14], range_rows, range_cols)
        ;view(matrices[1], range_rows, range_cols) view(matrices[5], range_rows, range_cols) view(matrices[9], range_rows, range_cols) view(matrices[13], range_rows, range_cols)]
    end
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


"""Returns the starting index of the first row in a cartesian cell as the actual tunnel row.
E.g. Assume a 10 x 10 tunnel divided into a 2 x 2 cartesian grid.

First index topleft (0, 1) corner = 10 ÷ 2 * (2 - 1 - 1) + 1 = 1

First index bottomleft (0, 0) corner = 10 ÷ 2 * (2 - 1 - 0) + 1 = 6

tunnel_rows ÷ cart_dims[rows] * (cart_dims[rows] - 1 - cart_coord[row]) + 1"""
tunnelStartRow(tunnel_rows, cart_dims, cart_coord) = tunnel_rows ÷ cart_dims[2] * (cart_dims[2] - 1 - cart_coord[2]) + 1


"""tunnel_cols ÷ cart_dims[columns] * cart_coord[column]"""
tunnelStartColumn(tunnel_cols, cart_dims, cart_coord) = tunnel_cols ÷ cart_dims[1] * cart_coord[1] + 1


"""Returns true if [r, c] is in own cartesian grid cell. Otherwise it returns false."""
inCartGrid(tsr, tsc, r, c, own_rows, own_cols) = (tsr <= r && r < tsr + own_rows && tsc <= c && c < tsc + own_cols) ? true : false


"""Pushes a particle from the ARGS to particles only if row and column are in the cartesian grid cell."""
function readFixedParticles(lb, ub, tsr, tsc, own_rows, own_cols, particles)
    for i in lb:ub
        r = parse(Int, ARGS[13 + i * 3])
        c = parse(Int, ARGS[14 + i * 3])
        if !inCartGrid(tsr, tsc, r, c, own_rows, own_cols)
            continue
        end
        # For fixed mass = 0.
        # TODO: Turn mass back to 0.
        mass = trunc(Int, PRECISION * (1 + 5  * rand(1)[1]))
        resistance = trunc(Int, parse(Float16, ARGS[15 + i * 3]) * PRECISION)
        speed = [0, 0]
        push!(particles, Particle([r, c] * PRECISION, mass, resistance, speed, 0))
    end
end


"""Update only the first (NOT the ghostcells) row of the ranks at the top of the topology. It generates new values in waves after STEPS iterations."""
function updateFan(iter, fan_pos, fan_size, flow, cart_dims, cart_coord, tunnel_cols, own_cols)
    if !(cart_dims[2] == cart_coord[2] + 1) || !(iter % STEPS == 1)
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


"""Function does a non-blocking send to distribute the flow to below. The directions are: S, SW, SE, W and E."""
function sendFlowToNeighbors(flow, cart_neighbors, cart_comm)
    reqs = MPI.Request[]
    for (k, v) in cart_neighbors
        if k == "S" && v >= 0
            req = MPI.Isend(view(flow, own_rows + 1, 2:own_cols + 1), cart_comm;dest=v, tag=0)
            push!(reqs, req)
        elseif k == "SW" && v >= 0
            req = MPI.Isend(view(flow, own_rows + 1, 2), cart_comm;dest=v, tag=0)
            push!(reqs, req)
        elseif k == "SE" && v >= 0
            req = MPI.Isend(view(flow, own_rows + 1, own_cols + 1), cart_comm;dest=v, tag=0)
            push!(reqs, req)
        elseif k == "W" && v >= 0
            req = MPI.Isend(view(flow, 2:own_rows + 1, 2), cart_comm;dest=v, tag=0)
            push!(reqs, req)
        elseif k == "E" && v >= 0
            req = MPI.Isend(view(flow, 2:own_rows + 1, own_cols + 1), cart_comm;dest=v, tag=0)
            push!(reqs, req)
        end
    end
    MPI.Waitall(reqs)
end


"""Function does a non-blocking receive to gather the flow from above. The directions are: N, NW, NE, W and E."""
function recvFlowFromNeighbors(flow, cart_neighbors, cart_comm)
    reqs = MPI.Request[]
    for (k, v) in cart_neighbors
        if k == "N" && v >= 0
            req = MPI.Irecv!(view(flow, 1, 2:own_cols + 1), cart_comm;source=v, tag=0)
            push!(reqs, req)
        elseif k == "NW" && v >= 0
            req = MPI.Irecv!(view(flow, 1, 1), cart_comm;source=v, tag=0)
            push!(reqs, req)
        elseif k == "NE" && v >= 0
            req = MPI.Irecv!(view(flow, 1, own_cols + 2), cart_comm;source=v, tag=0)
            push!(reqs, req)
        elseif k == "W" && v >= 0
            req = MPI.Irecv!(view(flow, 2:own_rows + 1, 1), cart_comm;source=v, tag=0)
            push!(reqs, req)
        elseif k == "E" && v >= 0
            req = MPI.Irecv!(view(flow, 2:own_rows + 1, own_cols + 2), cart_comm;source=v, tag=0)
            push!(reqs, req)
        end
    end
    MPI.Waitall(reqs)
end


"""Based on the flow, the function calculates the new position of a particle p.

Returns 1 if p is outside the boundaries of the rank. Otherwise returns 0."""
function moveParticle(flow, p, tsr, tsc, own_rows, own_cols, tunnel_rows, tunnel_cols, cart_comm, cart_neighbors)
    r = p.position[1] ÷ PRECISION
    c = p.position[2] ÷ PRECISION

    # Because flow is a (own_rows + 2) x (own_cols + 2) matrix.
    fr = r % own_rows == 0 ? own_rows + 1 : r % own_rows + 1
    fc = c % own_cols == 0 ? own_cols + 1 : c % own_cols + 1

    # Data dependencies received by recvFlowFromNeighbors.
    pressure = flow[fr - 1, fc]
    left = c == 1 ? 0 : pressure - flow[fr - 1, fc - 1]
    right = c == tunnel_cols ? 0 : pressure - flow[fr - 1, fc + 1]

    flow_row = trunc(Int, pressure ÷ p.mass * PRECISION)
    flow_col = trunc(Int, (right - left) ÷ p.mass * PRECISION)

    # Speed change.
    p.speed[1] = (p.speed[1] + flow_row) ÷ 2
    p.speed[2] = (p.speed[2] + flow_col) ÷ 2
    
    # Movement.
    p.position[1] = p.position[1] + p.speed[1] ÷ STEPS ÷ 2
    p.position[2] = p.position[2] + p.speed[2] ÷ STEPS ÷ 2        
    
    # Control limits. Particle must stay in the tunnel.
    if p.position[1] > tunnel_rows * PRECISION
        p.position[1] = tunnel_rows * PRECISION
    end
    if p.position[2] < 1
        p.position[2] = 1
    elseif p.position[2] > tunnel_cols * PRECISION
        p.position[2] = tunnel_cols * PRECISION
    end
    
    if !inCartGrid(tsr, tsc, r, c, own_rows, own_cols)
        return 1
    end

    return 0
end


"""Function first appends outgoing_particles in vectors of following directions: S, SW, SE, W and E. 
If there are no outgoing_particles, then the vectors contain one value of -1.

Second, the function does a non-blocking send to all directions mentioned above.

Last, after sending the outgoing_particles, all elements of outgoing_particles are deleted.
"""
function sendParticlesToNeighbors(outgoing_particles, cart_comm, cart_neighbors, tsr, tsc, own_rows, own_cols)
    
    reqs = MPI.Request[]
    send_SE::Vector{Int} = []
    send_SW::Vector{Int} = []
    send_S::Vector{Int} = []
    send_E::Vector{Int} = []
    send_W::Vector{Int} = []

    # Insert all outgoing_particles in the buffer corresponding to the neighbor.
    if length(outgoing_particles) > 0
        for p in outgoing_particles
            send_particle = [p.position[1], p.position[2], p.mass, p.resistance, p.speed[1], p.speed[2], p.old_flow]
            if tsr + own_rows - 1 < p.position[1]
                if tsc + own_cols - 1 < p.position[2]
                    append!(send_SE, send_particle)
                elseif p.position[2] < tsc
                    append!(send_SW, send_particle)
                else
                    append!(send_S, send_particle)
                end
            else
                if p.position[2] < tsc
                    append!(send_W, send_particle)
                else
                    append!(send_E, send_particle)
                end
            end
        end
    end

    # If there is something to send, then insert in correct request buffer.
    for (k, v) in cart_neighbors
        if k == "SE" && v >= 0
            if length(send_SE) == 0
                send_SE = [-1]
            end
            send_req = MPI.Isend(send_SE, cart_comm;dest=v, tag=0)
            push!(reqs, send_req)
        elseif k == "SW" && v >= 0
            if length(send_SW) == 0
                send_SW = [-1]
            end
            send_req = MPI.Isend(send_SW, cart_comm;dest=v, tag=0)
            push!(reqs, send_req)
        elseif k == "S" && v >= 0
            if length(send_S) == 0
                send_S = [-1]
            end
            send_req = MPI.Isend(send_S, cart_comm;dest=v, tag=0)
            push!(reqs, send_req)
        elseif k == "W" && v >= 0
            if length(send_W) == 0
                send_W = [-1]
            end
            send_req = MPI.Isend(send_W, cart_comm;dest=v, tag=0)
            push!(reqs, send_req)
        elseif k == "E" && v >= 0
            if length(send_E) == 0
                send_E = [-1]
            end
            send_req = MPI.Isend(send_E, cart_comm;dest=v, tag=0)
            push!(reqs, send_req)
        end
    end
    
    MPI.Waitall(reqs)

    while true
        length(outgoing_particles) == 0 ? break : deleteat!(outgoing_particles, 1)
    end
end


"""Receives from one of the following directions: N, NW, NE, W, E. 
Appends total_recv only if the size of received message > 1."""
function recvParticlesFromDirection(cart_comm, src, total_recv)
    status = MPI.Probe(cart_comm, MPI.Status;source=src, tag=0)
    num = MPI.Get_count(status, Int)
    recvbuf = Array{Int}(undef, num)
    MPI.Recv!(recvbuf, cart_comm;source=src, tag=0)
    if num > 1
        append!(total_recv, recvbuf)
    end
end


"""Function goes over all neighbors and receives particles from them. 
Received particles are insterted in incoming_particles."""
function recvParticlesFromNeighbors(incoming_particles, cart_neighbors, cart_comm)

    total_recv::Vector{Int} = []

    for (k, v) in cart_neighbors
    
        if k == "N" && v >= 0
            recvParticlesFromDirection(cart_comm, v, total_recv)
        elseif k == "NW" && v >= 0
            recvParticlesFromDirection(cart_comm, v, total_recv)
        elseif k == "NE" && v >= 0 
            recvParticlesFromDirection(cart_comm, v, total_recv)
        elseif k == "W" && v >= 0
            recvParticlesFromDirection(cart_comm, v, total_recv)
        elseif k == "E" && v >= 0
            recvParticlesFromDirection(cart_comm, v, total_recv)
        end

    end

    partitioned = collect(Iterators.partition(total_recv, 7))
    for i in partitioned
        push!(incoming_particles, Particle([i[1], i[2]], i[3], i[4], [i[5], i[6]], i[7]))
    end
end


"""Sets the values of p_locs to 0 when either iter <= tunnel_rows or on all tunnel_rows."""
function cleanParticleLocations(p_locs, own_rows, own_cols, iter, cart_dims, cart_coord)
    if iter <= tunnel_rows        
        ub = min(iter - own_rows * (cart_dims[2] - 1 - cart_coord[2]), own_rows) + 1
        for i in 2:ub , j in 2:own_cols + 1
            p_locs[i, j] = 0
        end
    else
        for i in 2:own_rows + 1, j in 2:own_cols + 1
            p_locs[i, j] = 0 
        end
    end
end


"""Updates the values of p_locs with the particles in own rank."""
function AnnotateParticleLocation(incoming_particles, own_rows, own_cols, p_locs)
    for p in incoming_particles
        r = p.position[1] ÷ PRECISION
        c = p.position[2] ÷ PRECISION
        # Because p_locs is a (own_rows + 2) x (own_cols + 2) matrix.
        pr = r % own_rows == 0 ? own_rows + 1 : r % own_rows + 1
        pc = c % own_cols == 0 ? own_cols + 1 : c % own_cols + 1
        p_locs[pr, pc] += 1
    end
end


"""
The outline of this function:
1. Clean particle locations.
2. Update flow by receving from neighbors.
3. Move particles with flow.
4. If odd first send then receive.
5. If even first receive then send.
6. Annotate particle locations.
"""
function particleMovements(iter, tsr, tsc, tunnel_rows, tunnel_cols, 
    p_locs, flow,
    incoming_particles, outgoing_particles, 
    cart_comm, cart_dims, cart_coord, cart_neighbors)
    if !(iter % STEPS == 1); return end

    # Clean particle locations in rank.
    cleanParticleLocations(p_locs, own_rows, own_cols, iter, cart_dims, cart_coord)
    
    # Move particles. First send/recv the flow of neighbors, then move the particles.    
    sendFlowToNeighbors(flow, cart_neighbors, cart_comm)
    recvFlowFromNeighbors(flow, cart_neighbors, cart_comm)

    i = 1
    del = 0
    while i <= length(incoming_particles)
        if incoming_particles[i].mass != 0
            for j in 1:STEPS
                del = moveParticle(flow, incoming_particles[i], tsr, tsc, own_rows, own_cols, tunnel_rows, tunnel_cols, cart_comm, cart_neighbors)
                if del == 1
                    break
                end
            end

            if del == 1
                push!(outgoing_particles, incoming_particles[i])
                deleteat!(incoming_particles, i)
                del = 0
            end
        end
        i += 1
    end

    sendParticlesToNeighbors(outgoing_particles, cart_comm, cart_neighbors, tsr, tsc, own_rows, own_cols)
    recvParticlesFromNeighbors(incoming_particles, cart_neighbors, cart_comm)
    
    # Annotate particle location in rank.
    AnnotateParticleLocation(incoming_particles, own_rows, own_cols, p_locs)
end


""""""
function updateFlow(flow, flow_copy, p_locs, r, c, fr, fc, tunnel_cols, skip_particles)
    # Skip update in particle positions.
    if (skip_particles && p_locs[r, c] != 0); return 0 end

    # Update if border left.
    if c == 1
        flow[fr, fc] = (flow_copy[fr, fc] + flow_copy[fr - 1, fc] * 2 + flow_copy[fr - 1, fc + 1]) ÷ 4
    end

    # Update if border right.
    if c == tunnel_cols
        flow[fr, fc] = (flow_copy[fr, fc] + flow_copy[fr - 1, fc] * 2 + flow_copy[fr - 1, fc - 1]) ÷ 4
    end

    # Update if central part.
    if c > 1 && c < tunnel_cols
        flow[fr, fc] = (flow_copy[fr, fc] + flow_copy[fr - 1, fc] * 2 + flow_copy[fr - 1, fc - 1] + flow_copy[fr - 1, fc + 1]) ÷ 5
    end

    # Return flow variation at this position.
    return abs(flow_copy[fr, fc] - flow[fr, fc])
end


""""""
function particleEffects(iter, incoming_particles, flow, flow_copy, p_locs, tunnel_cols, own_rows, own_cols)
    if !(iter % STEPS == 1); return end

    for p in incoming_particles
        r = p.position[1] ÷ PRECISION
        c = p.position[2] ÷ PRECISION
        # Because flow is a (own_rows + 2) x (own_cols + 2) matrix.
        fr = r % own_rows == 0 ? own_rows + 1 : r % own_rows + 1
        fc = c % own_cols == 0 ? own_cols + 1 : c % own_cols + 1
        updateFlow(flow, flow_copy, p_locs, r, c, fr, fc, tunnel_cols, false)
        p.old_flow = flow[fr, fc]
    end
    
    for p in incoming_particles
        r = p.position[1] ÷ PRECISION
        c = p.position[2] ÷ PRECISION
        # Because flow is a (own_rows + 2) x (own_cols + 2) matrix.
        fr = r % own_rows == 0 ? own_rows + 1 : r % own_rows + 1
        fc = c % own_cols == 0 ? own_cols + 1 : c % own_cols + 1
        back = trunc(Int, (p.old_flow * p.resistance ÷ PRECISION) ÷ p_locs[fr, fc])
        flow[fr, fc] -= back
        flow[fr - 1, fc] += back ÷ 2

        c > 1 ? (flow[fr - 1, fc - 1] += back ÷ 4) : (flow[fr - 1, fc] += back ÷ 4)

        c < tunnel_cols ? (flow[fr - 1, fc + 1] += back ÷ 4) : (flow[fr - 1, fc] += back ÷ 4)
    end
end


""""""
function propagateWaveFront(iter, tunnel_rows, tunnel_cols, tsr, tsc, own_rows, own_cols, flow, flow_copy, p_locs)

    wave_front = iter % STEPS

    if wave_front == 1
        max_var = 0
    end

    if wave_front == 0
        wave_front = STEPS
    end

    wave = wave_front
    own_wave_fronts::Vector{Int} = []
    for wave in 1:tunnel_rows
        if inCartGrid(tsr, tsc, wave, tsc, own_rows, own_cols)
            push!(own_wave_fronts, wave)
        end
    end

    for wave in own_wave_fronts

        if wave > iter
            break
        end
        
        for c in 1:tunnel_cols
            if inCartGrid(tsr, tsc, wave, c, own_rows, own_cols)
                # Because flow is a (own_rows + 2) x (own_cols + 2) matrix.
                fr = r % own_rows == 0 ? own_rows + 1 : r % own_rows + 1
                fc = c % own_cols == 0 ? own_cols + 1 : c % own_cols + 1
                var = updateFlow(flow, flow_copy, p_locs, wave, c, fr, fc, tunnel_cols, true)
                if var > max_var
                    max_var = var
                end
            end
        end
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
    cart_neighbors = getNeighborRanks(cart_dims, cart_coord, cart_comm)

    # Initialize own rows and cols.
    own_cols = tunnel_cols ÷ cart_dims[1]
    own_rows = tunnel_rows ÷ cart_dims[2]

    # Read particles form the arugments and only add them if the are within own cartesian grid cell.
    n_particles = (length(ARGS) - 15) ÷ 3
    incoming_particles::Vector{Particle} = []
    outgoing_particles::Vector{Particle} = []

    tsr = tunnelStartRow(tunnel_rows, cart_dims, cart_coord)
    tsc = tunnelStartColumn(tunnel_cols, cart_dims, cart_coord)
    readFixedParticles(1, n_particles, tsr, tsc, own_rows, own_cols, incoming_particles)

    # Initialization for parallelization.
    flow = zeros(Int, own_rows + 2, own_cols + 2)          # With ghostcell. Tunnel air flow.
    flow_copy = zeros(Int, own_rows + 2, own_cols + 2)          # With ghostcell. Tunnel air flow.
    p_locs = zeros(Int, own_rows + 2, own_cols + 2)        # With ghostcell. Quickly locate particle positions.
    
    # Initialize ghostcells with value -1.
    flow[1, :] = flow[end, :] .= -1; flow[:, 1] = flow[:, end] .= -1
    flow_copy[1, :] = flow_copy[end, :] .= -1; flow_copy[:, 1] = flow_copy[:, end] .= -1
    p_locs[1, :] = p_locs[end, :] .= -1; p_locs[:, 1] = p_locs[:, end] .= -1

    max_var = typemax(Int64)

    # Simulation
    for iter in 1:max_iter
        
        if max_var <= threshold
            break
        end

        # Change the fan values each STEP iterations.
        updateFan(iter, fan_pos + 1, fan_size, flow, cart_dims, cart_coord, tunnel_cols, own_cols)    
        
        # Particles movement each STEPS iterations.
        particleMovements(iter, tsr, tsc, tunnel_rows, tunnel_cols, p_locs, flow, incoming_particles, outgoing_particles, cart_comm, cart_dims, cart_coord, cart_neighbors)        
        
        # Effects due to particles each STEPS iterations.
        particleEffects(iter, incoming_particles, flow, flow_copy, p_locs, tunnel_cols, own_rows, own_cols)

        # Copy data in ancillary structure.
        copy!(flow_copy, flow)

        # Propagation stage.
        propagateWaveFront(iter, tunnel_rows, tunnel_cols, tsr, tsc, own_rows, own_cols, flow, flow_copy, p_locs)
    end

    # TODO: Stop global timer.
    MPI.Barrier(cart_comm)

    # TODO: Print results.
end

# Comment out call to main() whenever running tests.jl.
# main()