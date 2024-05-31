using MPI
using Random


PRECISION = 10000
STEPS = 8
TSR = 0
TSC = 0


"""Particle that has position represented as [tunnel_row, tunnel_column] and speed as [row, column]."""
mutable struct Particle
    position::Vector{Int}
    mass::Int
    resistance::Int
    speed::Vector{Int}
    old_flow::Int
end


"""Returns the start of the tunnel row of the current rank in the cartesian grid."""
tunnel_start_row(tunnel_rows, cart_dims, cart_coord) = tunnel_rows ÷ cart_dims[2] * (cart_dims[2] - 1 - cart_coord[2]) + 1


"""Returns the start of the tunnel col of the current rank in the cartesian grid."""
tunnel_start_col(tunnel_cols, cart_dims, cart_coord) = tunnel_cols ÷ cart_dims[1] * cart_coord[1] + 1


"""Returns true if [r, c] is inside the current rank of the cartesian grid. Otherwise the return value is false."""
in_cartesian_segment(r, c, own_rows, own_cols) = (TSR <= r && r < TSR + own_rows && TSC <= c && c < TSC + own_cols) ? true : false


"""Reads values from the arguments and pushes them as values of a particle in the correct rank of the cartesian grid."""
function read_fixed_particles(lb, ub, particles, own_rows, own_cols)
    for i in lb:ub
        r = parse(Int, ARGS[13 + i * 3])
        c = parse(Int, ARGS[14 + i * 3])
        if in_cartesian_segment(r, c, own_rows, own_cols)
            mass = 0
            resistance = trunc(Int, parse(Float16, ARGS[15 + i * 3]) * PRECISION)
            speed = [0, 0]
            push!(particles, Particle([r, c] * PRECISION, mass, resistance, speed, 0))
        end
    end
end


"""Generates fixed particles in the tunnel and pushes them in the correct rank of the cartesian grid."""
function generate_fixed_particles(lb, ub, particles, own_rows, own_cols, fbp, fbs, tc)
    for i in lb:ub
        r = trunc(Int, fbp + fbs * rand(1)[1])
        c = trunc(Int, tc * rand(1)[1])
        if in_cartesian_segment(r, c, own_rows, own_cols)
            mass = 0
            resistance = trunc(Int, PRECISION * rand(1)[1])
            speed = [0, 0]
            push!(particles, Particle([r, c] * PRECISION, mass, resistance, speed, 0))
        end
    end
end


"""Generates moving particles in the tunnel and pushes them in the correct rank of the cartesian grid."""
function generate_moving_particles(lb, ub, particles, own_rows, own_cols, mbp, mbs, tc)
    for i in lb:ub
        r = trunc(Int, mbp + mbs * rand(1)[1])
        c = trunc(Int, tc * rand(1)[1])
        if in_cartesian_segment(r, c, own_rows, own_cols)
            mass = trunc(Int, PRECISION * (1 + 5 * rand(1)[1]))
            resistance = trunc(Int, PRECISION * rand(1)[1])
            speed = [0, 0]
            push!(particles, Particle([r, c] * PRECISION, mass, resistance, speed, 0))
        end
    end
end


"""Returns cartesian dimensions as [columns, rows] given the number of ranks."""
function get_cart_dims(n_ranks)
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


"""Returns a Dict of the cartesian neighbors (N, E, S, W, NE, NW, SE, SW) based on the topology of get_cart_dims."""
function get_cart_neighbors(comm, dims, coord)
    north = coord[2] + 1 < dims[2] ? coord[2] + 1 : -1
    south = coord[2] - 1 >= dims[2] - 2 ? coord[2] - 1 : -1
    east = coord[1] + 1 < dims[1] ? coord[1] + 1 : -1
    west = coord[1] - 1 >= dims[1] - 2 ? coord[1] - 1 : -1

    rank_north = north == -1 ? -1 : MPI.Cart_rank(comm, [coord[1], north])
    rank_south = south == -1 ? -1 : MPI.Cart_rank(comm, [coord[1], south])
    rank_east = east == -1 ? -1 : MPI.Cart_rank(comm, [east, coord[2]])
    rank_west = west == -1 ? -1 : MPI.Cart_rank(comm, [west, coord[2]])

    rank_ne = (north != -1 && east != -1) ? MPI.Cart_rank(comm, [east, north]) : -1 
    rank_nw = (north != -1 && west != -1) ? MPI.Cart_rank(comm, [west, north]) : -1 
    rank_se = (south != -1 && east != -1) ? MPI.Cart_rank(comm, [east, south]) : -1 
    rank_sw = (south != -1 && west != -1) ? MPI.Cart_rank(comm, [west, south]) : -1 

    Dict([("N", rank_north), ("S", rank_south), ("E", rank_east), ("W", rank_west), ("NE", rank_ne), ("NW", rank_nw), ("SE", rank_se), ("SW", rank_sw)])
end


"""Updates the first row (NOT the ghostcells) of the ranks at the top of the cartesian grid."""
function update_fan(iter, dims, coord, own_cols, fp, fs, flow)
    if dims[2] != coord[2] + 1 || iter % STEPS != 1; return end
    
    for c in 2:own_cols + 1
        if (fp + 1 < TSC + c) && (TSC + c <= fp + 1 + fs)
            phase = iter ÷ STEPS * (π ÷ 4)
            phase_step = π ÷ 2 ÷ fs
            pressure = 9 + 2 * sin(phase + (c - fp + 1) * phase_step)
            noise = 0.5 - rand(1)[1]
            flow[2, c] = trunc(Int, PRECISION * (pressure + noise))
        end
    end
end


"""Cleans the particle locations in p_locs by setting their value to 0 when either iter <= tunnel rows or on all tunnel rows."""
function clean_p_locs(iter, dims, coord, tr, own_rows, own_cols, p_locs)
    if iter <= tr
        ub = min(iter - own_rows * (dims[2] - 1 - coord[2]), own_rows) + 1
        for r in 2:ub, c in 2:own_cols + 1
            p_locs[r, c] = 0
        end
    else
        for r in 2:own_rows + 1, c in 2:own_cols + 1
            p_locs[r, c] = 0
        end
    end
end


"""Update the values of p_locs with the items in particles."""
function annotate_p_locs(particles, p_locs, own_rows, own_cols)
    for p in particles
        r = (p.position[1] ÷ PRECISION) % own_rows
        c = (p.position[2] ÷ PRECISION) % own_cols
        pr = r == 0 ? own_rows + 1 : r + 1
        pc = c == 0 ? own_cols + 1 : c + 1
        p_locs[pr, pc] += 1
    end
end


"""Function does a non-blocking send to distribute the flow to neighbors (S, E, W, SE, SW)."""
function send_flow(comm, neighbors, flow, own_rows, own_cols)
    reqs = MPI.Request[]

    for (k, v) in neighbors
        if k == "S" && v >= 0
            req = MPI.Isend(view(flow, own_rows + 1, 2:own_cols + 1), comm;dest=v, tag=0)
            push!(reqs, req)
        elseif k == "SE" && v >= 0
            req = MPI.Isend(view(flow, own_rows + 1, own_cols + 1), comm;dest=v, tag=0)
            push!(reqs, req)
        elseif k == "SW" && v >= 0
            req = MPI.Isend(view(flow, own_rows + 1, 2), comm;dest=v, tag=0)
            push!(reqs, req)
        elseif k == "E" && v >= 0
            req = MPI.Isend(view(flow, 2:own_rows + 1, own_cols + 1), comm;dest=v, tag=0)
            push!(reqs, req)
        elseif k == "W" && v >= 0
            req = MPI.Isend(view(flow, 2:own_rows + 1, 2), comm;dest=v, tag=0)
            push!(reqs, req)
        end
    end

    MPI.Waitall(reqs)
end


"""Function does a non-blocking recv to gather the flow from neighbors (N, E, W, NE, NW)."""
function recv_flow(comm, neighbors, flow, own_rows, own_cols)
    reqs = MPI.Request[]

    for (k, v) in neighbors
        if k == "N" && v >= 0
            req = MPI.Irecv!(view(flow, 1, 2:own_cols + 1), comm;source=v, tag=0)
            push!(reqs, req)
        elseif k == "NE" && v >= 0
            req = MPI.Irecv!(view(flow, 1, own_cols + 2), comm;source=v, tag=0)
            push!(reqs, req)
        elseif k == "NW" && v >= 0
            req = MPI.Irecv!(view(flow, 1, 1), comm;source=v, tag=0)
            push!(reqs, req)
        elseif k == "E" && v >= 0
            req = MPI.Irecv!(view(flow, 2:own_rows + 1, own_cols + 2), comm;source=v, tag=0)
            push!(reqs, req)
        elseif k == "W" && v >= 0
            req = MPI.Irecv!(view(flow, 2:own_rows + 1, 1), comm;source=v, tag=0)
            push!(reqs, req)
        end
    end

    MPI.Waitall(reqs)
end


""""""
function p_move(flow, p, tr, tc, own_rows, own_cols)
    for step in 0:STEPS
        
        r = p.position[1] ÷ PRECISION
        c = p.position[2] ÷ PRECISION == 0 ? 1 : p.position[2] ÷ PRECISION
    
        fr = r % own_rows == 0 ? own_rows + 1 : r % own_rows + 1
        fc = c % own_cols == 0 ? own_cols + 1 : c % own_cols + 1

        pressure = flow[fr - 1, fc]

        c == 1 ? left = 0 : left = pressure - flow[fr - 1, fc - 1]
        c == tc ? right = 0 : right = pressure - flow[fr - 1, fc + 1]

        flow_r = trunc(Int, pressure ÷ p.mass * PRECISION)
        flow_c = trunc(Int, (right - left) ÷ p.mass * PRECISION)

        # Speed change.
        p.speed[1] = (p.speed[1] + flow_r) ÷ 2
        p.speed[2] = (p.speed[2] + flow_c) ÷ 2

        # Movement.
        p.position[1] = p.position[1] + p.speed[1] ÷ STEPS ÷ 2
        p.position[2] = p.position[2] + p.speed[2] ÷ STEPS ÷ 2

        # Boundaries.
        if p.position[1] >= PRECISION * tr
            p.position[1] = PRECISION * tr
        end
        if p.position[2] <= 1
            p.position[2] = 1
        elseif p.position[2] >= PRECISION * tc
            p.position[2] = PRECISION * tc
        end
    end
end


""""""
function p_send(comm, neighbors, particles, own_rows, own_cols)
    reqs = MPI.Request[]
    send_SE::Vector{Int} = []
    send_SW::Vector{Int} = []
    send_S::Vector{Int} = []
    send_E::Vector{Int} = []
    send_W::Vector{Int} = []

    if length(particles) > 0
        for p in particles
            sp = [p.position[1], p.position[2], p.mass, p.resistance, p.speed[1], p.speed[2], p.old_flow]
            if TSR + own_rows - 1 < p.position[1]
                if TSC + own_cols - 1 < p.position[2]
                    append!(send_SE, sp)
                elseif p.position[2] < TSC
                    append!(send_SW, sp)
                else
                    append!(send_S, sp)
                end
            else
                if p.position[2] < TSC
                    append!(send_W, sp)
                else
                    append!(send_E, sp)
                end
            end
        end
    end

    for (k, v) in neighbors
        if k == "S" && v >= 0
            if length(send_S) == 0; send_S = [-1] end
            sr = MPI.Isend(send_S, comm;dest=v, tag=0)
            push!(reqs, sr)
        elseif k == "SE" && v >= 0
            if length(send_SE) == 0; send_SE = [-1] end
            sr = MPI.Isend(send_SE, comm;dest=v, tag=0)
            push!(reqs, sr)
        elseif k == "SW" && v >= 0
            if length(send_SW) == 0; send_SW = [-1] end
            sr = MPI.Isend(send_SW, comm;dest=v, tag=0)
            push!(reqs, sr)
        elseif k == "E" && v >= 0
            if length(send_E) == 0; send_E = [-1] end
            sr = MPI.Isend(send_E, comm;dest=v, tag=0)
            push!(reqs, sr)
        elseif k == "W" && v >= 0
            if length(send_W) == 0; send_W = [-1] end
            sr = MPI.Isend(send_W, comm;dest=v, tag=0)
            push!(reqs, sr)
        end
    end

    MPI.Waitall(reqs)

    while true
        length(particles) == 0 ? break : deleteat!(particles, 1)
    end
end


""""""
function p_recv_neighbor(comm, src, total_recv)
    status = MPI.Probe(comm, MPI.Status;source=src, tag=0)
    num = MPI.Get_count(status, Int)
    rb = Array{Int}(undef, num)
    MPI.Recv!(rb, comm;source=src, tag=0)
    if num > 1
        append!(total_recv, rb)
    end
end


""""""
function p_recv(comm, neihgbors, particles)
    total_recv::Vector{Int} = []

    for (k, v) in neihgbors
        if k == "N" && v >= 0
            p_recv_neighbor(comm, v, total_recv)
        elseif k == "NE" && v >= 0
            p_recv_neighbor(comm, v, total_recv)
        elseif k == "NW" && v >= 0
            p_recv_neighbor(comm, v, total_recv)
        elseif k == "E" && v >= 0
            p_recv_neighbor(comm, v, total_recv)
        elseif k == "W" && v >= 0
            p_recv_neighbor(comm, v, total_recv)
        end
    end

    p_elements = 7
    partitioned = collect(Iterators.partition(total_recv, p_elements))
    for i in partitioned
        push!(particles, Particle([i[1], i[2]], i[3], i[4], [i[5], i[6]], i[7]))
    end
end


""""""
function p_movements(iter, comm, dims, coord, neighbors, tr, tc, own_rows, own_cols, flow, p_locs, incoming_particles, outgoing_particles)
    if iter % STEPS != 1; return end

    clean_p_locs(iter, dims, coord, tr, own_rows, own_cols, p_locs)

    send_flow(comm, neighbors, flow, own_rows, own_cols)
    recv_flow(comm, neighbors, flow, own_rows, own_cols)

    i = 1
    del = 0
    while i <= length(incoming_particles)
        if incoming_particles[i].mass != 0
            for j in 1:STEPS
                del = p_move(flow, incoming_particles[i], tr, tc, own_rows, own_cols)
                if del == 1; break end
            end

            if del == 1
                push!(outgoing_particles, incoming_particles[i])
                deleteat!(incoming_particles, i)
                del = 0
            end
        end
        i += 1
    end

    p_send(comm, neighbors, outgoing_particles, own_rows, own_cols)
    p_recv(comm, neighbors, incoming_particles)

    annotate_p_locs(incoming_particles, p_locs, own_rows, own_cols)
end


""""""
function update_flow(flow, flow_copy, p_locs, r, c, fr, fc, tc, skip)
    if skip && p_locs[fr, fc] != 0 || r == 1; return 0 end

    # update border left.
    if c == 1
        flow[fr, fc] = (flow_copy[fr, fc] + flow_copy[fr - 1, fc] + flow_copy[fr - 1, fc + 1]) ÷ 4
    end
    # Update border right.
    if c == tc
        flow[fr, fc] = (flow_copy[fr, fc] + flow_copy[fr - 1, fc] + flow_copy[fr - 1, fc - 1]) ÷ 4
    end
    # Update central part.
    if 1 < c && c < tc
        flow[fr, fc] = (flow_copy[fr, fc] + flow_copy[fr - 1, fc] + flow_copy[fr - 1, fc - 1] + flow_copy[fr - 1, fc + 1]) ÷ 5
    end
    
    return abs(flow_copy[fr, fc] - flow[fr, fc])
end


""""""
function p_effects(iter, particles, tc, own_rows, own_cols, flow, flow_copy, p_locs)
    if iter % STEPS != 1; return end

    for p in particles
        r = p.position[1] ÷ PRECISION
        c = p.position[2] ÷ PRECISION == 0 ? 1 : p.position[2] ÷ PRECISION
        fr = r % own_rows == 0 ? own_rows + 1 : r % own_rows + 1
        fc = c % own_cols == 0 ? own_cols + 1 : c % own_cols + 1
        update_flow(flow, flow_copy, p_locs, r, c, fr, fc, tc, false)
        p.old_flow = flow[fr, fc]
    end

    for p in particles
        r = p.position[1] ÷ PRECISION
        c = p.position[2] ÷ PRECISION == 0 ? 1 : p.position[2] ÷ PRECISION
        fr = r % own_rows == 0 ? own_rows + 1 : r % own_rows + 1
        fc = c % own_cols == 0 ? own_cols + 1 : c % own_cols + 1

        back = trunc(Int, (p.old_flow * p.resistance ÷ PRECISION) ÷ p_locs[fr, fc])
        flow[fr, fc] -= back
        flow[fr - 1, fc] += back ÷ 2

        c > 1 ? flow[fr - 1, fc - 1] += back ÷ 4 : flow[fr - 1, fc] += back ÷ 4
        c < tc ? flow[fr - 1, fc + 1] += back ÷ 4 : flow[fr - 1, fc] += back ÷ 4
    end    
end


""""""
function p_wavefront(iter, max_var, tr, tc, own_rows, own_cols, flow, flow_copy, p_locs)
    wavefront = iter % STEPS
    
    if wavefront == 1
        max_var = 0
    elseif wavefront == 0
        wavefront = STEPS
    end

    own_wavefronts::Vector{Int} = []
    for w in wavefront:tr
        if in_cartesian_segment(w, TSC, own_rows, own_cols)
            push!(own_wavefronts, w)
        end
    end

    for w in own_wavefronts
        if w > iter; break end

        for c in TSC:TSC + own_cols - 1
            fr = w % own_rows == 0 ? own_rows + 1 : w % own_rows + 1
            fc = c % own_cols == 0 ? own_cols + 1 : c % own_cols + 1
            var = update_flow(flow, flow_copy, p_locs, w, c, fr, fc, tc, true)
            if var > max_var; max_var = var end
        end
    end
end


# TODO: Only works for n_ranks 2 and 4!
"""Returns a matrix made out of matrices m ordered depending on the cartesian topology of get_cart_dims."""
function combine_matrices(n_ranks, m, range_r, range_c)
    if n_ranks == 1
        return [view(m[1], range_r, range_c)]
    elseif n_ranks == 2
        return [view(m[1], range_r, range_c) view(m[2], range_r, range_c)]
    elseif n_ranks == 4
        return [
            view(m[2], range_r, range_c) view(m[4], range_r, range_c);
            view(m[1], range_r, range_c) view(m[3], range_r, range_c)
        ]
    elseif n_ranks == 8
        return [
            view(m[4], range_r, range_c) view(m[8], range_r, range_c);
            view(m[3], range_r, range_c) view(m[7], range_r, range_c);
            view(m[2], range_r, range_c) view(m[6], range_r, range_c);
            view(m[1], range_r, range_c) view(m[5], range_r, range_c)
        ]
    elseif n_ranks == 16
        return [
            view(m[4], range_r, range_c) view(m[8], range_r, range_c) view(m[12], range_r, range_c) view(m[16], range_r, range_c);
            view(m[3], range_r, range_c) view(m[7], range_r, range_c) view(m[11], range_r, range_c) view(m[15], range_r, range_c);
            view(m[2], range_r, range_c) view(m[6], range_r, range_c) view(m[10], range_r, range_c) view(m[14], range_r, range_c);
            view(m[1], range_r, range_c) view(m[5], range_r, range_c) view(m[9], range_r, range_c) view(m[13], range_r, range_c)
        ]
    end
end


function main()
    if length(ARGS) < 15; println(" -> Not enough arugments!") end

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
    cart_dims = get_cart_dims(n_ranks)
    periodic = map(_ -> false, cart_dims)
    reorder = false
    cart_comm = MPI.Cart_create(comm, cart_dims; periodic, reorder)
    cart_coord = MPI.Cart_coords(cart_comm)
    cart_neighbors = get_cart_neighbors(cart_comm, cart_dims, cart_coord)

    # Initialize own rows and cols and TSR and TSC.
    own_cols = tunnel_cols ÷ cart_dims[1]
    own_rows = tunnel_rows ÷ cart_dims[2]
    global TSR = tunnel_start_row(tunnel_rows, cart_dims, cart_coord)
    global TSC = tunnel_start_col(tunnel_cols, cart_dims, cart_coord)

    # Read and generate particles.
    remaining_args = (length(ARGS) - 15) ÷ 3
    n_fixed_band = trunc(Int, fixed_band_size * tunnel_cols * fixed_density)
    n_moving_band = trunc(Int, moving_band_size * tunnel_cols * moving_density)
    n_particles = remaining_args + n_fixed_band + n_moving_band
    incoming_particles::Vector{Particle} = []
    outgoing_particles::Vector{Particle} = []

    read_fixed_particles(1, remaining_args, incoming_particles, own_rows, own_cols)
    ub = n_particles - remaining_args - n_moving_band
    generate_fixed_particles(1, ub, incoming_particles, own_rows, own_cols, fixed_band_pos, fixed_band_size, tunnel_cols)
    ub = n_particles - remaining_args - n_fixed_band
    generate_moving_particles(1, ub, incoming_particles, own_rows, own_cols, moving_band_pos, moving_band_size, tunnel_cols)

    max_var = typemax(Int64)

    flow = zeros(Int, own_rows + 2, own_cols + 2)
    flow[1, :] = flow[end, :] .= -1; flow[:, 1] = flow[:, end] .= -1
    flow_copy = copy(flow)
    p_locs = copy(flow)

    start_iter = 0

    for iter in start_iter:max_iter

        if max_var <= threshold
            break
        end

        update_fan(iter, cart_dims, cart_coord, own_cols, fan_pos, fan_size, flow)

        p_movements(iter, cart_comm, cart_dims, cart_coord, cart_neighbors, tunnel_rows, tunnel_cols, own_rows, own_cols, flow, p_locs, incoming_particles, outgoing_particles)

        p_effects(iter, incoming_particles, tunnel_cols, own_rows, own_cols, flow, flow_copy, p_locs)
    
        copy!(flow_copy, flow)

        p_wavefront(iter, max_var, tunnel_rows, tunnel_cols, own_rows, own_cols, flow, flow_copy, p_locs)
    end

    recv_flow = MPI.Gather(flow, cart_comm;root=0)
    recv_p_locs = MPI.Gather(p_locs, cart_comm;root=0)
    if rank == 0
        range_r = 2:own_rows + 1
        range_c = 2:own_cols + 1

        res_flow = zeros(Int, tunnel_rows, tunnel_cols)
        partitioned = collect(Iterators.partition(recv_flow, (own_rows + 2) * (own_cols + 2)))
        matrices = [reshape(partitioned[i], own_rows + 2, own_cols + 2) for i in 1:n_ranks]
        res_flow .= combine_matrices(n_ranks, matrices, range_r, range_c)

        res_p_locs = zeros(Int, tunnel_rows, tunnel_cols)
        partitioned = collect(Iterators.partition(recv_p_locs, (own_rows + 2) * (own_cols + 2)))
        matrices = [reshape(partitioned[i], own_rows + 2, own_cols + 2) for i in 1:n_ranks]
        res_p_locs .= combine_matrices(n_ranks, matrices, range_r, range_c)

        result = Matrix{String}(undef, tunnel_rows, tunnel_cols)
        for r in 1:tunnel_rows, c in 1:tunnel_cols
            if res_flow[r, c] >= 10 * PRECISION
                symbol = '*'
            elseif res_flow[r, c] >= 1 * PRECISION
                symbol = '0' + res_flow[r, c] ÷ PRECISION;
            elseif res_flow[r, c] >= 0.5 * PRECISION
                symbol = '+'
            elseif res_flow[r, c] > 0
                symbol = '.'
            else
                symbol = ' '
            end
        
            res_p_locs[r, c] > 0 ? result[r, c] = "[$symbol]" : result[r, c] = string(symbol)
        end

        display(result)

        res = ""
        for r in 1:tunnel_rows
            for c in 1:tunnel_cols
                length(result[r, c]) > 1 ? res *= " $(result[r,c]) " : res *= "  $(result[r,c])  "
            end
            res *= "\n"
        end
    
        write("restult-mpi.txt", res)
    end
end

main()