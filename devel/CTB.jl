# Generate Candidate Trajectory Bundles (CTB) from TOSs

module CTB
using Dates
using JuMP, Gurobi, LinearAlgebra, Luxor
using MAT
import ..TOS: Flight, ParsedFlight, add_trajectory_option!, parse_flight_info
import ..GA: MI_LXPM

function GetParam()
    nm2km = 1.812
    
    conflictRadius = 15
    speed = 300 # km/h
    sampleRate = 5 # min (for every sampleRate [min], a trajectory point is generated)

    wptPos = read(matopen("./Wpt_pos.mat"))
    wptPos = wptPos["test"]
    fabInfo = read(matopen("./FIR_coord.mat"))
    fabInfo = fabInfo["data"]
    fabKeys = collect(keys(fabInfo))

    conflictRadius = conflictRadius * nm2km
    return (; conflictRadius, speed, sampleRate, wptPos, fabInfo, fabKeys)
end

function haversine(lat1, lon1, lat2, lon2)
    R = 6371.0  # Earth radius in km
    dlat = deg2rad(lat2 - lat1)
    dlon = deg2rad(lon2 - lon1)
    a = sin(dlat / 2)^2 + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dlon / 2)^2
    c = 2 * atan(sqrt(a), sqrt(1 - a))
    return R * c
end

function compute_heading(lat1, lon1, lat2, lon2)
    dlon = deg2rad(lon2 - lon1)
    x = cos(deg2rad(lat2)) * sin(dlon)
    y = cos(deg2rad(lat1)) * sin(deg2rad(lat2)) - sin(deg2rad(lat1)) * cos(deg2rad(lat2)) * cos(dlon)
    heading = rad2deg(atan(x, y))
    return mod(heading, 360)  # Ensure heading is between 0 and 360 degrees
end

function local_FAB_cost(X,V)
    C = V*pinv(X)
    A = C[:,1:end-1]
    ei = real(eigen(A).values)
    ei = filter(x -> x < 0, ei)
    if isempty(ei)
        H = 0
    else
        H = maximum(abs.(ei))
    end
    return H
end

function is_in_fir(flight::ParsedFlight, option_idx::Int, fabName::String, t::Int64)
    fir_records = flight.trajectory_options[option_idx].firTime_s
    for (name, t_in, t_out) in fir_records
        if name == fabName && t_in ≤ t ≤ t_out
            return true
        end
    end
    return false
end

function FAB_cost(x, tos_list, fabIdx, param)
    r = param.conflictRadius  # Conflict radius (15nm)
    v = param.speed  # Aircraft speed (km/h)
    d = param.sampleRate  # Sampling rate (min)
    wptPos = param.wptPos  # Waypoint positions (lat, lon)
    fabName = param.fabKeys[fabIdx]
    num_flights = length(tos_list)

    tosIdx = zeros(Int, num_flights)  # Selected trajectory indices
    startTime = Array{Int64}(undef, num_flights)  # Start time for each flight

    # Step 1: Identify selected TOS for each flight
    for i = 1:num_flights
        localVar = x[5*(i-1)+1:5*i] # For GA
        tosIdx[i] = findmax(localVar)[2]  # Extract selected TOS index
        startTime[i] = tos_list[i].trajectory_options[tosIdx[i]].valid_start_time  # Store start time
    end

    # Step 2: Synchronize Start Time & Define Unified Timeline
    globalStartTime = minimum(startTime)  # Earliest start time
    globalEndTime = maximum([startTime[i] + d * 500 * 60 for i in 1:num_flights])  # Approximate max time
    timeline = collect(globalStartTime:d*60:globalEndTime)  # Time grid

    # Step 3: Generate & Filter Spatiotemporal Trajectory Samples
    traj_samples = Dict{Int, Dict{Int64, Tuple{Float64, Float64, Float64}}}()  # (time => (lat, lon, heading))

    for i = 1:num_flights
        traj_samples[i] = Dict{Int64, Tuple{Float64, Float64, Float64}}()
        
        # Get trajectory waypoints for selected TOS
        waypoints = tos_list[i].trajectory_options[tosIdx[i]].route
        waypoints = parse.(Int, split(waypoints[2:end-1]))
        dep_time = startTime[i]  # Flight start time

        # Find the closest timeline point after the flight's start time
        start_idx = findfirst(t -> t >= dep_time, timeline)

        # Generate trajectory points from the closest timeline point
        time = timeline[start_idx]  # Set time to the closest point after dep_time
        for j = 1:length(waypoints)-1
            wp1 = wptPos[waypoints[j],:]  # (lat1, lon1)
            wp2 = wptPos[waypoints[j+1],:]  # (lat2, lon2)
            
            dist = haversine(wp1[1], wp1[2], wp2[1], wp2[2])  # Distance in km
            flight_time = dist / v * 60  # Time to next waypoint (min)
            heading = compute_heading(wp1[1], wp1[2], wp2[1], wp2[2])  # Compute heading
            
            # Sample points along trajectory every `d` minutes
            num_samples = ceil(Int, flight_time / d)
            for k in 0:num_samples
                alpha = k / num_samples
                lat = (1 - alpha) * wp1[1] + alpha * wp2[1]
                lon = (1 - alpha) * wp1[2] + alpha * wp2[2]
                if time ∈ timeline && is_in_fir(tos_list[i], tosIdx[i], fabName, time)
                    traj_samples[i][time] = (lat, lon, heading)
                end
                time += d * 60  # Increment time by d minutes
            end
        end
    end

    # Step 4: Compute Global FAB Cost Based on Conflict Detection
    globalCost = 0
    for t in timeline
        points_at_t = []  # List of points at current time step
        X = Float64[]
        X = reshape(X,3,0)
        V = Float64[]
        V = reshape(V,2,0)

        # Collect all points at time `t` that are inside FAB
        for i = 1:num_flights
            if haskey(traj_samples[i], t)
                lat, lon, heading = traj_samples[i][t]
                push!(points_at_t, (i, lon, lat, heading))
            end
        end

        # Step 5: Detect Conflicts & Compute FAB Cost for Each Time Step
        num_points = length(points_at_t)
        if num_points == 1 || num_points == 0
            continue
        else
            for a = 1:num_points
                i, lon1, lat1, heading1 = points_at_t[a]
                X = hcat(X, [lon1;lat1;1])
                V = hcat(V, [cos(heading1*pi/180);sin(heading1*pi/180)])
                for b = 1:num_points
                    j, lon2, lat2, heading2 = points_at_t[b]
                    dist = haversine(lat1, lon1, lat2, lon2)
                    if dist ≤ r  # If within conflict radius
                        X = hcat(X, [lon2;lat2;1])
                        V = hcat(V, [cos(heading2*pi/180);sin(heading2*pi/180)])
                    end
                end
            end
            local_cost = local_FAB_cost(X, V)  # Compute local FAB cost with direction
            globalCost += local_cost  # Add to total cost
        end
    end

    return globalCost
end

function SetMiscData(tos_list, fabIdx, param, w, flight_to_airline, valid_options, num_flights, idxs)
    global MiscTosList = tos_list
    global MiscFabIdx = fabIdx
    global MiscParam = param
    global MiscWeights = w
    global MiscFlightToAirline = flight_to_airline
    global MiscValidOptions = valid_options
    global MiscNumFlights = num_flights
    global MiscTotIdxs = idxs
end

function GetMiscData()
    return (; MiscTosList, MiscFabIdx, MiscParam, MiscWeights, MiscFlightToAirline, MiscValidOptions, MiscNumFlights, MiscTotIdxs)
end

function cost_function(x, coordFactor)
    misc = GetMiscData()
    idxs = misc.MiscTotIdxs
    fabIdx = misc.MiscFabIdx
    cost = 0
    for i = 1:length(idxs)
        if fabIdx != idxs[i] && coordFactor[i] > 0
            cost += computeCost(x,idxs[i]) * coordFactor[i]
        elseif fabIdx == idxs[i]
            cost += computeCost(x,idxs[i]) * (1 + coordFactor[i])
            # cost += computeCost(x,idxs[i])
        end
    end
    return cost
end

function computeCost(x, fabIdx)
    misc = GetMiscData()
    w_airlines = misc.MiscWeights[1] * 0.0
    w_FAB = misc.MiscWeights[2]
    flight_to_airline = misc.MiscFlightToAirline
    valid_options = misc.MiscValidOptions
    tos_list = misc.MiscTosList
    param = misc.MiscParam
    num_flights = misc.MiscNumFlights

    cost = 0
    for i in 1:num_flights
        for j in valid_options[i]
            cost += w_airlines[flight_to_airline[i]] * tos_list[i].trajectory_options[j].relative_trajectory_cost * x[5*(i-1)+j]
        end
    end
    cost += w_FAB * FAB_cost(x, tos_list,fabIdx,param)
    return cost
end

function solveCtb(tos_list::Vector{ParsedFlight}, weights::Vector{Float64}, fabIdx, param, idxs, coordFactor)
    num_flights = length(tos_list)
    
    # Identify unique airlines by extracting first two letters of flight_id
    airline_codes = unique(first(tos_list[i].flight_id, 2) for i in 1:num_flights)
    num_airlines = length(airline_codes)
    
    # Map flights to their respective airline index
    flight_to_airline = Dict(i => findfirst(==(first(tos_list[i].flight_id, 2)), airline_codes) 
                   for i in 1:num_flights)
    
    # Dictionary for valid options per flight
    valid_options = Dict(i => 1:length(tos_list[i].trajectory_options) for i in 1:num_flights)
    
    # Weight variables for airlines and FAB
    w_airlines = weights[1:num_airlines]
    # w_FAB = weights[num_airlines+1]
    w_FAB = 1
    w = [w_airlines, w_FAB]

    # Transfer data for cost function computation
    SetMiscData(tos_list, fabIdx, param, w, flight_to_airline, valid_options, num_flights, idxs)

    # Set bound
    bounds = []
    for i = 1:sum(length(valid_options[i]) for i in 1:num_flights)
        push!(bounds,(0,1))
    end

    # Solve!
    best_solution, best_cost = MI_LXPM(cost_function, bounds, 20, 50, 7, 0.8, 0.1, 0.3, 2.0, coordFactor)
    
    return best_solution, best_cost
end

function generateCtb(tos_list::Vector{ParsedFlight}, weights::Vector{Float64}, fabIdx, param, idxs, coordFactor)
    # Generate CTBs
    # ctb = solveCtb2(tos_list, weights, fabIdx, param)
    ctb = solveCtb(tos_list, weights, fabIdx, param, idxs, coordFactor)
    return ctb
end

function generateCtbSet(tos_list::Vector{ParsedFlight}, fabIdx, param, idxs, coordFactor)
    ctbSet = []
    num_flights = length(tos_list)

    # Identify unique airlines by extracting first two letters of flight_id
    airline_codes = unique(first(tos_list[i].flight_id, 2) for i in 1:num_flights)
    num_airlines = length(airline_codes)

    # Generate CTBs
    # for i = 1:num_airlines + 1
    for i in 1:1
        weights = ones(num_airlines + 1)
        weights = weights / sum(weights)
        # Acquire a CTB
        ctb = generateCtb(tos_list, weights, fabIdx, param, idxs, coordFactor)
        push!(ctbSet, ctb)
    end

    # Return the set of CTBs
    return ctbSet
end

function mat_parser(matTOS)
    flightNum = matTOS["flightNum"]
    origin = string(matTOS["options"][1][1])
    destination = string(matTOS["options"][1][end])
    earliest_departure = Int.(matTOS["depTime"])
    flight = Flight(flightNum,origin,destination,earliest_departure)

    num_options = length(matTOS["options"])
    for i = 1:num_options
        option = matTOS["options"][i]
        route = string(Int.(option))
        altitude = Int(matTOS["altitude"])
        speed = 300
        RTC = Int(round(matTOS["RTK"][i]*60))
        add_trajectory_option!(flight,route,altitude,speed,RTC,earliest_departure,nothing,nothing)
    end

    return flight
end

function mat_ctb_generation(fabIdx)
    flightSet = Vector{ParsedFlight}(undef,0)

    # Load mat files
    out=read(matopen("./TOS.mat"))
    TOS_set = out["TOS_set"]
    num_flights = length(TOS_set)
    for i = 1:num_flights
        for j = 1:5
            TOS_set[i]["options"][j] = Int.(TOS_set[i]["options"][j])
        end
        TOS_set[i]["RTK"] = Int.(TOS_set[i]["RTK"])
    end

    # Parse TOS set
    for i = 1:num_flights
        matTOS = TOS_set[i]
        flight = mat_parser(matTOS)

        firTime_all = matTOS["firTime"]
        timeRef = matTOS["depTime"]

        for j = 1:length(flight.trajectory_options)
            firLog = firTime_all[j]
            firTime_s = Vector{Tuple{String, Int64, Int64}}()

            for r = 1:size(firLog, 1)
                firName = string(firLog[r,1])
                t_in = Int64(round(firLog[r,2]))
                t_out = Int64(round(firLog[r,3]))
                push!(firTime_s, (firName, t_in, t_out))
            end

            flight.trajectory_options[j].firTime_s = firTime_s
            flight.trajectory_options[j].timeRef = timeRef
        end

        parsedFlight = parse_flight_info(flight)
        push!(flightSet, parsedFlight)
    end

    return flightSet
end

function get_all_ctbs(idxs, coordFactor)
    selectedRoute = Vector{Any}(undef,0)
    rawRoute = Vector{Any}(undef,0)
    for i = idxs
        println("Generating CTB for sector $i")
        flightSet = mat_ctb_generation(i)
        param = GetParam()
        # @time begin
        ctbSet = generateCtbSet(flightSet, i, param, idxs, coordFactor)
        # end
        raw = Int.(ctbSet[1][1])
        selected = findall(x -> x == 1, Int.(ctbSet[1][1]))
        tos_per_flight = 5
        temp = [mod(i - 1, tos_per_flight) + 1 for i in selected]
        push!(selectedRoute, temp)
        push!(rawRoute,raw)
    end
    return (; selectedRoute, rawRoute)
end

function exportCTB(result, idxs)
    param = GetParam()
    fabNames = param.fabKeys
    raws = result.rawRoute
    n = length(result.selectedRoute)
    scores = Array{Any}(undef,n,n)
    
    names = fabNames[idxs]
    selectedRoute = result.selectedRoute
    for i = 1:n
        for j = 1:n
            localIdx = idxs[i]
            scores[i, j] = computeCost(raws[j],localIdx)
        end
    end

    matwrite("negoTable.mat", Dict(
        "names" => names,
        "selectedRoute" => selectedRoute,
        "scores" => scores
    ); version="v7.4")

    return (; scores, selectedRoute)

end

end
