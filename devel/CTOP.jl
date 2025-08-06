# Baseline CTOP algorithm.

module CTOP
using Dates
import ..TOS: Flight, ParsedFlight, add_trajectory_option!, parse_flight_info
using JuMP, Gurobi, LinearAlgebra, Luxor
using MAT
import ..CTB
import ..GA: MI_LXPM

function cost_function(x, coordFactor)
    misc = CTB.GetMiscData()
    idxs = misc.MiscTotIdxs
    totalCost = 0
    for i = 1:length(idxs)
        totalCost += CTB.computeCost(x,idxs[i])
    end
    return totalCost
end

function solveCtb(tos_list::Vector{ParsedFlight}, weights::Vector{Float64}, fabIdx, param, idxs)
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
    CTB.SetMiscData(tos_list, fabIdx, param, w, flight_to_airline, valid_options, num_flights, idxs)

    # Set bound
    bounds = []
    for i = 1:sum(length(valid_options[i]) for i in 1:num_flights)
        push!(bounds,(0,1))
    end
    coordFactor = zeros(length(idxs))

    # Solve!
    best_solution, best_cost = MI_LXPM(CTOP.cost_function, bounds, 20, 50, 7, 0.8, 0.1, 0.3, 2.0, coordFactor)
    
    return best_solution, best_cost
end

function RunCTOP(idxs)
    # 1. Load all flights (as ParsedFlight)
    flightSet = CTB.mat_ctb_generation(0)

    # 2. Prepare CTOP parameters
    param = CTB.GetParam()  # Shared utility
    num_flights = length(flightSet)

    # Extract all FAB indices (assume complete for centralized eval)
    fabIdx = 0  # Can choose 0 or arbitrary FAB to evaluate global cost

    # 3. Set uniform weights (or fixed fairness priorities)
    airline_codes = unique(first(flightSet[i].flight_id, 2) for i in 1:num_flights)
    num_airlines = length(airline_codes)
    weights = ones(num_airlines + 1)
    weights = weights / sum(weights)  # Normalize

    # 5. Solve the centralized CTB optimization
    raw, cost = CTOP.solveCtb(flightSet, weights, fabIdx, param, idxs)

    # Optional: convert x into TOS indices per flight
    selected_indices = []
    for i in 1:num_flights
        local_var = raw[5*(i-1)+1:5*i]
        push!(selected_indices, findmax(local_var)[2])
    end

    indCost = []
    for i = 1:length(idxs)
        push!(indCost,CTB.computeCost(raw,idxs[i]))
    end

    return (; raw, selected_indices, cost, indCost)
end


###

function fcfs_ctop_solver(tos_list::Vector{ParsedFlight}, idxs, param)
    num_flights = length(tos_list)
    tos_per_flight = length(tos_list[1].trajectory_options)
    x = zeros(Int, num_flights * tos_per_flight)  # binary vector

    # 정렬된 항공편 인덱스
    sorted_idx = sortperm([f.earliest_departure for f in tos_list])
    
    # Identify unique airlines by extracting first two letters of flight_id
    airline_codes = unique(first(tos_list[i].flight_id, 2) for i in 1:num_flights)
    num_airlines = length(airline_codes)
    
    # Map flights to their respective airline index
    flight_to_airline = Dict(i => findfirst(==(first(tos_list[i].flight_id, 2)), airline_codes) 
                   for i in 1:num_flights)
    
    # Dictionary for valid options per flight
    valid_options = Dict(i => 1:length(tos_list[i].trajectory_options) for i in 1:num_flights)
    
    # Weight variables for airlines and FAB
    weights = zeros(num_airlines+1,1)
    w_airlines = weights[1:num_airlines]
    # w_FAB = weights[num_airlines+1]
    w_FAB = 1
    w = [w_airlines, w_FAB]
    fabIdx = 0
    
    CTB.SetMiscData(tos_list, fabIdx, param, w, flight_to_airline, valid_options, num_flights, idxs)

    for count in 1:num_flights
        i = sorted_idx[count]
        best_cost = Inf
        best_tos = 0

        for j in 1:tos_per_flight
            temp_x = copy(x)
            temp_x[tos_per_flight*(i-1)+j] = 1

            cost = 0
            for fabIdx in idxs
                cost += CTB.computeCost(temp_x, fabIdx)
            end

            if cost < best_cost
                best_cost = cost
                best_tos = j
            end
        end

        # 최적 TOS를 확정
        x[tos_per_flight*(i-1)+best_tos] = 1
    end

    return x
end

function RunCTOPFCFS(idxs)
    # 1. Load all flights (as ParsedFlight)
    flightSet = CTB.mat_ctb_generation(0)

    # 2. Prepare CTOP parameters
    param = CTB.GetParam()  # Shared utility
    num_flights = length(flightSet)

    # Extract all FAB indices (assume complete for centralized eval)
    fabIdx = 0  # Can choose 0 or arbitrary FAB to evaluate global cost

    # 3. Set uniform weights (or fixed fairness priorities)
    airline_codes = unique(first(flightSet[i].flight_id, 2) for i in 1:num_flights)
    num_airlines = length(airline_codes)
    weights = ones(num_airlines + 1)
    weights = weights / sum(weights)  # Normalize

    #4. FCFS greedy algorithm
    raw = fcfs_ctop_solver(flightSet::Vector{ParsedFlight}, idxs, param)

    #5. Prepare output data
    selected_indices = []
    for i in 1:num_flights
        local_var = raw[5*(i-1)+1:5*i]
        push!(selected_indices, findmax(local_var)[2])
    end

    indCost = []
    for i = 1:length(idxs)
        push!(indCost,CTB.computeCost(raw,idxs[i]))
    end

    cost = sum(indCost)

    return (; raw, selected_indices, cost, indCost)
end

end # module CTOP
