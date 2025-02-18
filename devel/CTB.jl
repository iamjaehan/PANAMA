# Generate Candidate Trajectory Bundles (CTB) from TOSs

module CTB
using Dates
using JuMP, Gurobi
import ..TOS: Flight, ParsedFlight, add_trajectory_option!, parse_flight_info

function FAB_cost(x)
    return sum(x)
end

function solveCtb(tos_list::Vector{ParsedFlight}, weights::Vector{Float64})
model = Model(Gurobi.Optimizer)
set_silent(model)

num_flights = length(tos_list)

# Identify unique airlines by extracting first two letters of flight_id
airline_codes = unique(first(tos_list[i].flight_id, 2) for i in 1:num_flights)
num_airlines = length(airline_codes)

# Map flights to their respective airline index
flight_to_airline = Dict(i => findfirst(==(first(tos_list[i].flight_id, 2)), airline_codes) 
               for i in 1:num_flights)

# Dictionary for valid options per flight
valid_options = Dict(i => 1:length(tos_list[i].trajectory_options) for i in 1:num_flights)

# Decision variables: x[i, j] is binary for valid options
@variable(model, x[i=1:num_flights, j in valid_options[i]], Bin)

# Weight variables for airlines and FAB
w_airlines = weights[1:num_airlines]
w_FAB = weights[num_airlines+1]

# Objective function: Sum of airline and FAB costs
@objective(model, Min, 
sum(w_airlines[flight_to_airline[i]] * tos_list[i].trajectory_options[j].relative_trajectory_cost * x[i, j] 
for i in 1:num_flights, j in valid_options[i]) + w_FAB * FAB_cost(x))

# Constraint: Each flight selects exactly one option
@constraint(model, [i=1:num_flights], sum(x[i, j] for j in valid_options[i]) == 1)

optimize!(model)

return value.(x)
end

function generateCtb(tos_list::Vector{ParsedFlight}, weights::Vector{Float64})
    # Generate CTBs
    ctb = solveCtb(tos_list, weights)
    return ctb
end

function generateCtbSet(tos_list::Vector{ParsedFlight})
    ctbSet = []
    num_flights = length(tos_list)

    # Identify unique airlines by extracting first two letters of flight_id
    airline_codes = unique(first(tos_list[i].flight_id, 2) for i in 1:num_flights)
    num_airlines = length(airline_codes)

    # Generate CTBs
    for i = 1:num_airlines + 1
        weights = rand(num_airlines + 1)
        weights = weights / sum(weights)
        # Acquire a CTB
        ctb = generateCtb(tos_list, weights)
        push!(ctbSet, ctb)
    end

    # Return the set of CTBs
    return ctbSet
end

function example_ctb_generation()
    # Create sample flights
    flight1 = Flight("AA123", "ORD", "LGA", DateTime(2023, 4, 25, 13, 0))
    add_trajectory_option!(flight1, "ORD..ELX..JHW..RKA..LGA", 350, 480, 0, DateTime(2023, 4, 25, 13, 0), nothing, nothing)
    add_trajectory_option!(flight1, "ORD..TVC..RKA..IGN..LGA", 370, 480, 10, DateTime(2023, 4, 25, 14, 0), nothing, nothing)
    flight1 = parse_flight_info(flight1)

    flight2 = Flight("OZ456", "DFW", "JFK", DateTime(2023, 4, 25, 13, 10))
    add_trajectory_option!(flight2, "DFW..ATL..JFK", 360, 500, 5, DateTime(2023, 4, 25, 13, 10), nothing, nothing)
    add_trajectory_option!(flight2, "DFW..MEM..JFK", 340, 480, 2, DateTime(2023, 4, 25, 13, 30), nothing, nothing)
    flight2 = parse_flight_info(flight2)

    # Generate CTBs
    ctbSet = generateCtbSet([flight1, flight2])
    return ctbSet
end

end