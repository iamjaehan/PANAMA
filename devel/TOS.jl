# Generates Trajectory Options Set(TOS) for a flight

module TOS
using Dates

# Define the Flight structure
struct Flight
    flight_id::String
    origin::String
    destination::String
    earliest_departure::DateTime
    trajectory_options::Vector{Dict{String, Any}}
end

struct TrajectoryOption
    route::Any
    altitude::Any
    speed::Any
    relative_trajectory_cost::Any
    valid_start_time::Any
    valid_end_time::Any
    required_notice_time::Any
end

struct ParsedFlight
    flight_id::String
    origin::String
    destination::String
    earliest_departure::DateTime
    trajectory_options::Vector{TrajectoryOption}
end

# Constructor for Flight
function Flight(flight_id::String, origin::String, destination::String, earliest_departure::DateTime)
    return Flight(flight_id, origin, destination, earliest_departure, Vector{Dict{String, Any}}())
end

function parse_flight_info(flight::Flight)
    # Convert trajectory_options into TrajectoryOption structs
    trajectory_options = [TrajectoryOption(
        option["route"],
        option["altitude"],
        option["speed"],
        option["relative_trajectory_cost"],
        option["valid_start_time"],
        option["valid_end_time"],
        option["required_notice_time"]
    ) for option in flight.trajectory_options]
    
    return ParsedFlight(
        flight.flight_id,
        flight.origin,
        flight.destination,
        flight.earliest_departure,
        trajectory_options
    )
end

# Add trajectory option to a flight
function add_trajectory_option!(flight::Flight, route::String, altitude::Int, speed::Int, rtc::Int, valid_start::DateTime, valid_end::Union{DateTime, Nothing}, required_notice::Union{Int, Nothing})
    option = Dict(
        "route" => route,
        "altitude" => altitude,
        "speed" => speed,
        "relative_trajectory_cost" => rtc,
        "valid_start_time" => valid_start,
        "valid_end_time" => valid_end,
        "required_notice_time" => required_notice
    )
    push!(flight.trajectory_options, option)
end

# Simulate data input
function generate_flight_data()
    flight = Flight("FL123", "ORD", "LGA", DateTime(2023, 4, 25, 13, 0))

    # Example trajectory options
    add_trajectory_option!(flight, "ORD..ELX..JHW..RKA..LGA", 350, 480, 0, DateTime(2023, 4, 25, 13, 0), nothing, nothing)
    add_trajectory_option!(flight, "ORD..ELX..JHW..RKA..LGA", 370, 480, 10, DateTime(2023, 4, 25, 13, 0), DateTime(2023, 4, 25, 15, 0), nothing)
    add_trajectory_option!(flight, "ORD..TVC..RKA..IGN..LGA", 350, 480, 20, DateTime(2023, 4, 25, 13, 0), DateTime(2023, 4, 25, 16, 0), 45)

    # parse flight
    parsed_flight = parse_flight_info(flight)

    return flight, parsed_flight
end

# Evaluate costs for trajectory options
function evaluate_trajectory_costs!(flight::Flight)
    for option in flight.trajectory_options
        option["total_cost"] = option["relative_trajectory_cost"]  # Extend with other factors if needed
    end
end

# Generate and print TOS
function generate_tos(flight::Flight)
    println("Submitting TOS for flight: ", flight.flight_id)
    for option in flight.trajectory_options
        println(option)
    end
end

# Main workflow to run TOS generation
function tos_workflow()
    flight_data = generate_flight_data()
    evaluate_trajectory_costs!(flight_data)
    generate_tos(flight_data)
end

end # module TOS_devel
