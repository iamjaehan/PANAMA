# Generates Trajectory Options Set(TOS) for a flight

module TOS
using Dates

# Define the Flight structure
mutable struct TrajectoryOption
    route::Any
    altitude::Any
    speed::Any
    relative_trajectory_cost::Any
    valid_start_time::Any
    valid_end_time::Any
    required_notice_time::Any
    firTime_s::Vector{Tuple{String, Int64, Int64}}
    timeRef::Int64
end

struct Flight
    flight_id::String
    origin::String
    destination::String
    earliest_departure::Int64
    trajectory_options::Vector{TrajectoryOption}
end

struct ParsedFlight
    flight_id::String
    origin::String
    destination::String
    earliest_departure::Int64
    trajectory_options::Vector{TrajectoryOption}
end

# Constructor for Flight
function Flight(flight_id::String, origin::String, destination::String, earliest_departure::Int64)
    return Flight(flight_id, origin, destination, earliest_departure, Vector{TrajectoryOption}())
end

function parse_flight_info(flight::Flight)
    return ParsedFlight(
        flight.flight_id,
        flight.origin,
        flight.destination,
        flight.earliest_departure,
        flight.trajectory_options  # 이미 TrajectoryOption 타입임
    )
end

# Add trajectory option to a flight
function add_trajectory_option!(
    flight::Flight,
    route::Any,
    altitude::Any,
    speed::Any,
    rtc::Any,
    valid_start::Any,
    valid_end::Any,
    required_notice::Any
)
    option = TrajectoryOption(
        route,
        altitude,
        speed,
        rtc,
        valid_start,
        valid_end,
        required_notice,
        Vector{Tuple{String, Int64, Int64}}(),
        0
    )
    push!(flight.trajectory_options, option)
end

end # module TOS_devel
