# Baseline CTOP algorithm.

module CTOP
using Dates
import ..TOS: Flight, add_trajectory_option!

# Function to process TOS from multiple airlines and assign optimal trajectories
function process_ctop!(tos_list::Vector{Flight})
    # Group trajectories by flight_id directly from Flight structs
    grouped_trajectories = Dict{String, Vector{Dict{String, Any}}}()

    for flight in tos_list
        flight_id = flight.flight_id
        trajectory_options = flight.trajectory_options

        if haskey(grouped_trajectories, flight_id)
            append!(grouped_trajectories[flight_id], trajectory_options)
        else
            grouped_trajectories[flight_id] = trajectory_options
        end
    end

    # Assign trajectories using traditional CTOP (sort by earliest valid start time and lowest RTC)
    assigned_flights = []
    for (flight_id, options) in grouped_trajectories
        # Sorting based on 'valid_start_time' and 'relative_trajectory_cost'
        sorted_options = sort(options, by = x -> (x["valid_start_time"], x["relative_trajectory_cost"]))
        assigned_trajectory = first(sorted_options)
        
        push!(assigned_flights, Dict(
            "flight_id" => flight_id,
            "assigned_trajectory" => assigned_trajectory
        ))
    end

    # Display assigned trajectories
    for flight in assigned_flights
        println("Assigned Flight: ", flight["flight_id"])
        println("  Route: ", flight["assigned_trajectory"]["route"])
        println("  Valid Start Time: ", flight["assigned_trajectory"]["valid_start_time"])
        println("  Relative Trajectory Cost (RTC): ", flight["assigned_trajectory"]["relative_trajectory_cost"])
        println()
    end

    return assigned_flights
end

# Example usage with simulated TOS from multiple airlines
function example_tos_processing()
    # Create sample flights
    flight1 = Flight("FL123", "ORD", "LGA", DateTime(2023, 4, 25, 13, 0))
    add_trajectory_option!(flight1, "ORD..ELX..JHW..RKA..LGA", 350, 480, 0, DateTime(2023, 4, 25, 13, 0), nothing, nothing)
    add_trajectory_option!(flight1, "ORD..TVC..RKA..IGN..LGA", 370, 480, 10, DateTime(2023, 4, 25, 14, 0), nothing, nothing)

    flight2 = Flight("FL456", "DFW", "JFK", DateTime(2023, 4, 25, 13, 10))
    add_trajectory_option!(flight2, "DFW..ATL..JFK", 360, 500, 5, DateTime(2023, 4, 25, 13, 10), nothing, nothing)
    add_trajectory_option!(flight2, "DFW..MEM..JFK", 340, 480, 2, DateTime(2023, 4, 25, 13, 30), nothing, nothing)

    # Process and assign TOS
    process_ctop!([flight1, flight2])
end

end # module CTOP
