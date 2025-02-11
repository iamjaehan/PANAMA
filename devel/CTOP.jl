module CTOP

using Dates

# Function to process TOS from multiple airlines and assign optimal trajectories
function process_tos!(tos_list::Vector{Dict{String, Any}})
    # Flatten all trajectory options with flight metadata
    all_trajectories = []
    for tos in tos_list
        flight_id = tos["flight_id"]
        for option in tos["trajectory_options"]
            push!(all_trajectories, merge(option, Dict("flight_id" => flight_id)))
        end
    end

    # Group trajectories by flight_id
    grouped_trajectories = Dict{String, Vector{Dict{String, Any}}}()
    for traj in all_trajectories
        flight_id = traj["flight_id"]
        if haskey(grouped_trajectories, flight_id)
            push!(grouped_trajectories[flight_id], traj)
        else
            grouped_trajectories[flight_id] = [traj]
        end
    end

    # Assign trajectories using traditional CTOP (sort by earliest arrival and lowest RTC)
    assigned_flights = []
    for (flight_id, options) in grouped_trajectories
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
end

# Example usage with simulated TOS from multiple airlines
function example_tos_processing()
    tos_list = [
        Dict(
            "flight_id" => "FL123",
            "trajectory_options" => [
                Dict("route" => "ORD..ELX..JHW..RKA..LGA", "altitude" => 350, "speed" => 480, "relative_trajectory_cost" => 0, "valid_start_time" => DateTime(2023, 4, 25, 13, 0)),
                Dict("route" => "ORD..TVC..RKA..IGN..LGA", "altitude" => 370, "speed" => 480, "relative_trajectory_cost" => 10, "valid_start_time" => DateTime(2023, 4, 25, 14, 0))
            ]
        ),
        Dict(
            "flight_id" => "FL456",
            "trajectory_options" => [
                Dict("route" => "DFW..ATL..JFK", "altitude" => 360, "speed" => 500, "relative_trajectory_cost" => 5, "valid_start_time" => DateTime(2023, 4, 25, 13, 10)),
                Dict("route" => "DFW..MEM..JFK", "altitude" => 340, "speed" => 480, "relative_trajectory_cost" => 2, "valid_start_time" => DateTime(2023, 4, 25, 13, 30))
            ]
        ),
        Dict(
            "flight_id" => "FL789",
            "trajectory_options" => [
                Dict("route" => "LAX..DEN..BOS", "altitude" => 370, "speed" => 490, "relative_trajectory_cost" => 0, "valid_start_time" => DateTime(2023, 4, 25, 13, 5)),
                Dict("route" => "LAX..ORD..BOS", "altitude" => 350, "speed" => 480, "relative_trajectory_cost" => 3, "valid_start_time" => DateTime(2023, 4, 25, 13, 25))
            ]
        )
    ]

    process_tos!(tos_list)
end

end # module CTOP
