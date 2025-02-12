# Generate Candidate Trajectory Bundles (CTB) from TOSs

module CTB
using Dates
import ..TOS: Flight, add_trajectory_option!

function generateCtb(tos_list::Vector{Flight}, weights::Vector{Float64})
    # Solve MDO problem with weights
    # Return a CTB
end # module CTB


function generateCtbSet(tos_list::Vector{Flight})
    iterNum = 10
    ctbSet = []

    # Generate weights
    weights = rand(iterNum)
    weights = weights / sum(weights)

    # Generate CTBs
    for i = 1:iterNum
        # Generate a CTB
        ctb = generateCtb(tos_list, weights)
        push!(ctbSet, ctb)
    end

    # Return the set of CTBs
    return ctbSet
end
