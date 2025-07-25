module TACo
using StatsBase
using LinearAlgebra

function RunTACo(C::Matrix{Float64}, b::Vector{Float64}, d0::Float64, gamma::Float64, epsilon::Float64)
    n, m = size(C)
    O = zeros(n, m)             # Offer matrix
    P = zeros(n, m)             # Pay matrix
    d = d0                      # Trading unit
    isConverged = false         # Termination flag
    selections = zeros(Int, n) # Agent choices
    recordedStates = Set{Vector{Float64}}()
    agentOrder = collect(1:n)  # Sequential order
    step = 0
    warnFlag = false

    while !isConverged
        step += 1
        if step > 1e4 * 5
            println("Warning: TACo is taking too long. Terminating...")
            warnFlag = true
            break
        end
        i = agentOrder[mod1(step, n)]  # Julia uses 1-based indexing

        # Profit computation
        J = b[i] .* (O[i, :] .- P[i, :]) .- C[i, :]
        j_star = argmax(J)

        # Update
        P[i, j_star] += n * d
        O[:, j_star] .+= d
        selections[i] = j_star

        # Display formatted output
        # println(rpad("  $step", 6), "|", rpad("  $i", 6), "| ",
        #     "[", join(O[i, :], " "), "] | ",
        #     "[", join(P[i, :], " "), "] | ",
        #     "[", join(round.(J, digits=1), " "), "] | ",
        #     "[", join(selections, " "), "]")

        # Cycle detection
        key = vcat(vec(O - P), [i])
        key_found = any(x -> x == key, recordedStates)

        if all(x -> x == selections[1], selections)
            break
        end

        if key_found
            d *= gamma
            empty!(recordedStates)

            # Epsilon termination check
            allSatisfied = true
            for k in 1:n
                Jk = b[k] .* (O[k, :] .- P[k, :]) .- C[k, :]
                if maximum(Jk) - minimum(Jk) > epsilon
                    allSatisfied = false
                    break
                end
            end

            if allSatisfied
                isConverged = true
            end
        end

        # Record current state
        push!(recordedStates, key)
    end

    outcome = mode(selections)
    trade = O - P
    profit = diagm(b) * (O - P) - C

    return (;outcome, trade, profit, step, warnFlag)
end

function TestTACo()
    C = Matrix{Float64}([1 2 3; 8 5 4; 9 7 8])
    b = Vector{Float64}([2, 2, 2])
    d0 = Float64(1)
    gamma = Float64(0.9)
    epsilon = Float64(0.1)
    outcome, trade, profit = RunTACo(C, b, d0, gamma, epsilon)
    println("Outcome: $outcome")
    println("Trade: $trade")
end

end