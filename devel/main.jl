# Run simulation
using LinearAlgebra

function RunCTBGeneration(sectorIdxs, coordFactor)
    ctbSet = CTB.get_all_ctbs(sectorIdxs, coordFactor)
    out = CTB.exportCTB(ctbSet, sectorIdxs)
    return out
end

function RunNegotiation(out, assetReserve, taxParam)
    C = Matrix{Float64}(out.scores)
    d0 = Float64(1);
    gamma = Float64(0.9);
    epsilon = Float64(10);

    # Compute asset value
    b = taxParam./(assetReserve .+ 0.01);

    @time outcome, trade, profit, step, warnFlag = TACo.RunTACo(C, b, d0, gamma, epsilon)
    return (; outcome, trade, profit, step, warnFlag, C, b)
end

function ComputeShortFall(negoOut, assetReserve)
    trade = negoOut.trade
    n = size(trade)[1]
    shortFall = zeros(n)
    for i = 1:n
        shortFall[i] = max(0, trade[i] - assetReserve[i]) / assetReserve[i]
    end
    if sum(shortFall) == 0
        return (; rawShortFall = shortFall, shortFall = zeros(n))
    end
    return (; rawShortFall = shortFall, shortFall = normalize(shortFall,1))
end

function RunSimulation(assetReserve, taxParam)
    sectorIdxs = [12, 3, 13]
    roundLimit = 100

    n = length(sectorIdxs)
    coordFactor = zeros(n)
    negoOut = nothing
    shortFall = nothing
    debug = nothing
    rounds = nothing
    shortFall_history = Vector{Vector{Float64}}() # Store history of shortfalls
    negoOut_history = Vector{Any}() # Store history of negotiation outcomes

    for i = 1:roundLimit
        println("Round $i")
        # 1. Generate CTB
        out = RunCTBGeneration(sectorIdxs, coordFactor)
        println("CTB generation done")
        # 2. Negotiate
        negoOut = RunNegotiation(out, assetReserve, taxParam)
        println("Negotiation done => Steps: $(negoOut.step)")
        if negoOut.warnFlag
            debug = out
        end
        # 3. Compute shortfall
        shortFall_out = ComputeShortFall(negoOut, assetReserve)
        shortFall = shortFall_out.shortFall
        println("Shortfall: $(shortFall)")
        push!(shortFall_history, shortFall_out.rawShortFall)
        push!(negoOut_history, negoOut)
        if sum(shortFall) == 0
            println("Termination condition met. Terminating...")
            rounds = i
            break
        end
        # 4. Update coordination factor
        println("Updating coordination factor => Prev: $coordFactor, New: $(coordFactor + shortFall)")
        coordFactor = coordFactor + shortFall
        rounds = i
    end
    return (;negoOut, shortFall, debug, rounds, shortFall_history, negoOut_history)
end
