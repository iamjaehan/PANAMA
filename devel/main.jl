# Run simulation
using LinearAlgebra

function RunCTBGeneration(sectorIdxs, coordFactor)
    ctbSet = CTB.get_all_ctbs(sectorIdxs, coordFactor)
    out = CTB.exportCTB(ctbSet, sectorIdxs)
    return out
end

function RunNegotiation(out, assetReserve)
    taxParam = 100

    C = Matrix{Float64}(out.scores)
    d0 = Float64(1);
    gamma = Float64(0.9);
    epsilon = Float64(10);

    # Compute asset value
    b = taxParam/(assetReserve .+ 0.01); b = b';

    @time outcome, trade, profit, step = TACo.RunTACo(C, b, d0, gamma, epsilon)
    return (; outcome, trade, profit, step)
end

function ComputeShortFall(negoOut, assetReserve)
    trade = negoOut.trade
    n = size(trade)[1]
    shortFall = zeros(n)
    for i = 1:n
        shortFall[i] = max(0, trade[i] - assetReserve[i])
    end
    if sum(shortFall) == 0
        return zeros(n)
    end
    return normalize(shortFall,1)
end

function RunSimulation()
    sectorIdxs = [12, 3, 13]
    assetReserve = [20,20,20]
    roundLimit = 10

    n = length(sectorIdxs)
    coordFactor = zeros(n)
    negoOut = nothing
    shortFall = nothing
    debug = nothing

    for i = 1:roundLimit
        println("Round $i")
        # 1. Generate CTB
        out = RunCTBGeneration(sectorIdxs, coordFactor)
        println("CTB generation done")
        # 2. Negotiate
        debug = out
        negoOut = RunNegotiation(out, assetReserve)
        println("Negotiation done => Steps: $(negoOut.step)")
        # 3. Compute shortfall
        shortFall = ComputeShortFall(negoOut, assetReserve)
        println("Shortfall: $shortFall")
        if sum(shortFall) == 0
            println("Termination condition met. Terminating...")
            break
        end
        # 4. Update coordination factor
        println("Updating coordination factor => Prev: $coordFactor, New: $(coordFactor + shortFall)")
        coordFactor = coordFactor + shortFall
    end
    return (;negoOut, shortFall, debug)
end

out = RunSimulation()