using Random
using MAT
include("main.jl")

# Define test grid
assetReserveCaseNum = 10
taxparamCaseNum = 10

# assetReserve_grid = [rand(1:100, 3) for _ in 1:assetReserveCaseNum]
# taxParam_grid = [rand(1:20) for _ in 1:taxparamCaseNum]
# repeats = 5

results = []

# for assetReserve in assetReserve_grid
#     for taxParam in taxParam_grid
#         for rep in 1:repeats
#             println("Running: assetReserve=$(assetReserve), taxParam=$(taxParam), rep=$rep")
#             out = RunSimulation(assetReserve, taxParam)
#             push!(results, Dict(
#                 "assetReserve" => assetReserve,
#                 "taxParam" => taxParam,
#                 "repeat" => rep,
#                 "rounds" => out.rounds,
#                 "shortFall" => out.shortFall,
#                 "debug" => out.debug,
#                 "negoOut" => out.negoOut,
#                 "shortFall_history" => out.shortFall_history,
#                 "negoOut_history" => out.negoOut_history
#             ))
#         end
#     end
# end

testCaseNum = 250
for testCase = 1:testCaseNum
    assetReserve = rand(1:50, 3)
    taxParam = rand(1:10) / 10 * 2
    println("Running: assetReserve=$(assetReserve), taxParam=$(taxParam), testIdx=$testCase")
    out = RunSimulation(assetReserve, taxParam)
    push!(results, Dict(
        "assetReserve" => assetReserve,
        "taxParam" => taxParam,
        "rounds" => out.rounds,
        "shortFall" => out.shortFall,
        "debug" => out.debug,
        "negoOut" => out.negoOut,
        "shortFall_history" => out.shortFall_history,
        "negoOut_history" => out.negoOut_history
    ))
end

# Prepare for saving: convert results to a savable format
# (MAT.jl cannot save arbitrary Julia objects, so we extract basic fields)
mat_results = Dict()
mat_results["assetReserve"] = [r["assetReserve"] for r in results]
mat_results["taxParam"] = [r["taxParam"] for r in results]
mat_results["repeat"] = [r["repeat"] for r in results]
mat_results["rounds"] = [r["rounds"] for r in results]
mat_results["shortFall"] = [r["shortFall"] for r in results]
mat_results["shortFall_history"] = [r["shortFall_history"] for r in results]
mat_results["negoOut_history"] = [r["negoOut_history"] for r in results]

matwrite("MC_test_results_randomSampling_250.mat", mat_results; version="v7.4")