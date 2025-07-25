module GA
using Random, Distributions

# Define the MI-LXPM algorithm
function MI_LXPM(cost_function, bounds, num_generations, population_size, tournament_size, crossover_prob, mutation_prob, alpha, beta, coordFactor)
    num_variables = length(bounds)

    # Initialize population (real-coded values within bounds)
    population = [[rand(bounds[j][1]:0.01:bounds[j][2]) for j in 1:num_variables] for _ in 1:population_size]

    for generation in 1:num_generations
        fitness = Vector{Float64}(undef, length(population))
        Threads.@threads for i in eachindex(population)
            x = round.(population[i])
            enforce_onehot!(x)
            fitness[i] = cost_function(x, coordFactor)
        end

        # Apply tournament selection
        selected_parents = tournament_selection(population, fitness, tournament_size)

        # Generate offspring
        offspring = laplace_crossover(selected_parents, crossover_prob, alpha, bounds)
        offspring = power_mutation(offspring, mutation_prob, beta, bounds)

        # Enforce integer constraints + one-hot structure
        for i in eachindex(offspring)
            offspring[i] = round.(offspring[i])
            enforce_onehot!(offspring[i])
        end

        # Evaluate new population
        offspring_fitness = Vector{Float64}(undef, length(offspring))
        Threads.@threads for i in eachindex(offspring)
            offspring_fitness[i] = cost_function(offspring[i], coordFactor)
        end

        # Select the best individuals for the next generation
        population = survival_selection(population, fitness, offspring, offspring_fitness)
    end

    # Return the best solution found
    best_idx = argmin([cost_function(ind, coordFactor) for ind in population])
    x_best = round.(population[best_idx])
    enforce_onehot!(x_best)
    return x_best, cost_function(x_best, coordFactor)
end

function enforce_onehot!(x)
    for i in 0:(length(x) ÷ 5 - 1)
        group = x[5*i + 1 : 5*(i + 1)]
        idx = argmax(group)
        x[5*i + 1 : 5*(i + 1)] .= 0
        x[5*i + idx] = 1
    end
end

# Tournament Selection
function tournament_selection(population, fitness, tournament_size)
    selected = []
    for _ in 1:length(population)
        candidates = randperm(length(population))[1:tournament_size]
        best = argmin([fitness[c] for c in candidates])
        push!(selected, population[candidates[best]])
    end
    return selected
end

# Laplace Crossover
function laplace_crossover(parents, crossover_prob, alpha, bounds)
    offspring = []
    for i in 1:2:length(parents)-1
        if rand() < crossover_prob
            p1, p2 = parents[i], parents[i+1]
            b = [rand(Laplace(0, alpha)) for _ in 1:length(p1)]
            c1 = p1 + b .* (p1 - p2)
            c2 = p2 + b .* (p2 - p1)
            c1 = enforce_bounds(c1, bounds)
            c2 = enforce_bounds(c2, bounds)
            push!(offspring, c1, c2)
        else
            push!(offspring, parents[i], parents[i+1])
        end
    end
    return offspring
end

# Power Mutation
function power_mutation(offspring, mutation_prob, beta, bounds)
    for i in eachindex(offspring)
        if rand() < mutation_prob
            for j in eachindex(offspring[i])
                delta = rand()^beta  # Power mutation perturbation
                offspring[i][j] += delta * (rand(Bool) ? 1 : -1)
                offspring[i][j] = clamp(offspring[i][j], bounds[j][1], bounds[j][2])
            end
        end
    end
    return offspring
end

# Bound enforcement
function enforce_bounds(individual, bounds)
    return [clamp(individual[i], bounds[i][1], bounds[i][2]) for i in 1:length(bounds)]
end

# Survival Selection (Elitist Strategy)
function survival_selection(old_pop, old_fit, new_pop, new_fit)
    combined_pop = vcat(old_pop, new_pop)
    combined_fit = vcat(old_fit, new_fit)
    sorted_indices = sortperm(combined_fit)
    return [combined_pop[i] for i in sorted_indices[1:length(old_pop)]]
end

function test_cost_function(x)
    return (x[1] - 5)^2 + (x[2] - 10)^2 - x[3]  # Squared distance from (5,10)
    # return (maximum(x) - minimum(x))
end

function tester()
bounds = [(0, 1), (0, 1), (0,1)]  # Lower and upper bounds for each variable
best_solution, best_cost = MI_LXPM(
    test_cost_function, 
    bounds, 
    100, 
    20, 
    3, 
    0.8, 
    0.1, 
    0.3, 
    2.0
)

println("Best solution found: ", best_solution)
println("Best cost: ", best_cost)
end

end