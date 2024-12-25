using Statistics, CSV, DataFrames

include("../src/ISD.jl")
using .ISD

include("../src/STERN.jl")
using .STERN

include("../src/utils/tests_utils.jl")
using .TestsUtils

function evaluate_stern(H, s, t, p_list, l_list, niter, trials)
    results = []

    for p in p_list
        for l in l_list
            execution_times = []
            iters = []

            for trial in 1:trials
                start_time = time()
                decoding_result = decodeSTERN_for_analysis(H, s, t, p, l; niter=niter)
                execution_time = time() - start_time

                push!(execution_times, execution_time)
                push!(iters, decoding_result[2])
            end

            mean_execution_time = mean(execution_times)
            mean_iterations = mean(iters)

            push!(results, (p, l, mean_execution_time, mean_iterations))
        end
    end

    return results
end

function write_CSV(results_p, results_l)
    p_list = []
    l_list = []
    times = []
    iterations = []

    for result in results_p
        push!(p_list, result[1])
        push!(l_list, result[2])
        push!(times, result[3])
        push!(iterations, result[4])
    end

    for result in results_l
        push!(p_list, result[1])
        push!(l_list, result[2])
        push!(times, result[3])
        push!(iterations, result[4])
    end

    data = DataFrame(
        p = p_list,
        l = l_list,
        time = times,
        iterations = iterations
    )

    CSV.write("p_and_l_analysis.csv", data)
end

function p_and_l_analysis()
    println("Start creation of CSV file with results of p and l analysis")
    println()

    niter = 10000

    k, n = 32, 64

    t_list = [16]
    t = t_list[1]

    p_list = [1, 2, 4, 6]
    l_list = [1, 2, 4, 8, 16]


    println("Start finding params procedure")

    _, H, _, s = find_params(k, n, niter, p_list[1], l_list[1], t_list, true)

    println()

    niter = -1

    println("Start p analysis procedure")

    results_p = evaluate_stern(H, s, t, p_list, [l_list[1]], niter, 100)

    println("Stop p analysis")

    println()

    println("Start l analysis procedure")

    results_l = evaluate_stern(H, s, t, [p_list[1]], l_list, niter, 100)

    println("Stop l analysis")

    write_CSV(results_p, results_l)

    println()
    println("Results successfully saved to p_and_l_analysis.csv")
    println("----------------------------------------------------------------")
end

p_and_l_analysis()
