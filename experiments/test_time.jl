using CSV, DataFrames

include("../src/ISD.jl")
using .ISD

include("../src/STERN.jl")
using .STERN

include("../src/utils/tests_utils.jl")
using .TestsUtils

function find_execution_time(G, H, y, s, t, p, l, niter, iter_count)
	STERN_full_time = 0
	ISD_full_time = 0

	for i in 1:iter_count
		STERN_execution_time = @elapsed begin
			decoding_result = decodeSTERN(H, s, t, p, l; niter=niter)
		end

		STERN_full_time += STERN_execution_time

		ISD_execution_time = @elapsed begin
			decoding_result = decodeISD(G, y, t; niter=niter)
		end

		ISD_full_time += ISD_execution_time
	end

	return STERN_full_time / iter_count, ISD_full_time / iter_count
end

function test_time_per_iteration(k, n)
	println("Start test for finding time per iteration for k=", k, " and n=", n)

	niter = 10

	p = 1

	l = 2

	t_list = [4, 8, 16, 32]

	println("Start finding params procedure")

	G, H, y, s = find_params(k, n, niter, p, l, t_list, false)

	println("Start time execution")

	n_list = []
	k_list = []
	STERN_times = []
	ISD_times = []

	for t in t_list
		STERN_time, ISD_time = find_execution_time(G, H, y, s, t, p, l, niter, 100)

		push!(STERN_times, STERN_time / niter)
		push!(ISD_times, ISD_time / niter)
		push!(k_list, k)
		push!(n_list, n)
	end

	println("Stop test for finding time per iteration for k=", k, " and n=", n)
	println()

	return k_list, n_list, t_list, STERN_times, ISD_times
end

function test_find_solution_time(k, n)
	println("Start test for finding time of solution finding for k=", k, " and n=", n)

	niter = 10000

	p = 1

	l = 2

	t_list = [24, 28, 32, 36]

	println("Start finding params procedure")

	G, H, y, s = find_params(k, n, niter, p, l, t_list, true)

	niter = -1

	println("Start time execution")

	n_list = []
	k_list = []
	STERN_times = []
	ISD_times = []

	for t in t_list
		STERN_time, ISD_time = find_execution_time(G, H, y, s, t, p, l, niter, 100)

		push!(STERN_times, STERN_time)
		push!(ISD_times, ISD_time)
		push!(k_list, k)
		push!(n_list, n)
	end

	println("Stop test for finding time of solution finding for k=", k, " and n=", n)
	println()

	return k_list, n_list, t_list, STERN_times, ISD_times
end

function create_test_time_per_iteration_csv()
	println("Start creation of CSV file with results of tests for finding time per iteration")
	println()

    k_list = [16, 24, 32, 40, 48]
    n = 64

    all_k_list = []
    all_n_list = []
    all_STERN_times = []
    all_ISD_times = []
    all_t_list = []

    for k in k_list
		k_list, n_list, t_list, STERN_times, ISD_times = test_time_per_iteration(k, n)
		append!(all_k_list, k_list)
		append!(all_n_list, n_list)
		append!(all_t_list, t_list)
		append!(all_STERN_times, STERN_times)
		append!(all_ISD_times, ISD_times)
	end

	data = DataFrame(
	    k = all_k_list,
	    n = all_n_list,
	    t = all_t_list,
	    ISD_time = all_ISD_times,
	    STERN_time = all_STERN_times
	)

	CSV.write("test_time_per_iteration.csv", data)

	println("Results successfully saved to test_time_per_iteration.csv")
	println("----------------------------------------------------------------------------------------------")
end

function create_test_find_solution_time_csv()
	println("Start creation of CSV file with results of tests for finding time of solution finding")
	println()

    k_list = [48, 64, 80]
    n = 128

    all_k_list = []
    all_n_list = []
    all_STERN_times = []
    all_ISD_times = []
    all_t_list = []

    for k in k_list
		k_list, n_list, t_list, STERN_times, ISD_times = test_find_solution_time(k, n)
		append!(all_k_list, k_list)
		append!(all_n_list, n_list)
		append!(all_t_list, t_list)
		append!(all_STERN_times, STERN_times)
		append!(all_ISD_times, ISD_times)
	end

	data = DataFrame(
	    k = all_k_list,
	    n = all_n_list,
	    t = all_t_list,
	    ISD_time = all_ISD_times,
	    STERN_time = all_STERN_times
	)

	CSV.write("test_find_solution_time.csv", data)

	println("Results successfully saved to test_find_solution_time.csv")
	println("----------------------------------------------------------------------------------------------")
end

create_test_time_per_iteration_csv()
create_test_find_solution_time_csv()
