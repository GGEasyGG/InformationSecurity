using Test

include("../src/ISD.jl")
using .ISD

include("../src/STERN.jl")
using .STERN

include("../src/utils/tests_utils.jl")
using .TestsUtils

function generate_k_n()
    while true
        k, n = rand(16:32), rand(24:40)

        if n - k > 7
            return k, n
        end
    end
end

@testset "Testing decodeISD function" begin
    tests_count_with_solution = 0
    tests_count_without_solution = 0

    @testset "Testing common usage" begin
        for _ in 1:20
            k, n = generate_k_n()

            G, _ = generate_input_matrixes(k, n)

            y = permutedims(rand(0:1, n))

            t = rand(1:16)

            niter = rand(1000:10000)

            decoding_result = decodeISD(G, y, t; niter=niter)

            if length(decoding_result) != 2
                @test decoding_result == niter

                tests_count_without_solution += 1
            else
                m, e = decoding_result

                computed_y = mod.(mod.(m * G, 2) + e, 2)

                @test sum(e) == t
                @test computed_y == y

                tests_count_with_solution += 2
            end
        end
    end

    println("Testing decodeISD (2 tests per iteration): the number of tests with the found solution - ", tests_count_with_solution)
    println("Testing decodeISD (1 tests per iteration): the number of tests without a solution found - ", tests_count_without_solution)
end

@testset "Testing decodeSTERN function" begin
    tests_count_with_solution = 0
    tests_count_without_solution = 0

    @testset "Testing common usage" begin
        for _ in 1:20
            k, n = generate_k_n()

            G, H = generate_input_matrixes(k, n)

            y = rand(0:1, n)

            s = permutedims(mod.(H * y , 2))

            t = rand(4:16)

            p = rand(1:t//2-1)

            l = rand(1:n-k)

            niter = rand(1000:10000)

            decoding_result = decodeSTERN(H, s, t, p, l; niter=niter)

            if length(decoding_result) != n
                @test decoding_result == niter

                tests_count_without_solution += 1
            else
                e = decoding_result

                computed_s = permutedims(mod.(H * e, 2))

                @test sum(e) == t
                @test computed_s == s

                tests_count_with_solution += 2
            end
        end
    end

    println("Testing decodeSTERN (2 tests per iteration): the number of tests with the found solution - ", tests_count_with_solution)
    println("Testing decodeSTERN (1 tests per iteration): the number of tests without a solution found - ", tests_count_without_solution)
end
