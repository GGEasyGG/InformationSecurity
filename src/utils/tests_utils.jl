module TestsUtils
    using LinearAlgebra

    include("../ISD.jl")
    using .ISD

    include("../STERN.jl")
    using .STERN

    export generate_input_matrixes, find_params

    function row_reduction(A, k)
        m, n = size(A)

        for col in 1:n
            pivot_row = argmax(abs.(A[col:m, col])) + col - 1

            if A[pivot_row, col] == 0
                continue
            end

            A[[col, pivot_row], :] .= A[[pivot_row, col], :]
            A[col, :] .= mod.(A[col, :], 2)

            for row in col+1:m
                if A[row, col] == 1
                    A[row, :] .= mod.(A[row, :] - A[col, :], 2)
                end
            end

            if A[1:k, 1:k] == Matrix(I, k, k) && A[k+1:m, 1:k] == zeros(m - k, k)
                return A
            end
        end

        return A
    end

    function find_verification_matrix(G)
        k, n = size(G)

        F = hcat(G', Matrix(I, n, n))

        H = row_reduction(F, k)[k+1:n, k+1:k+n]

        if all(check_verification_matrix(H, G, n, k))
            return H
        else
            return -1
        end
    end

    function check_verification_matrix(H, G, n, k)
        return mod.(H * G', 2) .== zeros(Int, n-k, k)
    end

    function generate_G_matrix(k, n)
        G = Matrix{Int}(I, k, k)

        remaining_columns = n - k

        for i in 1:remaining_columns
            G = hcat(G, rand(Bool, k))
        end

        if rank(G) != k
            return -1
        end
        
        return G
    end

    function generate_input_matrixes(k ,n)
        while true
            G = generate_G_matrix(k, n)
            H = find_verification_matrix(G)

            if G != -1 && H != -1
                return G, H
            end
        end
    end

    function find_params(k, n, niter, p, l, t_list, has_solution)
        while true
            G, H = generate_input_matrixes(k, n)

            y = rand(0:1, n)

            s = permutedims(mod.(H * y , 2))

            y = permutedims(y)

            count = 0

            for t in t_list
                decoding_result_ISD = decodeISD(G, y, t; niter=niter)
                decoding_result_STERN = decodeSTERN(H, s, t, p, l; niter=niter)

                if has_solution
                    if length(decoding_result_ISD) == 2 && length(decoding_result_STERN) == n
                        count += 1
                    end
                else
                    if length(decoding_result_ISD) != 2 && length(decoding_result_STERN) != n
                        count += 1
                    end
                end
            end

            if count == length(t_list)
                println("Params successfully finded")
                return G, H, y, s
            end
        end
    end
end