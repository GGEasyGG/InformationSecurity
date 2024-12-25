module ISD
    using Random

    include("utils/linear_algebra_utils.jl")
    using .LinearAlgebraUtils

    export decodeISD

    function decodeISD(G, y, t; niter=-1)
        k, n = size(G)

        iter_count = 1

        while true
            if niter >= 0 && iter_count >= niter
                return niter
            end

            I = sort(shuffle(1:n)[1:k])

            G_I = G[:, I]

            decomposition_result = LUP_decomposition(G_I)

            if decomposition_result == -1
                iter_count += 1
                continue
            end

            L, U, P, swaps_count = decomposition_result

            if determinant_LUP(U, P, swaps_count) == 0.0
                iter_count += 1
                continue
            end

            G_I_inv = inverse_LUP(L, U, P)

            y_I = permutedims(y[I])

            m = mod.(y_I * G_I_inv, 2)

            a = mod.(m * G, 2)

            e = mod.(y - a, 2)

            if sum(e) == t
                return Int.(m), Int.(e)
            end

            iter_count += 1
        end
    end
end
