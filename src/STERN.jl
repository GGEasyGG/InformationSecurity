module STERN
    using Random

    include("utils/linear_algebra_utils.jl")
    using .LinearAlgebraUtils

    export decodeSTERN, decodeSTERN_for_analysis

    function find_information_set_and_H_parts(H, r, n)
        I = sort(shuffle(1:n)[1:r])

        H_L = H[:, I]
        H_R = H[:, setdiff(1:n, I)]

        decomposition_result = LUP_decomposition(H_L)

        if decomposition_result == -1
            return -1
        end

        L, U, P, swaps_count = decomposition_result

        if determinant_LUP(U, P, swaps_count) == 0.0
            return -1
        end

        return I, inverse_LUP(L, U, P), H_R
    end

    function partition_set(S)
        new_S = []

        for i in 1:length(S)
            push!(new_S, [i, S[i]])
        end

        shuffle(new_S)

        half_size = div(length(new_S), 2)
        X = new_S[1:half_size]
        Y = new_S[(half_size + 1):end]
        return X, Y
    end

    function compute_U(R, subset, p)
        U = Dict()

        for combination in combinations(subset, p)
            U[combination] = mod.(sum(R[:, elem[1]] for elem in combination), 2)
        end

        return U
    end

    function find_matching_projections(U_A_dict, U_B_dict, L, z)
        matching_pairs = []

        for (A, U_A) in U_A_dict
            for (B, U_B) in U_B_dict
                if mod.(U_A[L] + z, 2) == U_B[L]
                    push!(matching_pairs, (A, B))
                end
            end
        end

        return matching_pairs
    end

    function construct_e(A, B, U, I, y, r, n)
        e = zeros(Int, n)

        for i in A
            e[i[2]] = 1
        end

        for i in B
            e[i[2]] = 1
        end

        for i in 1:r
            e[I[i]] = mod(y[i] + U[i], 2)
        end

        return e
    end

    function combinations(S, k)
        n = length(S)

        result = []

        function combine(idx, current_comb)
            if length(current_comb) == k
                push!(result, current_comb)
                return
            end

            for i in idx:n
                combine(i + 1, push!(copy(current_comb), S[i]))
            end
        end

        combine(1, [])

        return shuffle(result)
    end

    function decodeSTERN(H, s, t, p, l; niter=-1)
        r, n = size(H)

        iteration = 0

        while niter == -1 || iteration < niter
            iteration += 1

            result = find_information_set_and_H_parts(H, r, n)

            if result == -1
                continue
            end

            I, H_L_inv, H_R = result

            R = mod.(H_L_inv * H_R, 2)

            y = mod.(H_L_inv * permutedims(s), 2)

            X, Y = partition_set(sort(setdiff(1:n, I)))

            L = sort(shuffle(collect(1:r))[1:l])

            z = s[L]

            U_A_dict = compute_U(R, X, p)
            U_B_dict = compute_U(R, Y, p)

            matching_pairs = find_matching_projections(U_A_dict, U_B_dict, L, z)

            for (A, B) in matching_pairs
                U = mod.(U_A_dict[A] + U_B_dict[B], 2)
                if sum(mod.(U + y, 2)) == t - 2 * p
                    return construct_e(A, B, U, I, y, r, n)
                end
            end
        end

        return niter
    end

    function decodeSTERN_for_analysis(H, s, t, p, l; niter=-1)
        r, n = size(H)

        iteration = 0

        while niter == -1 || iteration < niter
            iteration += 1

            result = find_information_set_and_H_parts(H, r, n)

            if result == -1
                continue
            end

            I, H_L_inv, H_R = result

            R = mod.(H_L_inv * H_R, 2)

            y = mod.(H_L_inv * permutedims(s), 2)

            X, Y = partition_set(sort(setdiff(1:n, I)))

            L = sort(shuffle(collect(1:r))[1:l])

            z = s[L]

            U_A_dict = compute_U(R, X, p)
            U_B_dict = compute_U(R, Y, p)

            matching_pairs = find_matching_projections(U_A_dict, U_B_dict, L, z)

            for (A, B) in matching_pairs
                U = mod.(U_A_dict[A] + U_B_dict[B], 2)
                if sum(mod.(U + y, 2)) == t - 2 * p
                    return construct_e(A, B, U, I, y, r, n), iteration
                end
            end
        end

        return niter
    end
end
