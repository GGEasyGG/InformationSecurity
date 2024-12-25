module LinearAlgebraUtils
    export LUP_decomposition, determinant_LUP, solve_LUP, inverse_LUP

    function LUP_decomposition(A)
        n = size(A, 1)

        U = copy(A)
        L = zeros(n, n)
        P = zeros(n, n)

        for i in 1:n
            L[i, i] = 1
            P[i, i] = 1
        end

        swaps_count = 0

        for i in 1:n-1
            i_max = argmax(abs.(U[i:n, i]))[1] + i - 1

            if i_max != i
                U[i, :], U[i_max, :] = U[i_max, :], U[i, :]
                P[i, :], P[i_max, :] = P[i_max, :], P[i, :]
                L[i, 1:i-1], L[i_max, 1:i-1] = L[i_max, 1:i-1], L[i, 1:i-1]

                swaps_count += 1
            end

            for j in i+1:n
                if U[i, i] == 0.0
                    return -1
                end

                L[j, i] = U[j, i] / U[i, i]
                U[j, i:n] = mod.(U[j, i:n] - L[j, i] * U[i, i:n], 2)
            end
        end

        return L, U, Int.(P), swaps_count
    end

    function determinant_LUP(U, P, swaps_count)
        det_P = (-1)^swaps_count
        
        U_diag = [U[i, i] for i in 1:size(U, 1)]

        return det_P * prod(U_diag)
    end

    function solve_LUP(L, U, P, b)
        n = length(b)
        y = zeros(n)
        x = zeros(n)

        b_perm = P * b
        
        for i in 1:n
            y[i] = mod.(b_perm[i] - sum(L[i, 1:i-1] .* y[1:i-1]), 2)
        end

        for i in n:-1:1
            x[i] = mod.(y[i] - sum(U[i, i+1:n] .* x[i+1:n]), 2) / U[i, i]
        end

        return x
    end

    function inverse_LUP(L, U, P)
        n = size(L, 1)

        inv_A = zeros(n, n)
        
        for i in 1:n
            b = zeros(n)
            b[i] = 1
            inv_A[:, i] = solve_LUP(L, U, P, b)
        end
        
        return inv_A
    end
end