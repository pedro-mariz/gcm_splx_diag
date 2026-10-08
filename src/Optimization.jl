using LinearAlgebra

# Build matrices for SCF
function build_matrices(a::Vector{Float64}, b::Vector{Float64}, grid, weights)
    # Número de funções Gumbel
    Ng = length(a)

    Sg = Matrix{Float64}(undef, N, N)
    hg = Matrix{Float64}(undef, N, N)
    Vee = Array{Float64, 4}(undef, N, N, N, N)

    # Building S and h blocks using diagonal symmetry
    for u in 1:N, v in 1:u
        Sg_val = calc_Sg_uv(a[u], b[u], a[v], b[v], grid, weights)
        hg_val = calc_hg_uv(a[u], b[u], a[v], b[v], grid, weights; Z=2.0)

        # Assing to both diagonals
        Sg[u, v], Sg[v, u] = Sg_val, Sg_val
        hg[u, v], hg[v, u] = hg_val, hg_val
    end # uv
    for u in 1:N, v in 1:u
        # Lienar index for (u,v), same as counting up all pairs (sum up to (u - 1) plus v)
        uv = (u * (u - 1)) ÷ 2 + v

        for m in 1:N, n in 1:m
            # Linear index for (m,n)
            mn = (m * (m + 1)) ÷ 2 + n

            # Calculate only if it is the upper diagonal block of tensor
            if uv >= mn
                Veeg_val = calc_Veeg_uvmn(
                    a[u], b[u], a[v], b[v],
                    a[m], b[m], a[n], b[n],
                    grid, weights
                )

                # Assing all symmetric elements
                Vee[u, v, m, n] = v_val
                # m|n permutation
                Vee[u, v, n, m] = v_val
                # u|v permutations
                Vee[v, u, m, n] = v_val
                Vee[v, u, n, m] = v_val
                # uv|mn permutations
                Vee[m, n, u, v] = v_val
                Vee[n, m, u, v] = v_val
                Vee[m, n, v, u] = v_val
                Vee[n, m, v, u] = v_val
            end # uv >= mn
        end # mn
    end # uv

    return Sg, hg, Veeg
end

function solve_lin_coef(Sg, h, Veeg, max_iters=100, tol=1e-8)
    N = size(h, 1)

    # Initial guess: Fg = hg
    Fg = copy(hg)

    # Solve the equation Fg*c = Sg*c*e
    E_val, C = eigen(Fg, Sg)
    
    # Get 1s coeficients and density matrix
    c = C[:, 1]
    P = 2.0 * (c * c')

    E_total_old = 0.0

    # Start SCF loop
    for iter in 1:max_iters
        # calculate new Fg
        Fg .= hg
        for u in 1:N, v in 1:N
            g_uv = 0.0
            for m in 1:N, n in 1:N
                g_uv += P[m, n] * (Veeg[u, v, m, n] - 0.5 * Veeg[u, m, v, n])
            end
            F[u, v] = += g_uv
        end

        # Solve de equation for new c
        E_val, C = eigen(Fg, Sg)
        c = C[:, 1]
        P = 2.0 * (c * c')

        # Calculate energy
        E_total = 0.5 * sum(P .* (hg + Fg))

        # Convergence check
        if abs(E_total - E_total_old) < tol
            return E_total, c
        end

        E_total_old = E_total
    end

    @warn "SCF did not converge after $max_iters iterations!"
    return E_total_old, c
end