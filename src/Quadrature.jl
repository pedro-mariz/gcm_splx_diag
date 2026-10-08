# Even-tempered quadrature
function even_tempered_grid(β_min, β_max, N)
    Δβ = (β_max - β_min) / (N - 1)

    grid = range(β_max, β_max, length=N)

    weights = fill(Δβ, N)

    return grid, weights
end

# Calculation of Sg matrix elements
function calc_Sg_uv(au, bu, av, bv, grid, weights)
    Sg_uv = 0.0
    for i in eachindex(grid), j in eachindex(grid)
        βi, βj = grid[i], grid[j]

        Giu = gumbel(βi, au, bu)
        Gjv = gumbel(βj, av, bv)

        αi, αj = exp(βi), exp(βj)

        Sg_uv += weights[i] * weights[j] * Giu * Gjv * S_int(αi, αj)
    end
    return Sg_uv
end

# Calculation of hg matrix elements
function calc_hg_uv(au, bu, av, bv, grid, weights; Z=2.0)
    hg_uv = 0.0
    for i in eachindex(grid), j in eachindex(grid)
        βi, βj = grid[i], grid[j]

        Giu = gumbel(βi, au, bu)
        Gjv = gumbel(βj, av, bv)

        αi, αj = exp(βi), exp(βj)

        hg_uv += weights[i] * weights[j] * Giu * Gjv * h_int(αi, αj; Z=Z)
    end
    return hg_uv
end

# Calculation of Veeg tensor elements
function calc_Veeg_uvmn(au, bu, av, bv, am, bm, an, bn, grid, weights)
    Veeg_uvmn = 0.0
    for i in eachindex(grid), j in eachindex(grid), k in eachindex(grid), l in eachindex(grid)
        βi, βj, βk, βl = grid[i], grid[j], grid[k], grid[l]

        Giu = gumbel(βi, au, bu)
        Gjv = gumbel(βj, av, bv)
        Gkm = gumbel(βk, am, bm)
        Gln = gumbel(βl, an, bn)

        αi, αj, αk, αl = exp(βi), exp(βj), exp(βk), exp(βl)

        Veeg_uvmn += weights[i] * weights[j] * weights[i] * weights[j] *
                    Giu * Gjv * Gkm * Gln * Vee_int(αi, αj, αk, αl)
    end
    return Veeg_uvmn
end