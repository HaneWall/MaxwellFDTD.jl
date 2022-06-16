abstract type Coeff end

struct GridCoefficients1D <: Coeff
    grid :: Grid1D
    Chyh :: Array{Float64, 1}
    Chye :: Array{Float64, 1}
    Cezh :: Array{Float64, 1}
    Ceze :: Array{Float64, 1}
    function GridCoefficients1D(g::Grid1D, m::Vector{T}) where T<:Medium
        Chyh = fill(1.0,(g.SizeX))
        Chye = fill(g.Δt/(μ_0 * g.Δx),(g.SizeX))
        Ceze = fill(1.0,(g.SizeX))
        Cezh = fill(g.Δt/(ϵ_0 * g.Δx),(g.SizeX))
        for medium in m
            Cezh[medium.location] .= Cezh[medium.location] ./ medium.ϵ_inf
        end
        new(g, Chyh, Chye, Cezh, Ceze)
    end
end