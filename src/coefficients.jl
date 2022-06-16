abstract type Coeff end


#=
 These are the coefficients in the one-dimensional Case.
=#
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

#=
 These are the coefficients in the two-dimensional Case.
=#

#=
 These are the coefficients in the three-dimensional Case.
=#

struct GridCoefficients3D <: Coeff
    grid :: Grid3D
    Chxh :: Array{Float64, 3}
    Chxe :: Array{Float64, 3}
    Chyh :: Array{Float64, 3}
    Chye :: Array{Float64, 3}
    Chzh :: Array{Float64, 3}
    Chze :: Array{Float64, 3}
    Cexh :: Array{Float64, 3}
    Cexe :: Array{Float64, 3}
    Ceyh :: Array{Float64, 3}
    Ceye :: Array{Float64, 3}
    Cezh :: Array{Float64, 3}
    Ceze :: Array{Float64, 3}
    function GridCoefficients3D(g::Grid3D, m::Vector{T}) where T<:Medium
        Chxh = fill(1.0,(g.SizeX, g.SizeY, g.SizeZ))
        Chxe = fill(g.Δt/(μ_0 * g.Δx),(g.SizeX, g.SizeY, g.SizeZ))
        Chyh = fill(1.0,(g.SizeX, g.SizeY, g.SizeZ))
        Chye = fill(g.Δt/(μ_0 * g.Δx),(g.SizeX, g.SizeY, g.SizeZ))
        Chzh = fill(1.0,(g.SizeX, g.SizeY, g.SizeZ))
        Chze = fill(g.Δt/(μ_0 * g.Δx),(g.SizeX, g.SizeY, g.SizeZ))
        Cexe = fill(1.0,(g.SizeX, g.SizeY, g.SizeZ))
        Cexh = fill(g.Δt/(ϵ_0 * g.Δx),(g.SizeX, g.SizeY, g.SizeZ))
        Ceye = fill(1.0,(g.SizeX, g.SizeY, g.SizeZ))
        Ceyh = fill(g.Δt/(ϵ_0 * g.Δx),(g.SizeX, g.SizeY, g.SizeZ))
        Ceze = fill(1.0,(g.SizeX, g.SizeY, g.SizeZ))
        Cezh = fill(g.Δt/(ϵ_0 * g.Δx),(g.SizeX, g.SizeY, g.SizeZ))
        for medium in m
            Cexh[medium.location] .= Cexh[medium.location] ./ medium.ϵ_inf
            Ceyh[medium.location] .= Ceyh[medium.location] ./ medium.ϵ_inf
            Cezh[medium.location] .= Cezh[medium.location] ./ medium.ϵ_inf
        end
        new(g, Chxh, Chxe, Chyh, Chye, Chzh, Chze, Cexh, Cexe, Ceyh, Ceye, Cezh, Ceze)
    end
end