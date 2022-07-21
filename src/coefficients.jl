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

mutable struct FieldIonizationCoefficients1D <: Coeff
    grid :: Grid1D
    gamma_au :: Vector{Float64}
    function FieldIonizationCoefficients1D(g::Grid1D)
        new(
            g, 
            zeros(Float64, g.SizeX)
        )
    end
end

#=
 These are the coefficients in the two-dimensional Case.
=#
struct GridCoefficients2D_w_CPML <: Coeff
    grid :: Grid2D
    Chxh :: Array{Float64, 2}
    Chxe :: Array{Float64, 2}
    Chyh :: Array{Float64, 2}
    Chye :: Array{Float64, 2}
    Cezh :: Array{Float64, 2}
    Ceze :: Array{Float64, 2}

    Den_Hx :: Array{Float64, 1}
    Den_Ex :: Array{Float64, 1}

    Den_Hy :: Array{Float64, 1}
    Den_Ey :: Array{Float64, 1}

    function GridCoefficients2D_w_CPML(g::Grid2D, m::Vector{T}, c::CPML_Parameters_2D) where T<:Medium
        Chxh = fill(1.0,(g.SizeX, g.SizeY))
        Chxe = fill(g.Δt/(μ_0),(g.SizeX, g.SizeY))
        Chyh = fill(1.0,(g.SizeX, g.SizeY))
        Chye = fill(g.Δt/(μ_0),(g.SizeX, g.SizeY))
        Ceze = fill(1.0,(g.SizeX, g.SizeY))
        Cezh = fill(g.Δt/(ϵ_0),(g.SizeX, g.SizeY))
        
        Den_Hx = c.denominator_H_x
        Den_Ex = c.denominator_E_x

        Den_Hy = c.denominator_H_y
        Den_Ey = c.denominator_E_y

        for medium in m
            Cezh[medium.location] .= Cezh[medium.location] ./ medium.ϵ_inf
        end
        
        new(g, Chxh, Chxe, Chyh, Chye, Cezh, Ceze, Den_Hx, Den_Ex, Den_Hy, Den_Ey)
    end
end

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

# have to write new coefficients, that respect the κ_profile values (only the coefficients that are )

struct GridCoefficients3D_WIP <: Coeff
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
    κ_E :: Array{Float64, 4}
    κ_H :: Array{Float64, 4}
    function GridCoefficients3D_WIP(g::Grid3D, m::Vector{T}, b::Vector) where T<:Medium
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
        κ_E = fill(1., (g.SizeX, g.SizeY, g.SizeZ, 3))
        κ_H = fill(1.0,(g.SizeX, g.SizeY, g.SizeZ, 3))
        for medium in m
            Cexh[medium.location] .= Cexh[medium.location] ./ medium.ϵ_inf
            Ceyh[medium.location] .= Ceyh[medium.location] ./ medium.ϵ_inf
            Cezh[medium.location] .= Cezh[medium.location] ./ medium.ϵ_inf
        end
        
        for bound in b
            κ_E[bound.location, :] .+= bound.κ_E .- 1. 
            κ_H[bound.location, :] .+= bound.κ_H .- 1.
        end
        new(g, Chxh, Chxe, Chyh, Chye, Chzh, Chze, Cexh, Cexe, Ceyh, Ceye, Cezh, Ceze, κ_E, κ_H)
    end
end

struct GridCoefficients3D_w_CPML <: Coeff
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

    Den_Hx :: Array{Float64, 1}
    Den_Ex :: Array{Float64, 1}

    Den_Hy :: Array{Float64, 1}
    Den_Ey :: Array{Float64, 1}

    Den_Hz :: Array{Float64, 1}
    Den_Ez :: Array{Float64, 1}

    function GridCoefficients3D_w_CPML(g::Grid3D, m::Vector{T}, c::CPML_Parameters_3D) where T<:Medium
        Chxh = fill(1.0,(g.SizeX, g.SizeY, g.SizeZ))
        Chxe = fill(g.Δt/(μ_0),(g.SizeX, g.SizeY, g.SizeZ))
        Chyh = fill(1.0,(g.SizeX, g.SizeY, g.SizeZ))
        Chye = fill(g.Δt/(μ_0),(g.SizeX, g.SizeY, g.SizeZ))
        Chzh = fill(1.0,(g.SizeX, g.SizeY, g.SizeZ))
        Chze = fill(g.Δt/(μ_0),(g.SizeX, g.SizeY, g.SizeZ))
        Cexe = fill(1.0,(g.SizeX, g.SizeY, g.SizeZ))
        Cexh = fill(g.Δt/(ϵ_0),(g.SizeX, g.SizeY, g.SizeZ))
        Ceye = fill(1.0,(g.SizeX, g.SizeY, g.SizeZ))
        Ceyh = fill(g.Δt/(ϵ_0),(g.SizeX, g.SizeY, g.SizeZ))
        Ceze = fill(1.0,(g.SizeX, g.SizeY, g.SizeZ))
        Cezh = fill(g.Δt/(ϵ_0),(g.SizeX, g.SizeY, g.SizeZ))
        
        Den_Hx = c.denominator_H_x
        Den_Ex = c.denominator_E_x

        Den_Hy = c.denominator_H_y
        Den_Ey = c.denominator_E_y

        Den_Hz = c.denominator_H_z
        Den_Ez = c.denominator_E_z

        for medium in m
            Cexh[medium.location] .= Cexh[medium.location] ./ medium.ϵ_inf
            Ceyh[medium.location] .= Ceyh[medium.location] ./ medium.ϵ_inf
            Cezh[medium.location] .= Cezh[medium.location] ./ medium.ϵ_inf
        end
        
        new(g, Chxh, Chxe, Chyh, Chye, Chzh, Chze, Cexh, Cexe, Ceyh, Ceye, Cezh, Ceze, Den_Hx, Den_Ex, Den_Hy, Den_Ey, Den_Hz, Den_Ez)
    end
end

struct GridCoefficients1D_w_CPML <: Coeff
    grid :: Grid1D
    Chyh :: Array{Float64, 1}
    Chye :: Array{Float64, 1}
    Cezh :: Array{Float64, 1}
    Ceze :: Array{Float64, 1}

    Den_Hx :: Array{Float64, 1}
    Den_Ex :: Array{Float64, 1}

    function GridCoefficients1D_w_CPML(g::Grid1D, m::Vector{T}, c::CPML_Parameters_1D) where T<:Medium
        Chyh = fill(1.0,(g.SizeX))
        Chye = fill(g.Δt/(μ_0),(g.SizeX))
        Ceze = fill(1.0,(g.SizeX))
        Cezh = fill(g.Δt/(ϵ_0),(g.SizeX))
        for medium in m
            Cezh[medium.location] .= Cezh[medium.location] ./ medium.ϵ_inf
        end

        Den_Hx = c.denominator_H_x
        Den_Ex = c.denominator_E_x

        new(g, Chyh, Chye, Cezh, Ceze, Den_Hx, Den_Ex)
    end
end