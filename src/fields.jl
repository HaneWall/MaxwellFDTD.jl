abstract type Field end

#=
 These are the Fields in  the one-dimesional Case.
=#
mutable struct Fields1D <: Field
    Hy :: Vector{Float64}
    Ez :: Vector{Float64}
    function Fields1D(g::Grid1D)
        new(zeros(Float64, (g.SizeX)), zeros(Float64, (g.SizeX)))
    end
end

mutable struct MaterialFields1D <:Field
    Jz :: Vector{Float64}
    Pz :: Vector{Float64}
    PzNl :: Vector{Float64}
    ρ_cb :: Vector{Float64}
    Γ_ADK :: Vector{Float64}
    function MaterialFields1D(g::Grid1D)
        new(zeros(Float64, g.SizeX),
            zeros(Float64, g.SizeX),
            zeros(Float64, g.SizeX),
            zeros(Float64, g.SizeX),
            zeros(Float64, g.SizeX)
            )
    end
end

mutable struct LorentzFields1D <: Field
    Jz :: Array{Float64}
    Pz :: Array{Float64}
    PzNl :: Array{Float64}
    function LorentzFields1D(m::LorentzMedium1D)
        new(
            zeros(Float64, (size(m.location)[1], m.oscillators)), 
            zeros(Float64, (size(m.location)[1], m.oscillators)),
            zeros(Float64, (size(m.location)[1], m.oscillators))
            )
    end
end

mutable struct DrudeFields1D <: Field
    Jz_free :: Array{Float64}
    function DrudeFields1D(m::DrudeMedium1D)
        new(
            zeros(Float64, (size(m.location)[1]))
        )
    end
end

mutable struct TunnelFields1D <: Field
    Jz_tunnel :: Vector{Float64}
    Pz_tunnel :: Vector{Float64}
    Jz_T_brunel :: Vector{Float64}
    dz_T :: Vector{Float64}
    function TunnelFields1D(m::TunnelMedium1D)
        new(
            zeros(Float64, (size(m.location)[1])), 
            zeros(Float64, (size(m.location)[1])),
            zeros(Float64, (size(m.location)[1])),
            zeros(Float64, (size(m.location)[1]))
            )
    end
end

#=
 These are the Fields in  the two-dimesional Case.
=#

#=
 These are the Fields in  the three-dimesional Case.
=#
mutable struct Fields3D <: Field
    Hx :: Array{Float64, 3}
    Hy :: Array{Float64, 3}
    Hz :: Array{Float64, 3}
    Ex :: Array{Float64, 3}
    Ey :: Array{Float64, 3}
    Ez :: Array{Float64, 3}
    function Fields1D(g::Grid3D)
        new(zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ)), 
            zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ)), 
            zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ)), 
            zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ)), 
            zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ)), 
            zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ))
            )
    end
end

mutable struct MaterialFields3D <:Field
    Jx :: Array{Float64, 3}
    Jy :: Array{Float64, 3}
    Jz :: Array{Float64, 3}
    Px :: Array{Float64, 3}
    Py :: Array{Float64, 3}
    Pz :: Array{Float64, 3}
    function MaterialFields3D(g::Grid3D)
        new(zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ)), 
            zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ)), 
            zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ)), 
            zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ)), 
            zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ)), 
            zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ))
            )
    end
end

mutable struct LorentzFields3D <: Field
    Jx :: Array{Float64}
    Jy :: Array{Float64}
    Jz :: Array{Float64}
    Px :: Array{Float64}
    Py :: Array{Float64}
    Pz :: Array{Float64}
    PxNl :: Array{Float64}
    PyNl :: Array{Float64}
    PzNl :: Array{Float64}
    function LorentzFields3D(m::LorentzMedium3D)
        new(
            zeros(Float64, (size(m.location)[1],size(m.location)[2],size(m.location)[3], m.oscillators)), 
            zeros(Float64, (size(m.location)[1],size(m.location)[2],size(m.location)[3], m.oscillators)),
            zeros(Float64, (size(m.location)[1],size(m.location)[2],size(m.location)[3], m.oscillators)),
            zeros(Float64, (size(m.location)[1],size(m.location)[2],size(m.location)[3], m.oscillators)),
            zeros(Float64, (size(m.location)[1],size(m.location)[2],size(m.location)[3], m.oscillators)),
            zeros(Float64, (size(m.location)[1],size(m.location)[2],size(m.location)[3], m.oscillators)),
            zeros(Float64, (size(m.location)[1],size(m.location)[2],size(m.location)[3], m.oscillators)),
            zeros(Float64, (size(m.location)[1],size(m.location)[2],size(m.location)[3], m.oscillators)),
            zeros(Float64, (size(m.location)[1],size(m.location)[2],size(m.location)[3], m.oscillators)),
            )
    end
end
