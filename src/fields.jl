abstract type Field end

mutable struct Fields1D <: Field
    Hy :: Vector{Float64}
    Ez :: Vector{Float64}
    function Fields1D(g::Grid1D)
        new(zeros(Float64, (g.SizeX)), zeros(Float64, (g.SizeX)))
    end
end

mutable struct MaterialFields1D <:Field
    Jz_old :: Vector{Float64}
    Jz :: Vector{Float64}
    Pz :: Vector{Float64}
    PzNl :: Vector{Float64}
    function MaterialFields1D(g::Grid1D)
        new(zeros(Float64, g.SizeX), 
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