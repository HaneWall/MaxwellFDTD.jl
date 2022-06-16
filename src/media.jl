abstract type Medium end

struct StaticMedium1D <: Medium
    grid :: Grid1D
    location :: CartesianIndices{1, Tuple{UnitRange{Int64}}}
    ϵ_inf :: Float64
end

struct LorentzMedium1D <: Medium
    grid :: Grid1D
    location :: CartesianIndices{1, Tuple{UnitRange{Int64}}}
    ϵ_inf :: Float64
    oscillators :: Int64
    γ :: Vector{Float64}
    Γ :: Vector{Float64}
    ω_0 :: Vector{Float64}
    χ_1 :: Vector{Float64}
    χ_2 :: Vector{Float64}
    χ_3 :: Vector{Float64}
    function LorentzMedium1D(g::Grid1D, location::CartesianIndices{1, Tuple{UnitRange{Int64}}}, ϵ_inf::Float64, γ::Vector{Float64}, ω_0::Vector{Float64}, χ_1::Vector{Float64}, χ_2::Vector{Float64}, χ_3::Vector{Float64})
        new(g, location, ϵ_inf, length(χ_1), γ, γ.*g.Δt/2, ω_0, χ_1, χ_2, χ_3)
    end
end


#=
 These are the medias in  the two-dimesional Case.
=#

#=
 These are the medias in  the three-dimesional Case.
=#

struct StaticMedium3D <: Medium
    grid :: Grid3D
    location :: CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}
    ϵ_inf :: Float64
end

struct LorentzMedium3D <: Medium
    grid :: Grid1D
    location :: CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}
    ϵ_inf :: Float64
    oscillators :: Int64
    γ :: Vector{Float64}
    Γ :: Vector{Float64}
    ω_0 :: Vector{Float64}
    χ_1 :: Vector{Float64}
    χ_2 :: Vector{Float64}
    χ_3 :: Vector{Float64}
    function LorentzMedium3D(g::Grid3D, location::CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}, ϵ_inf::Float64, γ::Vector{Float64}, ω_0::Vector{Float64}, χ_1::Vector{Float64}, χ_2::Vector{Float64}, χ_3::Vector{Float64})
        new(g, location, ϵ_inf, length(χ_1), γ, γ.*g.Δt/2, ω_0, χ_1, χ_2, χ_3)
    end
end