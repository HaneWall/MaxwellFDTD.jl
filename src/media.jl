abstract type Medium end

struct StaticMedium <: Medium
    grid :: Grid
    location :: CartesianIndices{1, Tuple{UnitRange{Int64}}}
    ϵ_inf :: Float64
end

struct LorentzMedium <: Medium
    grid :: Grid
    location :: CartesianIndices{1, Tuple{UnitRange{Int64}}}
    ϵ_inf :: Float64
    γ :: Array{Float64, 1}
    ω_0 :: Array{Float64, 1}
    χ_1 :: Array{Float64, 1}
    χ_2 :: Array{Float64, 1}
    χ_3 :: Array{Float64, 1}
end
