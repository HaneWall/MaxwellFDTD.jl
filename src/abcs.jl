abstract type ABC end
abstract type CPML end

function ABC!(F::Fields1D, g::Grid1D)
    F.Ez[1] = F.Ez[2]
    F.Ez[g.SizeX] = F.Ez[g.SizeX-1]
end

mutable struct LeftSideMurABC <: ABC
    location :: CartesianIndex
    courant :: Float64
    prev_array :: Array{Float64, 1}
    function LeftSideMurABC(g::Grid1D, p::CartesianIndex)
        prev_array = zeros(Float64, 2)
        new(p, g.S_c, prev_array)
    end
end

function saveFields!(B::LeftSideMurABC, F::Fields1D)
    B.prev_array[1] = F.Ez[B.location]
    B.prev_array[2] = F.Ez[B.location + CartesianIndex((1,))]
end

function stepABC!(F::Fields1D, B::LeftSideMurABC)
    F.Ez[B.location] = B.prev_array[2] + (B.courant - 1)/(B.courant + 1) * (F.Ez[B.location + CartesianIndex((1,))] - B.prev_array[1])
end

mutable struct RightSideMurABC <: ABC
    location :: CartesianIndex
    courant :: Float64
    prev_array :: Array{Float64, 1}
    function RightSideMurABC(g::Grid1D, p::CartesianIndex)
        prev_array = zeros(Float64, 2)
        new(p, g.S_c, prev_array)
    end
end

function saveFields!(B::RightSideMurABC, F::Fields1D)
    B.prev_array[1] = F.Ez[B.location - CartesianIndex((1,))]
    B.prev_array[2] = F.Ez[B.location]
end

function stepABC!(F::Fields1D, B::RightSideMurABC)
    F.Ez[B.location] = B.prev_array[1] + (B.courant - 1)/(B.courant + 1) * (F.Ez[B.location - CartesianIndex((1,))] - B.prev_array[2])
end


# # 3D - boundaries

# function σ_opt(m::Float64, Δx::Float64)
#     return 0.8*(m-1)/(376.730 * Δx)
# end

# function σ_profile(arr::Array{Float64, 1}, thickness::Int64, Δx::Float64; m::Float64=4.0)
#     return σ_opt(m, Δx) .* (arr ./ thickness).^m
# end

# function b_coeff(σ::Array{Float64}, κ_profile::Array{Float64}, α_profile::Array{Float64}, g::Grid3D)
#     return exp.(-(σ./(ϵ_0 .* κ_profile) .+ α_profile./ϵ_0) .* g.Δt)
# end

# function c_coeff(b::Array{Float64}, σ::Array{Float64}, κ::Array{Float64}, α::Array{Float64})
#     return (b .- 1).*σ ./ (σ .* κ .+ α .* κ.^2)
# end

# function κ_profile(arr::Array{Float64, 1}, thickness::Int64; m::Float64=4.0, κ_max::Float64=14.0)
#     return 1 .+ (κ_max .-1).*(arr./thickness).^m
# end

# function α_profile(arr::Array{Float64, 1}, thickness::Int64; m_a::Float64=1.0)
#     return 0.25.* ((thickness .- arr)./thickness).^m_a
# end

# mutable struct PMLXlow <: CPML
#     grid :: Grid3D
#     location :: CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}
#     thickness :: Int64

#     κ_E :: Array{Float64, 4}
#     κ_H :: Array{Float64, 4}

#     c_E :: Array{Float64, 4}
#     b_E :: Array{Float64, 4}
#     c_H :: Array{Float64, 4}
#     b_H :: Array{Float64, 4}

#     Ψ_Ex :: Array{Float64, 4}
#     Ψ_Ey :: Array{Float64, 4}
#     Ψ_Ez :: Array{Float64, 4}
#     Ψ_Hx :: Array{Float64, 4}
#     Ψ_Hy :: Array{Float64, 4}
#     Ψ_Hz :: Array{Float64, 4}


#     Φ_E :: Array{Float64, 4}
#     Φ_H :: Array{Float64, 4}
#     function PMLXlow(g::Grid3D, location::CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}})
#         σ_E_profile = σ_profile(Array(size(location)[1]-0.5:-1.:0.), size(location)[1], g.Δx) 
#         κ_E_profile = κ_profile(Array(size(location)[1]-0.5:-1.:0.), size(location)[1])
#         α_E_profile = α_profile(Array(size(location)[1]-0.5:-1.:0.), size(location)[1])

#         σ_E = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         κ_E = fill(1., size(location)[1], size(location)[2], size(location)[3], 3)
#         α_E = fill(1e-8, size(location)[1], size(location)[2], size(location)[3], 3)
#         na = [CartesianIndex()]

#         @. σ_E[:, :, :, 1] = σ_E_profile[:, na, na]
#         @. κ_E[:, :, :, 1] = κ_E_profile[:, na, na]
#         @. α_E[:, :, :, 1] = α_E_profile[:, na, na]

        
#         σ_H_profile = σ_profile(Array(size(location)[1]-1.:-1.:1.), size(location)[1], g.Δx) 
#         κ_H_profile = κ_profile(Array(size(location)[1]-1:-1.:1.), size(location)[1])
#         α_H_profile = α_profile(Array(size(location)[1]-1:-1.:1.), size(location)[1])
        

#         σ_H = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         κ_H = fill(1., size(location)[1], size(location)[2], size(location)[3], 3)
#         α_H = fill(1e-8, size(location)[1], size(location)[2], size(location)[3], 3)

#         @. σ_H[1:end-1, :, :, 1] = σ_H_profile[:, na, na]
#         @. κ_H[1:end-1, :, :, 1] = κ_H_profile[:, na, na]
#         @. α_H[1:end-1, :, :, 1] = α_H_profile[:, na, na]
        
#         b_E = b_coeff(σ_E, κ_E, α_E, g)
#         c_E = c_coeff(b_E, σ_E, κ_E, α_E)
#         b_H = b_coeff(σ_H, κ_H, α_H, g)
#         c_H = c_coeff(b_H, σ_H, κ_H, α_H)
        
#         Ψ_Ex = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Ψ_Ey = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Ψ_Ez = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3) 
#         Ψ_Hx = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3) 
#         Ψ_Hy = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Ψ_Hz = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)

#         Φ_E = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Φ_H = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)

#         new(
#             g, 
#             location, 
#             size(location)[1],
#             κ_E,
#             κ_H,
#             c_E,
#             b_E,
#             c_H,
#             b_H, 
#             Ψ_Ex,
#             Ψ_Ey,
#             Ψ_Ez,
#             Ψ_Hx,
#             Ψ_Hy,
#             Ψ_Hz,
#             Φ_E,
#             Φ_H 
#         )
#     end
# end

# # TODO: fix kappas and alphas 
# mutable struct PMLXhigh <: CPML
#     grid :: Grid3D
#     location :: CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}
#     thickness :: Int64

#     κ_E :: Array{Float64, 4}
#     κ_H :: Array{Float64, 4}

#     c_E :: Array{Float64, 4}
#     b_E :: Array{Float64, 4}
#     c_H :: Array{Float64, 4}
#     b_H :: Array{Float64, 4}

#     Ψ_Ex :: Array{Float64, 4}
#     Ψ_Ey :: Array{Float64, 4}
#     Ψ_Ez :: Array{Float64, 4}
#     Ψ_Hx :: Array{Float64, 4}
#     Ψ_Hy :: Array{Float64, 4}
#     Ψ_Hz :: Array{Float64, 4}


#     Φ_E :: Array{Float64, 4}
#     Φ_H :: Array{Float64, 4}
#     function PMLXhigh(g::Grid3D, location::CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}})
#         σ_E_profile = σ_profile(Array(0.5:1.:size(location)[1]), size(location)[1], g.Δx) 
#         κ_E_profile = κ_profile(Array(0.5:1.:size(location)[1]), size(location)[1])
#         α_E_profile = α_profile(Array(0.5:1.:size(location)[1]), size(location)[1])

#         σ_E = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         κ_E = fill(1., size(location)[1], size(location)[2], size(location)[3], 3)
#         α_E = fill(1e-8, size(location)[1], size(location)[2], size(location)[3], 3)
#         na = [CartesianIndex()]

#         @. σ_E[:, :, :, 1] = σ_E_profile[:, na, na]
#         @. κ_E[:, :, :, 1] = κ_E_profile[:, na, na]
#         @. α_E[:, :, :, 1] = α_E_profile[:, na, na]
        
#         σ_H_profile = σ_profile(Array(1.:1.:size(location)[1]-1), size(location)[1], g.Δx) 
#         κ_H_profile = κ_profile(Array(1.:1.:size(location)[1]-1), size(location)[1])
#         α_H_profile = α_profile(Array(1.:1.:size(location)[1]-1), size(location)[1])

#         σ_H = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         κ_H = fill(1., size(location)[1], size(location)[2], size(location)[3], 3)
#         α_H = fill(1e-8, size(location)[1], size(location)[2], size(location)[3], 3)

#         @. σ_H[1:end-1, :, :, 1] = σ_H_profile[:, na, na]
#         @. κ_H[1:end-1, :, :, 1] = κ_H_profile[:, na, na]
#         @. α_H[1:end-1, :, :, 1] = α_H_profile[:, na, na]
        
#         b_E = b_coeff(σ_E, κ_E, α_E, g)
#         c_E = c_coeff(b_E, σ_E, κ_E, α_E)
#         b_H = b_coeff(σ_H, κ_H, α_H, g)
#         c_H = c_coeff(b_H, σ_H, κ_H, α_H)
        
#         Ψ_Ex = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Ψ_Ey = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Ψ_Ez = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3) 
#         Ψ_Hx = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3) 
#         Ψ_Hy = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Ψ_Hz = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)

#         Φ_E = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Φ_H = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)

#         new(
#             g, 
#             location, 
#             size(location)[1],
#             κ_E,
#             κ_H,
#             c_E,
#             b_E,
#             c_H,
#             b_H, 
#             Ψ_Ex,
#             Ψ_Ey,
#             Ψ_Ez,
#             Ψ_Hx,
#             Ψ_Hy,
#             Ψ_Hz,
#             Φ_E,
#             Φ_H 
#         )
#     end
# end

# mutable struct PMLYlow <: CPML
#     grid :: Grid3D
#     location :: CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}
#     thickness :: Int64

#     κ_E :: Array{Float64, 4}
#     κ_H :: Array{Float64, 4}

#     c_E :: Array{Float64, 4}
#     b_E :: Array{Float64, 4}
#     c_H :: Array{Float64, 4}
#     b_H :: Array{Float64, 4}

#     Ψ_Ex :: Array{Float64, 4}
#     Ψ_Ey :: Array{Float64, 4}
#     Ψ_Ez :: Array{Float64, 4}
#     Ψ_Hx :: Array{Float64, 4}
#     Ψ_Hy :: Array{Float64, 4}
#     Ψ_Hz :: Array{Float64, 4}


#     Φ_E :: Array{Float64, 4}
#     Φ_H :: Array{Float64, 4}
#     function PMLYlow(g::Grid3D, location::CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}})
#         σ_E_profile = σ_profile(Array(size(location)[2]-0.5:-1.:0.), size(location)[2], g.Δx) 
#         κ_E_profile = κ_profile(Array(size(location)[2]-0.5:-1.:0.), size(location)[2])
#         α_E_profile = α_profile(Array(size(location)[2]-0.5:-1.:0.), size(location)[2])
        
        
        
#         σ_E = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         κ_E = fill(1., size(location)[1], size(location)[2], size(location)[3], 3)
#         α_E = fill(1e-8, size(location)[1], size(location)[2], size(location)[3], 3)
        
#         na = [CartesianIndex()]
#         @. σ_E[:, :, :, 2] = σ_E_profile[na, :, na]
#         @. κ_E[:, :, :, 2] = κ_E_profile[na, :, na]
#         @. α_E[:, :, :, 2] = α_E_profile[na, :, na]
        

#         σ_H_profile = σ_profile(Array(size(location)[2]-1.:-1.:1.), size(location)[2], g.Δx) 
#         κ_H_profile = κ_profile(Array(size(location)[2]-1:-1.:1.), size(location)[2])
#         α_H_profile = α_profile(Array(size(location)[2]-1:-1.:1.), size(location)[2])

#         σ_H = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         κ_H = fill(1., size(location)[1], size(location)[2], size(location)[3], 3)
#         α_H = fill(1e-8, size(location)[1], size(location)[2], size(location)[3], 3)

#         @. σ_H[:, 1:end-1, :, 2] = σ_H_profile[na, :, na]
#         @. κ_H[:, 1:end-1, :, 2] = κ_H_profile[na, :, na]
#         @. α_H[:, 1:end-1, :, 2] = α_H_profile[na, :, na]
        
#         b_E = b_coeff(σ_E, κ_E, α_E, g)
#         c_E = c_coeff(b_E, σ_E, κ_E, α_E)
#         b_H = b_coeff(σ_H, κ_H, α_H, g)
#         c_H = c_coeff(b_H, σ_H, κ_H, α_H)
        
#         Ψ_Ex = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Ψ_Ey = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Ψ_Ez = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3) 
#         Ψ_Hx = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3) 
#         Ψ_Hy = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Ψ_Hz = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)

#         Φ_E = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Φ_H = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)

#         new(
#             g, 
#             location, 
#             size(location)[2],
#             κ_E,
#             κ_H,
#             c_E,
#             b_E,
#             c_H,
#             b_H, 
#             Ψ_Ex,
#             Ψ_Ey,
#             Ψ_Ez,
#             Ψ_Hx,
#             Ψ_Hy,
#             Ψ_Hz,
#             Φ_E,
#             Φ_H 
#         )
#     end
# end

# mutable struct PMLYhigh <: CPML
#     grid :: Grid3D
#     location :: CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}
#     thickness :: Int64

#     κ_E :: Array{Float64, 4}
#     κ_H :: Array{Float64, 4}

#     c_E :: Array{Float64, 4}
#     b_E :: Array{Float64, 4}
#     c_H :: Array{Float64, 4}
#     b_H :: Array{Float64, 4}

#     Ψ_Ex :: Array{Float64, 4}
#     Ψ_Ey :: Array{Float64, 4}
#     Ψ_Ez :: Array{Float64, 4}
#     Ψ_Hx :: Array{Float64, 4}
#     Ψ_Hy :: Array{Float64, 4}
#     Ψ_Hz :: Array{Float64, 4}


#     Φ_E :: Array{Float64, 4}
#     Φ_H :: Array{Float64, 4}
#     function PMLYhigh(g::Grid3D, location::CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}})
#         σ_E_profile = σ_profile(Array(0.5:1.:size(location)[2]), size(location)[2], g.Δx) 
#         κ_E_profile = κ_profile(Array(0.5:1.:size(location)[2]), size(location)[2])
#         α_E_profile = α_profile(Array(0.5:1.:size(location)[2]), size(location)[2])
        
#         σ_E = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         κ_E = fill(1., size(location)[1], size(location)[2], size(location)[3], 3)
#         α_E = fill(1e-8, size(location)[1], size(location)[2], size(location)[3], 3)
        
#         na = [CartesianIndex()]
#         @. σ_E[:, :, :, 2] = σ_E_profile[na, :, na]
#         @. κ_E[:, :, :, 2] = κ_E_profile[na, :, na]
#         @. α_E[:, :, :, 2] = α_E_profile[na, :, na]
        
#         σ_H_profile = σ_profile(Array(1.:1.:size(location)[2]-1), size(location)[2], g.Δx) 
#         κ_H_profile = κ_profile(Array(1.:1.:size(location)[2]-1), size(location)[2])
#         α_H_profile = α_profile(Array(1.:1.:size(location)[2]-1), size(location)[2])


#         σ_H = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         κ_H = fill(1., size(location)[1], size(location)[2], size(location)[3], 3)
#         α_H = fill(1e-8, size(location)[1], size(location)[2], size(location)[3], 3)

#         @. σ_H[:, 1:end-1, :, 2] = σ_H_profile[na, :, na]
#         @. κ_H[:, 1:end-1, :, 2] = κ_H_profile[na, :, na]
#         @. α_H[:, 1:end-1, :, 2] = α_H_profile[na, :, na]
        
#         b_E = b_coeff(σ_E, κ_E, α_E, g)
#         c_E = c_coeff(b_E, σ_E, κ_E, α_E)
#         b_H = b_coeff(σ_H, κ_H, α_H, g)
#         c_H = c_coeff(b_H, σ_H, κ_H, α_H)
        
#         Ψ_Ex = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Ψ_Ey = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Ψ_Ez = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3) 
#         Ψ_Hx = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3) 
#         Ψ_Hy = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Ψ_Hz = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)

#         Φ_E = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Φ_H = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)

#         new(
#             g, 
#             location, 
#             size(location)[2],
#             κ_E,
#             κ_H,
#             c_E,
#             b_E,
#             c_H,
#             b_H, 
#             Ψ_Ex,
#             Ψ_Ey,
#             Ψ_Ez,
#             Ψ_Hx,
#             Ψ_Hy,
#             Ψ_Hz,
#             Φ_E,
#             Φ_H 
#         )
#     end
# end
# mutable struct PMLZlow <: CPML
#     grid :: Grid3D
#     location :: CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}
#     thickness :: Int64

#     κ_E :: Array{Float64, 4}
#     κ_H :: Array{Float64, 4}

#     c_E :: Array{Float64, 4}
#     b_E :: Array{Float64, 4}
#     c_H :: Array{Float64, 4}
#     b_H :: Array{Float64, 4}

#     Ψ_Ex :: Array{Float64, 4}
#     Ψ_Ey :: Array{Float64, 4}
#     Ψ_Ez :: Array{Float64, 4}
#     Ψ_Hx :: Array{Float64, 4}
#     Ψ_Hy :: Array{Float64, 4}
#     Ψ_Hz :: Array{Float64, 4}


#     Φ_E :: Array{Float64, 4}
#     Φ_H :: Array{Float64, 4}
#     function PMLZlow(g::Grid3D, location::CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}) 
#         σ_E_profile = σ_profile(Array(size(location)[3]-0.5:-1.:0.), size(location)[3], g.Δx) 
#         κ_E_profile = κ_profile(Array(size(location)[3]-0.5:-1.:0.), size(location)[3])
#         α_E_profile = α_profile(Array(size(location)[3]-0.5:-1.:0.), size(location)[3])

#         σ_E = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         κ_E = fill(1., size(location)[1], size(location)[2], size(location)[3], 3)
#         α_E = fill(1e-8, size(location)[1], size(location)[2], size(location)[3], 3)
        
#         na = [CartesianIndex()]
#         @. σ_E[:, :, :, 3] = σ_E_profile[na, na, :]
#         @. κ_E[:, :, :, 3] = κ_E_profile[na, na, :]
#         @. α_E[:, :, :, 3] = α_E_profile[na, na, :]

#         σ_H_profile = σ_profile(Array(size(location)[3]-1.:-1.:1.), size(location)[3], g.Δx) 
#         κ_H_profile = κ_profile(Array(size(location)[3]-1:-1.:1.), size(location)[3])
#         α_H_profile = α_profile(Array(size(location)[3]-1:-1.:1.), size(location)[3])

#         σ_H = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         κ_H = fill(1., size(location)[1], size(location)[2], size(location)[3], 3)
#         α_H = fill(1e-8, size(location)[1], size(location)[2], size(location)[3], 3)

#         @. σ_H[:, :, 1:end-1, 3] = σ_H_profile[na, na, :]
#         @. κ_H[:, :, 1:end-1, 3] = κ_H_profile[na, na, :]
#         @. α_H[:, :, 1:end-1, 3] = α_H_profile[na, na, :]
        
#         b_E = b_coeff(σ_E, κ_E, α_E, g)
#         c_E = c_coeff(b_E, σ_E, κ_E, α_E)
#         b_H = b_coeff(σ_H, κ_H, α_H, g)
#         c_H = c_coeff(b_H, σ_H, κ_H, α_H)
        
#         Ψ_Ex = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Ψ_Ey = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Ψ_Ez = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3) 
#         Ψ_Hx = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3) 
#         Ψ_Hy = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Ψ_Hz = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)

#         Φ_E = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Φ_H = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)

#         new(
#             g, 
#             location, 
#             size(location)[3],
#             κ_E,
#             κ_H,
#             c_E,
#             b_E,
#             c_H,
#             b_H, 
#             Ψ_Ex,
#             Ψ_Ey,
#             Ψ_Ez,
#             Ψ_Hx,
#             Ψ_Hy,
#             Ψ_Hz,
#             Φ_E,
#             Φ_H 
#         )
#     end
# end

# mutable struct PMLZhigh <: CPML
#     grid :: Grid3D
#     location :: CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}
#     thickness :: Int64

#     κ_E :: Array{Float64, 4}
#     κ_H :: Array{Float64, 4}

#     c_E :: Array{Float64, 4}
#     b_E :: Array{Float64, 4}
#     c_H :: Array{Float64, 4}
#     b_H :: Array{Float64, 4}

#     Ψ_Ex :: Array{Float64, 4}
#     Ψ_Ey :: Array{Float64, 4}
#     Ψ_Ez :: Array{Float64, 4}
#     Ψ_Hx :: Array{Float64, 4}
#     Ψ_Hy :: Array{Float64, 4}
#     Ψ_Hz :: Array{Float64, 4}


#     Φ_E :: Array{Float64, 4}
#     Φ_H :: Array{Float64, 4}
#     function PMLZhigh(g::Grid3D, location::CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}})
#         σ_E_profile = σ_profile(Array(0.5:1.:size(location)[3]), size(location)[3], g.Δx) 
#         κ_E_profile = κ_profile(Array(0.5:1.:size(location)[3]), size(location)[3])
#         α_E_profile = α_profile(Array(0.5:1.:size(location)[3]), size(location)[3])

#         σ_E = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         κ_E = fill(1., size(location)[1], size(location)[2], size(location)[3], 3)
#         α_E = fill(1e-8, size(location)[1], size(location)[2], size(location)[3], 3)
#         na = [CartesianIndex()]

#         @. σ_E[:, :, :, 3] = σ_E_profile[na, na, :]
#         @. κ_E[:, :, :, 3] = κ_E_profile[na, na, :]
#         @. α_E[:, :, :, 3] = α_E_profile[na, na, :]
        
#         σ_H_profile = σ_profile(Array(1.:1.:size(location)[3]-1), size(location)[3], g.Δx) 
#         κ_H_profile = κ_profile(Array(1.:1.:size(location)[3]-1), size(location)[3])
#         α_H_profile = α_profile(Array(1.:1.:size(location)[3]-1), size(location)[3])

#         σ_H = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         κ_H = fill(1., size(location)[1], size(location)[2], size(location)[3], 3)
#         α_H = fill(1e-8, size(location)[1], size(location)[2], size(location)[3], 3)
       
#         @. σ_H[:, :, 1:end-1, 3] = σ_H_profile[na, na, :]
#         @. κ_H[:, :, 1:end-1, 3] = κ_H_profile[na, na, :]
#         @. α_H[:, :, 1:end-1, 3] = α_H_profile[na, na, :]
        

#         b_E = b_coeff(σ_E, κ_E, α_E, g)
#         c_E = c_coeff(b_E, σ_E, κ_E, α_E)
#         b_H = b_coeff(σ_H, κ_H, α_H, g)
#         c_H = c_coeff(b_H, σ_H, κ_H, α_H)
        
#         Ψ_Ex = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Ψ_Ey = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Ψ_Ez = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3) 
#         Ψ_Hx = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3) 
#         Ψ_Hy = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Ψ_Hz = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)

#         Φ_E = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)
#         Φ_H = zeros(Float64, size(location)[1], size(location)[2], size(location)[3], 3)

#         new(
#             g, 
#             location, 
#             size(location)[3],
#             κ_E,
#             κ_H,
#             c_E,
#             b_E,
#             c_H,
#             b_H, 
#             Ψ_Ex,
#             Ψ_Ey,
#             Ψ_Ez,
#             Ψ_Hx,
#             Ψ_Hy,
#             Ψ_Hz,
#             Φ_E,
#             Φ_H 
#         )
#     end
# end


# function update_ϕ_H!(PML::CPML, F::Fields3D)
#     PML.Ψ_Hx .*= PML.b_H
#     PML.Ψ_Hy .*= PML.b_H
#     PML.Ψ_Hz .*= PML.b_H    

#     Ex = F.Ex[PML.location]
#     Ey = F.Ey[PML.location]
#     Ez = F.Ez[PML.location]
#     @inbounds begin
#         @. PML.Ψ_Hx[:, 1:end-1, :, 2] += (Ez[:, 2:end, :] - Ez[:, 1:end-1, :]) * PML.c_H[:, 1:end-1, :, 2]
#         @. PML.Ψ_Hx[:, :, 1:end-1, 3] += (Ey[:, :, 2:end] - Ey[:, :, 1:end-1]) * PML.c_H[:, :, 1:end-1, 3]

#         @. PML.Ψ_Hy[:, :, 1:end-1, 3] += (Ex[:, :, 2:end] - Ex[:, :, 1:end-1]) * PML.c_H[:, :, 1:end-1, 3]
#         @. PML.Ψ_Hy[1:end-1, :, :, 1] += (Ez[2:end, :, :] - Ez[1:end-1, :, :]) * PML.c_H[1:end-1, :, :, 1]

#         @. PML.Ψ_Hz[1:end-1, :, :, 1] += (Ey[2:end, :, :] - Ey[1:end-1, :, :]) * PML.c_H[1:end-1, :, :, 1]
#         @. PML.Ψ_Hz[:, 1:end-1, :, 2] += (Ex[:, 2:end, :] - Ex[:, 1:end-1, :]) * PML.c_H[:, 1:end-1, :, 2]

#         @. PML.Φ_H[:, :, :, 1] = PML.Ψ_Hx[:, :, :, 2] -  PML.Ψ_Hx[:, :, :, 3]
#         @. PML.Φ_H[:, :, :, 2] = PML.Ψ_Hy[:, :, :, 3] -  PML.Ψ_Hy[:, :, :, 1]
#         @. PML.Φ_H[:, :, :, 3] = PML.Ψ_Hz[:, :, :, 1] -  PML.Ψ_Hz[:, :, :, 2]
#     end
# end

# function update_ϕ_E!(PML::CPML, F::Fields3D)
#     PML.Ψ_Ex .*= PML.b_E
#     PML.Ψ_Ey .*= PML.b_E
#     PML.Ψ_Ez .*= PML.b_E    

#     Hx = F.Hx[PML.location]
#     Hy = F.Hy[PML.location]
#     Hz = F.Hz[PML.location]
#     @inbounds begin
#         @. PML.Ψ_Ex[:, 2:end, :, 2] += (Hz[:, 2:end, :] - Hz[:, 1:end-1, :]) * PML.c_E[:, 2:end, :, 2]
#         @. PML.Ψ_Ex[:, :, 2:end, 3] += (Hy[:, :, 2:end] - Hy[:, :, 1:end-1]) * PML.c_E[:, :, 2:end, 3]

#         @. PML.Ψ_Ey[:, :, 2:end, 3] += (Hx[:, :, 2:end] - Hx[:, :, 1:end-1]) * PML.c_E[:, :, 2:end, 3]
#         @. PML.Ψ_Ey[2:end, :, :, 1] += (Hz[2:end, :, :] - Hz[1:end-1, :, :]) * PML.c_E[2:end, :, :, 1]

#         @. PML.Ψ_Ez[2:end, :, :, 1] += (Hy[2:end, :, :] - Hy[1:end-1, :, :]) * PML.c_E[2:end, :, :, 1]
#         @. PML.Ψ_Ez[:, 2:end, :, 2] += (Hx[:, 2:end, :] - Hx[:, 1:end-1, :]) * PML.c_E[:, 2:end, :, 2]

#         @. PML.Φ_E[:, :, :, 1] = PML.Ψ_Ex[:, :, :, 2] -  PML.Ψ_Ex[:, :, :, 3]
#         @. PML.Φ_E[:, :, :, 2] = PML.Ψ_Ey[:, :, :, 3] -  PML.Ψ_Ey[:, :, :, 1]
#         @. PML.Φ_E[:, :, :, 3] = PML.Ψ_Ez[:, :, :, 1] -  PML.Ψ_Ez[:, :, :, 2]
#     end
# end

# function updateE!(F::Fields3D, g::Grid3D, PML::CPML)
#     @. F.Ex[PML.location] += g.Δt/g.Δx/ϵ_0 * PML.Φ_E[:, :, :, 1]
#     @. F.Ey[PML.location] += g.Δt/g.Δx/ϵ_0 * PML.Φ_E[:, :, :, 2]
#     @. F.Ez[PML.location] += g.Δt/g.Δx/ϵ_0 * PML.Φ_E[:, :, :, 3]
# end

# function updateH!(F::Fields3D, g::Grid3D, PML::CPML)
#     @. F.Hx[PML.location] -= g.Δt/g.Δx/μ_0 * PML.Φ_H[:, :, :, 1]
#     @. F.Hy[PML.location] -= g.Δt/g.Δx/μ_0 * PML.Φ_H[:, :, :, 2]
#     @. F.Hz[PML.location] -= g.Δt/g.Δx/μ_0 * PML.Φ_H[:, :, :, 3]
# end


# maybe order matters more than anticipatep, maybe first zlow, than ylow, xhigh, yhigh, xlow, zhigh

# function update_ψ!(PML::PMLYlow_new, F::Fields3D, g::Grid3D)
#     for k in 1:g.SizeZ-1
#         for j in 2:PML.thickness
#             for i in 1:g.SizeX-1
#                 PML.Ψ_Exy1[i, j, k] = be_y[j] * PML.Ψ_Exy1[i, j, k] + ce_y[j] * (F.Hz[i, j, k] - F.Hz[i, j - 1, k])/g.Δy
#                 PML.Ψ_Ezy1[i, j, k] = be_y[j] * PML.Ψ_Ezy1[i, j, k] + ce_y[j] * (F.Hx[i, j, k] - F.Hx[i, j - 1, k])/g.Δy
#     end;end;end
# end

# function update_ψ!(PML::PMLYhigh_new, F::Fields3D, g::Grids3D)
#     for k in 1:g.SizeZ-1
#         for j in (g.SizeY-PML.thickness):g.SizeY-1
#             for k in 1:g.SizeX-1
#                 PML.Ψ_Exy2 = [i, j, k] = be_y[j] * PML.Ψ_Exy2
#     end;end;end
# end