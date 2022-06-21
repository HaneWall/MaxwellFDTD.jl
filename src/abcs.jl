abstract type ABC end
abstract type CPML end

function ABC!(F::Fields1D, g::Grid1D)
    F.Ez[1] = F.Ez[2]
    F.Ez[g.SizeX] = F.Ez[g.SizeX-1]
end

struct LeftSideMurABC <: ABC
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

struct RightSideMurABC <: ABC
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


# 3D - boundaries

function σ_profile(arr::Array{Float64, 1}, thickness::Int64)
    return 40 .* arr.^3 ./ (thickness + 1).^4
end

function b_coeff(σ::Array{Float64}, g::Grid3D; κ::Float64=1.0, a::Float64=1e-8)
    return exp.(-(σ./κ .+ a) .* c_0.*g.Δt./g.Δx)
end

function c_coeff(b::Array{Float64}, σ::Array{Float64}; κ::Float64=1.0, a::Float64=1e-8)
    return (b .- 1).*σ ./ (σ .* κ .+ a .* κ.^2)
end



mutable struct PMLXlow <: CPML
    κ :: Float64
    a :: Float64
    thickness :: Int64

    grid :: Grid3D
    location :: CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}
    σ_E :: Array{Float64, 4}
    σ_H :: Array{Float64, 4} 

    c_E :: Array{Float64, 4}
    b_E :: Array{Float64, 4}
    c_H :: Array{Float64, 4}
    b_H :: Array{Float64, 4}

    Ψ_Ex :: Array{Float64, 4}
    Ψ_Ey :: Array{Float64, 4}
    Ψ_Ez :: Array{Float64, 4}
    Ψ_Hx :: Array{Float64, 4}
    Ψ_Hy :: Array{Float64, 4}
    Ψ_Hz :: Array{Float64, 4}


    Φ_E :: Array{Float64, 4}
    Φ_H :: Array{Float64, 4}
    function PMLXlow(g::Grid3D, location::CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}})
        new(
            1.0,
            1e-8,

        )
    end
end

mutable struct PMLXhigh <: CPML
    κ :: Float64
    a :: Float64
    thickness :: Int64

    grid :: Grid3D
    location :: CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}
    σ_E :: Array{Float64, 4}
    σ_H :: Array{Float64, 4} 

    c_E :: Array{Float64, 4}
    b_E :: Array{Float64, 4}
    c_H :: Array{Float64, 4}
    b_H :: Array{Float64, 4}

    Ψ_Ex :: Array{Float64, 4}
    Ψ_Ey :: Array{Float64, 4}
    Ψ_Ez :: Array{Float64, 4}
    Ψ_Hx :: Array{Float64, 4}
    Ψ_Hy :: Array{Float64, 4}
    Ψ_Hz :: Array{Float64, 4}


    Φ_E :: Array{Float64, 4}
    Φ_H :: Array{Float64, 4}
end

mutable struct PMLYlow <: CPML
    κ :: Float64
    a :: Float64
    thickness :: Int64

    grid :: Grid3D
    location :: CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}
    σ_E :: Array{Float64, 4}
    σ_H :: Array{Float64, 4} 

    c_E :: Array{Float64, 4}
    b_E :: Array{Float64, 4}
    c_H :: Array{Float64, 4}
    b_H :: Array{Float64, 4}

    Ψ_Ex :: Array{Float64, 4}
    Ψ_Ey :: Array{Float64, 4}
    Ψ_Ez :: Array{Float64, 4}
    Ψ_Hx :: Array{Float64, 4}
    Ψ_Hy :: Array{Float64, 4}
    Ψ_Hz :: Array{Float64, 4}


    Φ_E :: Array{Float64, 4}
    Φ_H :: Array{Float64, 4}
end

mutable struct PMLYhigh <: CPML
    κ :: Float64
    a :: Float64
    thickness :: Int64

    grid :: Grid3D
    location :: CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}
    σ_E :: Array{Float64, 4}
    σ_H :: Array{Float64, 4} 

    c_E :: Array{Float64, 4}
    b_E :: Array{Float64, 4}
    c_H :: Array{Float64, 4}
    b_H :: Array{Float64, 4}

    Ψ_Ex :: Array{Float64, 4}
    Ψ_Ey :: Array{Float64, 4}
    Ψ_Ez :: Array{Float64, 4}
    Ψ_Hx :: Array{Float64, 4}
    Ψ_Hy :: Array{Float64, 4}
    Ψ_Hz :: Array{Float64, 4}


    Φ_E :: Array{Float64, 4}
    Φ_H :: Array{Float64, 4}
end
mutable struct PMLZlow <: CPML
    κ :: Float64
    a :: Float64
    thickness :: Int64

    grid :: Grid3D
    location :: CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}
    σ_E :: Array{Float64, 4}
    σ_H :: Array{Float64, 4} 

    c_E :: Array{Float64, 4}
    b_E :: Array{Float64, 4}
    c_H :: Array{Float64, 4}
    b_H :: Array{Float64, 4}

    Ψ_Ex :: Array{Float64, 4}
    Ψ_Ey :: Array{Float64, 4}
    Ψ_Ez :: Array{Float64, 4}
    Ψ_Hx :: Array{Float64, 4}
    Ψ_Hy :: Array{Float64, 4}
    Ψ_Hz :: Array{Float64, 4}


    Φ_E :: Array{Float64, 4}
    Φ_H :: Array{Float64, 4}
end
mutable struct PMLZhigh <: CPML
    κ :: Float64
    a :: Float64
    thickness :: Int64

    grid :: Grid3D
    location :: CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}
    σ_E :: Array{Float64, 4}
    σ_H :: Array{Float64, 4} 

    c_E :: Array{Float64, 4}
    b_E :: Array{Float64, 4}
    c_H :: Array{Float64, 4}
    b_H :: Array{Float64, 4}

    Ψ_Ex :: Array{Float64, 4}
    Ψ_Ey :: Array{Float64, 4}
    Ψ_Ez :: Array{Float64, 4}
    Ψ_Hx :: Array{Float64, 4}
    Ψ_Hy :: Array{Float64, 4}
    Ψ_Hz :: Array{Float64, 4}


    Φ_E :: Array{Float64, 4}
    Φ_H :: Array{Float64, 4}
end


function update_ϕ_H!(PML::CPML, F::Fields3D)
        PML.Ψ_Hx .*= PML.b_H
        PML.Ψ_Hy .*= PML.b_H
        PML.Ψ_Hz .*= PML.b_H    

        Ex = F.Ex[PML.location]
        Ey = F.Ey[PML.location]
        Ez = F.Ez[PML.location]

        @. PML.Ψ_Hx[:, 1:end-1, :, 2] += (Ez[:, 2:end, :] - Ez[:, 1:end-1, :]) * PML.c_H[:, 1:end-1, :, 2]
        @. PML.Ψ_Hx[:, :, 1:end-1, 3] += (Ey[:, :, 2:end] - Ey[:, :, 1:end-1]) * PML.c_H[:, :, 1:end-1, 3]

        @. PML.Ψ_Hy[:, :, 1:end-1, 3] += (Ex[:, :, 2:end] - Ex[:, :, 1:end-1]) * PML.c_H[:, :, 1:end-1, 3]
        @. PML.Ψ_Hy[1:end-1, :, :, 1] += (Ez[2:end, :, :] - Ez[1:end-1, :, :]) * PML.c_H[1:end-1, :, :, 1]

        @. PML.Ψ_Hz[1:end-1, :, :, 1] += (Ey[2:end, :, :] - Ey[1:end-1, :, :]) * PML.c_H[1:end-1, :, :, 1]
        @. PML.Ψ_Hz[:, 1:end-1, :, 2] += (Ex[:, 2:end, :] - Ex[:, 1:end-1, :]) * PML.c_H[:, 1:end-1, :, 2]

        @. PML.Φ_H[:, :, :, 1] = PML.Ψ_Hx[:, :, :, 2] -  PML.Ψ_Hx[:, :, :, 3]
        @. PML.Φ_H[:, :, :, 2] = PML.Ψ_Hy[:, :, :, 3] -  PML.Ψ_Hy[:, :, :, 1]
        @. PML.Φ_H[:, :, :, 3] = PML.Ψ_Hz[:, :, :, 1] -  PML.Ψ_Hz[:, :, :, 2]
end

function update_ϕ_E!(PML::CPML, F::Fields3D)
    PML.Ψ_Ex .*= PML.b_E
    PML.Ψ_Ey .*= PML.b_E
    PML.Ψ_Ez .*= PML.b_E    

    Hx = F.Hx[PML.location]
    Hy = F.Hy[PML.location]
    Hz = F.Hz[PML.location]

    @. PML.Ψ_Ex[:, 2:end, :, 2] += (Hz[:, 2:end, :] - Hz[:, 1:end-1, :]) * PML.c_E[:, 2:end, :, 2]
    @. PML.Ψ_Ex[:, :, 2:end, 3] += (Ey[:, :, 2:end] - Hy[:, :, 1:end-1]) * PML.c_E[:, :, 2:end, 3]

    @. PML.Ψ_Ey[:, :, 2:end, 3] += (Hx[:, :, 2:end] - Hx[:, :, 1:end-1]) * PML.c_E[:, :, 2:end, 3]
    @. PML.Ψ_Ey[2:end, :, :, 1] += (Hz[2:end, :, :] - Hz[1:end-1, :, :]) * PML.c_E[2:end, :, :, 1]

    @. PML.Ψ_Ez[2:end, :, :, 1] += (Hy[2:end, :, :] - Hy[1:end-1, :, :]) * PML.c_E[2:end, :, :, 1]
    @. PML.Ψ_Ez[:, 2:end, :, 2] += (Hx[:, 2:end, :] - Hx[:, 1:end-1, :]) * PML.c_E[:, 2:end, :, 2]

    @. PML.Φ_E[:, :, :, 1] = PML.Ψ_Ex[:, :, :, 2] -  PML.Ψ_Ex[:, :, :, 3]
    @. PML.Φ_E[:, :, :, 2] = PML.Ψ_Ey[:, :, :, 3] -  PML.Ψ_Ey[:, :, :, 1]
    @. PML.Φ_E[:, :, :, 3] = PML.Ψ_Ez[:, :, :, 1] -  PML.Ψ_Ez[:, :, :, 2]
end

function updateE!(F::Fields3D, g::Grid3D, PML::CPML)
    @. F.Ex[PML.location] .+= g.Δt/g.Δx/ϵ_0 * PML.Φ_E[PML.location, 1]
    @. F.Ey[PML.location] .+= g.Δt/g.Δx/ϵ_0 * PML.Φ_E[PML.location, 2]
    @. F.Ez[PML.location] .+= g.Δt/g.Δx/ϵ_0 * PML.Φ_E[PML.location, 3]
end

function updateH!(F::Fields3D, g::Grid3D, PML::CPML)
    @. F.Hx[PML.location] .-= g.Δt/g.Δx/μ_0 * PML.Φ_H[PML.location, 1]
    @. F.Hy[PML.location] .-= g.Δt/g.Δx/μ_0 * PML.Φ_H[PML.location, 2]
    @. F.Hz[PML.location] .-= g.Δt/g.Δx/μ_0 * PML.Φ_H[PML.location, 3]
end