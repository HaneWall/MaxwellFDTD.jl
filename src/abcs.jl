abstract type ABC end # absorbing boundary condition 
abstract type CPML end # convolutional absorbing boundary condition 
abstract type PBC end # periodic boundary Condition 
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

mutable struct PeriodicBoundaryX <: PBC
    grid :: Grid
end

function updateE!(P::PeriodicBoundaryX, F::Fields1D)
    F.Ez[2] = F.Ez[end - 1]
end 

function updateE!(P::PeriodicBoundaryX, F::Fields2D)
    F.Ez[2, :] = F.Ez[end - 1, :]
end

function updateE!(P::PeriodicBoundaryX, F::Fields3D)
    F.Ex[2, :, :] = F.Ex[end - 1, :, :]
    F.Ez[2, :, :] = F.Ez[end - 1, :, :]
    F.Ey[2, :, :] = F.Ey[end - 1, :, :]
end

function updateH!(P::PeriodicBoundaryX, F::Fields1D)
    F.Hy[end - 1] = F.Hy[2]
end 

function updateH!(P::PeriodicBoundaryX, F::Fields2D)
    F.Hx[end - 1, :] = F.Hx[2, :]
    F.Hy[end - 1, :] = F.Hy[2, :]
end

function updateH!(P::PeriodicBoundaryX, F::Fields3D)
    F.Hx[end - 1, :, :] = F.Hx[2, :, :]
    F.Hy[end - 1, :, :] = F.Hy[2, :, :]
    F.Hz[end - 1, :, :] = F.Hz[2, :, :]
end

mutable struct PeriodicBoundaryY <: PBC
    grid :: Grid
end


function updateE!(P::PeriodicBoundaryY, F::Fields2D)
    F.Ez[:, 2] = F.Ez[:, end - 1]
end

function updateE!(P::PeriodicBoundaryY, F::Fields3D)
    F.Ex[:, 2, :] = F.Ex[:, end - 1, :]
    F.Ez[:, 2, :] = F.Ez[:, end - 1, :]
    F.Ey[:, 2, :] = F.Ey[:, end - 1, :]
end


function updateH!(P::PeriodicBoundaryY, F::Fields2D)
    F.Hx[:, end - 1] = F.Hx[:, 2]
    F.Hy[:, end - 1] = F.Hy[:, 2]
end

function updateH!(P::PeriodicBoundaryY, F::Fields3D)
    F.Hx[:, end - 1, :] = F.Hx[:, 2, :]
    F.Hy[:, end - 1, :] = F.Hy[:, 2, :]
    F.Hz[:, end - 1, :] = F.Hz[:, 2, :]
end


mutable struct PeriodicBoundaryZ <: PBC
    grid :: Grid
end

function updateH!(P::PeriodicBoundaryZ, F::Fields3D)
    F.Hx[:, :, end - 1] = F.Hx[:, :, 2]
    F.Hy[:, :, end - 1] = F.Hy[:, :, 2]
    F.Hz[:, :, end - 1] = F.Hz[:, :, 2]
end

function updateE!(P::PeriodicBoundaryZ, F::Fields3D)
    F.Ex[:, :, 2] = F.Ex[:, :, end - 1]
    F.Ez[:, :, 2] = F.Ez[:, :, end - 1]
    F.Ey[:, :, 2] = F.Ey[:, :, end - 1]
end

