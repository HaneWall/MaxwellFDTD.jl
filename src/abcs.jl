abstract type ABC end

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

#=
 These are the update equations in  the three-dimesional Case
=#


struct CPMLABC <: ABC
    thickness :: Int64
    
end