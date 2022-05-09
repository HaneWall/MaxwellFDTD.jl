abstract type Grid end

struct Grid1D <: Grid
    SizeX :: Int64
    S_c :: Float64
    Δx :: Float64
    MaxTime :: Int64
    Δt :: Float64
    function Grid1D(SizeX::Int64, S_c::Float64, Δx::Float64, MaxTime::Int64)
        new(SizeX, S_c, Δx, MaxTime, S_c*Δx/c_0)
    end
end

struct Grid2D <: Grid
    SizeX :: Int64
    SizeY :: Int64
    S_c :: Float64
    Δx :: Float64
    Δy :: Float64
    MaxTime :: Int64
    Δt :: Float64
    function Grid2D(SizeX::Int64, SizeY::Int64, S_c::Float64, Δx::Float64, Δy::Float64, MaxTime::Int64)
        new(SizeX, SizeY, S_c, Δx, Δy, MaxTime, S_c*Δx/c_0)
    end
end

struct Grid3D <: Grid
    SizeX :: Int64
    SizeY :: Int64
    SizeZ :: Int64
    S_c :: Float64
    Δx :: Float64
    Δy :: Float64
    Δz :: Float64
    MaxTime :: Int64
    Δt :: Float64
    function Grid3D(SizeX::Int64, SizeY::Int64, S_c::Float64, Δx::Float64, Δy::Float64, Δz::Float64, MaxTime::Int64)
        new(SizeX, SizeY, S_c, Δx, Δy, Δz, MaxTime, S_c*Δx/c_0)
    end
end