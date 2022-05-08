abstract type Grid end

struct Grid1D <: Grid
    SizeX :: Int64
    S_c :: Float64
    Δx :: Float64
    MaxTime :: Int64
end

struct Grid2D <: Grid
    SizeX :: Int64
    SizeY :: Int64
    S_c :: Float64
    Δx :: Float64
    Δy :: Float64
    MaxTime :: Int64
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
end