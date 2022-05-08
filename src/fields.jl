abstract type Field end

mutable struct Fields1D <: Field
    Hy :: Vector{Float64}
    Ez :: Vector{Float64}
    function Fields1D(SizeX::Int64)
        new(zeros(SizeX), zeros(SizeX))
    end
end
