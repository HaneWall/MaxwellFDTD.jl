abstract type Detector end

struct PointDetector <: Detector
    location :: CartesianIndex
    t_step_start :: Int64
    t_step_end :: Int64
    Ez :: Vector{Float64}
    function PointDetector(location::CartesianIndex, t_step_start::Int64, t_step_end::Int64)
        new(location, t_step_start, t_step_end, zeros(Float64, (t_step_end-t_step_start+1)))
    end
end

struct LineDetector <: Detector
    location :: CartesianIndices{1, Tuple{UnitRange{Int64}}}
    t_step_start :: Int64
    t_step_end :: Int64
    Ez :: Array{Float64, 2}
    function LineDetector(location::CartesianIndices, t_step_start::Int64, t_step_end::Int64)
        new(location, t_step_start, t_step_end, zeros(Float64, (t_step_end-t_step_start+1, length(location))))
    end
end

function safeE!(D::PointDetector, F::Fields1D, timestep::Int64)
    if timestep >= D.t_step_start && timestep <= D.t_step_end
        D.Ez[timestep-D.t_step_start+1] = F.Ez[D.location]
    end
end

function safeE!(D::LineDetector, F::Fields1D, timestep::Int64)
    if timestep >= D.t_step_start && timestep <= D.t_step_end
        D.Ez[timestep-D.t_step_start+1, :] = F.Ez[D.location]
    end
end