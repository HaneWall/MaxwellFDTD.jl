abstract type Detector end

struct PointDetector <: Detector
    location :: CartesianIndex
    t_step_start :: Int64
    t_step_end :: Int64
    Ez :: Vector{Float64}
    Pz :: Vector{Float64}
    PzNl :: Vector{Float64}
    Jz :: Vector{Float64}
    function PointDetector(location::CartesianIndex, t_step_start::Int64, t_step_end::Int64)
        new(location, t_step_start, t_step_end, 
            zeros(Float64, (t_step_end-t_step_start+1)),
            zeros(Float64, (t_step_end-t_step_start+1)),
            zeros(Float64, (t_step_end-t_step_start+1)),
            zeros(Float64, (t_step_end-t_step_start+1))
            )
    end
end

function safeE!(D::PointDetector, F::Fields1D, timestep::Int64)
    if timestep >= D.t_step_start && timestep <= D.t_step_end
        D.Ez[timestep-D.t_step_start+1] = F.Ez[D.location]
    end
end

function safeP!(D::PointDetector, MF::MaterialFields1D, timestep::Int64)
    if timestep >= D.t_step_start && timestep <= D.t_step_end
        D.Pz[timestep-D.t_step_start+1] = MF.Pz[D.location]
    end
end

function safeJ!(D::PointDetector, MF::MaterialFields1D, timestep::Int64)
    if timestep >= D.t_step_start && timestep <= D.t_step_end
        D.Jz[timestep-D.t_step_start+1] = MF.Jz[D.location]
    end
end

function safePNl!(D::PointDetector, MF::MaterialFields1D, timestep::Int64)
    if timestep >= D.t_step_start && timestep <= D.t_step_end
        D.PzNl[timestep-D.t_step_start+1] = MF.PzNl[D.location]
    end
end

struct LineDetector <: Detector
    location :: CartesianIndices{1, Tuple{UnitRange{Int64}}}
    t_step_start :: Int64
    t_step_end :: Int64
    Ez :: Array{Float64, 2}
    Pz :: Array{Float64, 2}
    Jz :: Array{Float64, 2}
    function LineDetector(location::CartesianIndices, t_step_start::Int64, t_step_end::Int64)
        new(location, t_step_start, t_step_end, 
            zeros(Float64, (t_step_end-t_step_start+1, length(location))), 
            zeros(Float64, (t_step_end-t_step_start+1, length(location))), 
            zeros(Float64, (t_step_end-t_step_start+1, length(location)))
            )
    end
end

function safeE!(D::LineDetector, F::Fields1D, timestep::Int64)
    if timestep >= D.t_step_start && timestep <= D.t_step_end
        D.Ez[timestep-D.t_step_start+1, :] = F.Ez[D.location]
    end
end

function safeP!(D::LineDetector, MF::MaterialFields1D, timestep::Int64)
    if timestep >= D.t_step_start && timestep <= D.t_step_end
        D.Pz[timestep-D.t_step_start+1, :] = MF.Pz[D.location]
    end
end

function safeJ!(D::LineDetector, MF::MaterialFields1D, timestep::Int64)
    if timestep >= D.t_step_start && timestep <= D.t_step_end
        D.Jz[timestep-D.t_step_start+1, :] = MF.Jz[D.location]
    end
end

function safePNl!(D::LineDetector, MF::MaterialFields1D, timestep::Int64)
end