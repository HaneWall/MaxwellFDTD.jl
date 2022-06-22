abstract type Detector end

struct PointDetector <: Detector
    location :: CartesianIndex
    t_step_start :: Int64
    t_step_end :: Int64
    Ez :: Vector{Float64}
    Pz :: Vector{Float64}
    PzNl :: Vector{Float64}
    Jz :: Vector{Float64}
    J_Tunnel :: Vector{Float64}
    J_Bound :: Vector{Float64}
    J_Free :: Vector{Float64}
    Γ_ADK :: Vector{Float64}
    function PointDetector(location::CartesianIndex, t_step_start::Int64, t_step_end::Int64)
        new(location, t_step_start, t_step_end, 
            zeros(Float64, (t_step_end-t_step_start+1)),
            zeros(Float64, (t_step_end-t_step_start+1)),
            zeros(Float64, (t_step_end-t_step_start+1)),
            zeros(Float64, (t_step_end-t_step_start+1)),
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

function safeJ_bound!(D::PointDetector, MF::MaterialFields1D, timestep::Int64)
    if timestep >= D.t_step_start && timestep <= D.t_step_end
        D.J_Bound[timestep-D.t_step_start+1] = MF.Jz_bound[D.location]
    end
end

function safeJ_free!(D::PointDetector, MF::MaterialFields1D, timestep::Int64)
    if timestep >= D.t_step_start && timestep <= D.t_step_end
        D.J_Free[timestep-D.t_step_start+1] = MF.Jz_free[D.location]
    end
end

function safeJ_tunnel!(D::PointDetector, MF::MaterialFields1D, timestep::Int64)
    if timestep >= D.t_step_start && timestep <= D.t_step_end
        D.J_Tunnel[timestep-D.t_step_start+1] = MF.Jz_tunnel[D.location]
    end
end

function safeΓ_ADK!(D::PointDetector, MF::MaterialFields1D, timestep::Int64)
    if timestep >= D.t_step_start && timestep <= D.t_step_end
        D.Γ_ADK[timestep-D.t_step_start+1] = MF.Γ_ADK[D.location]
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

function safeJ_tunnel!(D::LineDetector, MF::MaterialFields1D, timestep::Int64)
end

function safeJ_bound!(D::LineDetector, MF::MaterialFields1D, timestep::Int64)
end

function safeJ_free!(D::LineDetector, MF::MaterialFields1D, timestep::Int64)
end

function safeΓ_ADK!(D::LineDetector, MF::MaterialFields1D, timestep::Int64)
end

function safePNl!(D::LineDetector, MF::MaterialFields1D, timestep::Int64)
end


# ------------------------------------------------------------------------------------------------
# Detectors for 3D geometry
# ------------------------------------------------------------------------------------------------


struct ZSliceDetector <: Detector
    location :: CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}
    t_step_start :: Int64
    t_step_end :: Int64
    Ex :: Array{Float64, 3}
    Ey :: Array{Float64, 3}
    Ez :: Array{Float64, 3}
    function ZSliceDetector(location::CartesianIndices, t_step_start::Int64, t_step_end::Int64)
        new(location, t_step_start, t_step_end, 
            zeros(Float64, t_step_end-t_step_start+1, size(location)[1], size(location)[2]), 
            zeros(Float64, t_step_end-t_step_start+1, size(location)[1], size(location)[2]), 
            zeros(Float64, t_step_end-t_step_start+1, size(location)[1], size(location)[2])
            )
    end
end

function safeE!(D::ZSliceDetector, F::Fields3D, timestep)
    if timestep >= D.t_step_start && timestep <= D.t_step_end
        D.Ex[timestep-D.t_step_start+1, :, :] = F.Ex[D.location]
        D.Ey[timestep-D.t_step_start+1, :, :] = F.Ey[D.location]
        D.Ez[timestep-D.t_step_start+1, :, :] = F.Ez[D.location]
    end
end

struct YSliceDetector <: Detector
    location :: CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}
    t_step_start :: Int64
    t_step_end :: Int64
    Ex :: Array{Float64, 3}
    Ey :: Array{Float64, 3}
    Ez :: Array{Float64, 3}
    function YSliceDetector(location::CartesianIndices, t_step_start::Int64, t_step_end::Int64)
        new(location, t_step_start, t_step_end, 
            zeros(Float64, t_step_end-t_step_start+1, size(location)[1], size(location)[3]), 
            zeros(Float64, t_step_end-t_step_start+1, size(location)[1], size(location)[3]), 
            zeros(Float64, t_step_end-t_step_start+1, size(location)[1], size(location)[3])
            )
    end
end

function safeE!(D::YSliceDetector, F::Fields3D, timestep)
    if timestep >= D.t_step_start && timestep <= D.t_step_end
        D.Ex[timestep-D.t_step_start+1, :, :] = F.Ex[D.location]
        D.Ey[timestep-D.t_step_start+1, :, :] = F.Ey[D.location]
        D.Ez[timestep-D.t_step_start+1, :, :] = F.Ez[D.location]
    end
end

struct XSliceDetector <: Detector
    location :: CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}
    t_step_start :: Int64
    t_step_end :: Int64
    Ex :: Array{Float64, 3}
    Ey :: Array{Float64, 3}
    Ez :: Array{Float64, 3}
    function XSliceDetector(location::CartesianIndices, t_step_start::Int64, t_step_end::Int64)
        new(location, t_step_start, t_step_end, 
            zeros(Float64, t_step_end-t_step_start+1, size(location)[2], size(location)[3]), 
            zeros(Float64, t_step_end-t_step_start+1, size(location)[2], size(location)[3]), 
            zeros(Float64, t_step_end-t_step_start+1, size(location)[2], size(location)[3])
            )
    end
end

function safeE!(D::XSliceDetector, F::Fields3D, timestep)
    if timestep >= D.t_step_start && timestep <= D.t_step_end
        D.Ex[timestep-D.t_step_start+1, :, :] = F.Ex[D.location]
        D.Ey[timestep-D.t_step_start+1, :, :] = F.Ey[D.location]
        D.Ez[timestep-D.t_step_start+1, :, :] = F.Ez[D.location]
    end
end

struct BlockDetector <: Detector
    location :: CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}
    t_step_start :: Int64
    t_step_end :: Int64
    Ex :: Array{Float64, 4}
    Ey :: Array{Float64, 4}
    Ez :: Array{Float64, 4}
    function BlockDetector(location::CartesianIndices, t_step_start::Int64, t_step_end::Int64)
        new(location, t_step_start, t_step_end, 
            zeros(Float64, t_step_end-t_step_start+1, size(location)[1], size(location)[2], size(location)[3]), 
            zeros(Float64, t_step_end-t_step_start+1, size(location)[1], size(location)[2], size(location)[3]), 
            zeros(Float64, t_step_end-t_step_start+1, size(location)[1], size(location)[2], size(location)[3])
            )
    end
end

function safeE!(D::BlockDetector, F::Fields3D, timestep)
    if timestep >= D.t_step_start && timestep <= D.t_step_end
        D.Ex[timestep-D.t_step_start+1, :, :, :] = F.Ex[D.location]
        D.Ey[timestep-D.t_step_start+1, :, :, :] = F.Ey[D.location]
        D.Ez[timestep-D.t_step_start+1, :, :, :] = F.Ez[D.location]
    end
end