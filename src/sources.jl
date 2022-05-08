abstract type Source end 

struct GaussianPointSource <: Source 
    location :: CartesianIndex
    soft :: Bool
    sf_left :: Bool 
    sf_right :: Bool
    amplitude :: Float64
    t_step_peak :: Int64
    t_width :: Float64
end

function sourceE!(S::GaussianPointSource, F::Fields1D, timestep::Int64)
    if S.soft
        F.Ez[S.location] += S.amplitude * exp(-(timestep - S.t_step_peak)^2/S.t_width)
    else
        F.Ez[S.location] = S.amplitude * exp(-(timestep - S.t_step_peak)^2/S.t_width)
    end
end

function sourceH!(S::GaussianPointSource, F::Fields1D, timestep::Int64)
    if S.soft && S.sf_left
        F.Hy[S.location - CartesianIndex((1,))] -= S.amplitude * exp(-(timestep - (S.t_step_peak+1))^2/S.t_width)/377.
    elseif S.soft && S.sf_right
        F.Hy[S.location + CartesianIndex((1,))] -= S.amplitude * exp(-(timestep - (S.t_step_peak-1))^2/S.t_width)/377.
    end
end