abstract type Source end 


function sigmoid(z::Real)
    return 1.0 ./ (1.0 + exp(-z/8))
end;

function intensity2amplitude(intensity::Float64)
    amplitude = sqrt(intensity * 2/(c_0*ϵ_0))
    return amplitude
end 

function laserfluence(amplitude, τ)
    f = sqrt(π/2) * c_0 * ϵ_0 * τ / 2 * amplitude^2
    return f
end

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
        F.Hy[S.location + CartesianIndex((1,))] += S.amplitude * exp(-(timestep - (S.t_step_peak-1))^2/S.t_width)/377.
    end
end
struct SinusoidalPointSource <: Source
    grid :: Grid1D
    location :: CartesianIndex
    soft :: Bool
    sf_left :: Bool
    amplitude :: Float64
    ppw :: Float64
end

function sourceE!(S::SinusoidalPointSource, F::Fields1D, timestep::Int64)
    if S.soft
        F.Ez[S.location] += tanh(timestep/(100*S.ppw)) * S.amplitude * sin(2*π/S.ppw * (S.grid.S_c * timestep - S.location[1]))
    else
        F.Ez[S.location] = tanh(timestep/(100*S.ppw)) * S.amplitude * sin(2*π/S.ppw * (S.grid.S_c * timestep - S.location[1]))
    end
end

function sourceH!(S::SinusoidalPointSource, F::Fields1D, timestep::Int64)
    if S.soft && S.sf_left
        # TFSF doesn't work yet
        F.Hy[S.location - CartesianIndex((1,))] -= tanh(timestep+0.5/(10*S.ppw))*S.amplitude * sin(2*π/S.ppw * (S.grid.S_c * (timestep+0.5) - (S.location[1]-0.5)))/377.
    end
end


struct GaussianWavePointSource <: Source 
    grid:: Grid1D
    location :: CartesianIndex
    soft :: Bool
    sf_left :: Bool 
    sf_right :: Bool
    amplitude :: Float64
    t_step_peak :: Int64
    t_fwhm :: Float64
    ppw :: Float64
end

function sourceE!(S::GaussianWavePointSource, F::Fields1D, timestep::Int64)
    if S.soft
        F.Ez[S.location] += S.amplitude * exp(-2*log(2)*(((timestep + 0.5) - (S.t_step_peak))*S.grid.Δt/S.t_fwhm)^2) * sin(2*π/S.ppw * (S.grid.S_c * (timestep + 0.5) - (S.location[1] - 0.5)))
        #F.Ez[S.location] += S.amplitude * exp(-(timestep - S.t_step_peak)^2/S.t_width) * sin(2*π/S.ppw * (S.grid.S_c * timestep - S.location[1]))
    else
        F.Ez[S.location] = S.amplitude * exp(-2*log(2)*((timestep - S.t_step_peak)*S.grid.Δt/S.t_fwhm)^2) * sin(2*π/S.ppw * (S.grid.S_c * timestep - S.location[1]))
    end
end

function sourceH!(S::GaussianWavePointSource, F::Fields1D, timestep::Int64)
    if S.soft && S.sf_left
        F.Hy[S.location - CartesianIndex((1,))] -= S.grid.S_c * S.amplitude * exp(-2*log(2)*((timestep - (S.t_step_peak))*S.grid.Δt/S.t_fwhm)^2) * sin(2*π/S.ppw * (S.grid.S_c * (timestep) - S.location[1])) / 376.730
    elseif S.soft && S.sf_right
        F.Hy[S.location + CartesianIndex((1,))] += S.amplitude * exp(-(timestep - (S.t_step_peak-1))^2/S.t_width) * sin(2*π/S.ppw * (S.grid.S_c * timestep - S.location[1])) /377.
    end
end

struct RickerPointSource <: Source
    grid:: Grid1D
    location:: CartesianIndex
    soft:: Bool
    sf_left:: Bool
    amplitude:: Float64
    delay :: Int64
    ppw:: Float64
end

function sourceE!(S::RickerPointSource, F::Fields1D, timestep::Int64)
    arg = π * ((S.grid.S_c *(timestep - (S.delay +1)) - S.location[1]) / S.ppw -1)
    if S.soft && timestep >= S.delay
        F.Ez[S.location] += S.amplitude * (1.0-2.0 * arg^2) * exp(-(arg^2))
    elseif timestep >= S.delay
        F.Ez[S.location] = S.amplitude * (1.0-2.0 * arg^2) * exp(-(arg^2))
    end
end

function sourceH!(S::RickerPointSource, F::Fields1D, timestep::Int64)
    arg = π * ((S.grid.S_c *(timestep+1 - (S.delay +1)) - S.location[1]+1) / S.ppw -1)
    if S.soft && S.sf_left
        F.Hy[S.location - CartesianIndex((1,))] -= S.amplitude * (1.0-2.0 * arg^2) * exp(-(arg^2))/377.
    end
end