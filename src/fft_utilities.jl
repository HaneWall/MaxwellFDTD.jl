using DSP

function mytukey(n::Integer, a::Real, padding::Integer=0, zerophase::Bool=false)
    return DSP.tukey(n, a; padding, zerophase)
end


function analytic_signal(x::Array{Float64}, t::Array{Float64}, t_ref::Float64)
    n_fft = size(x, 1)
    n_half = ceil(Int64, n_fft / 2)
    idx_ref = argmin(abs.(t .- t_ref))
    shift = n_half - idx_ref
    x_shifted = circshift(x, shift)
    X = fft(ifftshift(x_shifted))
    # mask out negative frequencies, double positive frequency spectrum -> exact amplitudes
    X_shifted = fftshift(X)
    negative_frequency_domain = 1:n_half
    X_shifted[negative_frequency_domain] .= 0 + 0im
    X_shifted .*= 2
    # ifft 
    z = ifftshift(ifft(X_shifted))
    Re_z = real(z)
    Im_z = imag(z)
    circshift!(Re_z, shift)
    circshift!(Im_z, shift)
    z = Re_z + 1im*Im_z
    # X2 = X_shifted
    # noise_filter = abs.(X_shifted) .< maximum(abs.(X_shifted)/4)
    # X2[noise_filter] .= 0 + 0im
    ϕ = angle.(X_shifted)
    return [z, ϕ]
end

function isolate_around_freq(f::Array{Float64}, spec::Array{Complex{Float64}}, f_isolated::Float64, Δf::Float64, window::Bool=true)
    n_fft = size(f, 1)
    n_half = Int(ceil(n_fft/2))
    δf = abs(f[1] - f[2])
    # find freqeuncy in f array
    f_pos = f[n_half:end] 
    f_idx = argmin(abs.(f_pos .- f_isolated))
    length_window = Δf / δf
    f_idx_low = Int(f_idx - floor(length_window/2)) + n_half
    f_idx_high = Int(f_idx + floor(length_window/2)) + n_half
    window_spec = f_idx_low:f_idx_high
    spec_new = zeros(Complex{Float64}, n_fft)
    spec_new[window_spec] .= spec[window_spec]
    if window
        spec_new[window_spec] .*= hanning(size(window_spec, 1))
    end
    return spec_new
end

