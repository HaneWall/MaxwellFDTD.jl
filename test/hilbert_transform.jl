using MaxwellFDTD
using CairoMakie
using FFTW
CairoMakie.activate!(type = "svg")

function analytic_signal2(x::Array{Float64}, t::Array{Float64}, t_ref::Float64)
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

function isolate_around_freq2(f::Array{Float64}, spec::Array{Complex{Float64}}, f_isolated::Float64, Δf::Float64, window::Bool=true)
    n_fft = size(f, 1)
    n_half = ceil(Int64, n_fft/2)
    δf = abs(f[1] - f[2])
    # find freqeuncy in f array
    f_pos = f[n_half:end] 
    f_idx = argmin(abs.(f_pos .- f_isolated))
    length_window = Δf / δf
    f_idx_low_p = n_half + f_idx - floor(Int64, length_window/2) 
    f_idx_high_p = n_half + f_idx + floor(Int64, length_window/2) 
    f_idx_low_n = n_half - f_idx - floor(Int64, length_window/2) 
    f_idx_high_n = n_half - f_idx + floor(Int64, length_window/2)  

    window_spec_p = f_idx_low_p:f_idx_high_p
    window_spec_n = f_idx_low_n:f_idx_high_n
    spec_new = zeros(Complex{Float64}, n_fft)
    spec_new[window_spec_p] .= 2*spec[window_spec_p]
    #spec_new[window_spec_n] .= spec[window_spec_n]
    if window
        spec_new[window_spec_n] .*= hanning(size(window_spec_n, 1))
        spec_new[window_spec_p] .*= hanning(size(window_spec_p, 1))
    end
    return spec_new
end


function get_phase_spectrum(ω::Array{Float64}, spec::Array{Complex{Float64}}, harmonics::Array{Float64}, padded_L::Int64)
    L = size(harmonics, 1)
    ϕs = zeros(Float64, L)
    for (idx, k) in enumerate(harmonics)
        kth_harmonic_idx = argmin(abs.(ω .- k))
        # isolate spectra
        FT_k = isolate_around_freq2(ω, spec, k, 2., false)
        # ϕ = angle.(FT_k)
        # ϕs[idx] = ϕ[kth_harmonic_idx]
        k_t = ifft(ifftshift(FT_k))
        t = Array(1.:1.:padded_L)
        t_ref = Float64(argmax(abs.(k_t)))
        trash, ϕ = analytic_signal2(real(k_t), t, t_ref)
        ϕs[idx] = ϕ[kth_harmonic_idx]
    end        
    return ϕs
end


fc = 20
Δt = 0.00001
t = Array(0.:Δt:5.)
n_fft = length(t)
f = 1/(n_fft*Δt)*(-(floor((n_fft)/2)):1:(floor((n_fft)/2)))

t_peak = 0.8
t_peak_2 = 2.3
t_ref = 15.
width = 0.01
signal = zeros(Float64, n_fft)
@. signal =  (cos(2*pi*fc*(t-t_peak) + π) * exp(-((t-t_peak).^2)/width) + 0.5*cos(2*pi*3*fc *(t-t_peak) + π/4) * exp(-((t-t_peak_2).^2)/width))
FT_signal = fftshift(fft(signal))

n_fft_2 = 2000001
signal_2 = zeros(Float64, n_fft_2)
signal_2[1:500001] = signal
FT_signal_2 = fftshift(fft(signal_2))
t_2 = Array(0.:Δt:20.)
f_2 = 1/(n_fft_2*Δt)*(-(floor((n_fft_2/2))):1:(floor((n_fft_2)/2)))
z, phase = analytic_signal2(signal_2, t_2, 0.8)
#z, phase = analytic_signal(signal_2, t_2, 0.8)
ks = [1.0 + 2. * i for i in 0:1:1]
kth_harmonic_idx = [argmin(abs.(f_2./fc .- k)) for k in ks]
phase[kth_harmonic_idx]

firstharm = ifft(ifftshift(isolate_around_freq2(Array(f_2./fc), FT_signal_2, 1., 2., true)))
pha = get_phase_spectrum(Array(f_2./fc), FT_signal_2, [1., 3.], n_fft_2)

using GLMakie
GLMakie.activate!()
lines!(real.(firstharm))
lines!(signal_2)



# # fig = Figure(resolution = (800, 800), font="CMU Serif", fontsize = 13)
# # ax1 = Axis(fig[1, 1])
# # lines!(ax1, t_2, abs.(z))
# # lines!(ax1, t_2, signal_2)
# # ax2 = Axis(fig[2, 1])
# # lines!(ax2, f_2./fc, real.(phase))
# # xlims!(ax2, [0, 4])
# # fig


# ϕs = get_phase_spectrum(Array(f)./fc, FT_signal_2, ks, n_fft_2)
