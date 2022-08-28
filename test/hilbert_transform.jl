using MaxwellFDTD
using CairoMakie
using FFTW
CairoMakie.activate!(type = "svg")

fc = 20
Δt = 0.00001
t = Array(-1.:Δt:2.)
n_fft = length(t)
f = 1/(n_fft*Δt)*(-(floor((n_fft)/2)):1:(floor((n_fft)/2)))

t_peak = -0.2
t_ref = -0.2
width = 0.01
signal = zeros(Float64, n_fft)
@. signal =  (cos(2*pi*fc*(t-t_peak) + π/4 ) .* exp(-((t-t_peak).^2)/width))
z, phase = analytic_signal(signal, t, t_ref)

fig = Figure(resolution = (800, 800), font="CMU Serif", fontsize = 13)
ax1 = Axis(fig[1, 1])
lines!(ax1, t, abs.(z))
lines!(ax1, t, signal)
ax2 = Axis(fig[2, 1])
lines!(ax2, f./fc, real.(phase))
xlims!(ax2, [0, 2])
fig