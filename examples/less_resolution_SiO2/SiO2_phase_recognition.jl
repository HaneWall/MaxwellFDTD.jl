
using MaxwellFDTD
using CPUTime
using FFTW
#using GLMakie
using ProgressBars
using LaTeXStrings
using DSP
using CairoMakie
#GLMakie.activate!()
CairoMakie.activate!(type="svg")

CPUtic()
start = time()

# 1. define grid
SizeX = 13000
courant = 0.5
Δx = 10e-9
MaxTime = 50000

# 1. define grid
g = Grid1D(SizeX, courant, Δx, MaxTime)
t = g.Δt:g.Δt:g.Δt*MaxTime

# varin paramters Si02
ρ_mol_density = 2.2e28
# bound electrons
γ_lorentz = [0.0]
ω_0 = [2.75e16] # this might not work, use 2*π*/g.Δt instead (old varin paper/bachelor thesis) 2.75e16
χ_1 = [1.1025]
χ_2 = [0.0]
#χ_3 = [0.0]
χ_3 = [2.2e-22]

# drude parameters
γ_plasma = 1e15
#γ_plasma = 0.0

# tunnel parameters 
E_gap = 7.5

# source parameters
λ = 2100e-9
ω_central = 2 * π * c_0 / λ
ppw = λ / Δx
t_fwhm = 20e-15 # intensity FWHM
amplitude_pump = intensity2amplitude(12e16) # 12TWcm^-2

λ_probe = 800e-9
ω_probe = 2 * π * c_0 / λ_probe
ppw_probe = λ_probe / Δx
t_fwhm_probe = 45e-15 # intensity FWHM
amplitude_probe = intensity2amplitude(1.5e14)

# 2. define fields that exist everywhere
F = Fields1D(g)
MF = MaterialFields1D(g)

# PML Fields and parameters
PML_Thickness = [500]
F_PML = CPML_Ψ_Fields_1D(g, PML_Thickness)
c_PML = CPML_Parameters_1D(g, PML_Thickness)

# init the media (superposition of different effects, that act inside the medium)
m1 = LorentzMedium1D(g, CartesianIndices((8000:12999,)), 1.0, γ_lorentz, ω_0, χ_1, χ_2, χ_3)
m2 = DrudeMedium1D(g, CartesianIndices((8000:12999,)), γ_plasma, ρ_mol_density)
m3 = TunnelMedium1D(g, CartesianIndices((8000:12999,)), E_gap, ρ_mol_density)

bound_media = [m1]
drude_media = [m2]
tunnel_media = [m3]

bound_media = [m1]
# 4. define grid coefficients that respect ϵ_inf from the media 

c_grid = GridCoefficients1D_w_CPML(g, bound_media, c_PML)
#c_grid = GridCoefficients1D(g, bound_media)
f_grid = FieldIonizationCoefficients1D(g)

# 5. define fields inside the media
LF1 = LorentzFields1D(m1)
LF = [LF1]
DF1 = DrudeFields1D(m2)
DF = [DF1]
TF1 = TunnelFields1D(m3)
TF = [TF1]

# 6. place detectors 
#d1 = LineDetector(CartesianIndices((1:g.SizeX,)), 1, g.MaxTime)
d2 = PointDetector(CartesianIndex((8000,)), 1, g.MaxTime)
#d3 = PointDetector(CartesianIndex((7500,)), 1, g.MaxTime)
d4 = PointDetector(CartesianIndex((5505,)), 1, g.MaxTime)
detectors = [d2, d3, d4]

# 7. place sources 
s0 = GaussianWavePointSource(g, CartesianIndex((5508,)), true, false, false, amplitude_pump, ceil(500e-15 / g.Δt), t_fwhm, ppw)
s1 = GaussianWavePointSource(g, CartesianIndex((5508,)), true, true, false, amplitude_probe, ceil(500e-15 / g.Δt), t_fwhm_probe, ppw_probe)
sources = [s0]

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3), " seconds")
println("Init Complete")
CPUtic()
start = time()




Ê = amplitude_pump
Γ̂ = Γ_ADK(Ê, E_gap * q_0)
#a = 13.39
a = 13.00

for timestep in ProgressBar(1:g.MaxTime)

    # for (m_idx, m) in enumerate(tunnel_media)
    #     updatePlasma!(MF, TF[m_idx], f_grid, F, m)
    #     updateJtunnel!(MF, TF[m_idx], m)
    # end

    for (m_idx, m) in enumerate(tunnel_media)
        updatePlasmaTangent!(MF, TF[m_idx], f_grid, F, m, a, Γ̂, Ê)
        updateJtunnel!(MF, TF[m_idx], m)
    end

    for (m_idx, m) in enumerate(drude_media)
        updateJfree!(MF, DF[m_idx], F, m)
    end

    for (m_idx, m) in enumerate(bound_media)
        updatePNl!(MF, LF[m_idx], F, m)
        updateJbound!(MF, LF[m_idx], m, g)
        updatePbound!(MF, LF[m_idx], m, g)
    end

    updateJ!(MF)

    update_Ψ_H!(F_PML, F, g, c_PML)

    updateH!(F, g, c_grid)

    # for source in sources
    #     sourceH!(source, F, timestep)
    # end

    apply_Ψ_H!(F_PML, F, g, c_PML)

    update_Ψ_E!(F_PML, F, g, c_PML)

    updateE!(F, MF, g, c_grid)

    for source in sources
        sourceE!(source, F, timestep)
    end

    apply_Ψ_E!(F_PML, F, g, c_PML)

    for d in detectors
        #safeΓ_ADK!(d, MF, timestep)
        safeJ_tunnel!(d, MF, timestep)
        safeJ_free!(d, MF, timestep)
        safeJ_bound!(d, MF, timestep)
        safeE!(d, F, timestep)
        #safeP!(d, MF, timestep)
        safeJ!(d, MF, timestep)
        #safePNl!(d, MF, timestep)
    end
end

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3), " seconds")
println("Computation Complete")


S_C = courant
t_arr = Array(25000.0:1.0:35000.0)
t_peak = 500e-15 / g.Δt
t_real = t_arr .* g.Δt

L = length(t_real)
padded_L = 4*nextpow(2, L)
center_help_L_low = floor(Int64, (padded_L - L)/2)
center_help_L_high = floor(Int64, L + center_help_L_low - 1) 
δt = abs(t[2]-t[1])
δf = 1/δt
ω = 2*π*fftshift(fftfreq(padded_L, δf))./ ω_central
window = blackman(L)

n_real = n_bound(γ_lorentz[1], ω_0[1], χ_1[1], ω .* ω_central)[1] 
t_theory = t_s.(1., n_real, 0.0, 0.0)
r_theory = r_s.(1., n_real, 0.0, 0.0)
r_nl_theory = r_s_nl.(1., n_real)
tr_a = t_s(1., 1.45, 0.0, 0.0)
tr_b = t_s.(n_real, 1., 0.0, 0.0)
# j_inj = (χ_injection .* (tr_a .* E_arr).^13) 
# P = (χ_1 .* (tr_a .* E_arr) + χ_3.*(tr_a .* E_arr).^3)


window_idx = 30000:40000
E_measured_bulk = d2.Ez[window_idx]
E_amplitude_spectrum = zeros(Float64, padded_L)
E_amplitude_spectrum[window_idx] = window .* E_measured_bulk
FT_E_bulk = fftshift(fft(E_amplitude_spectrum))

j_brunel = d2.J_Free[window_idx]
j_inj = d2.J_Tunnel[window_idx]
j_kerr = d2.J_Bound[window_idx]
j_all = d2.Jz[window_idx]


signal_j_inj = zeros(Float64, padded_L)
signal_j_inj[window_idx] = window .* j_inj 
FT_J_inj = fftshift(fft(signal_j_inj))

signal_j_brunel = zeros(Float64, padded_L)
signal_j_brunel[window_idx] = window .* j_brunel
FT_J_brunel = fftshift(fft(signal_j_brunel))

signal_j_bound = zeros(Float64, padded_L)
signal_j_bound[window_idx] = window .* j_kerr
FT_J_bound= fftshift(fft(signal_j_bound))

signal_j_all = zeros(Float64, padded_L)
signal_j_all[window_idx] = window .* j_all
FT_J_all= fftshift(fft(signal_j_all)) 

window_ref = 35000:45000
E_measured_reflection = d4.Ez[window_ref]
E_amplitude_spectrum_refl = zeros(Float64, padded_L)
E_amplitude_spectrum_refl[window_ref] = window .* E_measured_reflection
FT_E_bulk_refl = fftshift(fft(E_amplitude_spectrum_refl))

function isolate_around_freq2(f::Array{Float64}, spec::Array{Complex{Float64}}, f_isolated::Float64, Δf::Float64, window::Bool=true)
    n_fft = size(f, 1)
    n_half = ceil(Int64, n_fft/2)
    δf = abs(f[1] - f[2])
    # find freqeuncy in f array
    f_pos = f[n_half:end] 
    f_idx = argmin(abs.(f_pos .- f_isolated))
    length_window = Δf / δf
    f_idx_low = f_idx - floor(Int64, length_window/2) + n_half
    f_idx_high = f_idx + floor(Int64, length_window/2) + n_half
    window_spec = f_idx_low:f_idx_high
    spec_new = zeros(Complex{Float64}, n_fft)
    spec_new[window_spec] .= spec[window_spec]
    if window
        spec_new[window_spec] .*= hanning(size(window_spec, 1))
    end
    return spec_new
end

ks = [1.0 + 2. * i for i in 0:1:6]


function get_phase_spectrum(ω::Array{Float64}, spec::Array{Complex{Float64}}, harmonics::Array{Float64}, padded_L::Int64)
    L = size(harmonics, 1)
    ϕs = zeros(Float64, L)
    for (idx, k) in enumerate(harmonics)
        kth_harmonic_idx = argmin(abs.(ω .- k))
        # isolate spectra
        FT_k = isolate_around_freq2(ω, spec, k, 1.7, true)
        k_t = ifft(ifftshift(FT_k))
        t = Array(1.:1.:padded_L)
        t_ref = Float64(argmax(abs.(k_t)))
        #t_ref = 0.
        trash, ϕ = analytic_signal(real(k_t), t, t_ref)
        ϕs[idx] = ϕ[kth_harmonic_idx]
    end        
    return ϕs
end

function get_phase_spectrum(ω::Array{Float64}, spec::Array{Complex{Float64}}, harmonics::Array{Float64}, padded_L::Int64, t_ref::Float64)
    L = size(harmonics, 1)
    ϕs = zeros(Float64, L)
    for (idx, k) in enumerate(harmonics)
        kth_harmonic_idx = argmin(abs.(ω .- k))
        # isolate spectra
        FT_k = isolate_around_freq2(ω, spec, k, 1.7, false)
        k_t = ifft(ifftshift(FT_k))
        t = Array(1.:1.:padded_L)
        trash, ϕ = analytic_signal(real(k_t), t, t_ref)
        ϕs[idx] = ϕ[kth_harmonic_idx]
    end        
    return ϕs
end

t_ref = 34980.
#t_ref = 2000.
t_ref_2 = 39970.

ϕ_brunel = get_phase_spectrum(Array(ω), FT_J_brunel, ks, padded_L) .+ π/2
ϕ_E_r = get_phase_spectrum(Array(ω), FT_E_bulk_refl, ks, padded_L) .+ π/2
#ϕ_E_r = get_phase_spectrum(Array(2*π*fftshift(fftfreq(10001, δf))./ ω_central), fftshift(fft(d4.Ez[window_ref])), ks, 10001)
ϕ_E_tr = get_phase_spectrum(Array(ω), FT_E_bulk, ks, padded_L).+ π/2
ϕ_Kerr = get_phase_spectrum(Array(ω), FT_J_bound, ks, padded_L).+ π/2
ϕ_Inj = get_phase_spectrum(Array(ω), FT_J_inj, ks, padded_L).+ π/2
ϕ_All = get_phase_spectrum(Array(ω), FT_J_all, ks, padded_L).+ π/2
# fig = Figure(resolution = (600, 600), font="CMU Serif")
# ax1 = Axis(fig[1, 1])
# scatterlines!(ks, ϕ_brunel, label=L"\phi_{Brunel}")
# scatterlines!(ks, ϕ_E_r,label=L"\phi_{Er}")
# scatterlines!(ks, ϕ_E_tr,label=L"\phi_{Etr}")
# scatterlines!(ks, ϕ_Kerr, label=L"\phi_{Kerr}")
# scatterlines!(ks, ϕ_Inj,label=L"\phi_{Inje}")
# scatterlines!(ks, ϕ_All, label=L"\phi_{All}")
# axislegend(ax1, position=:lt)
# ax1.yticks = -3π:π/4:3π
# #ax1.ytickformat = ys -> ["$(y/pi)π" for y in ys]
# xlims!(ax1, [0, 14])
# fig

k_phase_harms = angle.([exp(1im*(2*π/(ppw/n) * 2492)) for n in ks])
k_phase_harms_2 = angle.([exp(1im*(2*π/(ppw/n) * 2*2492)) for n in ks])

using PyPlot
#r = [1, 2, 3, 4, 5, 6, 7]
r = ks .+ 1
ϕ_new = ϕ_E_r .+ k_phase_harms
fig = plt.figure()
ax = fig.add_subplot(111, projection="polar")
l = ax.scatter(ϕ_Inj, r, label=L"J_{Injection}")
#d = ax.scatter(ϕ_new, r, label="k") 
e = ax.scatter(ϕ_brunel, r, label=L"J_{Brunel}")
f = ax.scatter(ϕ_E_r, r, label=L"E_{Reflection}")
p = ax.scatter(ϕ_E_tr, r, label=L"E_{Transmission}")
h = ax.scatter(ϕ_Kerr, r, label=L"J_{Kerr}")
ax.legend(loc="best", bbox_to_anchor=(0.05, 0.3))
ax[:grid](true)
savefig("Phasen_gamma_1e15_new_test_bigger_domainfg_morevac.pdf")

 
spectra_E_reflect = [isolate_around_freq2(Array(ω), FT_E_bulk_refl, k, 1.7) for k in ks]
time_E_reflect = [ifft(ifftshift(spec)) for spec in spectra_E_reflect]
fig_E_r = Makie.Figure(resolution = (800, 800), font="CMU Serif")
ax = Axis(fig_E_r[1, 1])
for (idx, k) in enumerate(1:2:9)
    lines!(ax, real.(time_E_reflect[idx])./maximum(abs.(time_E_reflect[idx])), label="$k")
    axislegend(ax, position=:rt)
end
xlims!(ax, [30000, 50100])
fig_E_r



#kth - Harmonics 
k = 5.
kth_harmonic_idx = Int(argmin(abs.(ω .- k)))
FT_E_k = isolate_around_freq2(Array(ω), FT_E_bulk, k, 1.7, false)
FT_Er_k = isolate_around_freq2(Array(ω), FT_E_bulk_refl, k, 1.7, false)

# ## test isolation of positive freq
# fig_test = Makie.Figure(resolution = (400, 400), font="CMU Serif")
# ax_test = Axis(fig_test[1, 1], yscale =log10, xlabel=L"\omega / \omega_{pump}", yticks = [0,0])
# lines!(ax_test, ω, abs.(FT_E_bulk_refl), linewidth=2.5, color=:black, label=L"E_{Refl}(\omega)")
# lines!(ax_test, ω, abs.(FT_Er_k), linestyle=:dash, linewidth=2.5, color=:red)
# xlims!(ax_test, [0., 13.])
# ylims!(ax_test, [1e4, 1e13])
# axislegend(ax_test, position=:rt)
# save("reflection_spectra_isolation.pdf", fig_test)
# fig_test



FT_j_harm_k = isolate_around_freq2(Array(ω), FT_J_all, k, 1.7, false)
FT_j_inj_k = isolate_around_freq2(Array(ω), FT_J_inj, k, 1.7, false)
FT_j_brunel_k = isolate_around_freq2(Array(ω), FT_J_brunel, k, 1.7, false)
FT_j_kerr_k = isolate_around_freq2(Array(ω), FT_J_bound, k, 1.7, false)

# fig_test = Makie.Figure(resolution = (800, 400), font="CMU Serif")
# ax_test_1 = Axis(fig_test[1, 1], yscale =log10, xlabel=L"\omega / \omega_{pump}", yticks = [0,0])
# ax_test_2 = Axis(fig_test[1, 2], yscale =log10, xlabel=L"\omega / \omega_{pump}", yticks = [0,0])

# lines!(ax_test_1, ω, abs.(FT_J_inj), linewidth=2.5, color=:black, label=L"j_{injection}(\omega)")
# lines!(ax_test_1, ω, abs.(FT_j_inj_k), linestyle=:dash, linewidth=2.5, color=:red)
# xlims!(ax_test_1, [0., 13.])
# ylims!(ax_test_1, [1e6, 1e16])
# axislegend(ax_test_1, position=:rt)
# lines!(ax_test_2, ω, abs.(FT_J_brunel), linewidth=2.5, color=:black, label=L"j_{brunel}(\omega)")
# lines!(ax_test_2, ω, abs.(FT_j_brunel_k), linestyle=:dash, linewidth=2.5, color=:red)
# xlims!(ax_test_2, [0., 13.])
# ylims!(ax_test_2, [1e6, 1e16])

# axislegend(ax_test_2, position=:rt)
# save("medium_spectra_isolation.pdf", fig_test)
# fig_test



#fith harmoncs in time
E_k_t = ifft(ifftshift(FT_E_k))
Er_k_t = ifft(ifftshift(FT_Er_k))
j_k_t = ifft(ifftshift(FT_j_harm_k))
j_inj_k_t = ifft(ifftshift(FT_j_inj_k))
j_brunel_k_t = ifft(ifftshift(FT_j_brunel_k))
j_kerr_k_t = ifft(ifftshift(FT_j_kerr_k))
# getting phaseinformation
t_help = Array(1.:1.:padded_L)
t_ref = 34980.
#t_ref = 2000.
t_ref_2 = 39970.
#t_ref_2 = 40100.
shift = 2*(padded_L/2 - t_ref_2)
z_e, ϕ_E = analytic_signal(real(E_k_t), t_help, t_ref)
z_er, ϕ_Er = analytic_signal(circshift(reverse(real(Er_k_t)), -shift), t_help, t_ref_2)
z_j, ϕ_j = analytic_signal(real(j_k_t), t_help, t_ref)
z_j_inj, ϕ_j_inj = analytic_signal(real(j_inj_k_t), t_help, t_ref)
z_j_brunel, ϕ_j_brunel = analytic_signal(real(j_brunel_k_t), t_help, t_ref)
z_j_kerr, ϕ_j_kerr = analytic_signal(real(j_kerr_k_t), t_help, t_ref)

fig = Makie.Figure(resolution = (800, 800), font="CMU Serif")
ax1 = Axis(fig[1, 1])
lines!(ax1, abs.(E_k_t)./maximum(abs.(E_k_t)), label=L"E_{Bulk}")
lines!(ax1, abs.(Er_k_t)./maximum(abs.(Er_k_t)), label=L"E_{Refl}")
lines!(ax1, abs.(j_k_t)./maximum(abs.(j_k_t)), label=L"J_{all}")
lines!(ax1, abs.(j_inj_k_t)./maximum(abs.(j_inj_k_t)), label=L"J_{J_inj}")
lines!(ax1, abs.(j_brunel_k_t)./maximum(abs.(j_brunel_k_t)), label=L"J_{Brunel}")
lines!(ax1, abs.(j_kerr_k_t)./maximum(abs.(j_kerr_k_t)), linestyle=:dash, label=L"J_{Kerr}")
axislegend(ax1, position=:rt)
xlims!(ax1, [30000, 45000])
ax2 = Axis(fig[2, 1])
lines!(ax2, ω, real(ϕ_E), label=L"E_{Bulk}")
lines!(ax2, ω, real(ϕ_Er), label=L"E_{Refl}")
lines!(ax2, ω, real(ϕ_j), label=L"J_{all}")
lines!(ax2, ω, real(ϕ_j_inj), label=L"J_{inj}")
lines!(ax2, ω, real(ϕ_j_brunel), label=L"J_{Brunel}")
lines!(ax2, ω, real(ϕ_j_kerr), label=L"J_{Kerr}")
#scatter!(k, real(ϕ_j_kerr)[kth_harmonic_idx])
axislegend(ax2, position=:rt)
ylims!(ax2, [-π, π])
ax2.yticks = -π:π/4:π
ax2.ytickformat = ys -> ["$(y/pi)π" for y in ys]
xlims!(ax2, [k-1, k+1])
fig


## plotting Bulk Currents

fig = Makie.Figure(resolution = (800, 400), font= "CMU Serif")
ax1 = Axis(fig[1, 1], xlabel=L"\Delta t", title="$k harmonic in medium",)
lines!(ax1, d2.Ez./maximum(d2.Ez), label=L"E_{Bulk}", color=:black, linewidth=1.5)
#lines!(ax1, real.(j_kerr_k_t)./maximum(abs.(j_kerr_k_t)), label=L"j_{\chi}")
lines!(ax1, real.(j_brunel_k_t)./maximum(abs.(j_brunel_k_t)), label=L"j_{Brunel}")
lines!(ax1, real.(j_inj_k_t)./maximum(abs.(j_inj_k_t)), label=L"j_{Injection}")
axislegend(ax1, position=:rt)
xlims!(ax1, [33500, 39500])
ylims!(ax1, [-1, 1])
fig




## plotting Refl with bulk currents 
## idee: ∫j dt sollte genau auf E_harm liegen 
## für höhere Harmonische nicht mehr der Fall da diese später im Medium enstehen --> anderer Offset 
CairoMakie.activate!(type="svg")
t_det = 39670:40270
t_med = 34680:35280
t_offset = Array(-300*g.Δt:g.Δt:300*g.Δt) .* 1e15
fig = Makie.Figure(resolution = (800, 400), font= "CMU Serif", fontsize=13)
ax1 = Axis(fig[1, 1], xlabel="t in fs", title="$k harmonic")
lines!(ax1, t_offset, d4.Ez[t_det]./maximum(d4.Ez[t_det]), label=L"E_{Refl}", color=:black, linewidth=2.5)
lines!(ax1, t_offset, d2.Ez[t_med]./maximum(d2.Ez[t_med]), label=L"E_{Bulk}", color=:gray, linewidth=2.5)
lines!(ax1, t_offset, real.(Er_k_t[t_det])./maximum(abs.(Er_k_t[t_det])), label=L"E_{Refl, harm}", color=:black, linewidth=1.5, linestyle=:dash)
lines!(ax1, t_offset, real.(E_k_t[t_med])./maximum(abs.(E_k_t[t_med])), label=L"E_{Bulk, harm}", color=:violet, linestyle=:dash)

ylims!(ax1, [-1, 1])
xlims!(ax1, [-300*g.Δt*1e15, 300*g.Δt*1e15])
#lines!(ax2, real.(E_k_t)./maximum(abs.(E_k_t)), label=L"E_{Bulk, Harm}", linestyle=:dash, color=:gray)
lines!(ax1, t_offset, real.(j_brunel_k_t[t_med])./maximum(abs.(j_brunel_k_t[t_med])), label=L"j_{Brunel}", color=:blue, linewidth=1.5)
lines!(ax1, t_offset, real.(j_inj_k_t[t_med])./maximum(abs.(j_inj_k_t[t_med])), label=L"j_{Injection}", color=:red, linewidth=1.5)
#lines!(ax1, t_offset, real.(j_kerr_k_t[t_med])./maximum(abs.(j_kerr_k_t[t_med])), label=L"j_{\chi}", color=:green, linewidth=1.5)
# #xlims!(ax2, [34000, 36000])
axislegend(ax1, position=:rt)
#save("seventh_harmonic_refl_with_currents_even_closer.pdf", fig)
fig






fig = Makie.Figure(resolution = (800, 400), font="CMU Serif")
ax1 = Axis(fig[1, 1], xlabel=L"\Delta t")
lines!(ax1, d2.Ez./maximum(d2.Ez), label=L"E_{Bulk}")
lines!(ax1, d4.Ez./maximum(d4.Ez[window_ref]), label=L"E_{Refl}")
lines!(ax1, real.(E_k_t)./maximum(abs.(E_k_t)), label=L"E_{Bulk, 3}")
#lines!(ax1, circshift(reverse(real.(Er_k_t)./maximum(abs.(Er_k_t))), -shift-4984), label=L"E_{Refl}")
lines!(ax1, real.(Er_k_t)./maximum(abs.(Er_k_t)), label=L"E_{Refl, 3}")
#lines!(ax1, E_amplitude_spectrum_refl./maximum(abs.(E_amplitude_spectrum_refl)), label="Ref", linestyle=:dash)
#lines!(ax1, d4.Ez./maximum(abs.(d4.Ez[35000:45000])), label=L"E_{refl_real}", linestyle=:dash)
#lines!(ax1, abs.(j_k_t)./maximum(abs.(j_k_t)), label=L"J_{all}")
#lines!(ax1, abs.(j_inj_k_t)./maximum(abs.(j_inj_k_t)), label=L"J_{J_inj}")
#lines!(ax1, abs.(j_brunel_k_t)./maximum(abs.(j_brunel_k_t)), label=L"J_{Brunel}")
#lines!(ax1, abs.(j_kerr_k_t)./maximum(abs.(j_kerr_k_t)), linestyle=:dash, label=L"J_{Kerr}")
axislegend(ax1, position=:rt)
ylims!(ax1, [-1, 1])
xlims!(ax1, [33000, 42000])
# ax2 = Axis(fig[2, 1], xlabel=L"\omega / \omega_{pump}")
# lines!(ax2, ω, real(ϕ_E), label=L"E_{Bulk, 3}")
# lines!(ax2, ω, real(ϕ_Er), label=L"E_{Refl}")
# #lines!(ax2, ω, real(ϕ_j), label=L"J_{all}")
# #lines!(ax2, ω, real(ϕ_j_inj), label=L"J_{inj}")
# #lines!(ax2, ω, real(ϕ_j_brunel), label=L"J_{Brunel}")
# #lines!(ax2, ω, real(ϕ_j_kerr), label=L"J_{Kerr}")
# #scatter!(k, real(ϕ_j_kerr)[kth_harmonic_idx])
# axislegend(ax2, position=:rt)
# ylims!(ax2, [-π, π])
# ax2.yticks = -π:π/4:π
# ax2.ytickformat = ys -> ["$(y/pi)π" for y in ys]
# xlims!(ax2, [k-1, k+1])
#save("antreibendes_feld.pdf", fig)
fig

