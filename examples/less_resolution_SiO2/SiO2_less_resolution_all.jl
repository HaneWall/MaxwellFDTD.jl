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
SizeX = 6000
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
χ_3 = [2.2e-22]

# drude parameters
γ_plasma = 1e15
#γ_plasma = 0.0

# tunnel parameters E_gap = 7.5

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
m1 = LorentzMedium1D(g, CartesianIndices((3000:5999,)), 1.0, γ_lorentz, ω_0, χ_1, χ_2, χ_3)
m2 = DrudeMedium1D(g, CartesianIndices((3000:5999,)), γ_plasma, ρ_mol_density)
m3 = TunnelMedium1D(g, CartesianIndices((3000:5999,)), E_gap, ρ_mol_density)

bound_media = [m1]
drude_media = [m2]
tunnel_media = [m3]

bound_media = [m1]
# 4. define grid coefficients that respect ϵ_inf from the media 
c_grid = GridCoefficients1D_w_CPML(g, bound_media, c_PML)
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
d2 = PointDetector(CartesianIndex((3000,)), 1, g.MaxTime)
d3 = PointDetector(CartesianIndex((5500,)), 1, g.MaxTime)
d4 = PointDetector(CartesianIndex((505,)), 1, g.MaxTime)
detectors = [d2, d3, d4]

# 7. place sources 
s0 = GaussianWavePointSource(g, CartesianIndex((508,)), true, false, false, amplitude_pump, ceil(500e-15 / g.Δt), t_fwhm, ppw)
s1 = GaussianWavePointSource(g, CartesianIndex((508,)), true, true, false, amplitude_probe, ceil(500e-15 / g.Δt), t_fwhm_probe, ppw_probe)
sources = [s0]

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3), " seconds")
println("Init Complete")
CPUtic()
start = time()




Ê = amplitude_pump
Γ̂ = Γ_ADK(Ê, E_gap * q_0)
a = 13.0

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

    for source in sources
        sourceH!(source, F, timestep)
    end

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

## same with only one cell --> nearly no cascade effects 
CPUtic()
start = time()

# 1. define grid
g_one_cell = Grid1D(SizeX, courant, Δx, MaxTime)

# 2. define fields that exist everywhere
F_one = Fields1D(g_one_cell)
MF_one = MaterialFields1D(g_one_cell)

# PML Fields and parameters
PML_Thickness = [500]
F_PML_one = CPML_Ψ_Fields_1D(g_one_cell, PML_Thickness)
c_PML_one = CPML_Parameters_1D(g_one_cell, PML_Thickness)

# init the media (superposition of different effects, that act inside the medium)
m1_ = LorentzMedium1D(g_one_cell, CartesianIndices((3000:3000,)), 1.0, γ_lorentz, ω_0, χ_1, χ_2, χ_3)
m2_ = DrudeMedium1D(g_one_cell, CartesianIndices((3000:3000,)), γ_plasma, ρ_mol_density)
m3_ = TunnelMedium1D(g_one_cell, CartesianIndices((3000:3000,)), E_gap, ρ_mol_density)

bound_media_ = [m1_]
drude_media_ = [m2_]
tunnel_media_ = [m3_]

# 4. define grid coefficients that respect ϵ_inf from the media 
c_grid_one = GridCoefficients1D_w_CPML(g_one_cell, bound_media_, c_PML_one)
f_grid_one = FieldIonizationCoefficients1D(g_one_cell)

# 5. define fields inside the media
LF1_ = LorentzFields1D(m1_)
LF_ = [LF1_]
DF1_ = DrudeFields1D(m2_)
DF_ = [DF1_]
TF1_ = TunnelFields1D(m3_)
TF_ = [TF1_]

# 6. place detectors 
#d1 = LineDetector(CartesianIndices((1:g.SizeX,)), 1, g.MaxTime)
d2_ = PointDetector(CartesianIndex((3000,)), 1, g_one_cell.MaxTime)
d3_ = PointDetector(CartesianIndex((5500,)), 1, g_one_cell.MaxTime)
d4_ = PointDetector(CartesianIndex((505,)), 1, g_one_cell.MaxTime)
detectors_ = [d2_, d3_, d4_]

# 7. place sources 
s0_ = GaussianWavePointSource(g_one_cell, CartesianIndex((508,)), true, false, false, amplitude_pump, ceil(500e-15 / g.Δt), t_fwhm, ppw)
s1_ = GaussianWavePointSource(g, CartesianIndex((508,)), true, true, false, amplitude_probe, ceil(500e-15 / g.Δt), t_fwhm_probe, ppw_probe)
sources_ = [s0_]

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3), " seconds")
println("Init Complete")
CPUtic()
start = time()


Ê = amplitude_pump
Γ̂ = Γ_ADK(Ê, E_gap * q_0)
a = 13.0

for timestep in ProgressBar(1:g_one_cell.MaxTime)

    # for (m_idx, m) in enumerate(tunnel_media_)
    #     updatePlasma!(MF_one, TF_[m_idx], f_grid_one, F_one, m)
    #     updateJtunnel!(MF_one, TF_[m_idx], m)
    # end

    for (m_idx, m) in enumerate(tunnel_media_)
        updatePlasmaTangent!(MF_one, TF_[m_idx], f_grid_one, F_one, m, a, Γ̂, Ê)
        updateJtunnel!(MF_one, TF_[m_idx], m)
    end

    for (m_idx, m) in enumerate(drude_media_)
        updateJfree!(MF_one, DF_[m_idx], F_one, m)
    end

    for (m_idx, m) in enumerate(bound_media_)
        updatePNl!(MF_one, LF_[m_idx], F_one, m)
        updateJbound!(MF_one, LF_[m_idx], m, g)
        updatePbound!(MF_one, LF_[m_idx], m, g)
    end

    updateJ!(MF_one)

    update_Ψ_H!(F_PML_one, F_one, g_one_cell, c_PML_one)

    updateH!(F_one, g_one_cell, c_grid_one)

    for source in sources_
        sourceH!(source, F_one, timestep)
    end

    apply_Ψ_H!(F_PML_one, F_one, g_one_cell, c_PML_one)

    update_Ψ_E!(F_PML_one, F_one, g_one_cell, c_PML_one)

    updateE!(F_one, MF_one, g_one_cell, c_grid_one)

    for source in sources_
        sourceE!(source, F_one, timestep)
    end

    apply_Ψ_E!(F_PML_one, F_one, g_one_cell, c_PML_one)

    for d in detectors_
        #safeΓ_ADK!(d, MF, timestep)
        safeJ_tunnel!(d, MF_one, timestep)
        safeJ_free!(d, MF_one, timestep)
        safeJ_bound!(d, MF_one, timestep)
        safeE!(d, F_one, timestep)
        #safeP!(d, MF, timestep)
        safeJ!(d, MF_one, timestep)
        #safePNl!(d, MF, timestep)
    end
end

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3), " seconds")
println("Computation Complete")

# ## 
# timeseries_plot(
#     d4.Ez,
#     "E_{z, refl}",
#     "detector_SF.pdf")

# ##
# timeseries_plot(
#     d2.Ez,
#     "E_{z, first cell medium}",
#     "detector_first_cell.pdf")

# ##
# timeseries_plot(
#     d3.Ez,
#     "E_{z, 7000 cells inside}",
#     "detector_mid.pdf")

##
# plot_log10_power_spectrum(
#     Array(30000*g.Δt:g.Δt:40000*g.Δt),
#     hcat(d2.J_Bound[30000:40000],
#         d2.J_Tunnel[30000:40000],
#         d2.J_Free[30000:40000],
#         d2.Jz[30000:40000]),
#     ω_central,
#     [0.0, 20.0],
#     [0.0, 25.0],
#     ["J_{Kerr}", "J_{Injection}", "J_{Brunel}", "J_{Overall}"],
#     true,
#     "ADK_injection_kerr_brunel_all.pdf")

# ##
# plot_log10_power_spectrum_current_and_E(
#     Array(30000*g.Δt:g.Δt:40000*g.Δt),
#     hcat(d2.Jz[30000:40000], d4.Ez[35000:45000], d2.J_Free[30000:40000], d2.J_Bound[30000:40000], d2.J_Tunnel[30000:40000]),
#     ω_central,
#     [0.0, 20.0],
#     [0.0, 1.0],
#     ["J_z", "E_{z, Refl}", "J_{Brunel}", "J_{Kerr}", "J_{Injection}"],
#     true,
#     "Reflection_a_13_less_res_ADK.pdf")

## Source for theoretical considerations 
function gaussian(amplitude::Float64, S_c::Float64, Δt::Float64, t::Array{Float64,1}, t_step_peak::Float64, t_fwhm::Float64, ppw::Float64)
    src = zeros(Float64, size(t)[1])
    for idx in 1:size(t)[1]
        src[idx] = amplitude * exp(-2 * log(2) * (((t[idx] - t_step_peak)) * Δt / t_fwhm)^2) * sin(2 * π / ppw * (S_c * t[idx]))
    end
    return src
end

function gaussian_envelope(amplitude::Float64, Δt::Float64, t::Array{Float64,1}, t_step_peak::Float64, t_fwhm::Float64)
    src = zeros(Float64, size(t)[1])
    for idx in 1:size(t)[1]
        src[idx] = amplitude * exp(-2 * log(2) * (((t[idx] - t_step_peak)) * Δt / t_fwhm)^2)
    end
    return src
end

S_C = courant
t_arr = Array(25000.0:1.0:35000.0)
t_peak = 500e-15 / g.Δt
t_real = t_arr .* g.Δt
E_arr = gaussian(amplitude_pump, S_C, g.Δt, t_arr, t_peak, t_fwhm, ppw)
E_arr_envelope = gaussian_envelope(amplitude_pump, g.Δt, t_arr, t_peak, t_fwhm)

## Theoretical Reflection coefficients - effective nl
eff_nl = effective_nonlinearity_m(amplitude_pump, E_gap * q_0)
χ_brunel = χ_brunel_stat(Γ̂, Ê, ρ_mol_density, eff_nl)
χ_injection = χ_injection_stat(Γ̂, Ê, E_gap * q_0, ρ_mol_density, eff_nl)

## measured E_reflection 
E_measured = d4.Ez[35000:45000]
E_measured_one_cell = d4_.Ez[35000:45000]
L = length(t_real)
padded_L = nextpow(2, L)
center_help_L_low = floor(Int64, (padded_L - L) / 2)
center_help_L_high = floor(Int64, L + center_help_L_low - 1)
δt = abs(t[2] - t[1])
δf = 1 / δt
ω = 2 * π * fftshift(fftfreq(padded_L, δf)) ./ ω_central
reflection_amplitude_spectrum = zeros(Float64, padded_L) 
reflection_amplitude_spectrum_one_cell = zeros(Float64, padded_L)
window = blackman(L)
reflection_amplitude_spectrum[center_help_L_low:center_help_L_high] = window .* E_measured
reflection_amplitude_spectrum_one_cell[center_help_L_low:center_help_L_high] = window .* E_measured_one_cell
FT_reflection = fftshift(fft(reflection_amplitude_spectrum))
FT_reflection_one_cell = fftshift(fft(reflection_amplitude_spectrum_one_cell))

j_measured = d2.Jz[30000:40000]
j_amplitude_spectrum = zeros(Float64, padded_L)
j_amplitude_spectrum[center_help_L_low:center_help_L_high] = window .* j_measured
FT_j = fftshift(fft(j_amplitude_spectrum))

j_measured_one_cell = d2_.Jz[30000:40000]
j_amplitude_spectrum_one_cell = zeros(Float64, padded_L)
j_amplitude_spectrum_one_cell[center_help_L_low:center_help_L_high] = window .* j_measured_one_cell
FT_j_one_cell = fftshift(fft(j_amplitude_spectrum_one_cell))


fig_refl = Figure(resolution = (800, 800), font="CMU Serif")
ax1 = Axis(fig_refl[1, 1], xlabel=L"\omega / \omega_{pump}", title="Reilection amplitudes", subtitle="Tangent ADK, a= $eff_nl")
lines!(ω, log10.((abs.(ω .*  FT_reflection)) ./ maximum(abs.( ω .* FT_reflection))), color=:blue, label=L"Bulk, first Cell, \log_{10}\mid \omega E_z \mid")
lines!(ω, log10.((abs.(ω .*  FT_reflection_one_cell)) ./ maximum(abs.(ω .* FT_reflection_one_cell))), color=:blue, linestyle=:dash, label=L"One Cell, \log_{10}\mid \omega E_z \mid")

lines!(ω, log10.((abs.(FT_j)) ./ maximum(abs.(FT_j))), color=:red, label=L"Bulk, first Cell, \log_{10} \mid J_z \mid")
lines!(ω, log10.((abs.(FT_j_one_cell)) ./ maximum(abs.(FT_j_one_cell))), color=:red, linestyle=:dash, label=L"One Cell, \log_{10}\mid J_z \mid")

ylims!(ax1, [-7.0, 0.0])
xlims!(ax1, [0.0, 15.0])
axislegend(ax1, position=:rt)
fig_refl
save("reflection_amplitude_Tangent_ADK.pdf", fig_refl)

# ## harmonics
# h_brunel = Array(1:2:13)
# h_injection = Array(1:2:11)
# degen_injection = binomial_degen.(eff_nl - 1, h_injection)
# degen_brunel = binomial_degen.(eff_nl + 1, h_brunel)
# normed_degen_injection = degen_injection ./ maximum(degen_injection)
# normed_degen_brunel = degen_brunel ./ degen_brunel[1] # first harmonic hits different

# ## amplitude_spectrum Brunel and Injection 
# L = length(t_real)
# padded_L = nextpow(2, L)
# center_help_L_low = floor(Int64, (padded_L - L) / 2)
# center_help_L_high = floor(Int64, L + center_help_L_low - 1)
# δt = abs(t[2] - t[1])
# δf = 1 / δt
# ω = 2 * π * fftshift(fftfreq(padded_L, δf)) ./ ω_central
# third_harmonic_index = argmin(abs.(ω .- 1))
# window = blackman(L)
# signal_injection = zeros(Float64, padded_L)
# signal_injection_one_cell = zeros(Float64, padded_L)
# signal_brunel = zeros(Float64, padded_L)
# signal_brunel_one_cell = zeros(Float64, padded_L)
# signal_injection[center_help_L_low:center_help_L_high] = window .* d2.J_Tunnel[30000:40000]
# signal_injection_one_cell[center_help_L_low:center_help_L_high] = window .* d2_.J_Tunnel[30000:40000]
# signal_brunel[center_help_L_low:center_help_L_high] = window .* d2.J_Free[30000:40000]
# signal_brunel_one_cell[center_help_L_low:center_help_L_high] = window .* d2_.J_Free[30000:40000]
# FT_injection = fftshift(fft(signal_injection))
# FT_injection_one_cell = fftshift(fft(signal_injection_one_cell))
# FT_brunel = fftshift(fft(signal_brunel))
# FT_brunel_one_cell = fftshift(fft(signal_brunel_one_cell))


# ## plot degens:
# fig_harm = Figure(resolution=(800, 800), font="CMU Serif")
# ax1 = Axis(fig_harm[1, 1], title="Normed Injection amplitude and Degens", subtitle="ADK, a = $eff_nl")
# scatter!(h_injection, (normed_degen_injection))
# lines!(ω, (abs.(FT_injection) ./ maximum(abs.(FT_injection))), color=:blue, label="Bulk, first Cell")
# lines!(ω, (abs.(FT_injection_one_cell) ./ maximum(abs.(FT_injection_one_cell))), color=:blue, linestyle=:dash, label="One Cell")
# ylims!(ax1, [0.0, 1.0])
# xlims!(ax1, [0.0, 15.0])
# ax2 = Axis(fig_harm[2, 1], title="log10 normed Injection amplitude and Degens", subtitle="ADK", xlabel=L"ω / ω_{pump}")
# lines!(ω, log10.(abs.(FT_injection) ./ maximum(abs.(FT_injection))), color=:red, label="Bulk, first Cell")
# lines!(ω, log10.(abs.(FT_injection_one_cell) ./ maximum(abs.(FT_injection_one_cell))), color=:red, linestyle=:dash, label="One Cell")
# scatter!(h_injection, log10.(normed_degen_injection), color=:red)
# ylims!(ax2, [-6.0, 0.0])
# xlims!(ax2, [0.0, 15.0])
# axislegend(ax2, position=:lb)
# ax3 = Axis(fig_harm[1,2], title="Normed Brunel amplitude and Degens", subtitle="ADK, a = $eff_nl")
# scatter!(h_brunel, (normed_degen_brunel))
# lines!(ω, (abs.(FT_brunel) ./ abs.(FT_brunel[third_harmonic_index])), color=:blue, label="Bulk, first Cell")
# lines!(ω, (abs.(FT_brunel_one_cell) ./ abs.(FT_brunel_one_cell[third_harmonic_index])), color=:blue, linestyle=:dash, label="One Cell")
# ylims!(ax3, [0.0, 1.0])
# xlims!(ax3, [0.0, 15.0])
# ax4 = Axis(fig_harm[2, 2], title="log10 normed Brunel amplitude and Degens", subtitle="ADK", xlabel=L"ω / ω_{pump}")
# lines!(ω, log10.(abs.(FT_brunel) ./ abs.(FT_brunel[third_harmonic_index])), color=:red, label="Bulk, first Cell")
# lines!(ω, log10.(abs.(FT_brunel_one_cell) ./ abs.(FT_brunel_one_cell[third_harmonic_index])), color=:red, linestyle=:dash, label="One Cell")
# scatter!(h_brunel, log10.(normed_degen_brunel), color=:red)
# ylims!(ax4, [-6.0, 0.0])
# xlims!(ax4, [0.0, 15.0])
# axislegend(ax4, position=:lb)
# fig_harm
# save("Degens_all_ADK.pdf", fig_harm)

# ## plot error degens:
# h_indices_injection = [argmin(abs.(ω .- i)) for i in 1:2:11]
# h_indices_brunel = [argmin(abs.(ω .- i)) for i in 1:2:13]
# fig_err = Figure(resolution=(800, 800), font="CMU Serif")
# ax1 = Axis(fig_err[1, 1])
# normed_peaks_injection = abs.(FT_injection[h_indices])./maximum(abs.(FT_injection))
# normed_peaks_injection_one_cell = abs.(FT_injection_one_cell[h_indices])./maximum(abs.(FT_injection_one_cell))
# normed_peaks_brunel = abs.(FT_injection[h_indices_brunel])./maximum(abs.(FT_brunel))
# normed_peaks_brunel_one_cell = abs.(FT_brunel_one_cell[h_indices_brunel])./maximum(abs.(FT_brunel_one_cell))
# scatter!(h_injection, log10.(abs.(normed_peaks_injection - normed_degen_injection)))
# scatter!(h_injection, log10.(abs.(normed_peaks_injection_one_cell - normed_degen_injection)))
# scatter!(h_brunel, log10.(abs.(normed_peaks_brunel - normed_degen_brunel)))
# scatter!(h_brunel, log10.(abs.(normed_peaks_brunel_one_cell - normed_degen_brunel)))
# fig_err

function n_bound(γ, ω_0, χ_1, ω)
    ϵ_complex = epsilon_complex(γ, ω_0, χ_1, ω)
    ϵ_real = real.(ϵ_complex)
    ϵ_imag = imag.(ϵ_complex)
    n_r =  zeros(Float64, size(ω)[1])
    n_i =  zeros(Float64, size(ω)[1])
    
    @. n_r = sqrt(1/2 *(sqrt(ϵ_real^2 + ϵ_imag^2) + ϵ_real))
    @. n_i = sqrt(1/2 *(sqrt(ϵ_real^2 + ϵ_imag^2) - ϵ_real))
    return  [n_r, n_i]
end

function epsilon_complex(γ, ω_0, χ_1, ω)
    eps_complex = zeros(ComplexF64, length(ω)) .+ 1.
    for idx = 1:length(γ)
        @. eps_complex += χ_1[idx] * (ω_0[idx]^2) / (ω_0[idx]^2 - ω^2 + 1im * γ[idx] * ω)
    end
    return eps_complex
end

# ## plot dispersion form bound electrons
# ϵ_complex = epsilon_complex(γ_lorentz, ω_0[1], χ_1[1], ω.*ω_central)
# ϵ_real = real.(ϵ_complex)
# ϵ_imag = imag.(ϵ_complex) 

# n_r =  zeros(Float64, size(ω)[1])
# n_i =  zeros(Float64, size(ω)[1])

# @. n_r = sqrt(1/2 *(sqrt(ϵ_real^2 + ϵ_imag^2) + ϵ_real))
# @. n_i = sqrt(1/2 *(sqrt(ϵ_real^2 + ϵ_imag^2) - ϵ_real))
# fig_disp = Figure(resolution = (800, 400), font = "CMU Serif")
# ax1 = Axis(fig_disp[1, 1], xlabel=L"\omega / \omega_{pump}", ylabel=L"n(\omega)")
# lines!(ω, n_r)
# xlims!(ax1, [0.0, 15.])
# ylims!(ax1, [1.4, 1.6])
# fig_disp
# save("bound_dispersion.pdf", fig_disp)

n_ior = 1.47




# E_brunel_theory = ϵ_0 .* abs.(E_arr .+ 0im).^(eff_nl) .* E_arr .* χ_brunel
# E_injection_theory = ϵ_0 .* abs.(E_arr .+ 0im).^(eff_nl) ./ E_arr .* χ_injection

#n_up = Array(7:1:13)
#n_down = Array(6:-1:0)
#degen_injection = binomial_degen.(eff_nl, n_up, n_down)

# function get_effective_reflection(χ_1::Float64, χ_a::Float64, E::Array{Float64,1})
#     E_L = E_reflection_for_χ(E, χ_1, 1., 1.0)
#     E_K = E_reflection_for_χ(E, χ_a, 13., 1.0)
#     E_ovr = E_L + E_K
#     E = hcat(E_L, E_K, E_ovr)
#     return E
# end

# E_refl_a = get_effective_reflection(χ_1[1], χ_13, E_arr)

# plot_log10_power_spectrum_current_and_E(
#     Array(80000*g.Δt:g.Δt:100000*g.Δt),
#     hcat(d2.Jz[80000:100000], E_refl_a[:, 2], d2.J_Free[80000:100000], d2.J_Bound[80000:100000], d2.J_Tunnel[80000:100000]),
#     ω_central,
#     [0.0, 20.0],
#     [0.0, 1.0],
#     ["J_z", "E_{z, Refl, a}", "J_{Brunel}", "J_{Kerr}", "J_{Injection}"],
#     true,
#     "Reflection_tangentadk_a13_theory.pdf")



