

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
MaxTime = 51000 

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

    # for (m_idx, m) in enumerate(drude_media)
    #     updateJfree!(MF, DF[m_idx], F, m)
    # end

    # for (m_idx, m) in enumerate(bound_media)
    #     updatePNl!(MF, LF[m_idx], F, m)
    #     updateJbound!(MF, LF[m_idx], m, g)
    #     updatePbound!(MF, LF[m_idx], m, g)
    # end

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
        safeP!(d, MF, timestep)
        safeJ!(d, MF, timestep)
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

# ##
#plot_log10_power_spectrum(
    #Array(30000*g.Δt:g.Δt:40000*g.Δt),
    #hcat(d2.J_Tunnel[30000:40000]),
    #ω_central,
    #[0.0, 20.0],
    #[0.0, 25.0],
    #["J_{Injection}"],
    #true,
    #"Tangent_only_injection.pdf")

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

    # for (m_idx, m) in enumerate(drude_media)
    #     updateJfree!(MF, DF[m_idx], F, m)
    # end

    # for (m_idx, m) in enumerate(bound_media)
    #     updatePNl!(MF, LF[m_idx], F, m)
    #     updateJbound!(MF, LF[m_idx], m, g)
    #     updatePbound!(MF, LF[m_idx], m, g)
    # end

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

# plot_log10_power_spectrum(
#     Array(30000*g.Δt:g.Δt:40000*g.Δt),
#     hcat(d2_.J_Tunnel[30000:40000]),
#     ω_central,
#     [0.0, 20.0],
#     [0.0, 25.0],
#     ["J_{Injection}"],
#     true,
#     "ADK_only_injection_one_cell.pdf")



## Source for theoretical considerations 
function gaussian(amplitude::Float64, S_c::Float64, Δt::Float64, t::Array{Float64,1}, t_step_peak::Float64, t_fwhm::Float64, ppw::Float64)
    src = zeros(Float64, size(t)[1])
    for idx in 1:size(t)[1]
        src[idx] = amplitude * exp(-2 * log(2) * (((t[idx] - t_step_peak)) * Δt / t_fwhm)^2) * sin(2 * π / ppw * (S_c * t[idx]))
    end
    return src
end

S_C = courant
t_arr = Array(25000.0:1.0:35000.0)
t_peak = 500e-15 / g.Δt
t_real = t_arr .* g.Δt
E_arr = gaussian(amplitude_pump, S_C, g.Δt, t_arr, t_peak, t_fwhm, ppw)

# ## Theoretical Reflection coefficients - effective nl
eff_nl = effective_nonlinearity_m(amplitude_pump, E_gap * q_0)
χ_brunel = χ_brunel_stat(Γ̂, Ê, ρ_mol_density, eff_nl)
χ_injection = χ_injection_stat(Γ̂, Ê, E_gap*q_0, ρ_mol_density, eff_nl)

## harmonics
h_brunel = Array(1:2:13)
h_injection = Array(1:2:11)
degen_injection = binomial_degen.(eff_nl - 1, h_injection)
degen_brunel = binomial_degen.(eff_nl + 1, h_brunel)
normed_degen_injection = degen_injection./maximum(degen_injection)
normed_degen_brunel = degen_brunel./maximum(degen_brunel)

## amplitude_spectrum Injection 
L = length(t_real)
padded_L = nextpow(2, L)
center_help_L_low = floor(Int64, (padded_L - L)/2)
center_help_L_high = floor(Int64, L + center_help_L_low - 1) 
δt = abs(t[2]-t[1])
δf = 1/δt
ω = 2*π*fftshift(fftfreq(padded_L, δf))./ ω_central
window = blackman(L)
signal_injection_theory = zeros(Float64, padded_L)
signal_injection = zeros(Float64, padded_L)
signal_injection_one = zeros(Float64, padded_L)
signal_injection_theory[center_help_L_low:center_help_L_high] = window .* χ_injection .* (E_arr.+0im).^13
signal_injection[center_help_L_low:center_help_L_high] = window.*d2.J_Tunnel[30000:40000]
signal_injection_one[center_help_L_low:center_help_L_high] = window.*d2_.J_Tunnel[30000:40000]

FT_injection_theory = fftshift(fft(signal_injection_theory))
FT_injection = fftshift(fft(signal_injection))
FT_injection_one = fftshift(fft(signal_injection_one)) 


## plot degens:
fig_harm = Figure(resolution = (800, 800), font = "CMU Serif")
ax = Axis(fig_harm[1, 1], title="Normed Injection amplitude and degeneracies", subtitle="ADK")
scatter!(h_injection, (normed_degen_injection))
lines!(ω, (abs.(FT_injection)./maximum(abs.(FT_injection))), label="Bulk, first Cell")
lines!(ω, (abs.(FT_injection_one)./maximum(abs.(FT_injection_one))), linestyle=:dash, color=:blue, label = "One Cell")
ylims!(ax, [0.0, 1.0])
xlims!(ax, [0.0, 15.0])
axislegend(ax, position=:rt)
ax2 = Axis(fig_harm[2, 1], title="log10 normed Injection amplitude and Degens, Tangent ADK", xlabel = L"ω / ω_{pump}")
lines!(ω, log10.(abs.(FT_injection)./maximum(abs.(FT_injection))), color=:red, label = "Bulk, first Cell")
lines!(ω, log10.((abs.(FT_injection_one)./maximum(abs.(FT_injection_one)))), color=:red, linestyle=:dash, label="One Cell")
#lines!(ω, log10.((abs.(FT_injection_theory))./maximum(abs.(FT_injection_theory))), color=:green, label="Theory, Reflection")
scatter!(h_injection, log10.(normed_degen_injection), color=:red)
ylims!(ax2, [-5.0, 0.0])
xlims!(ax2, [0.0, 15.0])
axislegend(ax2, position=:rt)
fig_harm
#save("Degens_injection_ADK.pdf", fig_harm)

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
ax1 = Axis(fig_refl[1, 1], xlabel=L"\omega / \omega_{pump}", title="Reflection amplitudes", subtitle="Tangent ADK, a= $eff_nl")
lines!(ω, log10.((abs.(ω .*  FT_reflection)) ./ maximum(abs.( ω .* FT_reflection))), color=:blue, label=L"Bulk, \log_{10}\mid \omega E_z \mid")
lines!(ω, log10.((abs.(ω .*  FT_reflection_one_cell)) ./ maximum(abs.(ω .* FT_reflection_one_cell))), color=:blue, linestyle=:dash, label=L"One Cell, \log_{10}\mid \omega E_z \mid")

lines!(ω, log10.((abs.(FT_j)) ./ maximum(abs.(FT_j))), color=:red, label=L"Bulk, first Cell, \log_{10} \mid J_z \mid")
lines!(ω, log10.((abs.(FT_j_one_cell)) ./ maximum(abs.(FT_j_one_cell))), color=:red, linestyle=:dash, label=L"One Cell, \log_{10}\mid J_z \mid")

ylims!(ax1, [-5.0, 0.0])
xlims!(ax1, [0.0, 15.0])
axislegend(ax1, position=:rt)
fig_refl
#save("reflection_amplitude_only_injection_Tangent_ADK.pdf", fig_refl)


## compare theory with simulation 
j_in =  χ_injection .* E_arr.^13 
j_real = χ_injection .* d2.Ez[30000:40000].^13
# P_diff = ϵ_0 .* (P_in[2:end] - P_in[1:end-1]) / g.Δt
fig_comp = Figure(resolution=(800, 800))
ax1 = Axis(fig_comp[1, 1], title = "Comparison")
lines!(log10.(abs.(j_in)), color=:red)
lines!(log10.(abs.(j_real)), color=:black, linestyle=:dash)
lines!(log10.(abs.(j_measured)), color=:green)
# lines!(log10.(abs.(P_diff)), color=:black)
ylims!(ax1, [0.0, 14.0])
xlims!(ax1, [3000, 7000])
ax2 = Axis(fig_comp[2, 1])
lines!(ω, log10.(abs.(FT_injection_theory)), color=:red)
lines!(ω, log10.(abs.(FT_injection)), color=:green)
ylims!(ax2, [10.0, 17.0])
xlims!(ax2, [0.0, 15.0])
fig_comp
# ##
# plot_log10_power_spectrum_current_and_E(
#     Array(80000*g.Δt:g.Δt:100000*g.Δt),
#     hcat(d2.Jz[80000:100000], d4.Ez[94500:114500], d2.J_Free[80000:100000], d2.J_Bound[80000:100000], d2.J_Tunnel[80000:100000]),
#     ω_central,
#     [0.0, 20.0],
#     [0.0, 1.0],
#     ["J_z", "E_{z, Refl}", "J_{Brunel}", "J_{Kerr}", "J_{Injection}"],
#     true,
#     "Reflection_tangentadk_a_13.pdf")

# ## Source for theoretical considerations 
# function gaussian(amplitude::Float64, S_c::Float64, Δt::Float64, t::Array{Float64,1}, t_step_peak::Float64, t_fwhm::Float64, ppw::Float64)
#     src = zeros(Float64, size(t)[1])
#     for idx in 1:size(t)[1]
#         src[idx] = amplitude * exp(-2 * log(2) * (((t[idx] - t_step_peak)) * Δt / t_fwhm)^2) * sin(2 * π / ppw * (S_c * t[idx]))
#     end
#     return src
# end

# S_C = courant
# t_arr = Array(65000.0:1.0:85000.0)
# t_peak = 500e-15 / g.Δt
# t_real = t_arr .* g.Δt
# E_arr = gaussian(amplitude_pump, S_C, g.Δt, t_arr, t_peak, t_fwhm, ppw)

# ## Theoretical Reflection coefficients - Only Bound Electrons - linear and kerr
# function get_kerr_reflection(χ_1::Float64, χ_3::Float64, E::Array{Float64,1})
#     E_L = E_reflection_for_χ(E, χ_1, 1.0, 1.0)
#     E_K = E_reflection_for_χ(E, χ_3, 3.0, 1.0)
#     E_ovr = E_L + E_K
#     E = hcat(E_L, E_K, E_ovr)
#     return E
# end

# E_refl_kerr = get_kerr_reflection(χ_1[1], χ_3[1], E_arr)

# plot_log10_power_spectrum_current_and_E(
#     Array(80000*g.Δt:g.Δt:100000*g.Δt),
#     hcat(d2.Jz[80000:100000], E_refl_kerr[:, end], d2.J_Free[80000:100000], d2.J_Bound[80000:100000], d2.J_Tunnel[80000:100000]),
#     ω_central,
#     [0.0, 20.0],
#     [0.0, 1.0],
#     ["J_z", "E_{z, Refl kerr}", "J_{Brunel}", "J_{Kerr}", "J_{Injection}"],
#     true,
#     "Reflection_tangentadk_a13_kerr.pdf")

# ## Theoretical Reflection coefficients - effective nl
# E_atom = 5.14e11
# χ_13 = χ_1[1] / E_atom^(a - 1)

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