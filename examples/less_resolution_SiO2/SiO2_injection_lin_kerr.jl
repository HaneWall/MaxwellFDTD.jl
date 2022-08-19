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
#amplitude_pump = intensity2amplitude(1.5e14)

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
a = 13.39

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

n_real = n_bound(γ_lorentz[1], ω_0[1], χ_1[1], ω .* ω_central)[1] 
t_theory = t_s.(1., n_real, 0.0, 0.0)
r_theory = r_s.(1., n_real, 0.0, 0.0)
r_nl_theory = r_s_nl.(1., n_real)
tr_a = t_s(1., 1.45, 0.0, 0.0)
j_eff = (χ_injection .* (tr_a .* E_arr).^13) 
P_1 = (χ_1 .* (tr_a .* E_arr) .+ χ_3 .* (tr_a .* E_arr).^3 )

E_measured_bulk = d2.Ez[30000:40000]
E_amplitude_spectrum = zeros(Float64, padded_L)
E_amplitude_spectrum[center_help_L_low:center_help_L_high] = window .* E_measured_bulk
FT_E_bulk = fftshift(fft(E_amplitude_spectrum))


signal_injection_theory = zeros(Float64, padded_L)
signal_injection_theory[center_help_L_low:center_help_L_high] = window .* P_1

signal_j_eff = zeros(Float64, padded_L)
signal_j_eff[center_help_L_low:center_help_L_high] = window .* j_eff

FT_P_1 = fftshift(fft(signal_injection_theory))
FT_J = fftshift(fft(signal_j_eff))

E_measured_reflection = d4.Ez[35000:45000]
E_amplitude_spectrum_refl = zeros(Float64, padded_L)
E_amplitude_spectrum_refl[center_help_L_low:center_help_L_high] = window .* E_measured_reflection
FT_E_bulk_refl = fftshift(fft(E_amplitude_spectrum_refl))

fig_relative_amplitudes = Figure(resolution = (800, 800), font="CMU Serif") 
ax1 = Axis(fig_relative_amplitudes[1, 1], xlabel=L"\omega / \omega_{pump}", title="log10 - Amplituden", yscale=log10) 
lines!(ω, abs.((FT_P_1)./(padded_L)) .+ 0.01, color=:blue, label=L"P(\omega) = FT[\chi^{(1)} t_{lin} E + \chi^{(3)}(t_{lin} E)^{3}](\omega)")
lines!(ω, abs.(r_theory .* (FT_P_1)./(padded_L)) .+ 0.01, color=:blue, label=L"r_{lin} P(\omega)", linestyle=:dash)
lines!(ω, abs.(r_nl_theory.* (FT_P_1)./(padded_L)) .+ 0.01, color=:seagreen, label=L"r_{nl} P(\omega)", linestyle=:dash)
lines!(ω, (29*10e-7 .* abs.((FT_J)./(padded_L)) .+ 0.01)./abs.(ω), color=:black, label=L"J_{eff} = C \cdot FT[\chi^{(Injection)} (t_{lin}E)^{13}](\omega) / \omega")
lines!(ω, r_nl_theory .* (29*10e-7 .* abs.((FT_J)./(padded_L)) .+ 0.01)./abs.(ω), color=:black, label=L"J_{eff}  = r_{nl} C \cdot FT[\chi^{(Injection)} (t_{lin}E)^{13}](\omega) / \omega", linestyle=:dash)

lines!(ω, (abs.(FT_E_bulk_refl)./padded_L .+ 0.01), color=:red, linestyle=:dash, label=L"E_{re}")
lines!(ω, (abs.(FT_E_bulk)./padded_L .+ 0.01), color=:red, label=L"E_{tr}")
xlims!(ax1, [0.0, 11.0])
ylims!(ax1, [10^2 , 10^(10)])
axislegend(ax1, position=:rt, nbanks=2)
fig_relative_amplitudes
save("kerr_linear_injection_reflection.pdf", fig_relative_amplitudes)