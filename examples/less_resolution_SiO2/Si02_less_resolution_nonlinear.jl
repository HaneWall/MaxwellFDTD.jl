

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

# tunnel parameters
E_gap = 7.5

# source parameters
λ = 2100e-9
ω_central = 2 * π * c_0 / λ
ppw = λ / Δx
t_fwhm = 20e-15 # intensity FWHM
#amplitude_pump = intensity2amplitude(12e16) # 12TWcm^-2
amplitude_pump = intensity2amplitude(1.5e14 ) # 12TWcm^-2

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

    # for (m_idx, m) in enumerate(tunnel_media)
    #     updatePlasmaTangent!(MF, TF[m_idx], f_grid, F, m, a, Γ̂, Ê)
    #     updateJtunnel!(MF, TF[m_idx], m)
    # end

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
m1_ = LorentzMedium1D(g_one_cell, CartesianIndices((3000:5999,)), 1.0, γ_lorentz, ω_0, χ_1, χ_2, χ_3)
m2_ = DrudeMedium1D(g_one_cell, CartesianIndices((3000:3000,)), γ_plasma, ρ_mol_density)
m3_ = TunnelMedium1D(g_one_cell, CartesianIndices((3000:5999,)), E_gap, ρ_mol_density)

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

    # for (m_idx, m) in enumerate(tunnel_media_)
    #     updatePlasmaTangent!(MF_one, TF_[m_idx], f_grid_one, F_one, m, a, Γ̂, Ê)
    #     updateJtunnel!(MF_one, TF_[m_idx], m)
    # end

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

## FFT Helper 
L = length(t_real)
padded_L = nextpow(2, 4*L)
center_help_L_low = floor(Int64, (padded_L - L)/2)
center_help_L_high = floor(Int64, L + center_help_L_low - 1) 
δt = abs(t[2]-t[1])
δf = 1/δt
ω = 2*π*fftshift(fftfreq(padded_L, δf))./ ω_central
window = blackman(L)
#window = ones(L)

## E-Amplitudenverhältnisse zwischen Reflektion und erster Zelle
#transmission
E_measured_bulk = d2.Ez[30000:40000]
E_amplitude_spectrum = zeros(Float64, padded_L)
E_amplitude_spectrum[center_help_L_low:center_help_L_high] = window .* E_measured_bulk
FT_E_bulk = fftshift(fft(E_amplitude_spectrum))

E_measured_vacuum = d2_.Ez[30000:40000]
E_amplitude_spectrum_vacuum = zeros(Float64, padded_L)
E_amplitude_spectrum_vacuum[center_help_L_low:center_help_L_high] = window .* E_measured_vacuum
FT_E_vacuum = fftshift(fft(E_amplitude_spectrum_vacuum))

#reflection 
E_measured_reflection = d4.Ez[35000:45000]
E_amplitude_spectrum_refl = zeros(Float64, padded_L)
E_amplitude_spectrum_refl[center_help_L_low:center_help_L_high] = window .* E_measured_reflection
FT_E_bulk_refl = fftshift(fft(E_amplitude_spectrum_refl))

E_measured_reflection_vacuum = d4_.Ez[35000:45000]
E_amplitude_spectrum_refl_vac = zeros(Float64, padded_L)
E_amplitude_spectrum_refl_vac[center_help_L_low:center_help_L_high] = window .* E_measured_reflection_vacuum
FT_E_vac_refl = fftshift(fft(E_amplitude_spectrum_refl_vac))

n_real = n_bound(γ_lorentz[1], ω_0[1], χ_1[1], ω .* ω_central)[1] 
t_theory = t_s.(1., n_real, 0.0, 0.0)
r_theory = r_s.(1., n_real, 0.0, 0.0)
r_nl_theory = r_s_nl.(1., n_real)

tr_a = t_s(1., 1.45, 0.0, 0.0)
P_theory = (χ_1 .* tr_a .* E_measured_vacuum .* (tr_a .* E_measured_vacuum).^3)
#P_theory = (χ_1 .* E_measured_vacuum .+ χ_3 .* (1. .* E_measured_vacuum).^3)

P_amplitude_spectrum = zeros(Float64, padded_L)
P_amplitude_spectrum[center_help_L_low:center_help_L_high] = window .* P_theory
FT_P = fftshift(fft(P_amplitude_spectrum))

P_theory_chi_3 = χ_3 .* tr_a .* E_measured_vacuum.^3

P_3_amplitude_spectrum = zeros(Float64, padded_L)
P_3_amplitude_spectrum[center_help_L_low:center_help_L_high] = window .* P_theory_chi_3
FT_P_3 = fftshift(fft(P_3_amplitude_spectrum))

fig_relative_amplitudes = Figure(resolution = (800, 800), font="CMU Serif") 
ax1 = Axis(fig_relative_amplitudes[1, 1], xlabel=L"\omega / \omega_{pump}", title="log10 - Amplituden", yscale=log10) 
lines!(ω, (abs.(FT_E_vac_refl)./padded_L) .+ 0.01, color=:black, linestyle=:dash, label=L"E_{re, vac}")
lines!(ω, abs.((FT_P)./(padded_L)) .+ 0.01, color=:blue, label=L"P(\omega) = FT[\chi^{(1)}t_{lin}E + \chi^{(3)}(t_{lin} E)^3](\omega)")
lines!(ω, (abs.(FT_E_bulk)./padded_L) .+0.01, color=:red, label=L"E_{tr}")
lines!(ω, (abs.(FT_E_vacuum)./padded_L) .+0.01, color=:black, label=L"E(\omega)")
lines!(ω, (abs.(FT_E_bulk_refl)./padded_L .+ 0.01), color=:red, linestyle=:dash, label=L"E_{re}")
lines!(ω, (abs.(FT_P_3)./(padded_L) .+ 0.01), color=:blue, label=L"P_{3}(\omega) = FT[\chi^{(3)}(t_{lin}E)^3](\omega)", linestyle=:dash)
lines!(ω, (abs.(FT_P .* r_nl_theory)./(padded_L) .+ 0.01), color=:green, label=L"r_{nl}(\omega) \cdot P(\omega)", linestyle=:dash)
lines!(ω, (abs.(FT_P_3 .* r_nl_theory)./(padded_L) .+ 0.01), color=:green, label=L"r_{nl}(\omega) \cdot P_{3}(\omega)", linestyle=:dashdot)
xlims!(ax1, [0.0, 5.0])
axislegend(ax1, position=:rt, nbanks=2)

ax2 =Axis(fig_relative_amplitudes[2, 1], xlabel=L"\omega / \omega_{pump}", title="Fresnelkoeffizienten - Amplitudenrelationen")
#lines!(ω, (abs.(FT_E_bulk))./(abs.(FT_P)./(ω.^2)), color=:red)
lines!(ω, (abs.(FT_E_bulk))./abs.(FT_E_vacuum), color=:red, label=L"E_{tr}")
lines!(ω, abs.(t_theory), color=:gray, linestyle=:dash, label=L"t_{lin}")
#lines!(ω, abs.(FT_E_bulk_refl)./(abs.(FT_P)./(ω.^2)), color=:blue)
lines!(ω, abs.(FT_E_bulk_refl)./(abs.(FT_E_vacuum)), color=:blue, label=L"E_{re}")
lines!(ω, abs.(r_theory), color=:black, linestyle=:dash, label=L"r_{lin}")
lines!(ω, abs.(r_nl_theory), color=:green, linestyle=:dash, label=L"r_{nl}")
xlims!(ax2, [0.5, 2.5])
ylims!(ax2, [0., 1.0])
axislegend(ax2, position=:lc)
fig_relative_amplitudes
save("kerr_reflection_fresnel_small_w.pdf", fig_relative_amplitudes)


fig_relative_amplitudes_3 = Figure(resolution = (800, 800), font="CMU Serif") 
ax1 =Axis(fig_relative_amplitudes_3[1, 1:2], xlabel=L"\omega / \omega_{pump}", title="log10 - Amplituden", yscale=log10) 
lines!(ω, (abs.(FT_E_vac_refl)./padded_L) .+ 0.01, color=:black, linestyle=:dash, label=L"E_{re, vac}")
lines!(ω, abs.((FT_P)./(padded_L)) .+ 0.01, color=:blue, label=L"P(\omega) = FT[\chi^{(1)}t_{lin}E + \chi^{(3)}(t_{lin}E)^3](\omega)")
#lines!(ω, abs.(amplitude_P_bulk) .+ 0.01, color=:yellow)
lines!(ω, (abs.(FT_E_bulk)./padded_L) .+0.01, color=:red, label=L"E_{tr}")
lines!(ω, (abs.(FT_E_vacuum)./padded_L) .+0.01, color=:black, label=L"E(\omega)")
lines!(ω, (abs.(FT_E_bulk_refl)./padded_L .+ 0.01), color=:red, linestyle=:dash, label=L"E_{re}")
lines!(ω, (abs.(FT_P_3)./(padded_L) .+ 0.01), color=:blue, label=L"P_3(\omega) = FT[\chi^{(3)}(t_{lin}E)^3](\omega)", linestyle=:dash)
lines!(ω, (abs.(FT_P .* r_nl_theory)./(padded_L )) .+ 0.01, color=:green, label=L"r_{nl}(\omega) \cdot P(\omega)", linestyle=:dash)
lines!(ω, (abs.(FT_P_3 .* r_nl_theory)./(padded_L) .+ 0.01), color=:green, label=L"r_{nl}(\omega) \cdot P_{3}(\omega)", linestyle=:dashdot)
xlims!(ax1, [0.0, 5.0])
axislegend(ax1, position=:rt, nbanks=2)

ax2 = Axis(fig_relative_amplitudes_3[2, 1], xlabel=L"\omega / \omega_{pump}", title="Fresnelähnlich")
lines!(ω, abs.(FT_E_bulk_refl) ./ (abs.(FT_P)), color=:blue, label=L"(E_{re} )/ (P(\omega))")
lines!(ω, abs.(r_theory), color=:black, linestyle=:dash, label=L"r_{lin}")
lines!(ω, abs.(r_nl_theory), color=:green, linestyle=:dash, label=L"r_{nl}")
#lines!(ω, abs.(FT_E_bulk) ./ (abs.(FT_E_bulk_refl)), color=:red, label=L"E_{tr}")
xlims!(ax2, [0., 5.])
ylims!(ax2, [0., 0.4])
axislegend(ax2, position=:ct)

ax3 = Axis(fig_relative_amplitudes_3[2, 2],  xlabel=L"\omega / \omega_{pump}", title="Amplitudenrelationen")
lines!(ω, abs.(FT_E_bulk_refl) ./ (abs.(r_nl_theory .* FT_P_3)), color=:red, label=L"(E_{re} )/ (r_{nl} \cdot P_3(\omega))")
lines!(ω, abs.(FT_E_bulk_refl) ./ (abs.(r_nl_theory .* FT_P)), color=:blue, label=L"(E_{re} )/ (r_{nl} \cdot P(\omega))")
axislegend(ax3, position=:ct)
xlims!(ax3, [1.5, 4.5])
ylims!(ax3, [0.5, 1.5])
fig_relative_amplitudes_3
save("kerr_reflection_fresnel_3_w.pdf", fig_relative_amplitudes_3)

## compare time amplitude
fig_time = Figure(resolution = (800, 800), font="CMU Serif")
ax1 = Axis(fig_time[1, 1], xlabel="Timestep", ylabel="Amplitude")
E_R_theory = r_s_nl(1., 1.45) .* (χ_1 .* E_arr .+ χ_3 .* (1. .* E_arr).^3)
#lines!(E_arr, label="source")
lines!(E_measured_vacuum, label="vacuum")
lines!(E_measured_bulk, label="bulk")
lines!(P_theory, label="Theory Bulk")
lines!(E_R_theory, label="reflection theory")
lines!(E_measured_reflection, label="reflection")
axislegend(ax1, position = :lb)
fig_time

