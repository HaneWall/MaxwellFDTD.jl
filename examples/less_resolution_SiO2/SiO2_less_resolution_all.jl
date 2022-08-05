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

    for (m_idx, m) in enumerate(tunnel_media)
        updatePlasma!(MF, TF[m_idx], f_grid, F, m)
        updateJtunnel!(MF, TF[m_idx], m)
    end

    # for (m_idx, m) in enumerate(tunnel_media)
    #     updatePlasmaTangent!(MF, TF[m_idx], f_grid, F, m, a, Γ̂, Ê)
    #     updateJtunnel!(MF, TF[m_idx], m)
    # end

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
plot_log10_power_spectrum(
    Array(30000*g.Δt:g.Δt:40000*g.Δt),
    hcat(d2.J_Bound[30000:40000],
        d2.J_Tunnel[30000:40000],
        d2.J_Free[30000:40000],
        d2.Jz[30000:40000]),
    ω_central,
    [0.0, 20.0],
    [0.0, 25.0],
    ["J_{Kerr}", "J_{Injection}", "J_{Brunel}", "J_{Overall}"],
    true,
    "ADK_injection_kerr_brunel_all.pdf")

# ##
plot_log10_power_spectrum_current_and_E(
    Array(30000*g.Δt:g.Δt:40000*g.Δt),
    hcat(d2.Jz[30000:40000], d4.Ez[35000:45000], d2.J_Free[30000:40000], d2.J_Bound[30000:40000], d2.J_Tunnel[30000:40000]),
    ω_central,
    [0.0, 20.0],
    [0.0, 1.0],
    ["J_z", "E_{z, Refl}", "J_{Brunel}", "J_{Kerr}", "J_{Injection}"],
    true,
    "Reflection_a_13_less_res_ADK.pdf")

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

## amplitude_spectrum Brunel and Injection 
L = length(t_real)
padded_L = nextpow(2, L)
center_help_L_low = floor(Int64, (padded_L - L)/2)
center_help_L_high = floor(Int64, L + center_help_L_low - 1) 
δt = abs(t[2]-t[1])
δf = 1/δt
ω = 2*π*fftshift(fftfreq(padded_L, δf))./ ω_central
window = blackman(L)
signal_injection = zeros(Float64, padded_L)
signal_brunel = zeros(Float64, padded_L)
signal_injection[center_help_L_low:center_help_L_high] = window.*d2.J_Tunnel[30000:40000]
signal_brunel[center_help_L_low:center_help_L_high] = window.*d2.J_Free[30000:40000]
FT_injection = fftshift(fft(signal_injection))
FT_brunel = fftshift(fft(signal_brunel)) 


## plot degens:
fig_harm = Figure(resolution = (800, 800), font = "CMU Serif")
ax = Axis(fig_harm[1, 1], title="Normed Injection amplitude and Degens")
scatter!(h_injection, (normed_degen_injection))
lines!(ω, (abs.(FT_injection)./maximum(abs.(FT_injection))))
ylims!(ax, [0.0, 1.0])
xlims!(ax, [0.0, 15.0])
ax2 = Axis(fig_harm[2, 1], title="Normed Brunel amplitude and Degens", xlabel = L"ω / ω_{pump}")
lines!(ω, (abs.(FT_brunel)./maximum(abs.(FT_brunel))), color=:red)
scatter!(h_brunel, (normed_degen_brunel), color=:red)
ylims!(ax2, [0.0, 1.0])
xlims!(ax2, [0.0, 15.0])
fig_harm


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