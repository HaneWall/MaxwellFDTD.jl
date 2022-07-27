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
SizeX = 30000
courant = 0.985
Δx = 2e-9
MaxTime = 140000

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
m1 = LorentzMedium1D(g, CartesianIndices((15000:29999,)), 1.0, γ_lorentz, ω_0, χ_1, χ_2, χ_3)
m2 = DrudeMedium1D(g, CartesianIndices((15000:29999,)), γ_plasma, ρ_mol_density)
m3 = TunnelMedium1D(g, CartesianIndices((15000:29999,)), E_gap, ρ_mol_density)

bound_media = [m1]
drude_media = [m2]
tunnel_media = [m3]

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
d2 = PointDetector(CartesianIndex((15000,)), 1, g.MaxTime)
d3 = PointDetector(CartesianIndex((22000,)), 1, g.MaxTime)
d4 = PointDetector(CartesianIndex((505,)), 1, g.MaxTime)
detectors = [d2, d3, d4]

# 7. place sources 
s0 = GaussianWavePointSource(g, CartesianIndex((508,)), true, true, false, amplitude_pump, ceil(500e-15 / g.Δt), t_fwhm, ppw)
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

## 
timeseries_plot(
    d4.Ez,
    "E_{z, refl}",
    "detector_SF.pdf")

##
timeseries_plot(
    d2.Ez,
    "E_{z, first cell medium}",
    "detector_first_cell.pdf")

##
timeseries_plot(
    d3.Ez,
    "E_{z, 7000 cells inside}",
    "detector_mid.pdf")

##
plot_log10_power_spectrum(
    Array(80000*g.Δt:g.Δt:100000*g.Δt),
    hcat(d2.J_Free[80000:100000], d2.J_Bound[80000:100000],
        d2.J_Tunnel[80000:100000]),
    ω_central,
    [0.0, 20.0],
    [0.0, 25.0],
    ["J_{Brunel}", "J_{Kerr}", "J_{Injection}"],
    true,
    "tangentadk_a_13.pdf")

##
plot_log10_power_spectrum_current_and_E(
    Array(80000*g.Δt:g.Δt:100000*g.Δt),
    hcat(d2.Jz[80000:100000], d4.Ez[94500:114500], d2.J_Free[80000:100000], d2.J_Bound[80000:100000], d2.J_Tunnel[80000:100000]),
    ω_central,
    [0.0, 20.0],
    [0.0, 1.0],
    ["J_z", "E_{z, Refl}", "J_{Brunel}", "J_{Kerr}", "J_{Injection}"],
    true,
    "Reflection_tangentadk_a_13.pdf")

## Source for theoretical considerations 
function gaussian(amplitude::Float64, S_c::Float64, Δt::Float64, t::Array{Float64,1}, t_step_peak::Float64, t_fwhm::Float64, ppw::Float64)
    src = zeros(Float64, size(t)[1])
    for idx in 1:size(t)[1]
        src[idx] = amplitude * exp(-2 * log(2) * (((t[idx] - t_step_peak)) * Δt / t_fwhm)^2) * sin(2 * π / ppw * (S_c * t[idx]))
    end
    return src
end

S_C = courant
t_arr = Array(65000.0:1.0:85000.0)
t_peak = 500e-15 / g.Δt
t_real = t_arr .* g.Δt
E_arr = gaussian(amplitude_pump, S_C, g.Δt, t_arr, t_peak, t_fwhm, ppw)

## Theoretical Reflection coefficients - Only Bound Electrons - linear and kerr
function get_kerr_reflection(χ_1::Float64, χ_3::Float64, E::Array{Float64,1})
    E_L = E_reflection_for_χ(E, χ_1, 1.0, 1.0)
    E_K = E_reflection_for_χ(E, χ_3, 3.0, 1.0)
    E_ovr = E_L + E_K
    E = hcat(E_L, E_K, E_ovr)
    return E
end

E_refl_kerr = get_kerr_reflection(χ_1[1], χ_3[1], E_arr)

plot_log10_power_spectrum_current_and_E(
    Array(80000*g.Δt:g.Δt:100000*g.Δt),
    hcat(d2.Jz[80000:100000], E_refl_kerr[:, end], d2.J_Free[80000:100000], d2.J_Bound[80000:100000], d2.J_Tunnel[80000:100000]),
    ω_central,
    [0.0, 20.0],
    [0.0, 1.0],
    ["J_z", "E_{z, Refl kerr}", "J_{Brunel}", "J_{Kerr}", "J_{Injection}"],
    true,
    "Reflection_tangentadk_a13_kerr.pdf")

## Theoretical Reflection coefficients - effective nl
E_atom = 5.14e11
χ_13 = χ_1[1] / E_atom^(a - 1)

function get_effective_reflection(χ_1::Float64, χ_a::Float64, E::Array{Float64,1})
    E_L = E_reflection_for_χ(E, χ_1, 1., 1.0)
    E_K = E_reflection_for_χ(E, χ_a, 13., 1.0)
    E_ovr = E_L + E_K
    E = hcat(E_L, E_K, E_ovr)
    return E
end

E_refl_a = get_effective_reflection(χ_1[1], χ_13, E_arr)

plot_log10_power_spectrum_current_and_E(
    Array(80000*g.Δt:g.Δt:100000*g.Δt),
    hcat(d2.Jz[80000:100000], E_refl_a[:, 2], d2.J_Free[80000:100000], d2.J_Bound[80000:100000], d2.J_Tunnel[80000:100000]),
    ω_central,
    [0.0, 20.0],
    [0.0, 1.0],
    ["J_z", "E_{z, Refl, a}", "J_{Brunel}", "J_{Kerr}", "J_{Injection}"],
    true,
    "Reflection_tangentadk_a13_theory.pdf")