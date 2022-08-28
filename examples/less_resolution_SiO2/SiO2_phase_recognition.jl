
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
m1 = LorentzMedium1D(g, CartesianIndices((2999:5999,)), 1.0, γ_lorentz, ω_0, χ_1, χ_2, χ_3)
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
h_kerr = Array(1:2:3)
degen_injection = binomial_degen.(a - 1, h_injection)
degen_brunel = binomial_degen.(a + 1, h_brunel)
degen_kerr = binomial_degen.(3., h_kerr)
normed_degen_injection = degen_injection./maximum(degen_injection)
normed_degen_brunel = degen_brunel./maximum(degen_brunel)
normed_degen_kerr = degen_kerr./maximum(degen_kerr)

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
tr_b = t_s.(n_real, 1., 0.0, 0.0)
# j_inj = (χ_injection .* (tr_a .* E_arr).^13) 
# P = (χ_1 .* (tr_a .* E_arr) + χ_3.*(tr_a .* E_arr).^3)


E_measured_bulk = d2.Ez[30000:40000]
E_amplitude_spectrum = zeros(Float64, padded_L)
E_amplitude_spectrum[1:size(E_measured_bulk, 1)] = window .* E_measured_bulk
FT_E_bulk = fftshift(fft(E_amplitude_spectrum))

j_brunel = d2.J_Free[30000:40000]
j_inj = d2.J_Tunnel[30000:40000]
j_kerr = d2.J_Bound[30000:40000]
j_all = d2.Jz[30000:40000]

signal_j_inj = zeros(Float64, padded_L)
signal_j_inj[1:size(j_inj, 1)] = window .* j_inj 
FT_J_inj = fftshift(fft(signal_j_inj))

signal_j_brunel = zeros(Float64, padded_L)
signal_j_brunel[1:size(j_brunel, 1)] = window .* j_brunel
FT_J_brunel = fftshift(fft(signal_j_brunel))

signal_j_bound = zeros(Float64, padded_L)
signal_j_bound[1:size(j_kerr, 1)] = window .* j_kerr
FT_J_bound= fftshift(fft(signal_j_bound))

signal_j_all = zeros(Float64, padded_L)
signal_j_all[1:size(j_all, 1)] = window .* j_all
FT_J_all= fftshift(fft(signal_j_all)) 

E_measured_reflection = d4.Ez[34984:44984]
E_amplitude_spectrum_refl = zeros(Float64, padded_L)
E_amplitude_spectrum_refl[1:size(E_measured_reflection, 1)] = window .* E_measured_reflection
FT_E_bulk_refl = fftshift(fft(E_amplitude_spectrum_refl))

function isolate_around_freq2(f::Array{Float64}, spec::Array{Complex{Float64}}, f_isolated::Float64, Δf::Float64, window::Bool=true)
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

#kth - Harmonics 
k = 1.
FT_E_k = isolate_around_freq2(Array(ω), FT_E_bulk, k, 2., false)
FT_Er_k = isolate_around_freq2(Array(ω), FT_E_bulk_refl, k, 2., false)
FT_j_harm_k = isolate_around_freq2(Array(ω), FT_J_all, k, 2., false)
FT_j_inj_k = isolate_around_freq2(Array(ω), FT_J_inj, k, 2., false)
FT_j_brunel_k = isolate_around_freq2(Array(ω), FT_J_brunel, k, 2., false)
FT_j_kerr_k = isolate_around_freq2(Array(ω), FT_J_bound, k, 2., false)
#fith harmoncs in time
E_k_t = ifft(ifftshift(FT_E_k))
Er_k_t = ifft(ifftshift(FT_Er_k))
j_k_t = ifft(ifftshift(FT_j_harm_k))
j_inj_k_t = ifft(ifftshift(FT_j_inj_k))
j_brunel_k_t = ifft(ifftshift(FT_j_brunel_k))
j_kerr_k_t = ifft(ifftshift(FT_j_kerr_k))
# getting phaseinformation
t = Array(1.:1.:padded_L)
t_ref = Float64(argmax(abs.(E_k_t)))
t_ref_2 = Float64(argmax(abs.(Er_k_t)))
z_e, ϕ_E = analytic_signal(real(E_k_t), t, t_ref)
z_er, ϕ_Er = analytic_signal(real(Er_k_t), t, t_ref_2)
z_j, ϕ_j = analytic_signal(real(j_k_t), t, t_ref)
z_j_inj, ϕ_j_inj = analytic_signal(real(j_inj_k_t), t, t_ref)
z_j_brunel, ϕ_j_brunel = analytic_signal(real(j_brunel_k_t), t, t_ref)
z_j_kerr, ϕ_j_kerr = analytic_signal(real(j_kerr_k_t), t, t_ref)

fig = Figure(resolution = (800, 800), font="CMU Serif")
ax1 = Axis(fig[1, 1])
lines!(ax1, abs.(E_k_t)./maximum(abs.(E_k_t)), label=L"E_{Bulk}")
lines!(ax1, abs.(Er_k_t)./maximum(abs.(Er_k_t)), label=L"E_{Refl}")
lines!(ax1, abs.(j_k_t)./maximum(abs.(j_k_t)), label=L"J_{all}")
lines!(ax1, abs.(j_inj_k_t)./maximum(abs.(j_inj_k_t)), label=L"J_{J_inj}")
lines!(ax1, abs.(j_brunel_k_t)./maximum(abs.(j_brunel_k_t)), label=L"J_{Brunel}")
lines!(ax1, abs.(j_kerr_k_t)./maximum(abs.(j_kerr_k_t)), linestyle=:dash, label=L"J_{Kerr}")
axislegend(ax1, position=:rt)
xlims!(ax1, [3000, 7000])
ax2 = Axis(fig[2, 1])
lines!(ax2, ω, real(ϕ_E), label=L"E_{Bulk}")
lines!(ax2, ω, real(ϕ_Er), label=L"E_{Refl}")
lines!(ax2, ω, real(ϕ_j), label=L"J_{all}")
lines!(ax2, ω, real(ϕ_j_inj), label=L"J_{inj}")
lines!(ax2, ω, real(ϕ_j_brunel), label=L"J_{Brunel}")
lines!(ax2, ω, real(ϕ_j_kerr), label=L"J_{Kerr}")
axislegend(ax2, position=:rt)
ylims!(ax2, [-π, π])
ax2.yticks = -π:π/4:π
ax2.ytickformat = ys -> ["$(y/pi)π" for y in ys]
xlims!(ax2, [k-1, k+1])
fig