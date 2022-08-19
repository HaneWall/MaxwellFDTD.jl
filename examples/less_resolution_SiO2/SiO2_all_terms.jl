
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
    j_inj = (χ_injection .* (tr_a .* E_arr).^13) 
    P = (χ_1 .* (tr_a .* E_arr) + χ_3.*(tr_a .* E_arr).^3)


    E_measured_bulk = d2.Ez[30000:40000]
    E_amplitude_spectrum = zeros(Float64, padded_L)
    E_amplitude_spectrum[center_help_L_low:center_help_L_high] = window .* E_measured_bulk
    FT_E_bulk = fftshift(fft(E_amplitude_spectrum))

    j_brunel = d2.J_Free[30000:40000]
    j_inj = d2.J_Tunnel[30000:40000]
    j_kerr = d2.J_Bound[30000:40000]
    j_all = d2.Jz[30000:40000]


    signal_injection_theory = zeros(Float64, padded_L)
    signal_injection_theory[center_help_L_low:center_help_L_high] = window .* P_1
    
    signal_j_inj = zeros(Float64, padded_L)
    signal_j_inj[center_help_L_low:center_help_L_high] = window .* j_inj 
    FT_J_inj = fftshift(fft(signal_j_inj))

    signal_j_brunel = zeros(Float64, padded_L)
    signal_j_brunel[center_help_L_low:center_help_L_high] = window .* j_brunel
    FT_J_brunel = fftshift(fft(signal_j_brunel))

    signal_j_bound = zeros(Float64, padded_L)
    signal_j_bound[center_help_L_low:center_help_L_high] = window .* j_kerr
    FT_J_bound= fftshift(fft(signal_j_bound))
    
    signal_j_all = zeros(Float64, padded_L)
    signal_j_all[center_help_L_low:center_help_L_high] = window .* j_all
    FT_J_all= fftshift(fft(signal_j_all)) 

    E_measured_reflection = d4.Ez[35000:45000]
    E_amplitude_spectrum_refl = zeros(Float64, padded_L)
    E_amplitude_spectrum_refl[center_help_L_low:center_help_L_high] = window .* E_measured_reflection
    FT_E_bulk_refl = fftshift(fft(E_amplitude_spectrum_refl))
    
    C = 0.25*10e-5
    idx_first_harm = argmin(abs.(ω .- 1))
    injection_amplitude = (C.*abs.(FT_J_inj)./(padded_L))[idx_first_harm]
    kerr_amplitude = (C.*abs.(FT_J_bound)./(padded_L))[idx_first_harm] 
    idx_fifth_harm = argmin(abs.(ω .- 5))
    injection_amplitude_fifth =((C.*abs.(FT_J_inj)./(padded_L))/5)[idx_fifth_harm] 
    E_re_amplitude_fifth = (abs.(FT_E_bulk_refl)./padded_L)[idx_fifth_harm]
    offset = injection_amplitude_fifth/E_re_amplitude_fifth
    
    fig_relative_amplitudes = Figure(resolution = (800, 800), font="CMU Serif") 
    ax1 = Axis(fig_relative_amplitudes[1, 1], xlabel=L"\omega / \omega_{pump}", title="log10 - Amplituden", yscale=log10) 
    lines!(ω, C.*(abs.(FT_J_bound)./(abs.(ω) .* padded_L) .+ 0.01), label=L"C \cdot (j_{lin} + j_{kerr}) / \omega")
    lines!(ω, C.*(abs.(FT_J_brunel)./(abs.(ω) .* padded_L) .+ 0.01), label=L"C \cdot j_{brunel} / \omega")
    lines!(ω, C.*(abs.(FT_J_inj)./(abs.(ω) .* padded_L) .+ 0.01), label=L"C \cdot j_{injection} / \omega")
    lines!(ω, C.*(abs.(FT_J_all)./(abs.(ω) .* padded_L) .+ 0.01), label=L"C \cdot j_{all} / \omega")
    scatter!(h_injection, (injection_amplitude*normed_degen_injection)./h_injection, color=:green, label=L"\Pi_{injection}")
    scatter!(h_injection, (injection_amplitude*normed_degen_injection)./h_injection .* (1/offset), color=:red, label=L"C_2 \cdot \Pi_{injection}")
    lines!(ω, (abs.(FT_E_bulk)./padded_L .+ 0.01), label=L"E_{tr}")
    lines!(ω, (abs.(FT_E_bulk_refl)./padded_L .+ 0.01), color=:red, linestyle=:dash, label=L"E_{re}")
    xlims!(ax1, [0.0, 13.0])
    ylims!(ax1, [10^(0) , 10^(10)])
    axislegend(ax1, position=:rt, nbanks=2)
    # ax2 = Axis(fig_relative_amplitudes[2, 1], xlabel=L"\omega / \omega_{pump}", title="Fresnelähnlich")
    # lines!(ω, (abs.(FT_E_bulk_refl)./padded_L .+ 0.01)./ (C.*(abs.(FT_J_inj)./(abs.(ω) .* padded_L) .+ 0.01)))
    # lines!(ω, (abs.(FT_E_bulk_refl)./padded_L .+ 0.01)./ (abs.(FT_E_bulk)./padded_L .+ 0.01))
    # lines!(ω, tr_b)
    # xlims!(ax2, [1., 11.])
    # ylims!(ax2, [0, 1.3])
    fig_relative_amplitudes
    save("reflection_all_with_injection_permutations_ADK.pdf", fig_relative_amplitudes)
    
    #save("amplituden_all_ref_tangent.pdf", fig_relative_amplitudes)

    ##relations trasnmission and refelction 
    fig_relations = Figure(resolution = (800, 800), font = "CMU Serif")
    ax1 = Axis(fig_relations[1, 1], xlabel= L"\omega / \omega_{pump}", title=L"E_{tr} / E_{re}")
    lines!(ω, abs.(FT_E_bulk_refl)./abs.(FT_E_bulk), color=:black)
    xlims!(ax1, [0.0, 11.0])
    ylims!(ax1, [0.9, 1.15])
    fig_relations
    ##comparing phases
    fig_time = Figure(resolution = (800, 800), font="CMU Serif", fontsize=13)
    ax1 = Axis(fig_time[1, 1], xlabel=L"arb. t", title="Amplituden") 
    lines!(t_arr, j_kerr./maximum(j_kerr), label=L"J_{Kerr} + J_{Lin}")
    lines!(t_arr, j_brunel./maximum(abs.(j_brunel)), label=L"J_{Brunel}")
    lines!(t_arr, j_inj./maximum(abs.(j_inj)), label=L"J_{Injection}")   
    lines!(t_arr, E_measured_bulk./maximum(E_measured_bulk), label=L"E_{bulk}")
    #lines!(t_arr, E_measured_reflection)  
    axislegend(ax1, position=:rt, nbanks=2)
    xlims!(ax1, [28000, 32000])
    fig_time
    #save("time_series_first_cell.pdf", fig_time)