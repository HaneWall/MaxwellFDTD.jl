using MaxwellFDTD
using CPUTime
using FFTW
using CairoMakie
using ProgressBars
using LaTeXStrings
using DSP
CairoMakie.activate!(type="svg")

CPUtic()
start = time()

# 1. define grid
SizeX = 400
courant = 0.985
Δx = 2e-9
MaxTime = 2^18

# 1. define grid
g = Grid1D(SizeX, courant, Δx, MaxTime)
t = g.Δt:g.Δt:g.Δt*MaxTime

# varin paramters Si02
ρ_mol_density = 2.2e28
# bound electrons
γ_lorentz = [0.]
ω_0 = [2.75e16] # this might not work, use 2*π*/g.Δt instead (old varin paper/bachelor thesis)
χ_1 = [1.1025]
χ_2 = [0.]
χ_3 = [2.2e-22]

# drude parameters
γ_plasma = 1e15

# tunnel parameters
E_gap = 7.5

# source parameters
λ = 2100e-9
ω_central = 2 * π * c_0 / λ
ppw = λ/Δx
t_fwhm = 140e-15 # intensity FWHM
amplitude_pump = intensity2amplitude(12e16) # 12TWcm^-2

λ_probe = 800e-9
ω_probe= 2 * π * c_0 / λ_probe
ppw_probe = λ_probe/Δx
t_fwhm_probe = 45e-15 # intensity FWHM
amplitude_probe = intensity2amplitude(1.5e14) 


# 2. define fields that exist everywhere
F = Fields1D(g)
MF = MaterialFields1D(g)

m1 = LorentzMedium1D(g, CartesianIndices((100:350,)), 1., γ_lorentz, ω_0, χ_1, χ_2, χ_3)
m2 = DrudeMedium1D(g, CartesianIndices((100:350,)), γ_plasma, ρ_mol_density)
m3 = TunnelMedium1D(g, CartesianIndices((100:350,)), E_gap, ρ_mol_density)

bound_media= [m1]
drude_media = [m2]
tunnel_media =  [m3]

# 4. define grid coefficients that respect ϵ_inf from the media 
c_grid = GridCoefficients1D(g, bound_media)
f_grid  = FieldIonizationCoefficients1D(g)

# 5. define fields inside the media
LF1 = LorentzFields1D(m1)
LF = [LF1]
DF1 = DrudeFields1D(m2)
DF = [DF1]
TF1 = TunnelFields1D(m3)
TF = [TF1]

# 6. place detectors 
d1 = LineDetector(CartesianIndices((1:g.SizeX,)), 1, g.MaxTime)
d2 = PointDetector(CartesianIndex((3,)), 1, g.MaxTime)
d3 = PointDetector(CartesianIndex((360,)), 1, g.MaxTime)
d4 = PointDetector(CartesianIndex((100,)), 1, g.MaxTime)
detectors = [d1, d2, d3, d4]

# 7. place sources 
s0 = GaussianWavePointSource(g, CartesianIndex((50,)),true, true, false, amplitude_pump, ceil(500e-15/g.Δt), t_fwhm, ppw)
s1 = GaussianWavePointSource(g, CartesianIndex((50,)),true, true, false, amplitude_probe, ceil(500e-15/g.Δt), t_fwhm_probe, ppw_probe)
sources = [s0]

# 8. place boundaries
b1 = LeftSideMurABC(g, CartesianIndex((1,)))
b2 = RightSideMurABC(g, CartesianIndex((SizeX,)))
boundaries = [b1, b2]

for timestep in ProgressBar(1:g.MaxTime)
    
    for (m_idx, m) in enumerate(tunnel_media)
        updatePlasma!(MF, TF[m_idx], f_grid, F, m)
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

    updateH!(F, g, c_grid)

    for source in sources
        sourceH!(source, F, timestep)
    end

    for b in boundaries
        saveFields!(b, F)
    end

    updateE!(F, MF, g, c_grid)

    for source in sources
        sourceE!(source, F, timestep)
    end

    for b in boundaries
        stepABC!(F, b)
    end

    for d in detectors 
        safeΓ_ADK!(d, MF, timestep)
        safeJ_tunnel!(d, MF, timestep)
        safeJ_free!(d, MF, timestep)
        safeJ_bound!(d, MF, timestep)
        safeE!(d, F, timestep)
        safeP!(d, MF, timestep)
        safeJ!(d, MF, timestep)
        safePNl!(d, MF, timestep)
    end
end

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3)," seconds")
println("Computation Complete")



function spectrum_plot()
    Δf = 1/g.Δt
    freqs = fftshift(fftfreq(MaxTime, Δf))
    harmonic_order = 2 * π * freqs./ω_central

    broadness = Int(ceil(6*t_fwhm/g.Δt))
    broad_idx_mean = ceil(broadness/2)
    signal_p_idx_max = argmax(abs.(d4.Jz))
    shift_p = Int(signal_p_idx_max - broad_idx_mean)
    window_p =  blackman(broadness; padding=length(t) - broadness)
    window_p = circshift(window_p, shift_p)
    #window_p = 1.

    spectrum_J_Tunnel = fftshift(fft(window_p .* d4.J_Tunnel))
    spectrum_J_Bound = fftshift(fft(window_p .* d4.J_Bound))
    spectrum_J_Free = fftshift(fft(window_p .* d4.J_Free))
    spectrum_J = fftshift(fft(window_p .* d4.Jz))
    spectrum_E_reflect = fftshift(fft(window_p .* d2.Ez))
    spectrum_E_trans = fftshift(fft(window_p .*d3.Ez))

    f = Figure(resolution = (1300, 1300))
    ax1 = Axis(f[1, 1],
            title = "First Cell Medium", 
            ylabel = L"$\log_{10}|F(P_z)|^2$", 
            xlabel = L"$\frac{\omega}{\omega_{central}}$", 
            ylabelsize = 18, 
            xlabelsize = 18, 
            xgridstyle = :dash, 
            ygridstyle = :dash, 
            xtickalign = 1,
            xticksize = 8, 
            ytickalign = 1, 
            yticksize = 8, 
            xlabelpadding = -8)
    
    lines!(ax1, harmonic_order, log10.(abs.(spectrum_J_Bound).^2 ./maximum(abs.(spectrum_J).^2)), label =L"j_{Kerr}")
    lines!(ax1, harmonic_order, log10.(abs.(spectrum_J_Tunnel).^2 ./maximum(abs.(spectrum_J).^2)), label = L"j_{Inject}")
    lines!(ax1, harmonic_order, log10.(abs.(spectrum_J_Free).^2 ./maximum(abs.(spectrum_J).^2)), label = L"j_{Brunel}")
    lines!(ax1, harmonic_order, log10.(abs.(spectrum_J).^2 ./maximum(abs.(spectrum_J).^2)),linestyle = :dash, label = L"j_{Ovr}")
    # lines!(harmonic_order, log10.(abs.(spectrum_J_Bound).^2 ./maximum(abs.(spectrum_J).^2)))
    # lines!(harmonic_order, log10.(abs.(spectrum_J_Tunnel).^2 ./maximum(abs.(spectrum_J).^2)))
    # lines!(harmonic_order, log10.(abs.(spectrum_J_Free).^2 ./maximum(abs.(spectrum_J).^2)))
    # lines!(harmonic_order, log10.(abs.(spectrum_J).^2 ./maximum(abs.(spectrum_J).^2)))
    xlims!(ax1, 0, 15)
    ylims!(ax1, -20, 0)
    axislegend(ax1, position=:lb, orientation=:horizontal)

    ax2 = Axis(f[1, 2],title = "Time Series P First Cell Medium", ylabel = L"J_z", xlabel = "t in fs")
    #lines!(ax2, t*10^15, (d4.J_Bound)./maximum(d4.J_Bound))
    #lines!(ax2, t*10^15, (d4.J_Tunnel)./maximum(d4.J_Tunnel))
    #lines!(ax2, t*10^15, (d4.J_Free)./maximum(d4.J_Free))
    lines!(ax2, t*10^15, (d4.Jz)./maximum(d4.Jz))
    
    lines!(ax2, t*10^15, window_p, color=:black)

    ax3 = Axis(f[2, 1],
                title = L"Reflection", 
                ylabel = L"$\log_{10}|F(E_z)|^2$", 
                xlabel = L"$\frac{\omega}{\omega_{central}}$",  
                ylabelsize = 18, 
                xlabelsize = 18, 
                xgridstyle = :dash, 
                ygridstyle = :dash, 
                xtickalign = 1,
                xticksize = 8, 
                ytickalign = 1, 
                yticksize = 8, 
                xlabelpadding = -8)
    lines!(ax3, harmonic_order, log10.(abs.(spectrum_E_reflect).^2 ./maximum(abs.(spectrum_E_reflect).^2)), label=L"E_{Reflection}")
    #lines!(ax3, harmonic_order, log10.(abs.(spectrum_J_Bound)./maximum(abs.(spectrum_J_Bound))))
    #lines!(ax3, harmonic_order, log10.(abs.(spectrum_J_Tunnel)./maximum(abs.(spectrum_J_Tunnel))))
    #lines!(ax3, harmonic_order, log10.(abs.(spectrum_J_Free)./maximum(abs.(spectrum_J_Free))))
    lines!(ax3, harmonic_order, log10.(abs.(spectrum_J).^2 ./maximum(abs.(spectrum_J).^2)), label=L"J_{Ovr}")
    axislegend(ax3, position=:lb, orientation=:horizontal)
    xlims!(ax3, 0, 15)
    ylims!(ax3, -20, 0)

    ax4 = Axis(f[2, 2],title = "Time Series E Reflection", ylabel = L"E_z", xlabel = "t in fs")
    lines!(ax4, t*10^15, d2.Ez)
    lines!(ax4, t*10^15, d2.Ez.* window_p, color=:black)

    ax5 = Axis(f[3, 1],
                title = L"Transmission", 
                ylabel = L"\log_{10}|F(E_z)|^2", 
                xlabel = L"$\frac{\omega}{\omega_{central}}$", 
                ylabelsize = 18, 
                xlabelsize = 18, 
                xgridstyle = :dash, 
                ygridstyle = :dash, 
                xtickalign = 1,
                xticksize = 8, 
                ytickalign = 1, 
                yticksize = 8, 
                xlabelpadding = -8)

    lines!(ax5, harmonic_order, log10.(abs.(spectrum_E_trans).^2 ./ maximum(abs.(spectrum_E_trans).^2)), label=L"E_{Transmission}")
    lines!(ax5, harmonic_order, log10.(abs.(spectrum_J).^2 ./maximum(abs.(spectrum_J).^2)), label=L"J_{Ovr}")
    axislegend(ax5, position=:lb, orientation=:horizontal)
    xlims!(ax5, 0, 15)
    ylims!(ax5, -20, 0)

    ax6 = Axis(f[3, 2],title = "Time Series E Transmission", ylabel = L"E_z", xlabel = "t in fs")
    lines!(ax6,t*10^15, d3.Ez)
    f
end


function waterfall_plot(ts_min, ts_max)
    f = Figure(resolution = (800, 1200))
    ax1 = Axis(f[1, 1],title = "Waterfall Plot", ylabel = "timestep%10", xlabel = L"E_z")
    for medium in bound_media
        vspan!(ax1, first(medium.location)[1], last(medium.location)[1], color=:gray90)
    end

    for (idx, t) in enumerate(ts_min:100:ts_max)
        lines!(ax1, first(d1.location)[1]:last(d1.location)[1], d1.Ez[t, :]./amplitude .+ idx, color=:black, linewidth=1.5)
    end

    ax2 = Axis(f[2, 1],title = "Waterfall Plot", ylabel = "timestep%10", xlabel = L"timestep")

    lines!(ax2, log10.(abs.(d2.Ez)), color=:black, linewidth=1.5)
    lines!(ax2, log10.(abs.(d3.Ez)), color=:red, linewidth=1.5)
    lines!(ax2, log10.(abs.(d4.J_Free)), color=:green, linewidth=1.5)
    lines!(ax2, log10.(abs.(d4.J_Tunnel)), color=:yellow, linewidth=1.5)
    lines!(ax2, log10.(abs.(d4.J_Bound)), color=:blue, linewidth=1.5)
    lines!(ax2, log10.(d4.Γ_ADK), color=:gray90, linewidth=1.5)
    xlims!(ax2, 1740, g.MaxTime)
    f
end