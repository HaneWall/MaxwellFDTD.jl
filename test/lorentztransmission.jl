using MaxwellFDTD
using CairoMakie
CairoMakie.activate!(type = "svg")
using CPUTime
using FFTW
using ProgressBars
using DSP

#testing code
CPUtic()
start = time()

# 1. define grid
SizeX = 500
courant = 0.975
Δx = 1/3 * 3.91e-07
MaxTime = 2^15


# Bachelor Parameters Lorentz Slab 
λ = 1.75e-5
ω_central = 2 * π * c_0 / λ
t_fwhm = 600e-15
ppw = λ / Δx
amplitude = 2.
γ = [8e12, 9e14]
ω_0 = [1.2566e14, 1.2e13]
χ_1 = [2.1, 2.4]
χ_2 = [30e-12, 0.]
χ_3 = [0., 0.]

function epsilon_complex(γ, ω_0, χ_1, ω)
    eps_complex = zeros(ComplexF64, length(ω)) .+ 1.
    for idx = 1:length(γ)
        @. eps_complex += χ_1[idx] * (ω_0[idx]^2) / (ω_0[idx]^2 - ω^2 + 1im * γ[idx] * ω)
    end
    return eps_complex
end

# 1. define grid
g = Grid1D(SizeX, courant, Δx, MaxTime)
t = g.Δt:g.Δt:g.Δt*MaxTime

# 2. define fields that exist everywhere
F = Fields1D(g)
MF = MaterialFields1D(g)

m1 = LorentzMedium1D(g, CartesianIndices((300:450,)), 1., γ, ω_0, χ_1, χ_2, χ_3)
media = [m1]

# 4. define grid coefficients that respect ϵ_inf from the media 
c_grid = GridCoefficients1D(g, media)

# 5. define fields inside the media
LF1 = LorentzFields1D(m1)
LF = [LF1]

# 6. place detectors 
d1 = LineDetector(CartesianIndices((1:g.SizeX,)), 1, g.MaxTime)
d2 = PointDetector(CartesianIndex((3,)), 1, g.MaxTime)
d3 = PointDetector(CartesianIndex((460,)), 1, g.MaxTime)
d4 = PointDetector(CartesianIndex((300,)), 1, g.MaxTime)
detectors = [d1, d2, d3, d4]

# 7. place sources 
s0 = GaussianWavePointSource(g, CartesianIndex((50,)),true, true, false, amplitude, 8500, 600e-15, ppw)
sources = [s0]

# 8. place boundaries
b1 = LeftSideMurABC(g, CartesianIndex((1,)))
b2 = RightSideMurABC(g, CartesianIndex((SizeX,)))
boundaries = [b1, b2]

for timestep in ProgressBar(1:g.MaxTime)
    
    for (m_idx, m) in enumerate(media)
        updatePNl!(MF, LF[m_idx], F, m)
        updateJbound!(MF, LF[m_idx], m, g)
        updatePbound!(MF, LF[m_idx], m, g)
    end

    updateH!(F, g, c_grid)

    for source in sources
        sourceH!(source, F, timestep)
    end

    for b in boundaries
        saveFields!(b, F)
    end
    #ABC!(F, g)

    updateE!(F, MF, g, c_grid)

    for source in sources
        sourceE!(source, F, timestep)
    end

    for b in boundaries
        stepABC!(F, b)
    end

    for d in detectors 
        safeE!(d, F, timestep)
        safeP!(d, MF, timestep)
        safeJ!(d, MF, timestep)
        safePNl!(d, MF, timestep)
    end
end

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3)," seconds")
println("Computation Complete")


function semilog(a::Float64)
    return sign(a) * log10(abs(a))
end

function semilog_timeseries()
    f = Figure(resolution = (800, 800))
    
    ax1 = Axis(f[1, 1], 
        title = "First Cell Medium", 
        ylabel = L"|(P_z)|", 
        xlabel = "t in ps", 
        ylabelsize = 18, 
        xlabelsize = 18, 
        yminorticksvisible=true,
        yminorticks=IntervalsBetween(9),
        yscale=log10,
        xgridstyle = :dash, 
        ygridstyle = :dash, 
        xtickalign = 1,
        xticksize = 8, 
        ytickalign = 1, 
        yticksize = 8, 
        xlabelpadding = -8)

    lines!(ax1,  t*10^12, abs.(d4.Pz) .+ 10e-12)
    f
end

function spectrum_plot()
    Δf = 1/g.Δt
    freqs = fftshift(fftfreq(MaxTime, Δf))
    harmonic_order = 2 * π * freqs./ω_central
    spectrum_P = fftshift(fft(d4.Pz))
    spectrum_E_reflect = fftshift(fft(d2.Ez))
    spectrum_E_trans = fftshift(fft(d3.Ez))

    broadness = Int(ceil(t_fwhm*13/g.Δt))
    broad_idx_mean = ceil(broadness/2)
    signal_p_idx_max = argmax(abs.(d4.Pz))
    shift_p = Int(signal_p_idx_max - broad_idx_mean)
    window_p =  blackman(broadness; padding=length(t) - broadness)
    window_p = circshift(window_p, shift_p)

    spectrum_E_reflect_window = fftshift(fft(d2.Ez .* window_p))

    ϵ_real = real.(epsilon_complex(γ, ω_0, χ_1, 2*π*freqs))
    ϵ_imag = imag.(epsilon_complex(γ, ω_0, χ_1, 2*π*freqs)) 
    
    n_r =  zeros(Float64,length(freqs))
    n_i =  zeros(Float64,length(freqs))


    @. n_r = sqrt(1/2 *(sqrt(ϵ_real^2 + ϵ_imag^2) + ϵ_real))
    @. n_i = sqrt(1/2 *(sqrt(ϵ_real^2 + ϵ_imag^2) - ϵ_real))


    f = Figure(resolution = (800, 800))
    
    ax1 = Axis(f[1, 1],
                title = "First Cell Medium", 
                ylabel = L"\log_{10}|F(P_z)|", 
                xlabel = L"$\omega / \omega_{central}$", 
                ylabelsize = 18, 
                xlabelsize = 18, 
                xgridstyle = :dash, 
                ygridstyle = :dash, 
                xtickalign = 1,
                xticksize = 8, 
                ytickalign = 1, 
                yticksize = 8, 
                xlabelpadding = -8)
    
    lines!(ax1, harmonic_order, log10.(abs.(spectrum_P)./maximum(abs.(spectrum_P))))
    xlims!(ax1, 0, 4)

    ax2 = Axis(f[1, 2],title = "Time Series P First Cell Medium", ylabel = L"P_z", xlabel = "t in ps")
    lines!(ax2, t*10^12, d4.Pz./maximum(d4.Pz))
    lines!(ax2, t*10^12, window_p)

    ax3 = Axis(f[2, 1],
                title = L"E_{Reflection}", 
                ylabel = L"\log_{10}|F(E_z)|", 
                xlabel = L"$\omega / \omega_{central}$",  
                ylabelsize = 18, 
                xlabelsize = 18, 
                xgridstyle = :dash, 
                ygridstyle = :dash, 
                xtickalign = 1,
                xticksize = 8, 
                ytickalign = 1, 
                yticksize = 8, 
                xlabelpadding = -8)
    lines!(ax3, harmonic_order, log10.(abs.(spectrum_E_reflect./MaxTime)))
    lines!(ax3, harmonic_order, log10.(abs.(spectrum_E_reflect_window./MaxTime)))
    xlims!(ax3, 0, 4)

    ax4 = Axis(f[2, 2],title = "Time Series E Reflection", ylabel = L"E_z", xlabel = "t in ps")
    lines!(ax4, t*10^12, d2.Ez)
    #lines!(ax4, t*10^12, d2.Ez.* window_p)

    ax5 = Axis(f[3, 1],
                title = L"E_{Transmission}", 
                ylabel = L"\log_{10}|F(E_z)|", 
                xlabel = L"$\omega / \omega_{central}$", 
                ylabelsize = 18, 
                xlabelsize = 18, 
                xgridstyle = :dash, 
                ygridstyle = :dash, 
                xtickalign = 1,
                xticksize = 8, 
                ytickalign = 1, 
                yticksize = 8, 
                xlabelpadding = -8)

    lines!(ax5, harmonic_order, log10.(abs.(spectrum_E_trans./MaxTime)))
    xlims!(ax5, 0, 4)

    ax6 = Axis(f[3, 2],title = "Time Series E Transmission", ylabel = L"E_z", xlabel = "t in ps")
    lines!(ax6,t*10^12, d3.Ez)

    ax7 = Axis(f[4, 1:2],title = "Refraction Index")
    lines!(ax7, harmonic_order, n_r, label=L"n_r")
    lines!(ax7, harmonic_order, n_i, label=L"\kappa")
    axislegend(L"n = n_r + i \kappa"; position = :rt, bgcolor = (:grey90, 0.25), nbanks = 2)
    xlims!(ax7, 0, 4)

    f
end

function waterfall_plot(ts_min, ts_max)
    f = Figure(resolution = (800, 1200))
    ax1 = Axis(f[1, 1],title = "Waterfall Plot", ylabel = "timestep%10", xlabel = L"E_z")
    for medium in media
        vspan!(ax1, first(medium.location)[1], last(medium.location)[1], color=:gray90)
    end

    for (idx, t) in enumerate(ts_min:100:ts_max)
        lines!(ax1, first(d1.location)[1]:last(d1.location)[1], d1.Ez[t, :]./amplitude .+ idx, color=:black, linewidth=1.5)
    end

    ax2 = Axis(f[2, 1],title = "Waterfall Plot", ylabel = "timestep%10", xlabel = L"timestep")

    lines!(ax2, log10.(abs.(d2.Ez)), color=:black, linewidth=1.5)
    lines!(ax2, log10.(abs.(d3.Ez)), color=:red, linewidth=1.5)
    lines!(ax2, log10.(abs.(d4.Pz)), color=:green, linewidth=1.5)
    lines!(ax2, log10.(abs.(d4.Jz)), color=:yellow, linewidth=1.5)
    lines!(ax2, log10.(abs.(d4.PzNl)), color=:blue, linewidth=1.5)
    xlims!(ax2, 1740, g.MaxTime)
    f
end