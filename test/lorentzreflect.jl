using MaxwellFDTD
using CairoMakie
CairoMakie.activate!(type = "svg")
using CPUTime
using FFTW

#testing code
CPUtic()
start = time()

# 1. define grid
SizeX = 1000
courant = 0.975
Δx = 4.e-9
MaxTime = 2^14

# 1. define grid
g = Grid1D(SizeX, courant, Δx, MaxTime)

# 2. define fields that exist everywhere
F = Fields1D(g)
MF = MaterialFields1D(g)

# 3. define and place media into the grid
# parameters from bachelor thesis
γ = [0.]
ω_0 = [1.5494e15]
χ_1 = [2.4272]
χ_2 = [30e-12]
χ_3 = [0.]
m1 = LorentzMedium1D(g, CartesianIndices((10:990,)), 1., γ, ω_0, χ_1, χ_2, χ_3)
media = [m1]

# 4. define grid coefficients that respect ϵ_inf from the media 
c_grid = GridCoefficients1D(g, media)

# 5. define fields inside the media
LF1 = LorentzFields1D(m1)
LF = [LF1]

# 6. place detectors 
d1 = LineDetector(CartesianIndices((1:g.SizeX,)), 1, g.MaxTime)
d2 = PointDetector(CartesianIndex((4,)), 1, g.MaxTime)
d3 = PointDetector(CartesianIndex((995,)), 1, g.MaxTime)
d4 = PointDetector(CartesianIndex((11,)), 1, g.MaxTime)
detectors = [d1, d2, d3, d4]

# 7. place sources 
amplitude = intensity2amplitude(10.0^12)

s0 = GaussianWavePointSource(g, CartesianIndex((7,)),false, false, false, amplitude, 8500, 20e-15, 266.)
sources = [s0]

# 8. place boundaries
b1 = LeftSideMurABC(g, CartesianIndex((1,)))
b2 = RightSideMurABC(g, CartesianIndex((SizeX,)))
boundaries = [b1, b2]

for timestep in 1:g.MaxTime

    for (m_idx, m) in enumerate(media)
        updatePNl!(MF, LF[m_idx], F, m)
        updateJ!(MF, LF[m_idx], m, g)
        updateP!(MF, LF[m_idx], m, g)
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
    end
end

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3)," seconds")
println("Computation Complete")

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

    # for (idx, t) in enumerate(11500:1:11540)
    #     lines!(ax2, d3.Jz[t], color=:black, linewidth=1.5)
    #     #ylims!(ax2, -19, -10)
    # end
    lines!(ax2, log10.(abs.(d2.Ez)), color=:black, linewidth=1.5)
    lines!(ax2, log10.(abs.(d3.Ez)), color=:red, linewidth=1.5)
    lines!(ax2, log10.(abs.(d4.Pz)), color=:green, linewidth=1.5)
    lines!(ax2, log10.(abs.(d4.Jz)), color=:yellow, linewidth=1.5)
    lines!(ax2, log10.(abs.(d4.PzNL)), color=:gray50, linewidth=1.5)
    xlims!(ax2, 1740, g.MaxTime)
    f
end

function spectrum_plot()
    Δf = 1/g.Δt
    freqs = fftshift(fftfreq(MaxTime, Δf))
    spectrum_P = fftshift(fft(d4.Pz./MaxTime))
    spectrum_E_reflect = fftshift(fft(d2.Ez./MaxTime))
    spectrum_E_trans = fftshift(fft(d3.Ez./MaxTime))

    f = Figure(resolution = (800, 800))
    ax1 = Axis(f[1, 1],title = "Spectrum P in Media", ylabel = L"|F(P_z)|", xlabel = L"$\omega$")
    lines!(ax1, 2*π*freqs, abs.(spectrum_P)./maximum(abs.(spectrum_P)))
    xlims!(ax1, 0, 5e15)

    ax2 = Axis(f[1, 2],title = "TimeSeries P in Media", ylabel = L"P_z", xlabel = L"$t$")
    lines!(ax2, d4.Pz)

    ax3 = Axis(f[2, 1],title = "Spectrum E-Reflection", ylabel = L"|$F(E_z)$|", xlabel = L"$\omega$")
    lines!(ax3, 2*π*freqs, abs.(spectrum_E_reflect)./maximum(abs.(spectrum_E_reflect)))
    xlims!(ax3, 0, 5e15)

    ax4 = Axis(f[2, 2],title = "TimeSeries E Reflection", ylabel = L"E_z", xlabel = L"$t$")
    lines!(ax4, d2.Ez)
    f
end