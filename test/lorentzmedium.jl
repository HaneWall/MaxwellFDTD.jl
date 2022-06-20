using MaxwellFDTD
using CairoMakie
CairoMakie.activate!(type = "svg")
using CPUTime
using FFTW

#testing code
CPUtic()
start = time()

# 1. define grid
SizeX = 2000
courant = 0.975
Δx = 4.e-9
MaxTime = 30000

# 1. define grid
g = Grid1D(SizeX, courant, Δx, MaxTime)

# 2. define fields that exist everywhere
F = Fields1D(g)
MF = MaterialFields1D(g)

# 3. define and place media into the grid
# parameters from bachelor thesis
γ = [0., 0., 0.]
ω_0 = [1.5494e16, 7.9514e15, 9.7766e13]
χ_1 = [2.4272, 1.4617, 9.6536]
χ_2 = [30e-12, 0., 0.]
χ_3 = [0., 0., 0.]
m1 = LorentzMedium1D(g, CartesianIndices((1700:1900,)), 1., γ, ω_0, χ_1, χ_2, χ_3)
media = [m1]

# 4. define grid coefficients that respect ϵ_inf from the media 
c_grid = GridCoefficients1D(g, media)

# 5. define fields inside the media
LF1 = LorentzFields1D(m1)
LF = [LF1]

# 6. place detectors 
d1 = LineDetector(CartesianIndices((1:g.SizeX,)), 1, g.MaxTime)
d2 = PointDetector(CartesianIndex((5,)), 1, g.MaxTime)
d3 = PointDetector(CartesianIndex((1751,)), 1, g.MaxTime)
detectors = [d1, d2, d3]

# 7. place sources 
s0 = GaussianWavePointSource(g, CartesianIndex((100,)),true, true, false, 1., 8500, 20e-15, 266.)
s1 = GaussianPointSource(CartesianIndex((20,)),true, true, false, 1., 90, 500.)
s2 = SinusoidalPointSource(g, CartesianIndex((5,)), true, false, 1., 266.)
s3 = RickerPointSource(g, CartesianIndex((15,)), true, false, 0.5, 1000, 266.)
sources = [s0]

# 8. place boundaries
b1 = LeftSideMurABC(g, CartesianIndex((1,)))
b2 = RightSideMurABC(g, CartesianIndex((SizeX,)))
boundaries = [b1, b2]

for timestep in 1:g.MaxTime

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
    end
end

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3)," seconds")
println("Computation Complete")


function epsilon_complex(γ, ω_0, χ_1, ω)
    eps_complex = zeros(ComplexF64, length(ω)) .+ 1.
    for idx = 1:length(γ)
        @. eps_complex += χ_1[idx] * (ω_0[idx]^2) / (ω_0[idx]^2 - ω^2 + 1im * γ[idx] * ω)
    end
    return eps_complex
end

function waterfall_plot(ts_min, ts_max)
    f = Figure(resolution = (800, 1200))
    ax1 = Axis(f[1, 1],title = "Waterfall Plot", ylabel = "timestep%10", xlabel = L"E_z")
    for medium in media
        vspan!(ax1, first(medium.location)[1], last(medium.location)[1], color=:gray90)
    end

    for (idx, t) in enumerate(ts_min:100:ts_max)
        lines!(ax1, first(d1.location)[1]:last(d1.location)[1], d1.Ez[t, :] .+ idx, color=:black, linewidth=1.5)
    end

    ax2 = Axis(f[2, 1],title = "Waterfall Plot", ylabel = "timestep%10", xlabel = L"timestep")

    lines!(ax2, log10.(abs.(d3.Pz)), color=:black, linewidth=1.5)
    lines!(ax2, log10.(abs.(d3.Jz)), color=:red, linewidth=1.5)
    lines!(ax2, log10.(abs.(d2.Ez)), color=:gray50, linewidth=1.5)
    lines!(ax2, log10.(abs.(d3.Ez)), color=:blue, linewidth=1.5)
    lines!(ax2, log10.(abs.(d3.PzNl)), color=:green, linewidth=1.5)
    xlims!(ax2, 1740, g.MaxTime)
    f
end

function spectrum_plot()
    Δf = 1/g.Δt
    spectrum = fftshift(fft(d2.Ez))
    spectrum_P = fftshift(fft(d3.Pz))
    spectrum_J = fftshift(fft(d3.Jz))
    freqs = fftshift(fftfreq(MaxTime, Δf))
    f1 = Figure(resolution = (800, 800))
    ax2 = Axis(f1[1, 1])
    lines!(ax2, 2*π*freqs, abs.(spectrum)./maximum(abs.(spectrum)))
    lines!(ax2, 2*π*freqs, abs.(spectrum_J./maximum(abs.(spectrum_J))))
    #lines!(ax2, 2*π*freqs, abs.(spectrum_P)./maximum(abs.(spectrum_P)))
    lines!(ax2, 2*π*freqs, real.(epsilon_complex(γ, ω_0, χ_1, 2*π*freqs)))
    #lines!(ax2, 2*π*freqs, imag.(epsilon_complex(γ, ω_0, χ_1, 2*π*freqs)))
    xlims!(ax2, 0, 17.6e15)
    ylims!(ax2, -10, 10)
    f1
end