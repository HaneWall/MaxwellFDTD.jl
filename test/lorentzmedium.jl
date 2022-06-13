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
courant = 1.
Δx = 4.e-9
MaxTime = 11900

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
χ_2 = [0., 0., 0.]
χ_3 = [0., 0., 0.]
m1 = LorentzMedium1D(g, CartesianIndices((1700:1900,)), 1.3, γ, ω_0, χ_1, χ_2, χ_3)
media = [m1]

# 4. define grid coefficients that respect ϵ_inf from the media 
c_grid = GridCoefficients1D(g, media)

# 5. define fields inside the media
LF1 = LorentzFields1D(m1)
LF = [LF1]

# 6. place detectors 
d1 = LineDetector(CartesianIndices((1:g.SizeX,)), 1, g.MaxTime)
d2 = PointDetector(CartesianIndex((5,)), 1, g.MaxTime)
detectors = [d1, d2]

# 7. place sources 
s0 = GaussianWavePointSource(g, CartesianIndex((5,)),true, false, false, 1., 10500, 5e-15, 266.)
s1 = GaussianPointSource(CartesianIndex((20,)),true, true, false, 1., 90, 500.)
s2 = SinusoidalPointSource(g, CartesianIndex((5,)), true, false, 1., 25.)
s3 = RickerPointSource(g, CartesianIndex((15,)), true, false, 0.5, 1000, 266.)
sources = [s0]

for timestep in 1:g.MaxTime

    for (m_idx, m) in enumerate(media)
        updateJ!(MF, LF[m_idx], m, g)
        updatePNl!(LF[m_idx], F, m)
        updateP!(MF, LF[m_idx], m, g)
    end

    updateH!(F, g, c_grid)

    for source in sources
        sourceH!(source, F, timestep)
    end

    ABC!(F, g)

    updateE!(F, MF, g, c_grid)

    for source in sources
        sourceE!(source, F, timestep)
    end

    for d in detectors 
        safeE!(d, F, timestep)
    end
end

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3)," seconds")
println("Computation Complete")

waterfall_plot=true
spectrum_plot=false

if waterfall_plot
    f = Figure(resolution = (800, 1200))
    ax1 = Axis(f[1, 1],title = "Waterfall Plot", ylabel = "timestep%10", xlabel = L"E_z")
    for medium in media
        vspan!(ax1, first(medium.location)[1], last(medium.location)[1], color=:gray90)
    end

    for (idx, t) in enumerate(11500:20:g.MaxTime)
        lines!(ax1, first(d1.location)[1]:last(d1.location)[1], d1.Ez[t, :] .+ idx, color=:black, linewidth=1.5)
    end
    f
end

if spectrum_plot 
    Δf = 1/g.Δt
    spectrum = fftshift(fft(d2.Ez))
    freqs = fftshift(fftfreq(MaxTime, Δf))
    f1 = Figure(resolution = (800, 400))
    ax = Axis(f1[1, 1])
    lines!(2*π*freqs, abs.(spectrum))
    #xlims!(ax, -5e-16, 5e-16)
    f1
end