using MaxwellFDTD
using CairoMakie
CairoMakie.activate!(type = "svg")
using CPUTime

#testing code
CPUtic()
start = time()

# 1. define grid
SizeX = 300
courant = 1.
Δx = 5e-07
MaxTime = 800
g = Grid1D(SizeX, courant, Δx, MaxTime)

# 2. define fields that exist everywhere
F = Fields1D(g.SizeX)

# 3. define and place media into the grid
m1 = StaticMedium1D(g, CartesianIndices((60:80,)), 5.)
m2 = StaticMedium1D(g, CartesianIndices((100:140,)), 9.)
media = [m1, m2]

# 4. define grid coefficients that respect ϵ_inf from the media 
c_grid = GridCoefficients1D(g, media)

# 5. define media fields that only have to be computed at the places where the media exist 


# 6. place detectors 
#d1 = PointDetector(CartesianIndex((3,)), 1, g.MaxTime)
d2 = LineDetector(CartesianIndices((1:g.SizeX,)), 1, g.MaxTime)
detectors = [d2]

# 7. place sources 
s1 = GaussianPointSource(CartesianIndex((20,)),true, true, false, 2., 70, 50.)
sources = [s1]

for timestep in 1:g.MaxTime
    updateH!(F, g, c_grid)

    for source in sources
        sourceH!(source, F, timestep)
    end

    ABC!(F, g)

    updateE!(F, g, c_grid)

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

if waterfall_plot
    f = Figure(resolution = (800, 1200))
    ax1 = Axis(f[1, 1],title = "Waterfall Plot", ylabel = "timestep%10", xlabel = L"E_z")
    for medium in media
        vspan!(ax1, first(medium.location)[1], last(medium.location)[1], color=:gray90)
    end

    for (idx, t) in enumerate(1:10:g.MaxTime)
        lines!(ax1, first(d2.location)[1]:last(d2.location)[1], d2.Ez[t, :] .+ idx, color=:black, linewidth=1.5)
    end
    f
end