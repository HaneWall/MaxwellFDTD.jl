using MaxwellFDTD
using GLMakie
GLMakie.activate!()
using CPUTime
using ProgressBars
using FFTW

#testing code
CPUtic()
start = time()

# 1. define grid
SizeX = 500
courant = 0.98
Δx = 10.e-9
MaxTime = 2000

# 1. define grid
g = Grid1D(SizeX, courant, Δx, MaxTime)

# 2. define fields that exist everywhere
F = Fields1D(g)
MF = MaterialFields1D(g)

function ricker(t,location)
    src = 0
    delay = 250
    ppw = 50
    if t > delay
        arg = π * ((1/sqrt(3.0) * (t - (delay+1)) - location) / ppw - 1)
        src = (1.0 - 2.0 * arg^2) * exp(-(arg^2))
    end
    return src
end

# 3. define and place media into the grid
# parameters from bachelor thesis
γ = [0.]
ω_0 = [1.5494e15]
χ_1 = [2.4272]
χ_2 = [30e-12]
χ_3 = [0.]
#m1 = LorentzMedium1D(g, CartesianIndices((10:990,)), 1., γ, ω_0, χ_1, χ_2, χ_3)
m1 = StaticMedium1D(g, CartesianIndices((4:6,)), 1.)
media = [m1]

# 4. define grid coefficients that respect ϵ_inf from the media 
c_grid = GridCoefficients1D(g, media)

# 6. place detectors 
d1 = LineDetector(CartesianIndices((1:g.SizeX,)), 1, g.MaxTime)
detectors = [d1]

# 7. place sources 
amplitude = intensity2amplitude(10.0^12)

s0 = GaussianWavePointSource(g, CartesianIndex((7,)),false, false, false, amplitude, 8500, 20e-15, 266.)
sources = []

# 8. place boundaries
b1 = PeriodicBoundaryX(g)
boundaries = [b1]

for timestep in ProgressBar(1:g.MaxTime)

    updateH!(F, g, c_grid)

    for source in sources
        sourceH!(source, F, timestep)
    end

    for b in boundaries
        updateH!(b, F)
    end
    #ABC!(F, g)

    updateE!(F, MF, g, c_grid)

    F.Ez[100] += ricker(timestep, 0)

    for b in boundaries
        updateE!(b, F)
    end

    for d in detectors 
        safeE!(d, F, timestep)
    end
end

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3)," seconds")
println("Computation Complete")