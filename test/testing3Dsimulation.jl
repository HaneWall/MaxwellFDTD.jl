using MaxwellFDTD
using FFTW
using GLMakie
GLMakie.activate!()
using ProgressBars
using CPUTime

CPUtic()
start = time()

const SizeX = 48
const SizeY = 48
const SizeZ = 48
const courant = 1.0/(sqrt(3.0))

const amplitude = 1.
const ppw = 25.

const Δx = Δy = Δz = 5e-07
const MaxTime = 120

function ricker(t,location)
    src = 0
    delay = 5
    ppw = 25
    if t > delay
        arg = π * ((1/sqrt(3.0) * (t - (delay+1)) - location) / ppw - 1)
        src = (1.0 - 2.0 * arg^2) * exp(-(arg^2))
    end
    return src
end

g = Grid3D(SizeX, SizeY, SizeZ, courant, Δx, Δy, Δz, MaxTime)

F = Fields3D(g)

const m_location = CartesianIndices((20:30, 20:30, 25:30))
m = StaticMedium3D(g, m_location, 1.)
media = [m]

c_grid = GridCoefficients3D(g, media)

const s_location = CartesianIndex((15,15,15))
s0 = RickerPointSource3D(g, s_location, true, false, amplitude, 0, ppw)
sources = [s0]


const loc_x_low = CartesianIndices((1:10, 1:SizeY, 1:SizeZ))
const loc_x_high = CartesianIndices((SizeX-9:SizeX, 1:SizeY, 1:SizeZ))

const loc_y_low = CartesianIndices((1:SizeX, 1:10, 1:g.SizeZ))
const loc_y_high = CartesianIndices((1:SizeX, SizeY-9:SizeY, 1:SizeZ))

const loc_z_low = CartesianIndices((1:SizeX, 1:SizeY, 1:10))
const loc_z_high = CartesianIndices((1:SizeX, 1:SizeY, SizeZ-9:SizeZ))

b_x_low = PMLXlow(g, loc_x_low)
b_x_high = PMLXhigh(g, loc_x_high)

b_y_low = PMLYlow(g, loc_y_low)
b_y_high = PMLYhigh(g, loc_y_high)

b_z_low = PMLZlow(g, loc_z_low)
b_z_high = PMLZhigh(g, loc_z_high)

boundaries = [b_x_low, b_x_high, b_y_low, b_y_high, b_z_low, b_z_high]


detector_plane = CartesianIndices((1:SizeX, 1:SizeY, 24:24))
d = ZSliceDetector(detector_plane, 1, MaxTime)
detectors = [d]

for timestep in ProgressBar(1:g.MaxTime)
    for bound in boundaries
        update_ϕ_H!(bound, F)
    end

    updateH!(F, g, c_grid)

    # for source in sources
    #     sourceH!(source, F, timestep)
    # end

    for bound in boundaries
        updateH!(F, g, bound)
    end

    for bound in boundaries
        update_ϕ_E!(bound, F)
    end

    updateE!(F, g, c_grid)

    # for source in sources
    #     sourceE!(source, F, timestep)
    # end

    F.Ex[24,24,24] += ricker(timestep, 0)

    for bound in boundaries
        updateE!(F, g, bound)
    end

    for det in detectors
        safeE!(d, F, timestep)
    end
end

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3)," seconds")
println("Computation Complete")