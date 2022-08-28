using MaxwellFDTD
using FFTW
using GLMakie
GLMakie.activate!()
using ProgressBars
using CPUTime

CPUtic()
start = time()

const SizeX = 100
const SizeY = 50
const SizeZ = 50
const courant = 1.0/(sqrt(3.0))

const amplitude = 1.
const ppw = 20.

const Δx = Δy = Δz = 2e-9
const MaxTime = 1000

function ricker(t,location)
    src = 0
    delay = 60
    ppw = 30
    if t > delay
        arg = π * ((1/sqrt(3.0) * (t - (delay+1)) - location) / ppw - 1)
        src = 30*(1.0 - 2.0 * arg^2) * exp(-(arg^2))
    end
    return src
end


PML_Thickness = [15, 15, 15]

g = Grid3D(SizeX, SizeY, SizeZ, courant, Δx, Δy, Δz, MaxTime)
F = Fields3D(g)
F_PML = CPML_Ψ_Fields_3D(g, PML_Thickness)
c_PML = CPML_Parameters_3D(g, PML_Thickness)

const m_location = CartesianIndices((20:30, 20:30, 20:30))
m = StaticMedium3D(g, m_location, 1.)
media = [m]

c_grid = GridCoefficients3D_w_CPML(g, media, c_PML)

s_location = CartesianIndices((17:17, 2:49, 2:49))
s1 = GaussianPlaneWaveSource3D_y(g, s_location, true, true, 30., 620, 3e-16, 40.)
sources = [s1]

plane_pos = CartesianIndices((1:SizeX, 1:SizeY, 1:SizeZ))
d1 = BlockDetector(plane_pos, 1, MaxTime)
detectors = [d1]


for timestep in ProgressBar(1:g.MaxTime)
    
    update_Ψ_H!(F_PML, F, g, c_PML)

    updateH!(F, g, c_grid)

    for source in sources
        sourceH!(source, F, timestep)
    end

    apply_Ψ_H!(F_PML, F, g, c_PML)
    
    update_Ψ_E!(F_PML, F, g, c_PML)

    updateE!(F, g, c_grid)

    for source in sources
        sourceE!(source, F, timestep)
    end

    #F.Ez[50, 50] += ricker(timestep, 0)

    apply_Ψ_E!(F_PML, F, g, c_PML)
    

    for det in detectors
        safeE!(det, F, timestep)
    end
end

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3)," seconds")
println("Computation Complete")