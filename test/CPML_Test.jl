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
const ppw = 20.

const Δx = Δy = Δz = 2e-9
const MaxTime = 280

function ricker(t,location)
    src = 0
    delay = 15
    ppw = 20
    if t > delay
        arg = π * ((1/sqrt(3.0) * (t - (delay+1)) - location) / ppw - 1)
        src = 30*(1.0 - 2.0 * arg^2) * exp(-(arg^2))
    end
    return src
end

PML_Thickness = [10, 10, 10]

g = Grid3D(SizeX, SizeY, SizeZ, courant, Δx, Δy, Δz, MaxTime)
F = Fields3D(g)
F_PML = CPML_Ψ_Fields_3D(g, PML_Thickness)
c_PML = CPML_Parameters_3D(g, PML_Thickness)

const m_location = CartesianIndices((20:30, 20:30, 25:30))
m = StaticMedium3D(g, m_location, 1.)
media = [m]

c_grid = GridCoefficients3D_WIP2(g, media, c_PML)

block_pos = CartesianIndices((1:SizeX,1:SizeY,1:SizeZ))
d3 = BlockDetector(block_pos, 1, MaxTime)
detectors = [d3]

for timestep in ProgressBar(1:g.MaxTime)
    
    update_Ψ_H!(F_PML, F, g, c_PML)

    updateH!(F, g, c_grid)

    # for source in sources
    #     sourceH!(source, F, timestep)
    # end

    apply_Ψ_H!(F_PML, F, g, c_PML)
    
    update_Ψ_E!(F_PML, F, g, c_PML)

    updateE!(F, g, c_grid)

    # for source in sources
    #     sourceE!(source, F, timestep)
    # end

    F.Ez[24,24,24] += ricker(timestep, 0)

    apply_Ψ_E!(F_PML, F, g, c_PML)
    

    for det in detectors
        safeE!(det, F, timestep)
    end
end

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3)," seconds")
println("Computation Complete")
