using MaxwellFDTD
using FFTW
using GLMakie
GLMakie.activate!()
using ProgressBars
using CPUTime

CPUtic()
start = time()

const SizeX = 800


const amplitude = 1.
const ppw = 50.

const Δx = 2e-9
const MaxTime = 2640

PML_Thickness = [200]
courant = 0.99

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


g = Grid1D(SizeX, courant, Δx, MaxTime)
F = Fields1D(g)

F_PML = CPML_Ψ_Fields_1D(g, PML_Thickness)
c_PML = CPML_Parameters_1D(g, PML_Thickness)

const m_location = CartesianIndices((20:30,))
m = StaticMedium1D(g, m_location, 1.)
media = [m]

c_grid = GridCoefficients1D_w_CPML(g, media, c_PML)

block_pos = CartesianIndices((1:SizeX,))
d3 = LineDetector(block_pos, 1, MaxTime)
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

    F.Ez[400] += ricker(timestep, 0)

    apply_Ψ_E!(F_PML, F, g, c_PML)
    

    for det in detectors
        safeE!(det, F, timestep)
    end
end

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3)," seconds")
println("Computation Complete")

