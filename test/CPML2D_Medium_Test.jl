using MaxwellFDTD
using FFTW
using GLMakie
GLMakie.activate!()
using ProgressBars
using CPUTime

CPUtic()
start = time()

const SizeX = 100
const SizeY = 100
const courant = 0.985 * 1.0/(sqrt(2.0))

const amplitude = 1.
const ppw = 20.

const Δx = Δy = 2e-9
const MaxTime = 5000

function ricker(t,location)
    src = 0
    delay = 700
    ppw = 30
    if t > delay
        arg = π * ((1/sqrt(3.0) * (t - (delay+1)) - location) / ppw - 1)
        src = 1000*(1.0 - 2.0 * arg^2) * exp(-(arg^2))
    end
    return src
end

PML_Thickness = [30,30]

g = Grid2D(SizeX, SizeY, courant, Δx, Δy, MaxTime)
F = Fields2D(g)
MF = MaterialFields2D(g)
F_PML = CPML_Ψ_Fields_2D(g, PML_Thickness)
c_PML = CPML_Parameters_2D(g, PML_Thickness)

const m_location = CartesianIndices((30:70, 2:40))
γ = [0., 0., 0.]
ω_0 = [1.5494e16, 7.9514e15, 9.7766e13]
χ_1 = [2.4272, 1.4617, 9.6536]
χ_2 = [30e-12, 0., 0.]
χ_3 = [0., 0., 0.]
m1 = LorentzMedium2D(g, m_location, 1., γ, ω_0, χ_1, χ_2, χ_3)
media = [m1]

LF1 = LorentzFields2D(m1)
LF = [LF1]

c_grid = GridCoefficients2D_w_CPML(g, media, c_PML)

plane_pos = CartesianIndices((1:SizeX, 1:SizeY))
d1 = PlaneDetector(plane_pos, 1, MaxTime)
detectors = [d1]

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3)," seconds")
println("Init Complete")
CPUtic()
start = time()

for timestep in ProgressBar(1:g.MaxTime)
    
    for (m_idx, m) in enumerate(media)
        updatePNl!(MF, LF[m_idx], F, m)
        updateJbound!(MF, LF[m_idx], m, g)
        updatePbound!(MF, LF[m_idx], m, g)
    end

    update_Ψ_H!(F_PML, F, g, c_PML)

    updateH!(F, g, c_grid)

    # for source in sources
    #     sourceH!(source, F, timestep)
    # end

    apply_Ψ_H!(F_PML, F, g, c_PML)
    
    update_Ψ_E!(F_PML, F, g, c_PML)

    updateE!(F, MF, g, c_grid)

    # for source in sources
    #     sourceE!(source, F, timestep)
    # end

    F.Ez[50, 50] += ricker(timestep, 0)

    apply_Ψ_E!(F_PML, F, g, c_PML)
    

    for det in detectors
        safeE!(det, F, timestep)
    end
end

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3)," seconds")
println("Computation Complete")