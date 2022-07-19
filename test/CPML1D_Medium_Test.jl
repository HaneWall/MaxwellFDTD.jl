using MaxwellFDTD
using FFTW
using GLMakie
GLMakie.activate!()
using ProgressBars
using CPUTime

CPUtic()
start = time()

SizeX = 2000
courant = 0.995
Δx = 4.e-9
MaxTime = 30000

PML_Thickness = [400]


g = Grid1D(SizeX, courant, Δx, MaxTime)
F = Fields1D(g)
MF = MaterialFields1D(g)

F_PML = CPML_Ψ_Fields_1D(g, PML_Thickness)
c_PML = CPML_Parameters_1D(g, PML_Thickness)

# Lorentz
γ = [0., 0., 0.]
ω_0 = [1.5494e16, 7.9514e15, 9.7766e13]
χ_1 = [2.4272, 1.4617, 9.6536]
χ_2 = [30e-12, 0., 0.]
χ_3 = [0., 0., 0.]
m1 = LorentzMedium1D(g, CartesianIndices((1000:2000,)), 1., γ, ω_0, χ_1, χ_2, χ_3)
media = [m1]

LF1 = LorentzFields1D(m1)
LF = [LF1]

s0 = GaussianWavePointSource(g, CartesianIndex((400,)),true, true, false, 1., 8500, 8e-15, 266.)
s1 = GaussianPointSource(CartesianIndex((400,)),true, true, false, 1., 90, 500.)
sources = [s0]

d1 = LineDetector(CartesianIndices((1:g.SizeX,)), 1, g.MaxTime)
detectors = [d1]

c_grid = GridCoefficients1D_w_CPML(g, media, c_PML)

for timestep in ProgressBar(1:g.MaxTime)

    for (m_idx, m) in enumerate(media)
        updatePNl!(MF, LF[m_idx], F, m)
        updateJbound!(MF, LF[m_idx], m, g)
        updatePbound!(MF, LF[m_idx], m, g)
    end

    update_Ψ_H!(F_PML, F, g, c_PML)

    updateH!(F, g, c_grid)

    for source in sources
        sourceH!(source, F, timestep-1)
    end

    apply_Ψ_H!(F_PML, F, g, c_PML)

    update_Ψ_E!(F_PML, F, g, c_PML)
    
    updateE!(F, MF, g, c_grid)

    for source in sources
        sourceE!(source, F, timestep)
    end

    apply_Ψ_E!(F_PML, F, g, c_PML)

    for d in detectors 
        safeE!(d, F, timestep)
    end
end

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3)," seconds")
println("Computation Complete")