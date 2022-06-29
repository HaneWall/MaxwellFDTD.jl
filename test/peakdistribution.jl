using MaxwellFDTD
using CPUTime
using CairoMakie
using ProgressBars
CairoMakie.activate!(type="svg")

CPUtic()
start = time()

Maxtime_steps = 150000
Δx = 2e-9
SizeX = 10
courant = 0.985


g = Grid1D(10, courant, Δx, Maxtime_steps)
F = Fields1D(g)
MF = MaterialFields1D(g)

t = Array(g.Δt:g.Δt:g.Δt*Maxtime_steps)

# material parameter
E_gap = 7.7 
ρ_mol_density = 2.2e28

# Drude parameter
γ_plasma = 1e15

# Lorentz parameter
γ_lorentz = [0.]
ω_0 = [2.75e16] # this might not work, use 2*π*/g.Δt instead (old varin paper/bachelor thesis) 2.75e16
χ_1 = [1.1025]
χ_2 = [0.]
χ_3 = [2.2e-22]

# Medium init 
m1 = LorentzMedium1D(g, CartesianIndices((6:6,)), 1., γ_lorentz, ω_0, χ_1, χ_2, χ_3)
m2 = DrudeMedium1D(g, CartesianIndices((6:6,)), γ_plasma, ρ_mol_density)
m3 = TunnelMedium1D(g, CartesianIndices((6:6,)), E_gap, ρ_mol_density)

bound_media = [m1]
drude_media = [m2]
tunnel_media = [m3]

c_grid =GridCoefficients1D(g, bound_media)
f_grid =FieldIonizationCoefficients1D(g)

LF1 = LorentzFields1D(m1)
LF = [LF1]
DF1 = DrudeFields1D(m2)
DF = [DF1]
TF1 = TunnelFields1D(m3)
TF = [TF1]

# Pump Pulse parameter 
λ_pump = 2100e-9
ω_central = 2 * π * c_0 / λ_pump
ppw_pump = λ_pump/g.Δx
I_pump = 11.5e16
A_pump = intensity2amplitude(I_pump)
τ_delay_pump = 500e-15
τ_FWHM_pump = 140e-15

# Probe Pulse parameter
λ_probe = 800e-9 
ppw_probe = λ_probe/g.Δx
I_probe = 0.005e16
A_probe = intensity2amplitude(I_probe)
τ_delay_probe = 500e-15
τ_FWHM_probe = 75e-15

# Source init
s1 = GaussianWavePointSource(g, CartesianIndex((4,)),true, true, false, A_pump, ceil(τ_delay_pump/g.Δt), τ_FWHM_pump, ppw_pump)
s2 = GaussianWavePointSource(g, CartesianIndex((4,)),true, true, false, A_probe, ceil(τ_delay_probe/g.Δt), τ_FWHM_probe, ppw_probe)
sources = [s1, s2]

# Detector init
d1 = PointDetector(CartesianIndex((6,)), 1, g.MaxTime)
detectors = [d1]

#Boundary init 
b1 = LeftSideMurABC(g, CartesianIndex((1,)))
b2 = RightSideMurABC(g, CartesianIndex((SizeX,)))
boundaries = [b1, b2]

for timestep in ProgressBar(1:g.MaxTime)
    
    for (m_idx, m) in enumerate(tunnel_media)
        updatePlasma!(MF, TF[m_idx], f_grid, F, m)
        updateJtunnel!(MF, TF[m_idx], m)
    end

    for (m_idx, m) in enumerate(drude_media)
        updateJfree!(MF, DF[m_idx], F, m)
    end

    for (m_idx, m) in enumerate(bound_media)
        updatePNl!(MF, LF[m_idx], F, m)
        updateJbound!(MF, LF[m_idx], m, g)
        updatePbound!(MF, LF[m_idx], m, g)
    end

    updateJ!(MF)

    updateH!(F, g, c_grid)

    for source in sources
        sourceH!(source, F, timestep)
    end

    for b in boundaries
        saveFields!(b, F)
    end

    updateE!(F, MF, g, c_grid)

    for source in sources
        sourceE!(source, F, timestep)
    end

    for b in boundaries
        stepABC!(F, b)
    end

    for d in detectors 
        safeJ_tunnel!(d, MF, timestep)
        safeJ_free!(d, MF, timestep)
        safeJ_bound!(d, MF, timestep)
        safeE!(d, F, timestep)
        safeJ!(d, MF, timestep)
    end
end

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3)," seconds")
println("Computation Complete")

