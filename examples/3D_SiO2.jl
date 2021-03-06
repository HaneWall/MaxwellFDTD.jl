
using MaxwellFDTD
using CPUTime
using FFTW
#using GLMakie
using ProgressBars
using LaTeXStrings
using DSP
using GLMakie
#GLMakie.activate!()
GLMakie.activate!()

CPUtic()
start = time()

# 1. define grid
SizeX = 50
SizeY = 50
SizeZ = 50
courant = 0.95 * 1/sqrt(3)
Δx = Δy = Δz = 2e-9
MaxTime = 100

# 1. define grid
g = Grid3D(SizeX, SizeY, SizeZ, courant, Δx, Δy, Δz, MaxTime)
t = g.Δt:g.Δt:g.Δt*MaxTime

# varin paramters Si02
ρ_mol_density = 2.2e28
# bound electrons
γ_lorentz = [0.]
ω_0 = [2.75e16] # this might not work, use 2*π*/g.Δt instead (old varin paper/bachelor thesis) 2.75e16
χ_1 = [1.1025]
χ_2 = [0.]
χ_3 = [2.2e-22]

# drude parameters
γ_plasma = 1e15

# tunnel parameters
E_gap = 7.5

# source parameters
λ = 2100e-9
ω_central = 2 * π * c_0 / λ
ppw = λ/Δx
t_fwhm = 20e-15 # intensity FWHM
amplitude_pump = intensity2amplitude(12e16) # 12TWcm^-2

λ_probe = 800e-9
ω_probe= 2 * π * c_0 / λ_probe
ppw_probe = λ_probe/Δx
t_fwhm_probe = 45e-15 # intensity FWHM
amplitude_probe = intensity2amplitude(1.5e14) 

F = Fields3D(g)
MF = MaterialFields3D(g)

PML_Thickness = [10, 10, 10]
F_PML = CPML_Ψ_Fields_3D(g, PML_Thickness)
c_PML = CPML_Parameters_3D(g, PML_Thickness)


# init the media (superposition of different effects, that act inside the medium)
m1 = LorentzMedium3D(g, CartesianIndices((30:50, 1:50, 1:50)), 1., γ_lorentz, ω_0, χ_1, χ_2, χ_3)
m2 = DrudeMedium3D(g, CartesianIndices((30:50, 1:50, 1:50)), γ_plasma, ρ_mol_density)
m3 = TunnelMedium3D(g, CartesianIndices((30:50, 1:50, 1:50)), E_gap, ρ_mol_density)

bound_media= [m1]
drude_media = [m2]
tunnel_media = [m3]

# 4. define grid coefficients that respect ϵ_inf from the media 
c_grid = GridCoefficients3D_w_CPML(g, bound_media, c_PML)
f_grid = FieldIonizationCoefficients3D(g)

# 5. define fields inside the media
LF1 = LorentzFields3D(m1)
LF = [LF1]
DF1 = DrudeFields3D(m2)
DF = [DF1]
TF1 = TunnelFields3D(m3)
TF = [TF1]

block_pos = CartesianIndices((1:SizeX, 1:SizeY, 1:SizeZ))
d1 = BlockDetector(block_pos, 6999, 7000)
detectors = []

s1 = GaussianWavePointSource3D(g, CartesianIndex((22, 25, 25)),false, true, false, amplitude_pump, ceil(50e-15/g.Δt), t_fwhm_probe, ppw_probe)
sources = [s1]

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3)," seconds")
println("Init Complete")
CPUtic()
start = time()

for timestep in ProgressBar(1:g.MaxTime)
    
    @inbounds for (m_idx, m) in enumerate(tunnel_media)
        updatePlasma!(MF, TF[m_idx], f_grid, F, m)
        updateJtunnel!(MF, TF[m_idx], m)
    end

    # for (m_idx, m) in enumerate(tunnel_media)
    #     updatePlasmaTangent!(MF, TF[m_idx], f_grid, F, m, a, Γ̂, Ê)
    #     updateJtunnel!(MF, TF[m_idx], m)
    # end

    @inbounds for (m_idx, m) in enumerate(drude_media)
        updateJfree!(MF, DF[m_idx], F, m)
    end

    @inbounds for (m_idx, m) in enumerate(bound_media)
        updatePNl!(MF, LF[m_idx], F, m)
        updateJbound!(MF, LF[m_idx], m, g)
        updatePbound!(MF, LF[m_idx], m, g)
    end

    updateJ!(MF)

    update_Ψ_H!(F_PML, F, g, c_PML)

    updateH!(F, g, c_grid)

   # for source in sources
   #     sourceH!(source, F, timestep)
   # end

    apply_Ψ_H!(F_PML, F, g, c_PML)
    
    update_Ψ_E!(F_PML, F, g, c_PML)

    updateE!(F, MF, g, c_grid)

    for source in sources
        sourceE!(source, F, timestep)
    end

    apply_Ψ_E!(F_PML, F, g, c_PML)

    for d in detectors 
        #safeΓ_ADK!(d, MF, timestep)
        #safeJ_tunnel!(d, MF, timestep)
        #safeJ_free!(d, MF, timestep)
        #safeJ_bound!(d, MF, timestep)
        safeE!(d, F, timestep)
        #safeP!(d, MF, timestep)
        #safeJ!(d, MF, timestep)
        #safePNl!(d, MF, timestep)
    end
end

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3)," seconds")
println("Computation Complete")

