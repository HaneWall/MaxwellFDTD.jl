using MaxwellFDTD
using FFTW
using GLMakie
GLMakie.activate!()
using ProgressBars
using CPUTime

CPUtic()
start = time()

SizeX = 48
SizeY = 47
SizeZ = 47
courant = 1.0/(sqrt(3.0))

amplitude = 1.
ppw = 15.

Δx = Δy = Δz = 5e-07
MaxTime = 150

g = Grid3D(SizeX, SizeY, SizeZ, courant, Δx, Δy, Δz, MaxTime)

F = Fields3D(g)

m_location = CartesianIndices((20:30, 20:30, 20:30))
m = StaticMedium3D(g, m_location, 1.1)
media = [m]

c_grid = GridCoefficients3D(g, media)

s_location = CartesianIndex((15,15,15))
s0 = RickerPointSource3D(g, s_location, true, false, amplitude, 100, ppw)
sources = [s0]


loc_x_low = CartesianIndices((1:10, 1:g.SizeY, 1:g.SizeZ))
loc_x_high = CartesianIndices((g.SizeX-9:g.SizeX, 1:g.SizeY, 1:g.SizeZ))

loc_y_low = CartesianIndices((1:g.SizeX, 1:10, 1:g.SizeZ))
loc_y_high = CartesianIndices((1:g.SizeX, g.SizeY-9:g.SizeY, 1:g.SizeZ))

loc_z_low = CartesianIndices((1:g.SizeX, 1:g.SizeY, 1:10))
loc_z_high = CartesianIndices((1:g.SizeX, 1:g.SizeY, g.SizeZ-9:g.SizeZ))

b_x_low = PMLXlow(g, loc_x_low)
b_x_high = PMLXhigh(g, loc_x_high)

b_y_low = PMLYlow(g, loc_y_low)
b_y_high = PMLYhigh(g, loc_y_high)

b_z_low = PMLZlow(g, loc_z_low)
b_z_high = PMLZhigh(g, loc_z_high)

boundaries = [b_x_low, b_x_high, b_y_low, b_y_high, b_z_low, b_z_high]

for timestep in ProgressBar(1:g.MaxTime)
    for bound in boundaries
        update_ϕ_H!(bound, F)
    end

    updateH!(F, g, c_grid)

    for source in sources
        sourceH!(source, F, timestep)
    end

    for bound in boundaries
        updateH!(F, g, bound)
    end

    for bound in boundaries
        update_ϕ_E!(bound, F)
    end

    updateE!(F, g, c_grid)

    for source in sources
        sourceE!(source, F, timestep)
    end

    for bound in boundaries
        updateE!(F, g, bound)
    end
end

CPUtoq()
println("elapsed real time: ", round(time() - start; digits=3)," seconds")
println("Computation Complete")