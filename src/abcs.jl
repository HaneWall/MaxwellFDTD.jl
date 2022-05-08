function ABC!(F::Fields1D, g::Grid1D)
    F.Ez[1] = F.Ez[2]
    F.Ez[g.SizeX] = F.Ez[g.SizeX-1]
end