# these are the update Functions for the one-dimesional Case
function updateE!(F::Fields1D, g::Grid1D, c::GridCoefficients1D)
    @inbounds for mm = 2:g.SizeX-1
        F.Ez[mm] = c.Ceze[mm] * F.Ez[mm] + c.Cezh[mm] * (F.Hy[mm] - F.Hy[mm-1])
    end
end

function updateH!(F::Fields1D, g::Grid1D, c::GridCoefficients1D)
    @inbounds for mm in 1:g.SizeX-1
        F.Hy[mm] = c.Chyh[mm] * F.Hy[mm] + c.Chye[mm] * (F.Ez[mm+1] - F.Ez[mm])
    end
end

function updateE!(F::Fields1D, M::LorentzFields1D, g::Grid1D, c::GridCoefficients1D)
    @inbounds for mm = 2:g.SizeX-1
        F.Ez[mm] = c.Ceze[mm] * F.Ez[mm] + c.Cezh[mm] * (F.Hy[mm] - F.Hy[mm-1] - g.Δx/2 * (M.Jz[mm]+M.Jz[mm-1]))
    end
end

function updateP!(MF::MaterialFields1D, LF::LorentzFields1D, M::LorentzMedium, g::Grid1D)
    @inbounds for osci in M.oscillators
        @. LF.Pz[M.location] = LF.Pz[M.location] + g.Δt/2 * (LF.Jz[M.location-CartesianIndex((1,))] + LF.Jz[M.location])
    end
end