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

function updateE!(F::Fields1D, MF::MaterialFields1D, g::Grid1D, c::GridCoefficients1D)
    @inbounds for mm = 2:g.SizeX-1
        F.Ez[mm] = c.Ceze[mm] * F.Ez[mm] + c.Cezh[mm] * (F.Hy[mm] - F.Hy[mm-1] - g.Δx *(MF.Jz[mm]))
    end
end

function updateJ!(MF::MaterialFields1D, LF::LorentzFields1D, M::LorentzMedium1D, g::Grid1D)
    @inbounds for osci in 1:M.oscillators
        for mm in 1:length(M.location)
            LF.Jz[mm, osci] = (1-M.Γ[osci])/(1+M.Γ[osci])*LF.Jz[mm, osci] + (g.Δt*(M.ω_0[osci])^2)/(1+M.Γ[osci]) * (LF.PzNl[mm, osci] - LF.Pz[mm, osci]) 
        end
    end
    MF.Jz_old[M.location] .= MF.Jz[M.location]
    MF.Jz[M.location] .= sum(LF.Jz, dims=2)[:]
end

function updateP!(MF::MaterialFields1D, LF::LorentzFields1D, M::LorentzMedium1D, g::Grid1D)
    @inbounds for osci in 1:M.oscillators
        for mm in 1:length(M.location)
            LF.Pz[mm, osci] = LF.Pz[mm, osci] + g.Δt * (LF.Jz[mm, osci])
        end
    end
    MF.Pz[M.location] .= sum(LF.Pz, dims=2)[:]
end

function updatePNl!(LF::LorentzFields1D, F::Fields1D, M::LorentzMedium1D)
    @inbounds for osci in 1:M.oscillators
        @. LF.PzNl[:, osci] = ϵ_0 * (M.χ_1[osci]*F.Ez[M.location] + M.χ_2[osci] * F.Ez[M.location]^2 + M.χ_3[osci] * F.Ez[M.location]^3)
    end
end