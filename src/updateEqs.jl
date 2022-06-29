using LinearAlgebra
#=
 These are the update equations in  the one-dimesional Case
=#

function updateE!(F::Fields1D, g::Grid1D, c::GridCoefficients1D)
    @inbounds for mm = 2:g.SizeX-1
        F.Ez[mm] = c.Ceze[mm] * F.Ez[mm] + c.Cezh[mm] * (F.Hy[mm] - F.Hy[mm-1])
    end
end

function updateE!(F::Fields1D, MF::MaterialFields1D, g::Grid1D, c::GridCoefficients1D)
    @inbounds for mm = 2:g.SizeX-1
        F.Ez[mm] = c.Ceze[mm] * F.Ez[mm] + c.Cezh[mm] * (F.Hy[mm] - F.Hy[mm-1] - g.Δx * MF.Jz[mm])
    end
end

function updateH!(F::Fields1D, g::Grid1D, c::GridCoefficients1D)
    @inbounds for mm in 1:g.SizeX-1
        F.Hy[mm] = c.Chyh[mm] * F.Hy[mm] + c.Chye[mm] * (F.Ez[mm+1] - F.Ez[mm])
    end
end

function updateJ!(MF::MaterialFields1D)
    @inbounds MF.Jz .= MF.Jz_bound .+ MF.Jz_free .+ MF.Jz_tunnel
end

function updateJbound!(MF::MaterialFields1D, LF::LorentzFields1D, M::LorentzMedium1D, g::Grid1D)
    @inbounds for osci in 1:M.oscillators
        for mm in 1:length(M.location)
            LF.Jz[mm, osci] = (1.0-M.Γ[osci])/(1.0+M.Γ[osci])*LF.Jz[mm, osci] + (g.Δt*(M.ω_0[osci])^2)/(1.0+M.Γ[osci]) * (LF.PzNl[mm, osci] - LF.Pz[mm, osci]) 
        end
    end
    @inbounds MF.Jz_bound[M.location] .= sum(LF.Jz, dims=2)[:]
end

function updatePbound!(MF::MaterialFields1D, LF::LorentzFields1D, M::LorentzMedium1D, g::Grid1D)
    @inbounds for osci in 1:M.oscillators
        for mm in 1:length(M.location)
            LF.Pz[mm, osci] = LF.Pz[mm, osci] + g.Δt * (LF.Jz[mm, osci])
        end
    end
    @inbounds MF.Pz[M.location] .= sum(LF.Pz, dims=2)[:]
end

function updatePNl!(MF::MaterialFields1D, LF::LorentzFields1D, F::Fields1D, M::LorentzMedium1D)
    @inbounds for osci in 1:M.oscillators
        @. LF.PzNl[:, osci] = ϵ_0 * (M.χ_1[osci]*F.Ez[M.location] + M.χ_2[osci] * F.Ez[M.location]^2 + M.χ_3[osci] * F.Ez[M.location]^3)
    end
    @inbounds MF.PzNl[M.location] .= sum(LF.PzNl, dims=2)[:]
end

function updateJfree!(MF::MaterialFields1D, DF::DrudeFields1D, F::Fields1D, M::DrudeMedium1D)
    @inbounds begin
        @. DF.Jz_free[:] = (1 - M.Γ)/(1 + M.Γ) * DF.Jz_free[:] + ϵ_0 * M.grid.Δt * ω_plasma(MF.ρ_cb[M.location] * M.ρ_mol_density)^2 * F.Ez[M.location]/(1 + M.Γ)
        MF.Jz_free[M.location] .= DF.Jz_free
    end
end

function updateJtunnel!(MF::MaterialFields1D, TF::TunnelFields1D, M::TunnelMedium1D)
    @inbounds begin 
        @. TF.Jz_tunnel[:] = q_0 * TF.dz_T[:]*(1 - MF.ρ_cb[M.location])*MF.Γ_ADK[M.location]
        MF.Jz_tunnel[M.location] .= TF.Jz_tunnel
    end
end

function updatePlasma!(MF::MaterialFields1D, TF::TunnelFields1D, FC::FieldIonizationCoefficients1D, F::Fields1D, M::TunnelMedium1D)
    update_Γ_ADK!(MF, F, FC, M)
    #update_cb_population!(MF, M)
    update_ρ!(MF, M)
    update_displacement!(TF, F, M)
end

function update_Γ_ADK!(MF::MaterialFields1D, F::Fields1D, FC::FieldIonizationCoefficients1D, M::TunnelMedium1D)
    @inbounds MF.Γ_ADK[M.location] = Γ_ADK(F.Ez[M.location], FC.gamma_au[M.location], M.E_gap * q_0, 1, 0, 0)
end

function update_displacement!(TF::TunnelFields1D, F::Fields1D, M::TunnelMedium1D)
    @inbounds begin 
        @. TF.dz_T[:] = M.ρ_mol_density * M.E_gap*q_0*F.Ez[M.location] / (q_0 * abs(F.Ez[M.location] + 1e-6)^2) 
    end
end

function update_cb_population!(MF::MaterialFields1D, M::TunnelMedium1D)
    # really simple Euler integration 
    @inbounds begin 
        @. MF.ρ_cb[M.location] +=  MF.Γ_ADK[M.location] *  M.grid.Δt * (1 - MF.ρ_cb[M.location])
    end
end


#=
 WIP equations that will be used in the Future ()
=#

function update_ρ!(MF::MaterialFields1D, M::TunnelMedium1D)
    @inbounds begin 
        @. MF.ρ_cb[M.location]  = (M.grid.Δt * MF.Γ_ADK[M.location]) / (1 + M.grid.Δt/2 * MF.Γ_ADK[M.location]) + (1 - M.grid.Δt/2 * MF.Γ_ADK[M.location])/(1 + M.grid.Δt/2 *MF.Γ_ADK[M.location]) * MF.ρ_cb[M.location]
    end
end

#=
 These are the update equations in  the two-dimesional Case
=#




#=
 These are the update equations in  the three-dimesional Case
=#


function updateH!(F::Fields3D, g::Grid3D, c::GridCoefficients3D)  
    @inbounds for pp = 1:g.SizeZ-1
        for nn = 1:g.SizeY-1
            for mm = 1:g.SizeX
                F.Hx[mm,nn,pp] = (c.Chxh[mm,nn,pp] * F.Hx[mm,nn,pp] + 
                                c.Chxe[mm,nn,pp] * ((F.Ey[mm,nn,pp+1] - F.Ey[mm,nn,pp]) - (F.Ez[mm,nn+1,pp] - F.Ez[mm,nn,pp])))
            end;end;end


    @inbounds for pp = 1:g.SizeZ-1
        for nn = 1:g.SizeY
            for mm = 1:g.SizeX-1
                F.Hy[mm,nn,pp] = (c.Chyh[mm,nn,pp] * F.Hy[mm,nn,pp] + 
                                c.Chye[mm,nn,pp] * ((F.Ez[mm+1,nn,pp] - F.Ez[mm,nn,pp]) - (F.Ex[mm,nn,pp+1] - F.Ex[mm,nn,pp])))
            end;end;end
    
    @inbounds for pp = 1:g.SizeZ
        for nn = 1:g.SizeY-1
            for mm = 1:g.SizeX-1
                F.Hz[mm,nn,pp] = (c.Chzh[mm,nn,pp] * F.Hz[mm,nn,pp] + 
                                c.Chze[mm,nn,pp] * ((F.Ex[mm,nn+1,pp] - F.Ex[mm,nn,pp]) - (F.Ey[mm+1,nn,pp] - F.Ey[mm,nn,pp])))
            end;end;end
end

function updateE!(F::Fields3D, g::Grid3D, c::GridCoefficients3D) 

    @inbounds for pp = 2:g.SizeZ-1
        for nn = 2:g.SizeY-1
            for mm = 1:g.SizeX-1
                F.Ex[mm,nn,pp] = (c.Cexe[mm,nn,pp] * F.Ex[mm,nn,pp] + 
                                c.Cexh[mm,nn,pp] * ((F.Hz[mm,nn,pp] - F.Hz[mm,nn-1,pp]) - (F.Hy[mm,nn,pp] - F.Hy[mm,nn,pp-1])))
            end;end;end

    @inbounds for pp = 2:g.SizeZ-1
        for nn = 1:g.SizeY-1
            for mm = 2:g.SizeX-1
                F.Ey[mm,nn,pp] = (c.Ceye[mm,nn,pp] * F.Ey[mm,nn,pp] + 
                                c.Ceyh[mm,nn,pp] * ((F.Hx[mm,nn,pp] - F.Hx[mm,nn,pp-1]) - (F.Hz[mm,nn,pp] - F.Hz[mm-1,nn,pp])))
            end;end;end
    
    @inbounds for pp = 1:g.SizeZ-1
        for nn = 2:g.SizeY-1
            for mm = 2:g.SizeX-1
                F.Ez[mm,nn,pp] = (c.Ceze[mm,nn,pp] * F.Ez[mm,nn,pp] + 
                                c.Cezh[mm,nn,pp] * ((F.Hy[mm,nn,pp] - F.Hy[mm-1,nn,pp]) - (F.Hx[mm,nn,pp] - F.Hx[mm,nn-1,pp])))
            end;end;end
end

function updateEWIP!(F::Fields3D, g::Grid3D, c::GridCoefficients3D_WIP) 

    @inbounds for pp = 2:g.SizeZ-1
        for nn = 2:g.SizeY-1
            for mm = 1:g.SizeX-1
                F.Ex[mm,nn,pp] = (c.Cexe[mm,nn,pp] * F.Ex[mm,nn,pp] + 
                                c.Cexh[mm,nn,pp] * ((F.Hz[mm,nn,pp] - F.Hz[mm,nn-1,pp])/c.κ_E[mm,nn,pp,2] - (F.Hy[mm,nn,pp] - F.Hy[mm,nn,pp-1])/c.κ_E[mm,nn,pp,3]))
            end;end;end

    @inbounds for pp = 2:g.SizeZ-1
        for nn = 1:g.SizeY-1
            for mm = 2:g.SizeX-1
                F.Ey[mm,nn,pp] = (c.Ceye[mm,nn,pp] * F.Ey[mm,nn,pp] + 
                                c.Ceyh[mm,nn,pp] * ((F.Hx[mm,nn,pp] - F.Hx[mm,nn,pp-1])/c.κ_E[mm,nn,pp,3] - (F.Hz[mm,nn,pp] - F.Hz[mm-1,nn,pp])/c.κ_E[mm,nn,pp,1]))
            end;end;end
    
    @inbounds for pp = 1:g.SizeZ-1
        for nn = 2:g.SizeY-1
            for mm = 2:g.SizeX-1
                F.Ez[mm,nn,pp] = (c.Ceze[mm,nn,pp] * F.Ez[mm,nn,pp] + 
                                c.Cezh[mm,nn,pp] * ((F.Hy[mm,nn,pp] - F.Hy[mm-1,nn,pp])/c.κ_E[mm,nn,pp,1] - (F.Hx[mm,nn,pp] - F.Hx[mm,nn-1,pp])/c.κ_E[mm,nn,pp,2]))
            end;end;end
end

function updateHWIP!(F::Fields3D, g::Grid3D, c::GridCoefficients3D_WIP)  
    @inbounds for pp = 1:g.SizeZ-1
        for nn = 1:g.SizeY-1
            for mm = 1:g.SizeX
                F.Hx[mm,nn,pp] = (c.Chxh[mm,nn,pp] * F.Hx[mm,nn,pp] + 
                                c.Chxe[mm,nn,pp] * ((F.Ey[mm,nn,pp+1] - F.Ey[mm,nn,pp])/c.κ_H[mm,nn,pp,3] - (F.Ez[mm,nn+1,pp] - F.Ez[mm,nn,pp])/c.κ_H[mm,nn,pp,2]))
            end;end;end


    @inbounds for pp = 1:g.SizeZ-1
        for nn = 1:g.SizeY
            for mm = 1:g.SizeX-1
                F.Hy[mm,nn,pp] = (c.Chyh[mm,nn,pp] * F.Hy[mm,nn,pp] + 
                                c.Chye[mm,nn,pp] * ((F.Ez[mm+1,nn,pp] - F.Ez[mm,nn,pp])/c.κ_H[mm,nn,pp,1] - (F.Ex[mm,nn,pp+1] - F.Ex[mm,nn,pp])/c.κ_H[mm,nn,pp,3]))
            end;end;end
    
    @inbounds for pp = 1:g.SizeZ
        for nn = 1:g.SizeY-1
            for mm = 1:g.SizeX-1
                F.Hz[mm,nn,pp] = (c.Chzh[mm,nn,pp] * F.Hz[mm,nn,pp] + 
                                c.Chze[mm,nn,pp] * ((F.Ex[mm,nn+1,pp] - F.Ex[mm,nn,pp])/c.κ_H[mm,nn,pp,2] - (F.Ey[mm+1,nn,pp] - F.Ey[mm,nn,pp])/c.κ_H[mm,nn,pp,1]))
            end;end;end
end

#assumes cubic grid
function updateE!(F::Fields3D, MF::MaterialFields3D, g::Grid3D, c::GridCoefficients3D) 

    @inbounds for pp = 2:g.SizeZ-1
        for nn = 2:g.SizeY-1
            for mm = 1:g.SizeX-1
                F.Ex[mm,nn,pp] = (c.Cexe[mm,nn,pp] * F.Ex[mm,nn,pp] + 
                                c.Cexh[mm,nn,pp] * ((F.Hz[mm,nn,pp] - F.Hz[mm,nn-1,pp]) - (F.Hy[mm,nn,pp] - F.Hy[mm,nn,pp-1])) - g.Δx * MF.Jx)
            end;end;end

    @inbounds for pp = 2:g.SizeZ-1
        for nn = 1:g.SizeY-1
            for mm = 2:g.SizeX-1
                F.Ey[mm,nn,pp] = (c.Ceye[mm,nn,pp] * F.Ey[mm,nn,pp] + 
                                c.Ceyh[mm,nn,pp] * ((F.Hx[mm,nn,pp] - F.Hx[mm,nn,pp-1]) - (F.Hz[mm,nn,pp] - F.Hz[mm-1,nn,pp])) - g.Δx * MF.Jy)
            end;end;end
    
    @inbounds for pp = 1:g.SizeZ-1
        for nn = 2:g.SizeY-1
            for mm = 2:g.SizeX-1
                F.Ez[mm,nn,pp] = (c.Ceze[mm,nn,pp] * F.Ez[mm,nn,pp] + 
                                c.Cezh[mm,nn,pp] * ((F.Hy[mm,nn,pp] - F.Hy[mm-1,nn,pp]) - (F.Hx[mm,nn,pp] - F.Hx[mm,nn-1,pp])) - g.Δx * MF.Jz)
            end;end;end
    
    return Ex,Ey,Ez
end

function updateJbound!(MF::MaterialFields3D, LF::LorentzFields3D, M::LorentzMedium3D, g::Grid3D)
    @inbounds for osci in 1:M.oscillators
        for pp in 1:size(M.location)[3]
            for nn in 1:size(M.location)[2]
                for mm in 1:size(M.location)[1]
                    LF.Jx[mm, nn, pp, osci] = (1.0-M.Γ[osci])/(1.0+M.Γ[osci])*LF.Jx[mm,nn,pp,osci] + (g.Δt*(M.ω_0[osci])^2)/(1.0+M.Γ[osci]) * (LF.PxNl[mm,nn,pp,osci] - LF.Px[mm,nn,pp,osci]) 
                    LF.Jy[mm, nn, pp, osci] = (1.0-M.Γ[osci])/(1.0+M.Γ[osci])*LF.Jy[mm,nn,pp,osci] + (g.Δt*(M.ω_0[osci])^2)/(1.0+M.Γ[osci]) * (LF.PyNl[mm,nn,pp,osci] - LF.Py[mm,nn,pp,osci]) 
                    LF.Jz[mm, nn, pp, osci] = (1.0-M.Γ[osci])/(1.0+M.Γ[osci])*LF.Jz[mm,nn,pp,osci] + (g.Δt*(M.ω_0[osci])^2)/(1.0+M.Γ[osci]) * (LF.PzNl[mm,nn,pp,osci] - LF.Pz[mm,nn,pp,osci]) 
        end;end;end
    end
    MF.Jx[M.location] .= sum(LF.Jx, dims=4)[:]
    MF.Jy[M.location] .= sum(LF.Jy, dims=4)[:]
    MF.Jz[M.location] .= sum(LF.Jz, dims=4)[:]
end

function updatePbound!(MF::MaterialFields3D, LF::LorentzFields3D, M::LorentzMedium3D, g::Grid1D)
    @inbounds for osci in 1:M.oscillators
        for pp in 1:size(M.location)[3]
            for nn in 1:size(M.location)[2]
                for mm in 1:size(M.location)[1]
                    LF.Px[mm, nn, pp, osci] = LF.Px[mm, nn, pp, osci] + g.Δt * (LF.Jx[mm,nn,pp, osci])
                    LF.Py[mm, nn, pp, osci] = LF.Py[mm, nn, pp, osci] + g.Δt * (LF.Jy[mm,nn,pp, osci])
                    LF.Pz[mm, nn, pp, osci] = LF.Pz[mm, nn, pp, osci] + g.Δt * (LF.Jz[mm,nn,pp, osci])
        end;end;end
    end
    MF.Pz[M.location] .= sum(LF.Pz, dims=4)[:]
end

function updatePNl!(MF::MaterialFields3D, LF::LorentzFields3D, F::Fields3D, M::LorentzMedium3D)
    @inbounds for osci in 1:M.oscillators
        @. LF.PxNl[:, :, :, osci] = ϵ_0 * (M.χ_1[osci]*F.Ex[M.location] + M.χ_2[osci] * F.Ex[M.location]^2 + M.χ_3[osci] * F.Ex[M.location]^3)
        @. LF.PyNl[:, :, :, osci] = ϵ_0 * (M.χ_1[osci]*F.Ey[M.location] + M.χ_2[osci] * F.Ey[M.location]^2 + M.χ_3[osci] * F.Ey[M.location]^3)
        @. LF.PzNl[:, :, :, osci] = ϵ_0 * (M.χ_1[osci]*F.Ez[M.location] + M.χ_2[osci] * F.Ez[M.location]^2 + M.χ_3[osci] * F.Ez[M.location]^3)
    end
    MF.PxNl[M.location] .= sum(LF.PxNl, dims=4)[:]
    MF.PyNl[M.location] .= sum(LF.PyNl, dims=4)[:]
    MF.PzNl[M.location] .= sum(LF.PzNl, dims=4)[:]
end
