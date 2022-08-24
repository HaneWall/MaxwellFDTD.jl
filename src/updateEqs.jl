using LinearAlgebra
#=
 These are the update equations in  the one-dimesional Case
=#

function updateE!(F::Fields1D, g::Grid1D, c::GridCoefficients1D)
    @inbounds for mm = 2:g.SizeX-1
        F.Ez[mm] = c.Ceze[mm] * F.Ez[mm] + c.Cezh[mm] * (F.Hy[mm] - F.Hy[mm-1])
    end
end

function updateE!(F::Fields1D, g::Grid1D, c::GridCoefficients1D_w_CPML)
    @inbounds for mm = 2:g.SizeX-1
        F.Ez[mm] = c.Ceze[mm] * F.Ez[mm] + c.Cezh[mm] * (F.Hy[mm] - F.Hy[mm-1]) * c.Den_Ex[mm]
    end
end

function updateE!(F::Fields1D, MF::MaterialFields1D, g::Grid1D, c::GridCoefficients1D)
    @inbounds for mm = 2:g.SizeX-1
        F.Ez[mm] = c.Ceze[mm] * F.Ez[mm] + c.Cezh[mm] * (F.Hy[mm] - F.Hy[mm-1] - g.Δx * MF.Jz[mm])
    end
end

function updateE!(F::Fields1D, MF::MaterialFields1D, g::Grid1D, c::GridCoefficients1D_w_CPML)
    @inbounds for mm = 2:g.SizeX-1
        F.Ez[mm] = c.Ceze[mm] * F.Ez[mm] + c.Cezh[mm]  * (c.Den_Ex[mm] * (F.Hy[mm] - F.Hy[mm-1]) -  MF.Jz[mm])
    end
end


function updateH!(F::Fields1D, g::Grid1D, c::GridCoefficients1D)
    @inbounds for mm in 1:g.SizeX-1
        F.Hy[mm] = c.Chyh[mm] * F.Hy[mm] + c.Chye[mm] * (F.Ez[mm+1] - F.Ez[mm])
    end
end

function updateH!(F::Fields1D, g::Grid1D, c::GridCoefficients1D_w_CPML)
    @inbounds for mm in 1:g.SizeX-1
        F.Hy[mm] = c.Chyh[mm] * F.Hy[mm] + c.Chye[mm] * (F.Ez[mm+1] - F.Ez[mm]) * c.Den_Hx[mm]
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

function updatePlasmaTangent!(MF::MaterialFields1D, TF::TunnelFields1D, FC::FieldIonizationCoefficients1D, F::Fields1D, M::TunnelMedium1D, a::Float64, Γ̂::Float64, Ê::Float64)
    update_Γ_Tangent!(MF, F, FC, M, a, Γ̂, Ê)
    update_ρ!(MF, M)
    update_displacement!(TF, F, M)
end

function update_Γ_ADK!(MF::MaterialFields1D, F::Fields1D, FC::FieldIonizationCoefficients1D, M::TunnelMedium1D)
    @inbounds MF.Γ_ADK[M.location] = Γ_ADK(F.Ez[M.location], FC.gamma_au[M.location], M.E_gap * q_0, 1, 0, 0)
end

function update_Γ_Tangent!(MF::MaterialFields1D, F::Fields1D, FC::FieldIonizationCoefficients1D, M::TunnelMedium1D, a::Float64, Γ̂::Float64, Ê::Float64)
    @inbounds MF.Γ_ADK[M.location] = Γ_Tangent(F.Ez[M.location], FC.gamma_au[M.location], a, Γ̂, Ê)
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

function update_ρ!(MF::MaterialFields1D, M::TunnelMedium1D)
    @inbounds begin 
        @. MF.ρ_cb[M.location]  = (M.grid.Δt * MF.Γ_ADK[M.location]) / (1 + M.grid.Δt/2 * MF.Γ_ADK[M.location]) + (1 - M.grid.Δt/2 * MF.Γ_ADK[M.location])/(1 + M.grid.Δt/2 *MF.Γ_ADK[M.location]) * MF.ρ_cb[M.location]
    end
end

#=
 These are the update equations in  the two-dimesional Case
=#
function updateE!(F::Fields2D, g::Grid2D, c::GridCoefficients2D)
    @inbounds for nn = 2:g.SizeY-1
        for mm = 2:g.SizeX-1
            F.Ez[mm,nn] = (c.Ceze[mm,nn] * F.Ez[mm,nn] + 
                            c.Cezh[mm,nn] * ((F.Hy[mm,nn] - F.Hy[mm-1,nn]) - (F.Hx[mm,nn] - F.Hx[mm,nn-1])))
        end;end
end

function updateE!(F::Fields2D, g::Grid2D, c::GridCoefficients2D_w_CPML)
        @inbounds for nn = 2:g.SizeY-1
            for mm = 2:g.SizeX-1
                F.Ez[mm,nn] = (c.Ceze[mm,nn] * F.Ez[mm,nn] + 
                                c.Cezh[mm,nn] * ((F.Hy[mm,nn] - F.Hy[mm-1,nn])*c.Den_Ex[mm] - (F.Hx[mm,nn] - F.Hx[mm,nn-1])*c.Den_Ey[nn]))
            end;end
end

function updateE!(F::Fields2D, MF::MaterialFields2D, g::Grid2D, c::GridCoefficients2D_w_CPML)
    @inbounds for nn = 2:g.SizeY-1
        for mm = 2:g.SizeX-1
            F.Ez[mm,nn] = (c.Ceze[mm,nn] * F.Ez[mm,nn] + 
                            c.Cezh[mm,nn] * ((F.Hy[mm,nn] - F.Hy[mm-1,nn])*c.Den_Ex[mm] - (F.Hx[mm,nn] - F.Hx[mm,nn-1])*c.Den_Ey[nn] - MF.Jz[mm,nn]))
        end;end
end

function updateH!(F::Fields2D, g::Grid2D, c::GridCoefficients2D)  
    @inbounds for nn = 1:g.SizeY-1
        for mm = 1:g.SizeX
            F.Hx[mm,nn] = (c.Chxh[mm,nn] * F.Hx[mm,nn] + 
                        c.Chxe[mm,nn] * ( - (F.Ez[mm,nn+1] - F.Ez[mm,nn])))
        end;end


    @inbounds for nn = 1:g.SizeY
        for mm = 1:g.SizeX-1
            F.Hy[mm,nn] = (c.Chyh[mm,nn] * F.Hy[mm,nn] + 
                            c.Chye[mm,nn] * ((F.Ez[mm+1,nn] - F.Ez[mm,nn])))
        end;end
end

function updateH!(F::Fields2D, g::Grid2D, c::GridCoefficients2D_w_CPML)  
    @inbounds for nn = 1:g.SizeY-1
        for mm = 1:g.SizeX
            F.Hx[mm,nn] = (c.Chxh[mm,nn] * F.Hx[mm,nn] + 
                        c.Chxe[mm,nn] * ( - (F.Ez[mm,nn+1] - F.Ez[mm,nn])*c.Den_Hy[nn]))
        end;end


    @inbounds for nn = 1:g.SizeY
        for mm = 1:g.SizeX-1
            F.Hy[mm,nn] = (c.Chyh[mm,nn] * F.Hy[mm,nn] + 
                            c.Chye[mm,nn] * ((F.Ez[mm+1,nn] - F.Ez[mm,nn])*c.Den_Hx[mm]))
        end;end
end

function update_ρ!(MF::MaterialFields2D, M::TunnelMedium2D)
    @inbounds begin 
        @. MF.ρ_cb[M.location]  = (M.grid.Δt * MF.Γ_ADK[M.location]) / (1 + M.grid.Δt/2 * MF.Γ_ADK[M.location]) + (1 - M.grid.Δt/2 * MF.Γ_ADK[M.location])/(1 + M.grid.Δt/2 *MF.Γ_ADK[M.location]) * MF.ρ_cb[M.location]
    end
end

function update_Γ_ADK!(MF::MaterialFields2D, F::Fields2D, FC::FieldIonizationCoefficients2D, M::TunnelMedium2D)
    @inbounds MF.Γ_ADK[M.location] = Γ_ADK(F.Ez[M.location], FC.gamma_au[M.location], M.E_gap * q_0, 1, 0, 0)
end

function update_Γ_Tangent!(MF::MaterialFields2D, F::Fields2D, FC::FieldIonizationCoefficients2D, M::TunnelMedium2D, a::Float64, Γ̂::Float64, Ê::Float64)
    @inbounds MF.Γ_ADK[M.location] = Γ_Tangent(F.Ez[M.location], FC.gamma_au[M.location], a, Γ̂, Ê)
end

function update_displacement!(TF::TunnelFields2D, F::Fields2D, M::TunnelMedium2D)
    @inbounds begin 
        @. TF.dz_T[:, :] = M.ρ_mol_density * M.E_gap*q_0*F.Ez[M.location] / (q_0 * abs(F.Ez[M.location] + 1e-6)^2) 
    end
end

function updateJ!(MF::MaterialFields2D)
    @inbounds MF.Jz .= MF.Jz_bound .+ MF.Jz_free .+ MF.Jz_tunnel
end

function updateJbound!(MF::MaterialFields2D, LF::LorentzFields2D, M::LorentzMedium2D, g::Grid2D)
    @inbounds for osci in 1:M.oscillators
        for nn in 1:size(M.location)[2]
            for mm in 1:size(M.location)[1]
                LF.Jz[mm, nn, osci] = (1.0-M.Γ[osci])/(1.0+M.Γ[osci])*LF.Jz[mm, osci] + (g.Δt*(M.ω_0[osci])^2)/(1.0+M.Γ[osci]) * (LF.PzNl[mm, osci] - LF.Pz[mm, osci]) 
            end
        end
    end
    @inbounds MF.Jz_bound[M.location] .= sum(LF.Jz, dims=3)
end

function updatePbound!(MF::MaterialFields2D, LF::LorentzFields2D, M::LorentzMedium2D, g::Grid2D)
    @inbounds for osci in 1:M.oscillators
        for nn in 1:size(M.location)[2]
            for mm in 1:size(M.location)[1]
                LF.Pz[mm, nn, osci] = LF.Pz[mm, nn, osci] + g.Δt * (LF.Jz[mm, nn, osci])
            end
        end
    end
    @inbounds MF.Pz[M.location] .= sum(LF.Pz, dims=3)
end

function updatePNl!(MF::MaterialFields2D, LF::LorentzFields2D, F::Fields2D, M::LorentzMedium2D)
    @inbounds for osci in 1:M.oscillators
        @. LF.PzNl[:, :, osci] = ϵ_0 * (M.χ_1[osci]*F.Ez[M.location] + M.χ_2[osci] * F.Ez[M.location]^2 + M.χ_3[osci] * F.Ez[M.location]^3)
    end
    @inbounds MF.PzNl[M.location] .= sum(LF.PzNl, dims=3)
end

function updateJfree!(MF::MaterialFields2D, DF::DrudeFields2D, F::Fields2D, M::DrudeMedium2D)
    @inbounds begin
        @. DF.Jz_free[:, :] = (1 - M.Γ)/(1 + M.Γ) * DF.Jz_free[:, :] + ϵ_0 * M.grid.Δt * ω_plasma(MF.ρ_cb[M.location] * M.ρ_mol_density)^2 * F.Ez[M.location]/(1 + M.Γ)
        MF.Jz_free[M.location] .= DF.Jz_free
    end
end

function updateJtunnel!(MF::MaterialFields2D, TF::TunnelFields2D, M::TunnelMedium2D)
    @inbounds begin 
        @. TF.Jz_tunnel[:, :] = q_0 * TF.dz_T[:, :]*(1 - MF.ρ_cb[M.location])*MF.Γ_ADK[M.location]
        MF.Jz_tunnel[M.location] .= TF.Jz_tunnel
    end
end

function updatePlasma!(MF::MaterialFields2D, TF::TunnelFields2D, FC::FieldIonizationCoefficients2D, F::Fields2D, M::TunnelMedium2D)
    update_Γ_ADK!(MF, F, FC, M)
    #update_cb_population!(MF, M)
    update_ρ!(MF, M)
    update_displacement!(TF, F, M)
end

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

function updateE!(F::Fields3D, g::Grid3D, c::GridCoefficients3D_w_CPML)
    @inbounds for pp = 2:g.SizeZ-1
        for nn = 2:g.SizeY-1
            for mm = 1:g.SizeX-1
                F.Ex[mm,nn,pp] = (c.Cexe[mm,nn,pp] * F.Ex[mm,nn,pp] + 
                                c.Cexh[mm,nn,pp] * ((F.Hz[mm,nn,pp] - F.Hz[mm,nn-1,pp])*c.Den_Ey[nn] - (F.Hy[mm,nn,pp] - F.Hy[mm,nn,pp-1])*c.Den_Ez[pp]))
            end;end;end

    @inbounds for pp = 2:g.SizeZ-1
        for nn = 1:g.SizeY-1
            for mm = 2:g.SizeX-1
                F.Ey[mm,nn,pp] = (c.Ceye[mm,nn,pp] * F.Ey[mm,nn,pp] + 
                                c.Ceyh[mm,nn,pp] * ((F.Hx[mm,nn,pp] - F.Hx[mm,nn,pp-1])*c.Den_Ez[pp] - (F.Hz[mm,nn,pp] - F.Hz[mm-1,nn,pp])*c.Den_Ex[mm]))
            end;end;end
    
    @inbounds for pp = 1:g.SizeZ-1
        for nn = 2:g.SizeY-1
            for mm = 2:g.SizeX-1
                F.Ez[mm,nn,pp] = (c.Ceze[mm,nn,pp] * F.Ez[mm,nn,pp] + 
                                c.Cezh[mm,nn,pp] * ((F.Hy[mm,nn,pp] - F.Hy[mm-1,nn,pp])*c.Den_Ex[mm] - (F.Hx[mm,nn,pp] - F.Hx[mm,nn-1,pp])*c.Den_Ey[nn]))
            end;end;end
end

function updateE!(F::Fields3D, MF::MaterialFields3D, g::Grid3D, c::GridCoefficients3D_w_CPML)
    @inbounds for pp = 2:g.SizeZ-1
        for nn = 2:g.SizeY-1
            for mm = 1:g.SizeX-1
                F.Ex[mm,nn,pp] = (c.Cexe[mm,nn,pp] * F.Ex[mm,nn,pp] + 
                                c.Cexh[mm,nn,pp] * ((F.Hz[mm,nn,pp] - F.Hz[mm,nn-1,pp])*c.Den_Ey[nn] - (F.Hy[mm,nn,pp] - F.Hy[mm,nn,pp-1])*c.Den_Ez[pp] - MF.Jx[mm,nn,pp]))
            end;end;end

    @inbounds for pp = 2:g.SizeZ-1
        for nn = 1:g.SizeY-1
            for mm = 2:g.SizeX-1
                F.Ey[mm,nn,pp] = (c.Ceye[mm,nn,pp] * F.Ey[mm,nn,pp] + 
                                c.Ceyh[mm,nn,pp] * ((F.Hx[mm,nn,pp] - F.Hx[mm,nn,pp-1])*c.Den_Ez[pp] - (F.Hz[mm,nn,pp] - F.Hz[mm-1,nn,pp])*c.Den_Ex[mm] - MF.Jy[mm,nn,pp]))
            end;end;end
    
    @inbounds for pp = 1:g.SizeZ-1
        for nn = 2:g.SizeY-1
            for mm = 2:g.SizeX-1
                F.Ez[mm,nn,pp] = (c.Ceze[mm,nn,pp] * F.Ez[mm,nn,pp] + 
                                c.Cezh[mm,nn,pp] * ((F.Hy[mm,nn,pp] - F.Hy[mm-1,nn,pp])*c.Den_Ex[mm] - (F.Hx[mm,nn,pp] - F.Hx[mm,nn-1,pp])*c.Den_Ey[nn] - MF.Jz[mm,nn,pp]))
            end;end;end
end
#            F.Ez[mm,nn] = (c.Ceze[mm,nn] * F.Ez[mm,nn] + c.Cezh[mm,nn] * ((F.Hy[mm,nn] - F.Hy[mm-1,nn])*c.Den_Ex[mm] - (F.Hx[mm,nn] - F.Hx[mm,nn-1])*c.Den_Ey[nn] - MF.Jz[mm,nn]))

function updateH!(F::Fields3D, g::Grid3D, c::GridCoefficients3D_w_CPML)  
    @inbounds for pp = 1:g.SizeZ-1
        for nn = 1:g.SizeY-1
            for mm = 1:g.SizeX
                F.Hx[mm,nn,pp] = (c.Chxh[mm,nn,pp] * F.Hx[mm,nn,pp] + 
                                c.Chxe[mm,nn,pp] * ((F.Ey[mm,nn,pp+1] - F.Ey[mm,nn,pp])*c.Den_Hz[pp] - (F.Ez[mm,nn+1,pp] - F.Ez[mm,nn,pp])*c.Den_Hy[nn]))
            end;end;end


    @inbounds for pp = 1:g.SizeZ-1
        for nn = 1:g.SizeY
            for mm = 1:g.SizeX-1
                F.Hy[mm,nn,pp] = (c.Chyh[mm,nn,pp] * F.Hy[mm,nn,pp] + 
                                c.Chye[mm,nn,pp] * ((F.Ez[mm+1,nn,pp] - F.Ez[mm,nn,pp])*c.Den_Hx[mm] - (F.Ex[mm,nn,pp+1] - F.Ex[mm,nn,pp])*c.Den_Hz[pp]))
            end;end;end
    
    @inbounds for pp = 1:g.SizeZ
        for nn = 1:g.SizeY-1
            for mm = 1:g.SizeX-1
                F.Hz[mm,nn,pp] = (c.Chzh[mm,nn,pp] * F.Hz[mm,nn,pp] + 
                                c.Chze[mm,nn,pp] * ((F.Ex[mm,nn+1,pp] - F.Ex[mm,nn,pp])*c.Den_Hy[nn] - (F.Ey[mm+1,nn,pp] - F.Ey[mm,nn,pp])*c.Den_Hx[mm]))
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

function update_ρ!(MF::MaterialFields3D, M::TunnelMedium3D)
    @inbounds begin 
        @. MF.ρ_cb[M.location]  = (M.grid.Δt * MF.Γ_ADK[M.location]) / (1 + M.grid.Δt/2 * MF.Γ_ADK[M.location]) + (1 - M.grid.Δt/2 * MF.Γ_ADK[M.location])/(1 + M.grid.Δt/2 *MF.Γ_ADK[M.location]) * MF.ρ_cb[M.location]
    end
end

function update_Γ_ADK!(MF::MaterialFields3D, F::Fields3D, FC::FieldIonizationCoefficients3D, M::TunnelMedium3D)
    @inbounds MF.Γ_ADK[M.location] = Γ_ADK(F.Ez[M.location], FC.gamma_au[M.location], M.E_gap * q_0, 1, 0, 0)
end

function update_Γ_Tangent!(MF::MaterialFields3D, F::Fields3D, FC::FieldIonizationCoefficients3D, M::TunnelMedium3D, a::Float64, Γ̂::Float64, Ê::Float64)
    @inbounds MF.Γ_ADK[M.location] = Γ_Tangent(F.Ez[M.location], FC.gamma_au[M.location], a, Γ̂, Ê)
end

function update_displacement!(TF::TunnelFields3D, F::Fields3D, M::TunnelMedium3D)
    @inbounds begin 
        @. TF.dz_T[:, :, :] = M.ρ_mol_density * M.E_gap*q_0*F.Ez[M.location] / (q_0 * abs(F.Ez[M.location] + 1e-6)^2) 
    end
end

function updateJ!(MF::MaterialFields3D)
    @inbounds MF.Jz .= MF.Jz_bound .+ MF.Jz_free .+ MF.Jz_tunnel
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
    @inbounds MF.Jx[M.location] .= sum(LF.Jx, dims=4)
    @inbounds MF.Jy[M.location] .= sum(LF.Jy, dims=4)
    @inbounds MF.Jz[M.location] .= sum(LF.Jz, dims=4)
end

function updatePbound!(MF::MaterialFields3D, LF::LorentzFields3D, M::LorentzMedium3D, g::Grid3D)
    @inbounds for osci in 1:M.oscillators
        for pp in 1:size(M.location)[3]
            for nn in 1:size(M.location)[2]
                for mm in 1:size(M.location)[1]
                    LF.Px[mm, nn, pp, osci] = LF.Px[mm, nn, pp, osci] + g.Δt * (LF.Jx[mm,nn,pp, osci])
                    LF.Py[mm, nn, pp, osci] = LF.Py[mm, nn, pp, osci] + g.Δt * (LF.Jy[mm,nn,pp, osci])
                    LF.Pz[mm, nn, pp, osci] = LF.Pz[mm, nn, pp, osci] + g.Δt * (LF.Jz[mm,nn,pp, osci])
        end;end;end
    end
    @inbounds MF.Pz[M.location] .= sum(LF.Pz, dims=4)
end

function updatePNl!(MF::MaterialFields3D, LF::LorentzFields3D, F::Fields3D, M::LorentzMedium3D)
    @inbounds for osci in 1:M.oscillators
        @. LF.PxNl[:, :, :, osci] = ϵ_0 * (M.χ_1[osci]*F.Ex[M.location] + M.χ_2[osci] * (F.Ex[M.location]^2 + F.Ey[M.location]^2 + F.Ez[M.location]^2) + M.χ_3[osci] * (F.Ex[M.location]^2 + F.Ey[M.location]^2 + F.Ez[M.location]^2) * F.Ex[M.location])
        @. LF.PyNl[:, :, :, osci] = ϵ_0 * (M.χ_1[osci]*F.Ey[M.location] + M.χ_2[osci] * (F.Ex[M.location]^2 + F.Ey[M.location]^2 + F.Ez[M.location]^2) + M.χ_3[osci] * (F.Ex[M.location]^2 + F.Ey[M.location]^2 + F.Ez[M.location]^2) * F.Ey[M.location])
        @. LF.PzNl[:, :, :, osci] = ϵ_0 * (M.χ_1[osci]*F.Ez[M.location] + M.χ_2[osci] * (F.Ex[M.location]^2 + F.Ey[M.location]^2 + F.Ez[M.location]^2) + M.χ_3[osci] * (F.Ex[M.location]^2 + F.Ey[M.location]^2 + F.Ez[M.location]^2) * F.Ez[M.location])
    end
    @inbounds MF.PxNl[M.location] .= sum(LF.PxNl, dims=4)
    @inbounds MF.PyNl[M.location] .= sum(LF.PyNl, dims=4)
    @inbounds MF.PzNl[M.location] .= sum(LF.PzNl, dims=4)
end

function updateJfree!(MF::MaterialFields3D, DF::DrudeFields3D, F::Fields3D, M::DrudeMedium3D)
    @inbounds begin
        @. DF.Jz_free[:, :, :] = (1 - M.Γ)/(1 + M.Γ) * DF.Jz_free[:, :, :] + ϵ_0 * M.grid.Δt * ω_plasma(MF.ρ_cb[M.location] * M.ρ_mol_density)^2 * F.Ez[M.location]/(1 + M.Γ)
        MF.Jz_free[M.location] .= DF.Jz_free
    end
end

function updateJtunnel!(MF::MaterialFields3D, TF::TunnelFields3D, M::TunnelMedium3D)
    @inbounds begin 
        @. TF.Jz_tunnel[:, :, :] = q_0 * TF.dz_T[:, :, :]*(1 - MF.ρ_cb[M.location])*MF.Γ_ADK[M.location]
        MF.Jz_tunnel[M.location] .= TF.Jz_tunnel
    end
end

function updatePlasma!(MF::MaterialFields3D, TF::TunnelFields3D, FC::FieldIonizationCoefficients3D, F::Fields3D, M::TunnelMedium3D)
    update_Γ_ADK!(MF, F, FC, M)
    #update_cb_population!(MF, M)
    update_ρ!(MF, M)
    update_displacement!(TF, F, M)
end

function update_Ψ_E!(PML::CPML_Ψ_Fields_3D, F::Fields3D, g::Grid3D, c::CPML_Parameters_3D)
    
    #X-Boundary
    @inbounds for pp in 2:g.SizeZ-1
        for nn in 1:g.SizeY-1
            # X-Bot
            for mm in 2:c.PML_Thickness[1]
                PML.Ψ_Eyx[mm,nn,pp,1] = c.b_E_x[mm, 1] * PML.Ψ_Eyx[mm,nn,pp,1] + c.c_E_x[mm, 1] * 1/g.Δx *(F.Hz[mm-1,nn,pp] - F.Hz[mm,nn,pp])
            end
            # X-Top
            m_rev = c.PML_Thickness[1]
            for mm in (g.SizeX+1-c.PML_Thickness[1]):g.SizeX-1
                PML.Ψ_Eyx[m_rev,nn,pp,2] = c.b_E_x[m_rev, 2] * PML.Ψ_Eyx[m_rev,nn,pp,2] + c.c_E_x[m_rev, 2] * 1/g.Δx * (F.Hz[mm-1,nn,pp] - F.Hz[mm,nn,pp])
                m_rev = m_rev - 1
            end
    end;end

    #@. PML.Ψ_Eyx[2:c.PML_Thickness[1], 1:end-1, 2:end-1, 1] = c.b_E_x[2:c.PML_Thickness[1], 1] * PML.Ψ_Eyx[2:c.PML_Thickness[1], 1:end-1, 2:end-1, 1] + c.b_E_x[2:c.PML_Thickness[1], 1] * 1/g.Δx *(F.Hz[1:c.PML_Thickness[1]-1, 1:end-1, 2:end-1] - F.Hz[2:c.PML_Thickness[1], 1:end-1, 2:end-1])
    #@. PML.Ψ_Eyx[end:-1:begin+1, 1:end-1, 2:end-1, 2] = c.b_E_x[end:-1:begin+1, 2] * PML.Ψ_Eyx[end:-1:begin+1, 1:end-1, 2:end-1, 2] + c.b_E_x[end:-1:begin+1,2]* 1/g.Δx *(F.Hz[(g.SizeX-c.PML_Thickness[1]):g.SizeX-2, 1:end-1, 2:end-1] - F.Hz[(g.SizeX+1-c.PML_Thickness[1]):g.SizeX-1, 1:end-1, 2:end-1])
    

    @inbounds for pp in 1:g.SizeZ-1
        for nn in 2:g.SizeY-1
            # X-Bot
            for mm in 2:c.PML_Thickness[1]
                PML.Ψ_Ezx[mm,nn,pp,1] = c.b_E_x[mm, 1] * PML.Ψ_Ezx[mm,nn,pp,1] + c.c_E_x[mm, 1] * 1/g.Δx *(F.Hy[mm,nn,pp] - F.Hy[mm-1,nn,pp])
            end
            # X-Top
            m_rev = c.PML_Thickness[1]
            for mm in (g.SizeX + 1 - c.PML_Thickness[1]):g.SizeX-1
                PML.Ψ_Ezx[m_rev,nn,pp,2] = c.b_E_x[m_rev, 2] * PML.Ψ_Ezx[m_rev,nn,pp,2] + c.c_E_x[m_rev, 2] * 1/g.Δx *(F.Hy[mm,nn,pp] - F.Hy[mm-1,nn,pp])
                m_rev = m_rev - 1  
            end
    end;end

    #Y-Boundary 
    @inbounds for pp in 2:g.SizeZ-1
        for mm in 1:g.SizeX-1
            # Y-Bot
            for nn in 2:c.PML_Thickness[2]
                PML.Ψ_Exy[mm,nn,pp,1] = c.b_E_y[nn, 1] * PML.Ψ_Exy[mm,nn,pp,1] + c.c_E_y[nn, 1] * 1/g.Δy *(F.Hz[mm,nn,pp] - F.Hz[mm,nn-1,pp])
            end
            # Y-Top
            n_rev = c.PML_Thickness[2]
            for nn in (g.SizeY+1-c.PML_Thickness[2]):g.SizeY-1
                PML.Ψ_Exy[mm,n_rev,pp,2] = c.b_E_y[n_rev, 2] * PML.Ψ_Exy[mm,n_rev,pp,2] + c.c_E_y[n_rev, 2] * 1/g.Δy * (F.Hz[mm,nn,pp] - F.Hz[mm,nn-1,pp])
                n_rev = n_rev - 1
            end
    end;end

    @inbounds for pp in 1:g.SizeZ-1
        for mm in 2:g.SizeX-1
            # Y-Bot
            for nn in 2:c.PML_Thickness[2]
                PML.Ψ_Ezy[mm,nn,pp,1] = c.b_E_y[nn, 1] * PML.Ψ_Ezy[mm,nn,pp,1] + c.c_E_y[nn, 1] * 1/g.Δy *(F.Hx[mm,nn-1,pp] - F.Hx[mm,nn,pp])
            end
            # Y-Top
            n_rev = c.PML_Thickness[2]
            for nn in (g.SizeY + 1 - c.PML_Thickness[2]):g.SizeY-1
                PML.Ψ_Ezy[mm,n_rev,pp,2] = c.b_E_y[n_rev, 2] * PML.Ψ_Ezy[mm,n_rev,pp,2] + c.c_E_y[n_rev, 2] * 1/g.Δy * (F.Hx[mm,nn-1,pp] - F.Hx[mm,nn,pp])
                n_rev = n_rev - 1  
            end
    end;end

    #Z-Boundary 
    @inbounds for mm in 1:g.SizeX-1
        for nn in 2:g.SizeY-1
            # Z-Bot
            for pp in 2:c.PML_Thickness[3]
                PML.Ψ_Exz[mm,nn,pp,1] = c.b_E_z[pp, 1] * PML.Ψ_Exz[mm,nn,pp,1] + c.c_E_z[pp, 1] * 1/g.Δz *(F.Hy[mm,nn,pp - 1] - F.Hy[mm,nn,pp])
            end
            # Z-Top
            p_rev = c.PML_Thickness[3]
            for pp in (g.SizeZ+1-c.PML_Thickness[3]):g.SizeZ-1
                PML.Ψ_Exz[mm,nn,p_rev,2] = c.b_E_z[p_rev, 2] * PML.Ψ_Exz[mm,nn,p_rev,2] + c.c_E_z[p_rev, 2] * 1/g.Δz * (F.Hy[mm,nn,pp-1] - F.Hy[mm,nn,pp])
                p_rev = p_rev - 1
            end
    end;end

    @inbounds for mm in 2:g.SizeX-1
        for nn in 1:g.SizeY-1
            # Z-Bot
            for pp in 2:c.PML_Thickness[3]
                PML.Ψ_Eyz[mm,nn,pp,1] = c.b_E_z[pp, 1] * PML.Ψ_Eyz[mm,nn,pp,1] + c.c_E_z[pp, 1] * 1/g.Δz *(F.Hx[mm,nn,pp] - F.Hx[mm,nn,pp-1])
            end
            # Z-Top
            p_rev = c.PML_Thickness[3]
            for pp in (g.SizeZ + 1 - c.PML_Thickness[3]):g.SizeZ-1
                PML.Ψ_Eyz[mm,nn,p_rev, 2] = c.b_E_z[p_rev, 2] * PML.Ψ_Eyz[mm,nn,p_rev,2] + c.c_E_z[p_rev, 2] * 1/g.Δz *(F.Hx[mm,nn,pp] - F.Hx[mm,nn,pp-1])
                p_rev = p_rev - 1  
            end
    end;end

end

function update_Ψ_E!(PML::CPML_Ψ_Fields_2D, F::Fields2D, g::Grid2D, c::CPML_Parameters_2D)
    #X-Boundary,  notice that there does not exist Ey in the TMz config 
   # @inbounds for nn in 1:g.SizeY-1
   #     # X-Bot
   #     for mm in 2:c.PML_Thickness[1]
   #         PML.Ψ_Eyx[mm,nn,1] = c.b_E_x[mm, 1] * PML.Ψ_Eyx[mm,nn,1] + c.c_E_x[mm, 1] * 1/g.Δx *(F.Hz[mm-1,nn] - F.Hz[mm,nn])
   #     end
   #     # X-Top
   #     m_rev = c.PML_Thickness[1]
   #     for mm in (g.SizeX+1-c.PML_Thickness[1]):g.SizeX-1
   #         PML.Ψ_Eyx[m_rev,nn,2] = c.b_E_x[m_rev, 2] * PML.Ψ_Eyx[m_rev,nn,2] + c.c_E_x[m_rev, 2] * 1/g.Δx * (F.Hz[mm-1,nn] - F.Hz[mm,nn])
   #         m_rev = m_rev - 1
   #     end
   # end

    @inbounds for nn in 2:g.SizeY-1
        # X-Bot
        for mm in 2:c.PML_Thickness[1]
            PML.Ψ_Ezx[mm,nn,1] = c.b_E_x[mm, 1] * PML.Ψ_Ezx[mm,nn,1] + c.c_E_x[mm, 1] * 1/g.Δx *(F.Hy[mm,nn] - F.Hy[mm-1,nn])
        end
        # X-Top
        m_rev = c.PML_Thickness[1]
        for mm in (g.SizeX + 1 - c.PML_Thickness[1]):g.SizeX-1
            PML.Ψ_Ezx[m_rev,nn,2] = c.b_E_x[m_rev, 2] * PML.Ψ_Ezx[m_rev,nn,2] + c.c_E_x[m_rev, 2] * 1/g.Δx *(F.Hy[mm,nn] - F.Hy[mm-1,nn])
            m_rev = m_rev - 1  
        end
    end

    #Y-Boundary, notice that there does not exist Ex in the TMz config  
    # @inbounds for mm in 1:g.SizeX-1
    #     # Y-Bot
    #     for nn in 2:c.PML_Thickness[2]
    #         PML.Ψ_Exy[mm,nn,1] = c.b_E_y[nn, 1] * PML.Ψ_Exy[mm,nn,1] + c.c_E_y[nn, 1] * 1/g.Δy *(F.Hz[mm,nn] - F.Hz[mm,nn-1])
    #     end
    #     # Y-Top
    #     n_rev = c.PML_Thickness[2]
    #     for nn in (g.SizeY+1-c.PML_Thickness[2]):g.SizeY-1
    #         PML.Ψ_Exy[mm,n_rev,2] = c.b_E_y[n_rev, 2] * PML.Ψ_Exy[mm,n_rev,2] + c.c_E_y[n_rev, 2] * 1/g.Δy * (F.Hz[mm,nn] - F.Hz[mm,nn-1])
    #         n_rev = n_rev - 1
    #     end
    # end

    @inbounds for mm in 2:g.SizeX-1
        # Y-Bot
        for nn in 2:c.PML_Thickness[2]
            PML.Ψ_Ezy[mm,nn,1] = c.b_E_y[nn, 1] * PML.Ψ_Ezy[mm,nn,1] + c.c_E_y[nn, 1] * 1/g.Δy *(F.Hx[mm,nn-1] - F.Hx[mm,nn])
        end
        # Y-Top
        n_rev = c.PML_Thickness[2]
        for nn in (g.SizeY + 1 - c.PML_Thickness[2]):g.SizeY-1
            PML.Ψ_Ezy[mm,n_rev,2] = c.b_E_y[n_rev, 2] * PML.Ψ_Ezy[mm,n_rev,2] + c.c_E_y[n_rev, 2] * 1/g.Δy * (F.Hx[mm,nn-1] - F.Hx[mm,nn])
            n_rev = n_rev - 1  
        end
    end
end


function update_Ψ_E!(PML::CPML_Ψ_Fields_1D, F::Fields1D, g::Grid1D, c::CPML_Parameters_1D)
    # this is the reason why in 1d it doesnt work properly, symmetry issues 
    #X-Boundary

    # # X-Bot
    # @inbounds for mm in 2:c.PML_Thickness[1]
    #     PML.Ψ_Eyx[mm,1] = c.b_E_x[mm, 1] * PML.Ψ_Eyx[mm,1] + c.c_E_x[mm, 1] * 1/g.Δx *(F.Hz[mm-1] - F.Hz[mm])
    # end
    
    # # X-Top
    # m_rev = c.PML_Thickness[1]
    # @inbounds for mm in (g.SizeX+1-c.PML_Thickness[1]):g.SizeX-1
    #     PML.Ψ_Eyx[m_rev, 2] = c.b_E_x[m_rev, 2] * PML.Ψ_Eyx[m_rev, 2] + c.c_E_x[m_rev, 2] * 1/g.Δx * (F.Hz[mm-1] - F.Hz[mm])
    #     m_rev = m_rev - 1
    # end


    # X-Top
    @inbounds for mm in 2:c.PML_Thickness[1]
        PML.Ψ_Ezx[mm,1] = c.b_E_x[mm, 1] * PML.Ψ_Ezx[mm,1] + c.c_E_x[mm, 1] * 1/g.Δx *(F.Hy[mm] - F.Hy[mm-1])
    end
    # X-Top
    m_rev = c.PML_Thickness[1]
    @inbounds for mm in (g.SizeX + 1 - c.PML_Thickness[1]):g.SizeX-1
        PML.Ψ_Ezx[m_rev,2] = c.b_E_x[m_rev, 2] * PML.Ψ_Ezx[m_rev, 2] + c.c_E_x[m_rev, 2] * 1/g.Δx *(F.Hy[mm] - F.Hy[mm-1])
        m_rev = m_rev - 1  
    end
end

function update_Ψ_H!(PML::CPML_Ψ_Fields_3D, F::Fields3D, g::Grid3D, c::CPML_Parameters_3D)
    #X-Boundary
    @inbounds for pp in 1:g.SizeZ-1
        for nn in 1:g.SizeY-1
            # X-Bot
            for mm in 1:c.PML_Thickness[1]-1
                PML.Ψ_Hyx[mm,nn,pp,1] = c.b_H_x[mm,1] * PML.Ψ_Hyx[mm,nn,pp,1] + c.c_H_x[mm,1]*1/g.Δx*(F.Ez[mm+1,nn,pp] - F.Ez[mm,nn,pp])
            end
            # X-top
            m_rev = c.PML_Thickness[1]-1
            for mm in (g.SizeX+1-c.PML_Thickness[1]):g.SizeX-1
                PML.Ψ_Hyx[m_rev,nn,pp,2] = c.b_H_x[m_rev,2] * PML.Ψ_Hyx[m_rev,nn,pp,2] + c.c_H_x[m_rev,2]*1/g.Δx*(F.Ez[mm+1,nn,pp] - F.Ez[mm,nn,pp])
                m_rev = m_rev - 1
            end
    end;end

    @inbounds for pp in 2:g.SizeZ-1
        for nn in 1:g.SizeY-1
            # X-Bot
            for mm in 1:c.PML_Thickness[1]-1
                PML.Ψ_Hzx[mm,nn,pp,1] = c.b_H_x[mm,1] * PML.Ψ_Hzx[mm,nn,pp,1] + c.c_H_x[mm, 1] * 1/g.Δx *(F.Ey[mm,nn,pp] - F.Ey[mm+1,nn,pp])
            end
            # X-Top
            m_rev = c.PML_Thickness[1]-1
            for mm in (g.SizeX+1-c.PML_Thickness[1]):g.SizeX-1
                PML.Ψ_Hzx[m_rev,nn,pp,2] = c.b_H_x[m_rev, 2] * PML.Ψ_Hzx[m_rev,nn,pp,2] +c.c_H_x[m_rev, 2] * 1/g.Δx *(F.Ey[mm,nn,pp] - F.Ey[mm+1,nn,pp])
                m_rev = m_rev - 1
            end
    end;end

    #Y-Boundary
    @inbounds for pp in 1:g.SizeZ-1
        for mm in 1:g.SizeX-1
            # Y-Bot
            for nn in 1:c.PML_Thickness[2]-1
                PML.Ψ_Hxy[mm,nn,pp,1] = c.b_H_y[nn,1] * PML.Ψ_Hxy[mm,nn,pp,1] + c.c_H_y[nn,1]*1/g.Δy*(F.Ez[mm,nn,pp] - F.Ez[mm,nn+1,pp])
            end
            # Y-top
            n_rev = c.PML_Thickness[2]-1
            for nn in (g.SizeY+1-c.PML_Thickness[2]):g.SizeY-1
                PML.Ψ_Hxy[mm,n_rev,pp,2] = c.b_H_y[n_rev,2] * PML.Ψ_Hxy[mm,n_rev,pp,2] + c.c_H_y[n_rev,2]*1/g.Δy*(F.Ez[mm,nn,pp] - F.Ez[mm,nn+1,pp])
                n_rev = n_rev - 1
            end
    end;end

    @inbounds for pp in 2:g.SizeZ-1
        for mm in 1:g.SizeX-1
            # Y-Bot
            for nn in 1:c.PML_Thickness[2]-1
                PML.Ψ_Hzy[mm,nn,pp,1] = c.b_H_y[nn,1] * PML.Ψ_Hzy[mm,nn,pp,1] + c.c_H_y[nn, 1] * 1/g.Δy * (F.Ex[mm,nn+1,pp] - F.Ex[mm,nn,pp])
            end
            # Y-Top
            n_rev = c.PML_Thickness[2]-1
            for nn in (g.SizeY+1-c.PML_Thickness[2]):g.SizeY-1
                PML.Ψ_Hzy[mm,n_rev,pp,2] = c.b_H_y[n_rev, 2] * PML.Ψ_Hzy[mm,n_rev,pp,2] +c.c_H_y[n_rev, 2] * 1/g.Δy * (F.Ex[mm,nn+1,pp] - F.Ex[mm,nn,pp])
                n_rev = n_rev - 1
            end
    end;end

    #Z-Boundary
    @inbounds for mm in 1:g.SizeX-1
        for nn in 1:g.SizeY-1
            # Z-Bot
            for pp in 1:c.PML_Thickness[3]-1
                PML.Ψ_Hyz[mm,nn,pp,1] = c.b_H_z[pp,1] * PML.Ψ_Hyz[mm,nn,pp,1] + c.c_H_z[pp,1] * 1/g.Δz*(F.Ex[mm,nn,pp] - F.Ex[mm,nn,pp+1])
            end
            # Z-top
            p_rev = c.PML_Thickness[3]-1
            for pp in (g.SizeZ+1-c.PML_Thickness[3]):g.SizeZ-1
                PML.Ψ_Hyz[mm,nn,p_rev,2] = c.b_H_z[p_rev,2] * PML.Ψ_Hyz[mm,nn,p_rev,2] + c.c_H_z[p_rev,2]*1/g.Δz*(F.Ex[mm,nn,pp] - F.Ex[mm,nn,pp+1])
                p_rev = p_rev - 1
            end
    end;end

    @inbounds for mm in 1:g.SizeX-1
        for nn in 1:g.SizeY-1
            # Z-Bot
            for pp in 1:c.PML_Thickness[3]-1
                PML.Ψ_Hxz[mm,nn,pp,1] = c.b_H_z[pp,1] * PML.Ψ_Hxz[mm,nn,pp,1] + c.c_H_z[pp, 1] * 1/g.Δz *(F.Ey[mm,nn,pp+1] - F.Ey[mm,nn,pp])
            end
            # Z-Top
            p_rev = c.PML_Thickness[3]-1
            for pp in (g.SizeZ+1-c.PML_Thickness[1]):g.SizeZ-1
                PML.Ψ_Hxz[mm,nn,p_rev,2] = c.b_H_z[p_rev, 2] * PML.Ψ_Hzx[mm,nn,p_rev, 2] +c.c_H_z[p_rev, 2] * 1/g.Δz *(F.Ey[mm,nn,pp+1] - F.Ey[mm,nn,pp])
                p_rev = p_rev - 1
            end
    end;end
end

# two dimensional update 

function update_Ψ_H!(PML::CPML_Ψ_Fields_2D, F::Fields2D, g::Grid2D, c::CPML_Parameters_2D)
    #X-Boundary
    @inbounds for nn in 1:g.SizeY-1
        # X-Bot
        for mm in 1:c.PML_Thickness[1]-1
            PML.Ψ_Hyx[mm,nn,1] = c.b_H_x[mm,1] * PML.Ψ_Hyx[mm,nn,1] + c.c_H_x[mm,1]*1/g.Δx*(F.Ez[mm+1,nn] - F.Ez[mm,nn])
        end
        # X-top
        m_rev = c.PML_Thickness[1]-1
        for mm in (g.SizeX+1-c.PML_Thickness[1]):g.SizeX-1
            PML.Ψ_Hyx[m_rev,nn,2] = c.b_H_x[m_rev,2] * PML.Ψ_Hyx[m_rev,nn,2] + c.c_H_x[m_rev,2]*1/g.Δx*(F.Ez[mm+1,nn] - F.Ez[mm,nn])
            m_rev = m_rev - 1
        end
    end

   # @inbounds for nn in 1:g.SizeY-1
   #     # X-Bot
   #     for mm in 1:c.PML_Thickness[1]-1
   #         PML.Ψ_Hzx[mm,nn,1] = c.b_H_x[mm,1] * PML.Ψ_Hzx[mm,nn,1] + c.c_H_x[mm, 1] * 1/g.Δx *(F.Ey[mm,nn] - F.Ey[mm+1,nn])
   #     end
   #     # X-Top
   #     m_rev = c.PML_Thickness[1]-1
   #     for mm in (g.SizeX+1-c.PML_Thickness[1]):g.SizeX-1
   #         PML.Ψ_Hzx[m_rev,nn,2] = c.b_H_x[m_rev, 2] * PML.Ψ_Hzx[m_rev,nn,2] +c.c_H_x[m_rev, 2] * 1/g.Δx *(F.Ey[mm,nn] - F.Ey[mm+1,nn])
   #         m_rev = m_rev - 1
   #     end
   # end

    #Y-Boundary
    @inbounds for mm in 1:g.SizeX-1
        # Y-Bot
        for nn in 1:c.PML_Thickness[2]-1
            PML.Ψ_Hxy[mm,nn,1] = c.b_H_y[nn,1] * PML.Ψ_Hxy[mm,nn,1] + c.c_H_y[nn,1]*1/g.Δy*(F.Ez[mm,nn] - F.Ez[mm,nn+1])
        end
        # Y-top
        n_rev = c.PML_Thickness[2]-1
        for nn in (g.SizeY+1-c.PML_Thickness[2]):g.SizeY-1
            PML.Ψ_Hxy[mm,n_rev,2] = c.b_H_y[n_rev,2] * PML.Ψ_Hxy[mm,n_rev,2] + c.c_H_y[n_rev,2]*1/g.Δy*(F.Ez[mm,nn] - F.Ez[mm,nn+1])
            n_rev = n_rev - 1
        end
    end

   # @inbounds for mm in 1:g.SizeX-1
   #     # Y-Bot
   #     for nn in 1:c.PML_Thickness[2]-1
   #         PML.Ψ_Hzy[mm,nn,1] = c.b_H_y[nn,1] * PML.Ψ_Hzy[mm,nn,1] + c.c_H_y[nn, 1] * 1/g.Δy * (F.Ex[mm,nn+1] - F.Ex[mm,nn])
   #     end
   #     # Y-Top
   #     n_rev = c.PML_Thickness[2]-1
   #     for nn in (g.SizeY+1-c.PML_Thickness[2]):g.SizeY-1
   #         PML.Ψ_Hzy[mm,n_rev,2] = c.b_H_y[n_rev, 2] * PML.Ψ_Hzy[mm,n_rev,2] +c.c_H_y[n_rev, 2] * 1/g.Δy * (F.Ex[mm,nn+1] - F.Ex[mm,nn])
   #         n_rev = n_rev - 1
   #     end
   # end

end


# one dimensional update
function update_Ψ_H!(PML::CPML_Ψ_Fields_1D, F::Fields1D, g::Grid1D, c::CPML_Parameters_1D)
    #X-Boundary
    
    # X-Bot
    @inbounds for mm in 1:c.PML_Thickness[1]-1
        PML.Ψ_Hyx[mm,1] = c.b_H_x[mm,1] * PML.Ψ_Hyx[mm,1] + c.c_H_x[mm,1]*1/g.Δx*(F.Ez[mm+1] - F.Ez[mm])
    end
    # X-top
    m_rev = c.PML_Thickness[1]-1
    @inbounds for mm in (g.SizeX+1-c.PML_Thickness[1]):g.SizeX-1
        PML.Ψ_Hyx[m_rev,2] = c.b_H_x[m_rev,2] * PML.Ψ_Hyx[m_rev,2] + c.c_H_x[m_rev,2]*1/g.Δx*(F.Ez[mm+1] - F.Ez[mm])
        m_rev = m_rev - 1
    end


    # # X-Bot
    # @inbounds for mm in 1:c.PML_Thickness[1]-1
    #     PML.Ψ_Hzx[mm,1] = c.b_H_x[mm,1] * PML.Ψ_Hzx[mm,1] + c.c_H_x[mm, 1] * 1/g.Δx *(F.Ey[mm] - F.Ey[mm+1])
    # end
    # # X-Top
    # m_rev = c.PML_Thickness[1]-1
    # @inbounds for mm in (g.SizeX+1-c.PML_Thickness[1]):g.SizeX-1
    #     PML.Ψ_Hzx[m_rev,2] = c.b_H_x[m_rev, 2] * PML.Ψ_Hzx[m_rev,2] +c.c_H_x[m_rev, 2] * 1/g.Δx *(F.Ey[mm] - F.Ey[mm+1])
    #     m_rev = m_rev - 1
    # end
end

function apply_Ψ_E!(PML::CPML_Ψ_Fields_3D, F::Fields3D, g::Grid3D, c::CPML_Parameters_3D)
    CB = (g.Δt/ϵ_0)
    CA = 1.

    @. F.Ex[1:g.SizeX-1, 2:c.PML_Thickness[2], 2:g.SizeZ-1] += CB * PML.Ψ_Exy[1:g.SizeX-1, 2:c.PML_Thickness[2], 2:g.SizeZ-1, 1]
    @. F.Ex[1:g.SizeX-1, (g.SizeY+1-c.PML_Thickness[2]):(g.SizeY-1), 2:g.SizeZ-1] += CB * PML.Ψ_Exy[1:g.SizeX-1,(c.PML_Thickness[2]):-1:begin+1, 2:g.SizeZ-1, 2]
    
    @. F.Ex[1:g.SizeX-1, 2:g.SizeY-1, 2:c.PML_Thickness[3]] += CB * PML.Ψ_Exz[1:g.SizeX-1, 2:g.SizeY-1, 2:c.PML_Thickness[3], 1]
    @. F.Ex[1:g.SizeX-1, 2:g.SizeY-1, (g.SizeZ+1-c.PML_Thickness[3]):(g.SizeZ-1)] += CB * PML.Ψ_Exz[1:g.SizeX-1, 2:g.SizeY-1,(c.PML_Thickness[3]):-1:begin+1, 2]

    @. F.Ey[2:c.PML_Thickness[1], 1:g.SizeY-1, 2:g.SizeZ-1] += CB * PML.Ψ_Eyx[2:c.PML_Thickness[1], 1:g.SizeY-1, 2:g.SizeZ-1, 1]
    @. F.Ey[(g.SizeX+1-c.PML_Thickness[1]):(g.SizeX-1), 1:g.SizeY-1, 2:g.SizeZ-1] += CB * PML.Ψ_Eyx[(c.PML_Thickness[1]):-1:begin+1,1:g.SizeY-1, 2:g.SizeZ-1, 2]

    @. F.Ey[2:g.SizeX-1, 1:g.SizeY-1, 2:c.PML_Thickness[3]] += CB * PML.Ψ_Eyz[2:g.SizeX-1, 1:g.SizeY-1, 2:c.PML_Thickness[3], 1]
    @. F.Ey[2:g.SizeX-1, 1:g.SizeY-1, (g.SizeZ+1-c.PML_Thickness[3]):(g.SizeZ-1)] += CB * PML.Ψ_Eyz[2:g.SizeX-1, 1:g.SizeY-1, (c.PML_Thickness[3]):-1:begin+1, 2]

    @. F.Ez[2:c.PML_Thickness[1], 2:g.SizeY-1, 1:g.SizeZ-1] += CB * PML.Ψ_Ezx[2:c.PML_Thickness[1], 2:g.SizeY-1, 1:g.SizeZ-1, 1]
    @. F.Ez[(g.SizeX+1-c.PML_Thickness[1]):(g.SizeX-1), 2:g.SizeY-1, 1:g.SizeZ-1] += CB * PML.Ψ_Ezx[(c.PML_Thickness[1]):-1:begin+1, 2:g.SizeY-1, 1:g.SizeZ-1, 2]

    @. F.Ez[2:g.SizeX-1, 2:c.PML_Thickness[2], 1:g.SizeZ-1] += CB * PML.Ψ_Ezy[2:g.SizeX-1, 2:c.PML_Thickness[2], 1:g.SizeZ-1, 1]
    @. F.Ez[2:g.SizeX-1, (g.SizeY+1-c.PML_Thickness[2]):(g.SizeY-1), 1:g.SizeZ-1] += CB * PML.Ψ_Ezy[2:g.SizeX-1, (c.PML_Thickness[2]):-1:begin+1, 1:g.SizeZ-1, 2]
end

function apply_Ψ_E!(PML::CPML_Ψ_Fields_2D, F::Fields2D, g::Grid2D, c::CPML_Parameters_2D)
    CB = (g.Δt/ϵ_0)
    CA = 1.
    @. F.Ez[2:c.PML_Thickness[1], 2:g.SizeY-1] += CB * PML.Ψ_Ezx[2:c.PML_Thickness[1], 2:g.SizeY-1, 1]
    @. F.Ez[(g.SizeX+1-c.PML_Thickness[1]):(g.SizeX-1), 2:g.SizeY-1] += CB * PML.Ψ_Ezx[(c.PML_Thickness[1]):-1:begin+1, 2:g.SizeY-1, 2]

    @. F.Ez[2:g.SizeX-1, 2:c.PML_Thickness[2]] += CB * PML.Ψ_Ezy[2:g.SizeX-1, 2:c.PML_Thickness[2], 1]
    @. F.Ez[2:g.SizeX-1, (g.SizeY+1-c.PML_Thickness[2]):(g.SizeY-1)] += CB * PML.Ψ_Ezy[2:g.SizeX-1, (c.PML_Thickness[2]):-1:begin+1, 2]
end


function apply_Ψ_E!(PML::CPML_Ψ_Fields_1D, F::Fields1D, g::Grid1D, c::CPML_Parameters_1D)
    CB = (g.Δt/ϵ_0)
    CA = 1.

    @. F.Ez[2:c.PML_Thickness[1]] += CB * PML.Ψ_Ezx[2:c.PML_Thickness[1], 1]
    @. F.Ez[(g.SizeX+1-c.PML_Thickness[1]):(g.SizeX-1)] += CB * PML.Ψ_Ezx[(c.PML_Thickness[1]):-1:begin+1, 2]

end

function apply_Ψ_H!(PML::CPML_Ψ_Fields_3D, F::Fields3D, g::Grid3D, c::CPML_Parameters_3D)
    DB = g.Δt/μ_0

    @. F.Hx[1:g.SizeX-1, 1:c.PML_Thickness[2]-1, 1:g.SizeZ-1] += DB * PML.Ψ_Hxy[1:g.SizeX-1, 1:c.PML_Thickness[2]-1, 1:g.SizeZ-1, 1]
    @. F.Hx[1:g.SizeX-1, (g.SizeY+1-c.PML_Thickness[2]):(g.SizeY-1), 1:g.SizeZ-1] += DB * PML.Ψ_Hxy[1:g.SizeX-1, (c.PML_Thickness[2]-1):-1:begin, 1:g.SizeZ-1, 2]
    
    @. F.Hx[1:g.SizeX-1, 1:g.SizeY-1, 1:c.PML_Thickness[3]-1] += DB * PML.Ψ_Hxz[1:g.SizeX-1, 1:g.SizeY-1, 1:c.PML_Thickness[3]-1, 1]
    @. F.Hx[1:g.SizeX-1, 1:g.SizeY-1, (g.SizeZ+1-c.PML_Thickness[3]):(g.SizeZ-1)] += DB * PML.Ψ_Hxz[1:g.SizeX-1, 1:g.SizeY-1, (c.PML_Thickness[3]-1):-1:begin, 2] # recheck

    @. F.Hy[1:c.PML_Thickness[1]-1, 1:g.SizeY-1, 1:g.SizeZ-1] += DB * PML.Ψ_Hyx[1:c.PML_Thickness[1]-1, 1:g.SizeY-1, 1:g.SizeZ-1, 1]
    @. F.Hy[(g.SizeX+1-c.PML_Thickness[1]):(g.SizeX-1), 1:g.SizeY-1, 1:g.SizeZ-1] += DB * PML.Ψ_Hyx[(c.PML_Thickness[1]-1):-1:begin, 1:g.SizeY-1, 1:g.SizeZ-1, 2]

    @. F.Hy[1:g.SizeX-1, 1:g.SizeY-1, 1:c.PML_Thickness[3]-1] += DB * PML.Ψ_Hyz[1:g.SizeX-1, 1:g.SizeY-1, 1:c.PML_Thickness[3]-1, 1]
    @. F.Hy[1:g.SizeX-1, 1:g.SizeY-1, (g.SizeZ+1-c.PML_Thickness[3]):(g.SizeZ-1)] += DB * PML.Ψ_Hyz[1:g.SizeX-1, 1:g.SizeY-1, (c.PML_Thickness[3]-1):-1:begin, 2]

    @. F.Hz[1:c.PML_Thickness[1]-1, 1:g.SizeY-1, 2:g.SizeZ-1] += DB * PML.Ψ_Hzx[1:c.PML_Thickness[1]-1, 1:g.SizeY-1, 2:g.SizeZ-1, 1]
    @. F.Hz[(g.SizeX+1-c.PML_Thickness[1]):(g.SizeX-1), 1:g.SizeY-1, 2:g.SizeZ-1] += DB * PML.Ψ_Hzx[(c.PML_Thickness[1]-1):-1:begin, 1:g.SizeY-1, 2:g.SizeZ-1, 2]

    @. F.Hz[1:g.SizeX-1, 1:c.PML_Thickness[2]-1, 2:g.SizeZ-1] += DB * PML.Ψ_Hzy[1:g.SizeX-1, 1:c.PML_Thickness[2]-1, 2:g.SizeZ-1, 1]
    @. F.Hz[1:g.SizeX-1, (g.SizeY+1-c.PML_Thickness[2]):(g.SizeY-1), 2:g.SizeZ-1] += DB * PML.Ψ_Hzy[1:g.SizeX-1, (c.PML_Thickness[2]-1):-1:begin, 2:g.SizeZ-1, 2]
end


function apply_Ψ_H!(PML::CPML_Ψ_Fields_2D, F::Fields2D, g::Grid2D, c::CPML_Parameters_2D)
    DB = g.Δt/μ_0

    @. F.Hx[1:g.SizeX-1, 1:c.PML_Thickness[2]-1] += DB * PML.Ψ_Hxy[1:g.SizeX-1, 1:c.PML_Thickness[2]-1, 1]
    @. F.Hx[1:g.SizeX-1, (g.SizeY+1-c.PML_Thickness[2]):(g.SizeY-1)] += DB * PML.Ψ_Hxy[1:g.SizeX-1, (c.PML_Thickness[2]-1):-1:begin, 2]
    
    @. F.Hy[1:c.PML_Thickness[1]-1, 1:g.SizeY-1] += DB * PML.Ψ_Hyx[1:c.PML_Thickness[1]-1, 1:g.SizeY-1, 1]
    @. F.Hy[(g.SizeX+1-c.PML_Thickness[1]):(g.SizeX-1), 1:g.SizeY-1] += DB * PML.Ψ_Hyx[(c.PML_Thickness[1]-1):-1:begin, 1:g.SizeY-1, 2]

end

function apply_Ψ_H!(PML::CPML_Ψ_Fields_1D, F::Fields1D, g::Grid1D, c::CPML_Parameters_1D)
    DB = g.Δt/μ_0

    @. F.Hy[1:c.PML_Thickness[1]-1] += DB * PML.Ψ_Hyx[1:c.PML_Thickness[1]-1, 1]
    @. F.Hy[(g.SizeX+1-c.PML_Thickness[1]):(g.SizeX-1)] += DB * PML.Ψ_Hyx[(c.PML_Thickness[1]-1):-1:begin, 2]

end
