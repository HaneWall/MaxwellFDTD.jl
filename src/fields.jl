abstract type Field end

#=
 These are the Fields in  the one-dimesional Case.
=#
mutable struct Fields1D <: Field
    Hy :: Vector{Float64}
    Ez :: Vector{Float64}
    function Fields1D(g::Grid1D)
        new(zeros(Float64, (g.SizeX)), zeros(Float64, (g.SizeX)))
    end
end

mutable struct MaterialFields1D <:Field
    Jz :: Vector{Float64}
    Pz :: Vector{Float64}
    PzNl :: Vector{Float64}

    Jz_bound :: Vector{Float64}
    Jz_free :: Vector{Float64}
    Jz_tunnel :: Vector{Float64}

    ρ_cb :: Vector{Float64}
    Γ_ADK :: Vector{Float64}
    function MaterialFields1D(g::Grid1D)
        new(zeros(Float64, g.SizeX),
            zeros(Float64, g.SizeX),
            zeros(Float64, g.SizeX),
            zeros(Float64, g.SizeX),
            zeros(Float64, g.SizeX),
            zeros(Float64, g.SizeX),
            zeros(Float64, g.SizeX),
            zeros(Float64, g.SizeX)
            )
    end
end

mutable struct LorentzFields1D <: Field
    Jz :: Array{Float64}
    Pz :: Array{Float64}
    PzNl :: Array{Float64}
    function LorentzFields1D(m::LorentzMedium1D)
        new(
            zeros(Float64, (size(m.location)[1], m.oscillators)), 
            zeros(Float64, (size(m.location)[1], m.oscillators)),
            zeros(Float64, (size(m.location)[1], m.oscillators))
            )
    end
end

mutable struct DrudeFields1D <: Field
    Jz_free :: Array{Float64}
    function DrudeFields1D(m::DrudeMedium1D)
        new(
            zeros(Float64, (size(m.location)[1]))
        )
    end
end

mutable struct TunnelFields1D <: Field
    Jz_tunnel :: Vector{Float64}
    dz_T :: Vector{Float64}
    function TunnelFields1D(m::TunnelMedium1D)
        new(
            zeros(Float64, (size(m.location)[1])), 
            zeros(Float64, (size(m.location)[1]))
            )
    end
end

#=
 These are the Fields in  the two-dimesional Case.
=#

#=
 These are the Fields in  the three-dimesional Case.
=#
mutable struct Fields3D <: Field
    Hx :: Array{Float64, 3}
    Hy :: Array{Float64, 3}
    Hz :: Array{Float64, 3}
    Ex :: Array{Float64, 3}
    Ey :: Array{Float64, 3}
    Ez :: Array{Float64, 3}
    function Fields3D(g::Grid3D)
        new(zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ)), 
            zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ)), 
            zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ)), 
            zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ)), 
            zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ)), 
            zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ))
            )
    end
end

mutable struct MaterialFields3D <:Field
    Jx :: Array{Float64, 3}
    Jy :: Array{Float64, 3}
    Jz :: Array{Float64, 3}
    Px :: Array{Float64, 3}
    Py :: Array{Float64, 3}
    Pz :: Array{Float64, 3}
    function MaterialFields3D(g::Grid3D)
        new(zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ)), 
            zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ)), 
            zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ)), 
            zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ)), 
            zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ)), 
            zeros(Float64, (g.SizeX, g.SizeY, g.SizeZ))
            )
    end
end

mutable struct LorentzFields3D <: Field
    Jx :: Array{Float64}
    Jy :: Array{Float64}
    Jz :: Array{Float64}
    Px :: Array{Float64}
    Py :: Array{Float64}
    Pz :: Array{Float64}
    PxNl :: Array{Float64}
    PyNl :: Array{Float64}
    PzNl :: Array{Float64}
    function LorentzFields3D(m::LorentzMedium3D)
        new(
            zeros(Float64, (size(m.location)[1],size(m.location)[2],size(m.location)[3], m.oscillators)), 
            zeros(Float64, (size(m.location)[1],size(m.location)[2],size(m.location)[3], m.oscillators)),
            zeros(Float64, (size(m.location)[1],size(m.location)[2],size(m.location)[3], m.oscillators)),
            zeros(Float64, (size(m.location)[1],size(m.location)[2],size(m.location)[3], m.oscillators)),
            zeros(Float64, (size(m.location)[1],size(m.location)[2],size(m.location)[3], m.oscillators)),
            zeros(Float64, (size(m.location)[1],size(m.location)[2],size(m.location)[3], m.oscillators)),
            zeros(Float64, (size(m.location)[1],size(m.location)[2],size(m.location)[3], m.oscillators)),
            zeros(Float64, (size(m.location)[1],size(m.location)[2],size(m.location)[3], m.oscillators)),
            zeros(Float64, (size(m.location)[1],size(m.location)[2],size(m.location)[3], m.oscillators)),
            )
    end
end

mutable struct CPML_Ψ_Fields_3D
    #=
        the last dimension is necessary to prevent creating 2 Array for top and bottom PML (--> symmetrical PML)
            example: Bottom: Ψ_Exy1 := Ψ_Exy[:, :, :, 1] , Top: Ψ_Exy2 := Ψ_Exy[:, :, :, 2]
    =#
    Ψ_Exy :: Array{Float64, 4}
    Ψ_Exz :: Array{Float64, 4}
   
    Ψ_Eyx :: Array{Float64, 4}
    Ψ_Eyz :: Array{Float64, 4}
    
    Ψ_Ezx :: Array{Float64, 4}
    Ψ_Ezy :: Array{Float64, 4}
    
    Ψ_Hxy :: Array{Float64, 4}
    Ψ_Hxz :: Array{Float64, 4}

    Ψ_Hyx :: Array{Float64, 4}
    Ψ_Hyz :: Array{Float64, 4}

    Ψ_Hzx :: Array{Float64, 4}
    Ψ_Hzy :: Array{Float64, 4}

    # inner constructor function:
    function CPML_Ψ_Fields_3D(g::Grid3D, PML_Thickness::Vector{Int64})
        
        # x-Boundaries
        Ψ_Ezx = zeros(Float64, PML_Thickness[1], g.SizeY, g.SizeZ, 2)
        Ψ_Eyx = zeros(Float64, PML_Thickness[1], g.SizeY - 1, g.SizeZ, 2)
        
        Ψ_Hyx = zeros(Float64, PML_Thickness[1] - 1, g.SizeY, g.SizeZ, 2)
        Ψ_Hzx = zeros(Float64, PML_Thickness[1] - 1 , g.SizeY - 1, g.SizeZ, 2)
        
        # y-Boundaries
        Ψ_Ezy = zeros(Float64, g.SizeX, PML_Thickness[2], g.SizeZ, 2)
        Ψ_Exy = zeros(Float64, g.SizeX - 1, PML_Thickness[2], g.SizeZ, 2)

        Ψ_Hxy = zeros(Float64, g.SizeX, PML_Thickness[2] - 1 , g.SizeZ, 2)
        Ψ_Hzy = zeros(Float64, g.SizeX - 1, PML_Thickness[2] - 1 , g.SizeZ, 2) 

        # z-Boundaries
        Ψ_Exz = zeros(Float64, g.SizeX - 1, g.SizeY, PML_Thickness[3], 2)
        Ψ_Eyz = zeros(Float64, g.SizeX, g.SizeY - 1, PML_Thickness[3], 2)
        
        Ψ_Hxz = zeros(Float64, g.SizeX, g.SizeY - 1, PML_Thickness[3] - 1, 2)
        Ψ_Hyz = zeros(Float64, g.SizeX - 1, g.SizeY, PML_Thickness[3] - 1, 2)
        
        new(
            Ψ_Exy, Ψ_Exz, Ψ_Eyx, Ψ_Eyz, Ψ_Ezx, Ψ_Ezy, Ψ_Hxy, Ψ_Hxz, Ψ_Hyx, Ψ_Hyz, Ψ_Hzx, Ψ_Hzy
            )
    end
end


mutable struct CPML_Ψ_Fields_1D
    #=
        the last dimension is necessary to prevent creating 2 Array for top and bottom PML (--> symmetrical PML)
            example: Bottom: Ψ_Exy1 := Ψ_Exy[:, :, :, 1] , Top: Ψ_Exy2 := Ψ_Exy[:, :, :, 2]
    =#
    Ψ_Ezx :: Array{Float64, 2}
    Ψ_Hyx :: Array{Float64, 2}

    # inner constructor function:
    function CPML_Ψ_Fields_1D(g::Grid1D, PML_Thickness::Vector{Int64})
        
        # x-Boundaries
        Ψ_Ezx = zeros(Float64, PML_Thickness[1], 2)
    
        Ψ_Hyx = zeros(Float64, PML_Thickness[1] - 1, 2)
        
        new(
            Ψ_Ezx, Ψ_Hyx
            )
    end
end


function σ_opt(m::Float64, Δx::Float64)
    return 0.8*(m+1)/(376.730 * Δx)
end

function σ_profile(arr::Array{Float64, 1}, thickness::Int64, Δx::Float64; m::Float64=4.0)
    return σ_opt(m, Δx) .* ((thickness .- arr) ./ (thickness - 1)).^m
end

function σ_profile(idx::Int64, thickness::Int64, Δx::Float64; m::Float64=4.0)
    return σ_opt(m, Δx) * (idx - 1 / (thickness - 1)).^m
end

function b_coeff(σ::Array{Float64}, κ_profile::Array{Float64}, α_profile::Array{Float64}, g::T) where T<:Grid
    return exp.(-(σ./(ϵ_0 .* κ_profile) .+ α_profile./ϵ_0) .* g.Δt)
end

function c_coeff(b::Array{Float64}, σ::Array{Float64}, κ::Array{Float64}, α::Array{Float64})
    return (b .- 1).*σ ./ (σ .* κ .+ α .* κ.^2)
end

function κ_profile(arr::Array{Float64, 1}, thickness::Int64; m::Float64=4.0, κ_max::Float64=14.)
    return 1 .+ (κ_max .-1).*((thickness .- arr)./(thickness - 1.)).^m
end

function α_profile(arr::Array{Float64, 1}, thickness::Int64; m_a::Float64=1.0)
    #return 0.1
    return 0.25*((arr)./(thickness .- 1)).^m_a
end

function α_profile(idx::Int64, thickness::Int64; m_a::Float64=1.0)
    return 0.25* ((thickness - idx)./thickness)^m_a
end
struct CPML_Parameters_3D
    # location_x_bot :: CartesianIndices(Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}) 
    # location_x_top :: CartesianIndices(Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}})
   
    # location_y_bot :: CartesianIndices(Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}})
    # location_y_top :: CartesianIndices(Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}})
    
    # location_z_bot :: CartesianIndices(Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}})
    # location_z_top :: CartesianIndices(Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}})
    PML_Thickness :: Vector{Int64}

    α_E_x :: Array{Float64, 2}
    α_E_y :: Array{Float64, 2}
    α_E_z :: Array{Float64, 2}
    κ_E_x :: Array{Float64, 2}
    κ_E_y :: Array{Float64, 2}
    κ_E_z :: Array{Float64, 2}
    σ_E_x :: Array{Float64, 2}
    σ_E_y :: Array{Float64, 2}
    σ_E_z :: Array{Float64, 2}

    b_E_x :: Array{Float64, 2}
    b_E_y :: Array{Float64, 2}
    b_E_z :: Array{Float64, 2}    

    c_E_x :: Array{Float64, 2}
    c_E_y :: Array{Float64, 2}
    c_E_z :: Array{Float64, 2}


    α_H_x :: Array{Float64, 2}
    α_H_y :: Array{Float64, 2}
    α_H_z :: Array{Float64, 2}
    κ_H_x :: Array{Float64, 2}
    κ_H_y :: Array{Float64, 2}
    κ_H_z :: Array{Float64, 2}
    σ_H_x :: Array{Float64, 2}
    σ_H_y :: Array{Float64, 2}
    σ_H_z :: Array{Float64, 2}

    b_H_x :: Array{Float64, 2}
    b_H_y :: Array{Float64, 2}
    b_H_z :: Array{Float64, 2} 

    c_H_x :: Array{Float64, 2}
    c_H_y :: Array{Float64, 2}
    c_H_z :: Array{Float64, 2} 

    denominator_E_x :: Array{Float64, 1}
    denominator_H_x :: Array{Float64, 1}

    denominator_E_y :: Array{Float64, 1}
    denominator_H_y :: Array{Float64, 1}
    
    denominator_E_z :: Array{Float64, 1}
    denominator_H_z :: Array{Float64, 1}

    # inner constructor
    function CPML_Parameters_3D(g::Grid3D, PML_Thickness::Vector{Int64})
        b_E_x, c_E_x, α_E_x, σ_E_x, κ_E_x = [zeros(Float64, PML_Thickness[1], 2) for _ in 1:5]
        b_H_x, c_H_x, α_H_x, σ_H_x, κ_H_x = [zeros(Float64, PML_Thickness[1] - 1, 2) for _ in 1:5]
        
        b_E_y, c_E_y, α_E_y, σ_E_y, κ_E_y = [zeros(Float64, PML_Thickness[2], 2) for _ in 1:5]
        b_H_y, c_H_y, α_H_y, σ_H_y, κ_H_y = [zeros(Float64, PML_Thickness[2] - 1, 2) for _ in 1:5]

        b_E_z, c_E_z, α_E_z, σ_E_z, κ_E_z = [zeros(Float64, PML_Thickness[3], 2) for _ in 1:5]
        b_H_z, c_H_z, α_H_z, σ_H_z, κ_H_z = [zeros(Float64, PML_Thickness[3] - 1, 2) for _ in 1:5]

        # Bot X
        arr_E_x = Array(1.:1.:PML_Thickness[1])
        #arr_E_x = Array(PML_Thickness[1]:-1.:1)
        σ_E_x[:, 1] = σ_profile(arr_E_x, PML_Thickness[1], g.Δx)
        α_E_x[:, 1] = α_profile(arr_E_x .- 1, PML_Thickness[1])
        κ_E_x[:, 1] = κ_profile(arr_E_x, PML_Thickness[1])
        b_E_x[:, 1] = b_coeff(σ_E_x[:, 1], κ_E_x[:, 1], α_E_x[:, 1], g)
        c_E_x[2:end-1, 1] = c_coeff(b_E_x[2:end-1, 1], σ_E_x[2:end-1, 1], κ_E_x[2:end-1, 1], α_E_x[2:end-1, 1])

        arr_H_x = Array(0.5:1.:PML_Thickness[1]-1.5)
        #arr_H_x = Array(PML_Thickness[1]-1.5:-1.:0.5)
        σ_H_x[:, 1] = σ_profile(arr_H_x, PML_Thickness[1], g.Δx)
        α_H_x[:, 1] = α_profile(arr_H_x .- 1, PML_Thickness[1])
        κ_H_x[:, 1] = κ_profile(arr_H_x, PML_Thickness[1])
        b_H_x[:, 1] = b_coeff(σ_H_x[:, 1], κ_H_x[:, 1], α_H_x[:, 1], g)
        c_H_x[:, 1] = c_coeff(b_H_x[:,1], σ_H_x[:,1], κ_H_x[:,1], α_H_x[:,1])

        # Top X 
        σ_E_x[:, 2] = σ_E_x[:, 1]
        α_E_x[:, 2] = α_E_x[:, 1]
        κ_E_x[:, 2] = κ_E_x[:, 1]
        b_E_x[:, 2] = b_E_x[:, 1]
        c_E_x[:, 2] = c_E_x[:, 1]

        σ_H_x[:, 2] = σ_H_x[:, 1]
        α_H_x[:, 2] = α_H_x[:, 1]
        κ_H_x[:, 2] = κ_H_x[:, 1]
        b_H_x[:, 2] = b_H_x[:, 1]
        c_H_x[:, 2] = c_H_x[:, 1]

        # Bot Y 
        arr_E_y = Array(1.:1.:PML_Thickness[2])
        #arr_E_y = Array(PML_Thickness[2]:-1.:1)

        α_E_y[:, 1] = α_profile(arr_E_y .- 1, PML_Thickness[2])
        σ_E_y[:, 1] = σ_profile(arr_E_y, PML_Thickness[2], g.Δy)
        κ_E_y[:, 1] = κ_profile(arr_E_y, PML_Thickness[2])
        b_E_y[:, 1] = b_coeff(σ_E_y[:,1], κ_E_y[:,1], α_E_y[:,1], g)
        c_E_y[2:end-1, 1] = c_coeff(b_E_y[2:end-1,1], σ_E_y[2:end-1,1], κ_E_y[2:end-1,1], α_E_y[2:end-1,1])

        arr_H_y = Array(0.5:1.:PML_Thickness[1]-1.5)
        #arr_H_y = Array(PML_Thickness[2]-1.5:-1.:0.5)

        σ_H_y[:, 1] = σ_profile(arr_H_y, PML_Thickness[2], g.Δy)
        α_H_y[:, 1] = α_profile(arr_H_y .- 1, PML_Thickness[2])
        κ_H_y[:, 1] = κ_profile(arr_H_y, PML_Thickness[2])
        b_H_y[:, 1] = b_coeff(σ_H_y[:,1], κ_H_y[:,1], α_H_y[:,1], g)
        c_H_y[:, 1] = c_coeff(b_H_y[:,1], σ_H_y[:,1], κ_H_y[:,1], α_H_y[:,1])

        # Top Y 
        σ_E_y[:, 2] = σ_E_y[:, 1]
        α_E_y[:, 2] = α_E_y[:, 1]
        κ_E_y[:, 2] = κ_E_y[:, 1]
        b_E_y[:, 2] = b_E_y[:, 1]
        c_E_y[:, 2] = c_E_y[:, 1]

        σ_H_y[:, 2] = σ_H_y[:, 1]
        α_H_y[:, 2] = α_H_y[:, 1]
        κ_H_y[:, 2] = κ_H_y[:, 1]
        b_H_y[:, 2] = b_H_y[:, 1]
        c_H_y[:, 2] = c_H_y[:, 1]

         # Bot Z
         arr_E_z = Array(1.:1.:PML_Thickness[3])
         #arr_E_z = Array(PML_Thickness[3]:-1.:1)

         α_E_z[:, 1] = α_profile(arr_E_z .- 1, PML_Thickness[3])
         σ_E_z[:, 1] = σ_profile(arr_E_z, PML_Thickness[3], g.Δz)
         κ_E_z[:, 1] = κ_profile(arr_E_z, PML_Thickness[3])
         b_E_z[:, 1] = b_coeff(σ_E_z[:, 1], κ_E_z[:, 1], α_E_z[:, 1], g)
         c_E_z[2:end-1, 1] = c_coeff(b_E_z[2:end-1, 1], σ_E_z[2:end-1, 1], κ_E_z[2:end-1, 1], α_E_z[2:end-1, 1])
 
         arr_H_z = Array(0.5:1.:PML_Thickness[3]-1.5)
         #arr_H_z = Array(PML_Thickness[3]-1.5:-1.:0.5)
 
         σ_H_z[:, 1] = σ_profile(arr_H_z, PML_Thickness[3], g.Δy)
         α_H_z[:, 1] = α_profile(arr_H_z .- 1, PML_Thickness[3])
         κ_H_z[:, 1] = κ_profile(arr_H_z, PML_Thickness[3])
         b_H_z[:, 1] = b_coeff(σ_H_z[:, 1], κ_H_z[:, 1], α_H_z[:, 1], g)
         c_H_z[:, 1] = c_coeff(b_H_z[:, 1], σ_H_z[:, 1], κ_H_z[:, 1], α_H_z[:, 1])
 
         # Top Z 
         σ_E_z[:, 2] = σ_E_z[:, 1]
         α_E_z[:, 2] = α_E_z[:, 1]
         κ_E_z[:, 2] = κ_E_z[:, 1]
         b_E_z[:, 2] = b_E_z[:, 1]
         c_E_z[:, 2] = c_E_z[:, 1]
         σ_H_z[:, 2] = σ_H_z[:, 1]
         α_H_z[:, 2] = α_H_z[:, 1]
         κ_H_z[:, 2] = κ_H_z[:, 1]
         b_H_z[:, 2] = b_H_z[:, 1]
         c_H_z[:, 2] = c_H_z[:, 1]

         denominator_E_x, denominator_H_x = [zeros(Float64, g.SizeX-1) for _ in 1:2]
         denominator_E_y, denominator_H_y = [zeros(Float64, g.SizeY-1) for _ in 1:2]
         denominator_E_z, denominator_H_z = [zeros(Float64, g.SizeZ-1) for _ in 1:2]

        ii_E = PML_Thickness[1]
        @inbounds for mm in 1:g.SizeX-1
            if mm <= PML_Thickness[1]
                denominator_E_x[mm] = 1/(κ_E_x[mm, 1] * g.Δx)
            elseif mm >= g.SizeX + 1 - PML_Thickness[1]
                denominator_E_x[mm] = 1/(κ_E_x[ii_E, 2] * g.Δx) 
                ii_E = ii_E - 1
            else 
                denominator_E_x[mm] = 1/g.Δx
            end   
        end

        ii_H = PML_Thickness[1] - 1
        @inbounds for mm in 1:g.SizeX-1
            if mm <= PML_Thickness[1]-1
                denominator_H_x[mm] = 1/(κ_H_x[mm, 1] * g.Δx) 
            elseif mm >= g.SizeX + 1 - PML_Thickness[1]
                denominator_H_x[mm] = 1/(κ_H_x[ii_H, 2] * g.Δx) 
                ii_H = ii_H - 1
            else 
                denominator_H_x[mm] = 1/g.Δx
            end   
        end


        jj_E = PML_Thickness[2]
        @inbounds for mm in 1:g.SizeY-1
            if mm <= PML_Thickness[2]
                denominator_E_y[mm] = 1/(κ_E_y[mm, 1] * g.Δy)
            elseif mm >= g.SizeY + 1 - PML_Thickness[2]
                denominator_E_y[mm] = 1/(κ_E_y[jj_E, 2] * g.Δy)
                jj_E = jj_E - 1
            else 
                denominator_E_y[mm] = 1/g.Δy
            end   
        end

        jj_H = PML_Thickness[2] - 1
        @inbounds for mm in 1:g.SizeY-1
            if mm <= PML_Thickness[2] - 1 
                denominator_H_y[mm] = 1/(κ_H_y[mm, 1] * g.Δy) 
            elseif mm >= g.SizeY + 1 - PML_Thickness[2]
                denominator_H_y[mm] = 1/(κ_H_y[jj_H, 2] * g.Δy) 
                jj_H = jj_H - 1
            else 
                denominator_H_y[mm] = 1/g.Δy
            end   
        end


        kk_E = PML_Thickness[3]
        @inbounds for mm in 1:g.SizeZ-1
            if mm <= PML_Thickness[3]
                denominator_E_z[mm] = 1/(κ_E_z[mm, 1] * g.Δz)
            elseif mm >= g.SizeZ + 1 - PML_Thickness[3]
                denominator_E_z[mm] = 1/(κ_E_z[kk_E, 2] * g.Δz)
                kk_E = kk_E - 1
            else 
                denominator_E_z[mm] = 1/g.Δz
            end   
        end

        kk_H = PML_Thickness[3] - 1 
        @inbounds for mm in 1:g.SizeZ-1
            if mm <= PML_Thickness[3] - 1 
                denominator_H_z[mm] = 1/(κ_H_z[mm, 1] * g.Δz) 
            elseif mm >= g.SizeZ + 1 - PML_Thickness[3]
                denominator_H_z[mm] = 1/(κ_H_z[kk_H, 2] * g.Δz) 
                kk_H = kk_H - 1
            else 
                denominator_H_z[mm] = 1/g.Δz
            end   
        end

        new(PML_Thickness,
            α_E_x,
            α_E_y,
            α_E_z,
            κ_E_x,
            κ_E_y,
            κ_E_z,
            σ_E_x,
            σ_E_y,
            σ_E_z,

            b_E_x,
            b_E_y,
            b_E_z,    

            c_E_x,
            c_E_y,
            c_E_z,

            α_H_x,
            α_H_y,
            α_H_z,
            κ_H_x,
            κ_H_y,
            κ_H_z,
            σ_H_x,
            σ_H_y,
            σ_H_z,

            b_H_x,
            b_H_y,
            b_H_z, 

            c_H_x,
            c_H_y,
            c_H_z,

            denominator_E_x, 
            denominator_H_x, 
            denominator_E_y, 
            denominator_H_y, 
            denominator_E_z, 
            denominator_H_z 
        )
    end
end

struct CPML_Parameters_1D
    # location_x_bot :: CartesianIndices(Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}) 
    # location_x_top :: CartesianIndices(Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}})
   
    # location_y_bot :: CartesianIndices(Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}})
    # location_y_top :: CartesianIndices(Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}})
    
    # location_z_bot :: CartesianIndices(Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}})
    # location_z_top :: CartesianIndices(Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}})
    PML_Thickness :: Vector{Int64}

    α_E_x :: Array{Float64, 2}
    κ_E_x :: Array{Float64, 2}
    σ_E_x :: Array{Float64, 2}

    b_E_x :: Array{Float64, 2}

    c_E_x :: Array{Float64, 2}

    α_H_x :: Array{Float64, 2}
    κ_H_x :: Array{Float64, 2}
    σ_H_x :: Array{Float64, 2}

    b_H_x :: Array{Float64, 2}

    c_H_x :: Array{Float64, 2}

    denominator_E_x :: Array{Float64, 1}
    denominator_H_x :: Array{Float64, 1}

    # inner constructor
    function CPML_Parameters_1D(g::Grid1D, PML_Thickness::Vector{Int64})
        b_E_x, c_E_x, α_E_x, σ_E_x, κ_E_x = [zeros(Float64, PML_Thickness[1], 2) for _ in 1:5]
        b_H_x, c_H_x, α_H_x, σ_H_x, κ_H_x = [zeros(Float64, PML_Thickness[1] - 1, 2) for _ in 1:5]

        # Bot X
        arr_E_x = Array(1.:1.:PML_Thickness[1])
        #arr_E_x = Array(PML_Thickness[1]:-1.:1)
        σ_E_x[:, 1] = σ_profile(arr_E_x, PML_Thickness[1], g.Δx)
        α_E_x[:, 1] = α_profile(arr_E_x .- 1, PML_Thickness[1])
        κ_E_x[:, 1] = κ_profile(arr_E_x, PML_Thickness[1])
        b_E_x[:, 1] = b_coeff(σ_E_x[:, 1], κ_E_x[:, 1], α_E_x[:, 1], g)
        c_E_x[2:end-1, 1] = c_coeff(b_E_x[2:end-1, 1], σ_E_x[2:end-1, 1], κ_E_x[2:end-1, 1], α_E_x[2:end-1, 1])

        arr_H_x = Array(0.5:1.:PML_Thickness[1]-1.5)
        #arr_H_x = Array(PML_Thickness[1]-1.5:-1.:0.5)
        σ_H_x[:, 1] = σ_profile(arr_H_x, PML_Thickness[1], g.Δx)
        α_H_x[:, 1] = α_profile(arr_H_x, PML_Thickness[1])
        κ_H_x[:, 1] = κ_profile(arr_H_x, PML_Thickness[1])
        b_H_x[:, 1] = b_coeff(σ_H_x[:, 1], κ_H_x[:, 1], α_H_x[:, 1], g)
        c_H_x[:, 1] = c_coeff(b_H_x[:,1], σ_H_x[:,1], κ_H_x[:,1], α_H_x[:,1])

        # Top X 
        σ_E_x[:, 2] = σ_E_x[:, 1]
        α_E_x[:, 2] = α_E_x[:, 1]
        κ_E_x[:, 2] = κ_E_x[:, 1]
        b_E_x[:, 2] = b_E_x[:, 1]
        c_E_x[:, 2] = c_E_x[:, 1]

        σ_H_x[:, 2] = σ_H_x[:, 1]
        α_H_x[:, 2] = α_H_x[:, 1]
        κ_H_x[:, 2] = κ_H_x[:, 1]
        b_H_x[:, 2] = b_H_x[:, 1]
        c_H_x[:, 2] = c_H_x[:, 1]


         denominator_E_x, denominator_H_x = [zeros(Float64, g.SizeX-1) for _ in 1:2]

        ii_E = PML_Thickness[1]
        @inbounds for mm in 1:g.SizeX-1
            if mm <= PML_Thickness[1]
                denominator_E_x[mm] = 1/(κ_E_x[mm, 1] * g.Δx)
            elseif mm >= g.SizeX + 1 - PML_Thickness[1]
                denominator_E_x[mm] = 1/(κ_E_x[ii_E, 2] * g.Δx) 
                ii_E = ii_E - 1
            else 
                denominator_E_x[mm] = 1/g.Δx
            end   
        end

        ii_H = PML_Thickness[1] - 1
        @inbounds for mm in 1:g.SizeX-1
            if mm <= PML_Thickness[1]-1
                denominator_H_x[mm] = 1/(κ_H_x[mm, 1] * g.Δx) 
            elseif mm >= g.SizeX + 1 - PML_Thickness[1]
                denominator_H_x[mm] = 1/(κ_H_x[ii_H, 2] * g.Δx) 
                ii_H = ii_H - 1
            else 
                denominator_H_x[mm] = 1/g.Δx
            end   
        end


        new(PML_Thickness,
            α_E_x,
            κ_E_x,
            σ_E_x,

            b_E_x,  

            c_E_x,

            α_H_x,
            κ_H_x,
            σ_H_x,

            b_H_x,

            c_H_x,

            denominator_E_x, 
            denominator_H_x, 
        )
    end
end