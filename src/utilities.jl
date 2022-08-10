function linear_predictor(arr_old::Array{Float64}, arr_current::Array{Float64})
    return 3/2 .* arr_old .- 1/2 .* arr_current
end

function polarization_rotation_matrix_e(Ψ::Float64, ϕ::Float64, θ::Float64)
    px = cos(Ψ)sin(ϕ) - sin(Ψ)cos(θ)cos(ϕ)
    py = -cos(Ψ)cos(ϕ) - sin(Ψ)cos(θ)sin(ϕ)
    pz = sin(Ψ)sin(θ)
    return [px, py, pz]
end

function polarization_rotation_matrix_h(Ψ::Float64, ϕ::Float64, θ::Float64)
    px = sin(Ψ)sin(ϕ) - cos(Ψ)cos(θ)cos(ϕ)
    py = -sin(Ψ)cos(ϕ) - cos(Ψ)cos(θ)sin(ϕ)
    pz = -cos(Ψ)sin(θ)
    P = vcat(px, py, pz)
    return P
end

function map_ijk_to_ir(xyz::Array{Int64, 1}, mnp::Array{Int64, 1})
    i_r =ceil(Int64, mnp[1]*xyz[1] + mnp[2]*xyz[2] + mnp[3]*xyz[3] )
    return i_r
end

function azimuthal_angle_2p(px::Int64, py::Int64)
    # notice that p and m are in general the same measurents for a uniform grid, notice that in this case Δr and Δxyz are equivalent
    ϕ = atan2(py/px)
    return ϕ
end

function polar_angle_3p(px::Int64, py::Int64, pz::Int64)
    Θ = atan2(sqrt(px^2 + py^2)/pz)    
    return Θ
end

function E_reflection_for_χ(E::Array{Float64, 1}, χ::Float64, order::Float64, n_r::Float64)
    ## if order bigger than 3 (i.e. effective nonlinearity), one should use the perturbation approximation (χ via E_{atom}) 
    Pz =  χ .* (E+0im).^order 
    return Pz./(2*n_r*(1+ n_r))
end

function χ_injection_stat(Γ̂::Float64, Ê::Float64, E_gap::Float64, n_0::Float64, eff_nl::Float64)
    return E_gap * q_0 * n_0 * Γ̂ / Ê^(eff_nl - 1)
end

function χ_brunel_stat(Γ̂::Float64, Ê::Float64, n_0::Float64, eff_nl::Float64)
    return n_0 * q_0^2 * Γ̂ / (m_e * Ê^(eff_nl))
end


function n_bound(γ, ω_0, χ_1, ω)
    ϵ_complex = epsilon_complex(γ, ω_0, χ_1, ω)
    ϵ_real = real.(ϵ_complex)
    ϵ_imag = imag.(ϵ_complex)
    n_r =  zeros(Float64, size(ω)[1])
    n_i =  zeros(Float64, size(ω)[1])
    
    @. n_r = sqrt(1/2 *(sqrt(ϵ_real^2 + ϵ_imag^2) + ϵ_real))
    @. n_i = sqrt(1/2 *(sqrt(ϵ_real^2 + ϵ_imag^2) - ϵ_real))
    return  [n_r, n_i]
end

function epsilon_complex(γ, ω_0, χ_1, ω)
    eps_complex = zeros(ComplexF64, length(ω)) .+ 1.
    for idx = 1:length(γ)
        @. eps_complex += χ_1[idx] * (ω_0[idx]^2) / (ω_0[idx]^2 - ω^2 + 1im * γ[idx] * ω)
    end
    return eps_complex
end

