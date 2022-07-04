using SpecialFunctions

function Γ_ADK(E::Vector{Float64}, Gamma_arr::Vector{Float64}, I_p::Float64, Z::Int64, l::Int64, m::Int64)
    I_p_au=I_p/q_0/27.21139609 #au in eV
    κ=sqrt(2*I_p_au)
    ns=Z/κ
    ls=ns-1
    F=abs.(E .+ 1e-6)./5.14220652e11    #E field in au

    A_lm=((2*l+1)*factorial(l+abs(m)))/(2^(abs(m))*factorial(abs(m))*factorial(l-abs(m)))

    C_nsls_abs_squared=2^(2*ns)/(ns*gamma(ns+ls+1)*gamma(ns-ls))
    for mm in 1:length(E)
        Gamma_arr[mm] = C_nsls_abs_squared*A_lm*I_p_au*(2*κ^3/F[mm])^(2*ns-abs(m)-1)*exp(-(2*κ^3/(3*F[mm])))
    end
    Γ=Gamma_arr./24.18884336e-18
    Γ[isnan.(Γ)].=0
    Γ[isinf.(Γ)].=0
    return Γ
end

function Γ_ADK(E_max::Float64, I_p::Float64; Z::Int64=1, l::Int64=0, m::Int64=0)
    I_p_au=I_p/q_0/27.21139609 #au in eV
    κ=sqrt(2*I_p_au)
    ns=Z/κ
    ls=ns-1
    F=abs(E_max .+ 1e-6)/5.14220652e11    #E field in au

    A_lm=((2*l+1)*factorial(l+abs(m)))/(2^(abs(m))*factorial(abs(m))*factorial(l-abs(m)))

    C_nsls_abs_squared=2^(2*ns)/(ns*gamma(ns+ls+1)*gamma(ns-ls))
    Γ_au = C_nsls_abs_squared*A_lm*I_p_au*(2*κ^3/F)^(2*ns-abs(m)-1)*exp(-(2*κ^3/(3*F)))
    Γ=Γ_au/24.18884336e-18
    return Γ
end

function Γ_Tangent(E::Vector{Float64}, Γ::Vector{Float64}, a::Float64, Γ̂::Float64, Ê::Float64)
    for mm in 1:length(E)
        Γ[mm] = Γ̂ * abs(E[mm]/Ê)^a ;
    end
    Γ[isnan.(Γ)].=0
    Γ[isinf.(Γ)].=0
    return Γ
end

function effective_nonlinearity_m(Ê::Float64, I_p::Float64)
    h = 1e-3;
    Γ̂ = Γ_ADK(Ê, I_p)
    m = Ê/Γ̂ * (Γ_ADK(Ê+h, I_p) - Γ_ADK(Ê, I_p))./h ;
    return m 
end

function ω_plasma(ρ::Float64)
    return abs(q_0^2 * ρ / (ϵ_0*m_e))^(1/2)
end

function linear_predictor(arr_old::Array{Float64}, arr_current::Array{Float64})
    return 3/2 .* arr_old .- 1/2 .* arr_current
end

function multinomial_degen(m::Integer, n_pump::Integer, n_probe::Integer)
    degeneracy = factorial(m)/(factorial(n_probe) * factorial((m + n_pump - n_probe)/2)* factorial((m - n_pump - n_probe)/2))
    return degeneracy
end

function multinomial_degen(m::Float64, n_pump::Float64, n_probe::Float64)
    degeneracy = gamma(m + 1)/(gamma(n_probe + 1) * gamma((m + n_pump - n_probe)/2 + 1)* gamma((m - n_pump - n_probe)/2 + 1))
    return degeneracy
end