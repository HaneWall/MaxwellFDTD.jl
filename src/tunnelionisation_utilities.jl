using SpecialFunctions

function Γ_ADK(I_p::Float64, E::Float64; Z::Integer=1, l::Integer=0, m::Integer=0)
    I_p_au=I_p/q_0./27.21139609 #au in eV
    κ=sqrt(2*I_p_au)
    ns=Z./κ
    ls=ns-1
    F=E/5.14220652e11;    #E field in au

    A_lm=((2*l+1)*factorial(l+abs(m)))/(2^(abs(m))*factorial(abs(m))*factorial(l-abs(m)))

    C_nsls_abs_squared=2^(2*ns)/(ns*gamma(ns+ls+1)*gamma(ns-ls))

    Gamma_au=C_nsls_abs_squared.*A_lm.*I_p_au*(2*κ^3/F)^(2*ns-abs(m)-1).*exp(-(2*κ^3/(3*F)));

    Γ=Gamma_au/24.18884336e-18
    Γ[isnan(Γ)]=0
    Γ[isinf(Gamma)]=0;
    return Γ
end

function ω_plasma(ρ::Float64)
    return abs((q_0^2 * ρ / ϵ_0* m_0))^(1/2)
end