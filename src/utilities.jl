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
    ϕ = atan2(sqrt(px^2 + py^2)/pz)    
    return ϕ
end

function reflection_for_χ(θ::Float64, χ::Float64, order::Float64, n_r::Float64)
    
end