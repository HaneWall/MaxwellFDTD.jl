using FFTW
using DSP
using GLMakie
GLMakie.activate!()

function slide_arr_over_time(arr::Array{Float64, 2})
    # norm_array 
    n_arr = arr./maximum(arr)
    t_index = Observable(1)
    t_slice = @lift(n_arr[$t_index, :])
    f = Figure()
    ax1 = Axis(f[1, 1])
    lines!(f[1,1], t_slice, color=:black)
    ylims!(ax1, -1, 1)
    sl = Slider(f[1, 2], horizontal = false, range = 1:size(arr, 1))
    connect!(t_index, sl.value)
    f
end

function slide_arr_over_time(arr::Array{Float64, 3})
    # norm_array 
    n_arr = arr./maximum(arr)
    t_index = Observable(1)
    t_slice = @lift(n_arr[$t_index, :, :])
    f = Figure()
    x = 1:1:size(arr, 2)
    y = 1:1:size(arr, 3)
    heatmap(f[1,1], x, y, t_slice, colormap=:turbo)
    sl = Slider(f[1, 2], horizontal = false, range = 1:size(arr, 1))
    connect!(t_index, sl.value)
    f
end


function slide_arr_over_time(arr::Array{Float64, 4})
    f = Figure()
    int_sl = IntervalSlider(
            f[3, 1],
            range = LinRange(0:0.01:maximum(arr)),
            startvalues = (0, maximum(arr)),
            halign = :center,
            linewidth = 10.0
            )

    # norm_array 
    t_index = Observable(1)
    t_slice = @lift(abs.(arr[$t_index, :, :, :]))
    x = 1:1:size(arr, 2)
    y = 1:1:size(arr, 3)
    z = 1:1:size(arr, 4)
    #limits_ref = tuple(@lift(minimum(arr[$t_index, :, :, :])), @lift(maximum(arr[$t_index, :, :, :])))
    volume(f[1,1], x, y, z, t_slice; colormap=:turbo, transparency=true, colorrange=int_sl.interval)
    Colorbar(f[2, 1], colormap = :turbo, limits=int_sl.interval, vertical = false)
    sl = Slider(f[1, 2], horizontal = false, range = 1:size(arr, 1))
    connect!(t_index, sl.value)
    #connect!(interval, int_sl.interval)
    f
end

function plot_amplitude_spectrum(t::Array{Float64,1}, signal::Array{Float64,1})
    L = length(t)
    δt = abs(t[2]-t[1])
    δf = 1/δt
    ω = 2*π*fftshift(fftfreq(L, δf))
    FT = fftshift(fft(signal))
    fig = Figure(resolution = (800, 800))
    ax = Axis(fig[1, 1], title="Amplitude Spectrum")
    lines!(ax, ω, abs.(FT)./L)
    fig
end

function plot_amplitude_spectrum(t::Array{Float64,1}, signal::Array{Float64,1}, padding::Bool)
    L = length(t)
    padded_L = nextpow(2, L)
    center_help_L_low = floor(Int64, (padded_L - L)/2)
    center_help_L_high = floor(Int64, L + center_help_L_low - 1) 
    δt = abs(t[2]-t[1])
    sigpad = zeros(eltype(signal), padded_L)
    window = blackman(L)
    sigpad[center_help_L_low:center_help_L_high] = window.*signal[:]
    δf = 1/δt
    ω = 2*π*fftshift(fftfreq(padded_L, δf))
    FT = fftshift(fft(sigpad))
    fig = Figure(resolution = (800, 800))
    ax = Axis(fig[1, 1], title="Amplitude Spectrum")
    lines!(ax, ω, abs.(FT)./padded_L)
    fig
end

function plot_log10_amplitude_spectrum(t::Array{Float64,1}, signal::Array{Float64,1}, padding::Bool)
    L = length(t)
    padded_L = nextpow(2, L)
    center_help_L_low = floor(Int64, (padded_L - L)/2)
    center_help_L_high = floor(Int64, L + center_help_L_low - 1) 
    δt = abs(t[2]-t[1])
    sigpad = zeros(eltype(signal), padded_L)
    window = blackman(L)
    sigpad[center_help_L_low:center_help_L_high] = window.*signal[:]
    δf = 1/δt
    ω = 2*π*fftshift(fftfreq(padded_L, δf))
    FT = fftshift(fft(sigpad))
    fig = Figure(resolution = (800, 800))
    ax = Axis(fig[1, 1], title="Amplitude Spectrum")
    lines!(ax, ω, log10.(abs.(FT)./padded_L))
    fig
end

function plot_log10_amplitude_spectrum(t::Array{Float64,1}, signal::Matrix, padding::Bool)
    fig = Figure(resolution = (800, 800), font = "CMU Serif")
    ax = Axis(fig[1, 1], title="Amplitude Spectrum")
    L = length(t)
    padded_L = nextpow(2, L)
    center_help_L_low = floor(Int64, (padded_L - L)/2)
    center_help_L_high = floor(Int64, L + center_help_L_low - 1) 
    δt = abs(t[2]-t[1])
    δf = 1/δt
    ω = 2*π*fftshift(fftfreq(padded_L, δf))
    window = blackman(L)
    for idx in 1:size(signal)[2]
        sigpad = zeros(eltype(signal[:,idx]), padded_L)
        sigpad[center_help_L_low:center_help_L_high] = window.*signal[:, idx]
        FT = fftshift(fft(sigpad))
        lines!(ax, ω, log10.(abs.(FT)./padded_L))
    end
    fig
end

function plot_log10_power_spectrum(t::Array{Float64,1}, signal::Matrix, ω_central::Float64, xlims::Vector{Float64}, ylims::Vector{Float64}, legend::Vector{String}, padding::Bool, name::String)
    fig = Figure(resolution = (800, 400), font = "CMU Serif")
    ax = Axis(fig[1, 1], xlabel=L"$\omega$ / $ \omega_{Pump}$", ylabel=L"\log_{10}|FFT(S(t)|^2")
    L = length(t)
    padded_L = nextpow(2, L)
    center_help_L_low = floor(Int64, (padded_L - L)/2) 
    center_help_L_high = floor(Int64, L + center_help_L_low - 1) 
    δt = abs(t[2]-t[1])
    δf = 1/δt
    ω = 2*π*fftshift(fftfreq(padded_L, δf))./ω_central
    window = blackman(L)
    for idx in 1:size(signal)[2]
        sigpad = zeros(eltype(signal[:,idx]), padded_L)
        sigpad[center_help_L_low:center_help_L_high] = window.*signal[:, idx]
        FT = fftshift(fft(sigpad))
        lines!(ax, ω, log10.((abs.(FT)./padded_L).^2), label=L"%$(legend[idx])")
    end
    axislegend(ax, position=:rt, orientation=:horizontal)
    xlims!(ax, xlims[1], xlims[2])
    ylims!(ax, ylims[1], ylims[2])
    save(name, fig)
end

function permutation_plot(t::Array{Float64,1}, signal::Matrix, ω_central::Float64, ω_harm::Float64, xlims::Vector{Float64}, ylims::Vector{Float64}, legend::Vector{String}, padding::Bool, name::String)
    fig = Figure(resolution = (800, 400), font = "CMU Serif")
    ax = Axis(fig[1, 1], xlabel=L"$\omega$ / $ \omega_{Pump}$", ylabel=L"\log_{10}|FFT(S(t)|^2")
    L = length(t)
    padded_L = nextpow(2, L)
    center_help_L_low = floor(Int64, (padded_L - L)/2) 
    center_help_L_high = floor(Int64, L + center_help_L_low - 1) 
    δt = abs(t[2]-t[1])
    δf = 1/δt
    ω = 2*π*fftshift(fftfreq(padded_L, δf))./ω_central
    idx_first_harm = argmin(abs.(ω - ω_harm))
    window = blackman(L)
    for idx in 1:size(signal)[2]
        sigpad = zeros(eltype(signal[:,idx]), padded_L)
        sigpad[center_help_L_low:center_help_L_high] = window.*signal[:, idx]
        FT = fftshift(fft(sigpad))
        lines!(ax, ω, log10.((abs.(FT)./padded_L).^2)./log10((abs.(FT[idx_first_harm])./padded_L).^2), label=L"%$(legend[idx])")
    end
    axislegend(ax, position=:rt, orientation=:horizontal)
    xlims!(ax, xlims[1], xlims[2])
    ylims!(ax, ylims[1], ylims[2])
    save(name, fig)
end