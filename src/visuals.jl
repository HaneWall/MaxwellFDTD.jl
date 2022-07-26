using FFTW
using DSP
using GLMakie
using LaTeXStrings
GLMakie.activate!()


function slide_arr_over_time(arr::Array{Float64, 1})
    # norm_array 
    n_arr = arr./maximum(arr)
    t_index = Observable(1)
    t_slice = @lift(n_arr[$t_index])
    f = Figure()
    ax1 = Axis(f[1, 1])
    scatter!(f[1,1], 1, t_slice, color=:black)
    ylims!(ax1, -1, 1)
    sl = Slider(f[1, 2], horizontal = false, range = 1:size(arr, 1))
    connect!(t_index, sl.value)
    f
end

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
    n_arr = arr
    t_index = Observable(1)
    t_slice = @lift(abs.(n_arr[$t_index, :, :]))
    f = Figure()
    int_sl = IntervalSlider(
            f[3, 1],
            range = LinRange(0:0.01:maximum(arr)),
            startvalues = (0, maximum(arr)),
            halign = :center,
            linewidth = 10.0
            )
    x = 1:1:size(arr, 2)
    y = 1:1:size(arr, 3)
    heatmap(f[1,1], x, y, t_slice, colormap=:turbo, colorrange=int_sl.interval)
    Colorbar(f[2, 1], colormap = :turbo, limits=int_sl.interval, vertical = false)
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

function record_arr_over_time(arr::Array{Float64, 3}, filename::String)
    f = Figure(resolution = (1000, 1000))
    ax = Axis(f[1, 1])
    t_index = Observable(1)
    t_slice = @lift(abs.(arr[$t_index, :, :]))
    clims = @lift 1.1 .* extrema(abs.(arr[$t_index, :, :]))
    frames = 1:size(arr, 1)
    x = 1:1:size(arr, 2)
    y = 1:1:size(arr, 3)
    vol = heatmap!(ax, x, y, t_slice; colormap=:turbo, transparency=true, colorrange=clims)
    #vol = volume!(ax, x, y, z, t_slice; colormap=:turbo, transparency=true)
    #Colorbar(f[1, 2], vol)
    record(f, filename * ".mp4", frames, framerate=8) do i
        msg = string("Plotting frame ", i, " of ", frames[end])
        print(msg * " \r")
        t_index[] = i
    end
end

function record_arr_over_time(arr::Array{Float64, 4}, filename::String)
    f = Figure(resolution = (1000, 1000))
    ax = Axis3(f[1,1], perspectiveness=0.5)
    t_index = Observable(1)
    t_slice = @lift(abs.(arr[$t_index, :, :, :]))
    clims = @lift 1.1 .* extrema(abs.(arr[$t_index, :, :, :]))
    frames = 1:size(arr, 1)
    x = 1:1:size(arr, 2)
    y = 1:1:size(arr, 3)
    z = 1:1:size(arr, 4)
    vol = volume!(ax, x, y, z, t_slice; colormap=:turbo, transparency=true, colorrange=clims)
    #vol = volume!(ax, x, y, z, t_slice; colormap=:turbo, transparency=true)
    #Colorbar(f[1, 2], vol)
    record(f, filename * ".mp4", frames, framerate=8) do i
        msg = string("Plotting frame ", i, " of ", frames[end])
        print(msg * " \r")
        t_index[] = i
    end
end

function timeseries_plot(arr::Array{Float64,1})
    timesteps = length(arr)
    fig = Figure(resolution = (800, 1200))
    ax = Axis(fig[1, 1], title="timeseries")
    lines!(ax, 1:1:timesteps, log10.(abs.(arr)./maximum(arr)))
    fig
end

function timeseries_plot(arr::Array{Float64,1}, legend::String, name::String)
    timesteps = length(arr)
    fig = Figure(resolution = (800, 800))
    ax = Axis(fig[1, 1], title="Log10-Timeseries")
    lines!(ax, 1:1:timesteps, log10.(abs.(arr)./maximum(arr)), label = L"%$legend")
    ax2 = Axis(fig[2, 1], title="Timeseries")
    lines!(ax2, 1:1:timesteps, arr, label = L"%$legend")
    axislegend(ax, position=:rb)
    save(name, fig)
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

function plot_amplitude_spectrum(t::Array{Float64,1}, signal::Matrix, padding::Bool)
    fig = Figure(resolution = (800, 800))
    ax = Axis(fig[1, 1], title="Amplitude Spectrum")
    L = length(t)
    padded_L = nextpow(2, L)
    center_help_L_low = floor(Int64, (padded_L - L)/2)
    center_help_L_high = floor(Int64, L + center_help_L_low - 1) 
    δt = abs(t[2]-t[1])
    window = blackman(L)
    δf = 1/δt
    ω = 2*π*fftshift(fftfreq(padded_L, δf))
    for idx in 1:size(signal)[2]
        sigpad = zeros(eltype(signal[:, idx]), padded_L)
        sigpad[center_help_L_low:center_help_L_high] = window.*signal[:, idx]
        FT = fftshift(fft(sigpad))
        lines!(ax, ω, abs.(FT)./padded_L)
    end
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
    fig = Figure(resolution = (400, 400), font = "CMU Serif")
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
    axislegend(ax, position=:rt)
    xlims!(ax, xlims[1], xlims[2])
    ylims!(ax, ylims[1], ylims[2])
    save(name, fig)
end

function plot_log10_power_spectrum_current_and_E(t::Array{Float64,1}, signal::Matrix, ω_central::Float64, xlims::Vector{Float64}, ylims::Vector{Float64}, legend::Vector{String}, padding::Bool, name::String)
    fig = Figure(resolution = (400, 400), font = "CMU Serif")
    ax = Axis(fig[1, 1], xlabel=L"$\omega$ / $ \omega_{Pump}$", ylabel=L"\log_{10}|FFT(S(t)|^2")
    L = length(t)
    padded_L = nextpow(2, L)
    center_help_L_low = floor(Int64, (padded_L - L)/2) 
    center_help_L_high = floor(Int64, L + center_help_L_low - 1) 
    δt = abs(t[2]-t[1])
    δf = 1/δt
    ω = 2*π*fftshift(fftfreq(padded_L, δf))./ω_central
    window = blackman(L)
    max_arr = 0.
    for idx in 1:size(signal)[2]
        sigpad = zeros(eltype(signal[:,idx]), padded_L)
        sigpad[center_help_L_low:center_help_L_high] = window.*signal[:, idx]
        FT = fftshift(fft(sigpad))
        if idx == 1
            jz = log10.(((abs.(FT)./padded_L).^2)) 
            max_arr = maximum(jz)
            lines!(ax, ω, jz./max_arr, label=L"%$(legend[idx])")
        elseif idx == 2
            res = log10.(((abs.(ω .* FT)./padded_L).^2))
            shift = 1- maximum(res./max_arr)
            lines!(ax, ω,  res./max_arr .+ shift, label=L"%$(legend[idx])", linestyle = :dash)
        else
            j = log10.(((abs.(FT)./padded_L).^2)) 
            lines!(ax, ω, j./max_arr, label=L"%$(legend[idx])")
        end
    end
    axislegend(ax, position=:rt)
    xlims!(ax, xlims[1], xlims[2])
    ylims!(ax, ylims[1], ylims[2])
    save(name, fig)
end


function plot_log10_power_spectrum(t::Array{Float64,1}, signal::Array{Float64,1}, ω_central::Float64, xlims::Vector{Float64}, ylims::Vector{Float64}, legend::String, padding::Bool, name::String)
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
    
    sigpad = zeros(eltype(signal[:]), padded_L)
    sigpad[center_help_L_low:center_help_L_high] = window.*signal[:]
    FT = fftshift(fft(sigpad))
    lines!(ax, ω, log10.((abs.(FT)./padded_L).^2), label=L"%$(legend)")
    
    axislegend(ax, position=:rt, orientation=:horizontal)
    xlims!(ax, xlims[1], xlims[2])
    ylims!(ax, ylims[1], ylims[2])
    save(name, fig)
end


function plot_Reflection_spectrum(t::Array{Float64,1}, signal::Matrix, ω_central::Float64, xlims::Vector{Float64}, ylims::Vector{Float64}, legend::LaTeXString, padding::Bool, name::String)
    fig = Figure(resolution = (800, 800), font = "CMU Serif")
    ax = Axis(fig[1, 1], xlabel=L"$\omega$ / $ \omega_{Pump}$", ylabel=L"Ratio")
    ax2 = Axis(fig[2, 1], xlabel=L"$\omega$ / $ \omega_{Pump}$", ylabel=L"Ratio")
    L = length(t)
    padded_L = nextpow(2, L)
    center_help_L_low = floor(Int64, (padded_L - L)/2) 
    center_help_L_high = floor(Int64, L + center_help_L_low - 1) 
    δt = abs(t[2]-t[1])
    δf = 1/δt
    ω = 2*π*fftshift(fftfreq(padded_L, δf))./ω_central
    window = blackman(L)
    idx_omega = argmin(abs.(ω .- 1)) - ceil(Int64,length(ω)/2)
    idx_omega_set = [idx_omega + 2*i*idx_omega for i in 0:10] .+ ceil(Int64,length(ω)/2)
    sigpad_1 = zeros(eltype(signal[:, 1]), padded_L)
    sigpad_1[center_help_L_low:center_help_L_high] = window.*signal[:, 1]
    sigpad_2 = zeros(eltype(signal[:, 2]), padded_L)
    sigpad_2[center_help_L_low:center_help_L_high] = window.*signal[:, 2]
    FT_1 = fftshift(fft(sigpad_1[:]))
    FT_2 = fftshift(fft(sigpad_2[:]))
    ratio = ((abs.(FT_1)./padded_L).^2)./((abs.(FT_2)./padded_L).^2)
    scatter_refl = [ratio[idx] for idx in idx_omega_set]
    lines!(ax, ω, ratio, label=L"%$(legend)")
    vlines!(ax, [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21], color = :red, linestyle=:dash)
    axislegend(ax, position=:rt, orientation=:horizontal)
    xlims!(ax, xlims[1], xlims[2])
    ylims!(ax, ylims[1], ylims[2])
    scatter!(ax2, [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21], scatter_refl, color=:red)
    xlims!(ax2, xlims[1], xlims[2])
    ylims!(ax2, 1.65, 2.1)
    save(name, fig)
end

function permutation_plot(t::Array{Float64,1}, E_gap::Float64, signal::Matrix, ω_central::Float64, ω_probe::Float64, xlims::Vector{Float64}, ylims::Vector{Float64}, legend::Vector{String}, padding::Bool, name::String)
    fig = Figure(resolution = (400, 400), font = "CMU Serif")
    ax = Axis(fig[1, 1], xlabel=L"$\omega$ / $ \omega_{Pump}$", ylabel=L"\log_{10}|FFT(S(t)|^2")
    L = length(t)
    padded_L = nextpow(2, L)
    center_help_L_low = floor(Int64, (padded_L - L)/2) 
    center_help_L_high = floor(Int64, L + center_help_L_low - 1) 
    δt = abs(t[2]-t[1])
    δf = 1/δt
    ω = 2*π*fftshift(fftfreq(padded_L, δf))./ω_central
    ω_harm = [n*2*ω_central + ω_probe for n in 1:10]./ω_central
    idx_first_harm = argmin(abs.(ω .- ω_harm[1]))
    window = blackman(L)
    for idx in 1:size(signal)[2]
        sigpad = zeros(eltype(signal[:,idx]), padded_L)
        sigpad[center_help_L_low:center_help_L_high] = window.*signal[:, idx]
        FT = fftshift(fft(sigpad))
        int = abs.(FT).^2
        harm_int = int[idx_first_harm]
        lines!(ax, ω, log10.(int./harm_int), label=L"%$(legend[idx])")
        
        #m_brunel = effective_nonlinearity_m(maximum(signal[:,idx]), E_gap * q_0) + 2
        m_brunel = 12.54 + 2
        #m_injection = effective_nonlinearity_m(maximum(signal[:,idx]), E_gap * q_0) 
        m_injection = 12.54
        
        n_brunel = Array(1.:1.:((ceil(Int64, m_brunel-1))/2))
        n_injection = Array(1.:1.:((ceil(Int64, m_injection-1))/2))
        permutations_brunel = zeros(Float64, length(n_brunel))
        permutations_injection = zeros(Float64, length(n_injection))
        n_pump_injection = 2*n_injection
        n_pump_brunel = 2*n_brunel
        n_probe = 1.

        for j in 1:length(n_injection)
            permutations_injection[j] = multinomial_degen(m_injection, n_pump_injection[j], n_probe)
        end

        for j in 1:length(n_brunel)
            permutations_brunel[j] = multinomial_degen(m_brunel, n_pump_brunel[j], n_probe)
        end

        scatter!(ax, n_pump_brunel .+ ω_probe/ω_central, log10.(permutations_brunel.^2 ./permutations_brunel[1]^2), color=:black, marker=:+, label="Brunel-Bi")
        scatter!(ax, n_pump_injection .+ ω_probe/ω_central, log10.(permutations_injection.^2 ./permutations_injection[1]^2), color=:red, marker=:+, label="Injection-Bi")
    end
    axislegend(ax, position=:rt)
    xlims!(ax, xlims[1], xlims[2])
    ylims!(ax, ylims[1], ylims[2])
    save(name, fig)
end