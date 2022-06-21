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
    t_slice = @lift(n_arr[$t_index, :])
    f = Figure()
    x = 1:1:size(arr, 2)
    y = 1:1:size(arr, 3)
    heatmap(f[1,1], x, y, t_slice, colormap=:turbo, colorrange=(-1, 1))
    sl = Slider(f[1, 2], horizontal = false, range = 1:size(arr, 1))
    connect!(t_index, sl.value)
    f
end

