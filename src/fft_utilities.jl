using DSP

function mytukey(n::Integer, a::Real, padding::Integer=0, zerophase::Bool=false)
    return DSP.tukey(n, a; padding, zerophase)
end

