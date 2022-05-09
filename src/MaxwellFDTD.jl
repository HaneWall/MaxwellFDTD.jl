module MaxwellFDTD
    include("constants.jl")
    include("grids.jl")
    include("fields.jl")
    include("media.jl")
    include("coefficients.jl")
    include("updateEqs.jl")
    include("detectors.jl")
    include("abcs.jl")
    include("sources.jl")
    
    export c_0, μ_0, ϵ_0
    export GaussianPointSource, sourceE!, sourceH!
    export PointDetector, LineDetector, safeE!
    export Fields1D, MaterialFields1D, LorentzFields1D
    export Grid1D, Grid2D, Grid3D
    export GridCoefficients1D
    export updateE!, updateH!
    export ABC!
    export LorentzMedium1D, StaticMedium1D
end