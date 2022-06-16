module MaxwellFDTD
    include("constants.jl")
    include("grids.jl")
    include("media.jl")
    include("fields.jl")
    include("coefficients.jl")
    include("updateEqs.jl")
    include("detectors.jl")
    include("abcs.jl")
    include("sources.jl")
    
    export c_0, μ_0, ϵ_0
    export GaussianPointSource, GaussianWavePointSource, SinusoidalPointSource, RickerPointSource, sourceE!, sourceH!, intensity2amplitude
    export PointDetector, LineDetector, safeE!, safeP!, safeJ!, safePNl!
    export Fields1D, MaterialFields1D, LorentzFields1D
    export Grid1D, Grid2D, Grid3D
    export GridCoefficients1D
    export updateE!, updateH!, updateJ!, updateP!, updatePNl!
    export ABC!, LeftSideMurABC, RightSideMurABC, saveFields!, stepABC!
    export LorentzMedium1D, StaticMedium1D
end