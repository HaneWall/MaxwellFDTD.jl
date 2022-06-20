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
    include("fft_utilities.jl")
    include("tunnelionisation_utilities.jl")
    
    export c_0, μ_0, ϵ_0, q_0, m_e
    export GaussianPointSource, GaussianWavePointSource, SinusoidalPointSource, RickerPointSource, sourceE!, sourceH!, intensity2amplitude, laserfluence
    export PointDetector, LineDetector, safeE!, safeP!, safeJ!, safePNl!
    export Fields1D, MaterialFields1D, LorentzFields1D, DrudeFields1D, TunnelFields1D
    export Grid1D, Grid2D, Grid3D
    export GridCoefficients1D
    export updateE!, updateH!, updateJbound!, updatePbound!, updatePNl!, updatePlasma!, updateJfree!, updateJtunnel!
    export ABC!, LeftSideMurABC, RightSideMurABC, saveFields!, stepABC!
    export LorentzMedium1D, StaticMedium1D, TunnelMedium1D, DrudeMedium1D
    export mytukey
    export Γ_ADK, ω_plasma
end