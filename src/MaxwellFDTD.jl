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
    include("visuals.jl")
    
    export c_0, μ_0, ϵ_0, q_0, m_e
    export GaussianPointSource, GaussianWavePointSource, SinusoidalPointSource, RickerPointSource, RickerPointSource3D, sourceE!, sourceH!, intensity2amplitude, laserfluence
    export Fields1D, MaterialFields1D, LorentzFields1D, DrudeFields1D, TunnelFields1D, Fields3D
    export PointDetector, LineDetector, ZSliceDetector, safeE!, safeP!, safeJ!, safePNl!, safeΓ_ADK!, safeJ_tunnel!, safeJ_bound!, safeJ_free!
    export Grid1D, Grid2D, Grid3D
    export GridCoefficients1D, GridCoefficients3D, FieldIonizationCoefficients1D
    export updateE!, updateH!, updateJ!, updateJbound!, updatePbound!, updatePNl!, updatePlasma!, updateJfree!, updateJtunnel!
    export ABC!, LeftSideMurABC, RightSideMurABC, saveFields!, stepABC!, PMLXlow, PMLXhigh, PMLYlow, PMLYhigh, PMLZlow, PMLZhigh, update_ϕ_E!, update_ϕ_H!, updateE!, updateH!
    export LorentzMedium1D, StaticMedium1D, TunnelMedium1D, DrudeMedium1D, StaticMedium3D
    export mytukey
    export Γ_ADK, ω_plasma
    export slide_arr_over_time
end