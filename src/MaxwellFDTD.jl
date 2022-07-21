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
    export Fields1D, Fields2D, MaterialFields1D, LorentzFields1D, DrudeFields1D, TunnelFields1D, Fields3D, MaterialFields2D, TunnelFields2D, DrudeFields2D, LorentzFields2D
    export PointDetector, LineDetector, ZSliceDetector, YSliceDetector, XSliceDetector, BlockDetector, PlaneDetector, safeE!, safeP!, safeJ!, safePNl!, safeΓ_ADK!, safeJ_tunnel!, safeJ_bound!, safeJ_free!
    export Grid1D, Grid2D, Grid3D
    export CPML_Ψ_Fields_3D, CPML_Ψ_Fields_2D, CPML_Ψ_Fields_1D, CPML_Parameters_3D, CPML_Parameters_2D, CPML_Parameters_1D
    export GridCoefficients1D, GridCoefficients3D,  GridCoefficients1D_w_CPML, GridCoefficients2D_w_CPML, GridCoefficients3D_WIP, GridCoefficients3D_w_CPML, FieldIonizationCoefficients1D, FieldIonizationCoefficients2D
    export updateE!, updateH!, updateEWIP!, updateHWIP!, updateJ!, updateJbound!, updatePbound!, updatePNl!, updatePlasma!, updatePlasmaTangent!, updateJfree!, updateJtunnel!, update_Ψ_E!, update_Ψ_H!, apply_Ψ_E!, apply_Ψ_H!
    export ABC!, LeftSideMurABC, RightSideMurABC, saveFields!, stepABC!
    export LorentzMedium1D, StaticMedium1D, TunnelMedium1D, DrudeMedium1D, StaticMedium2D, TunnelFields2D, LorentzMedium2D, DrudeMedium2D, StaticMedium3D
    export mytukey
    export Γ_ADK, Γ_Tangent, ω_plasma, linear_predictor, multinomial_degen, effective_nonlinearity_m
    export record_arr_over_time, slide_arr_over_time, plot_amplitude_spectrum, plot_log10_amplitude_spectrum, plot_log10_power_spectrum, plot_Reflection_spectrum, permutation_plot, timeseries_plot
end