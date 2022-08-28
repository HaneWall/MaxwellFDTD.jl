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
    include("utilities.jl")
    
    export c_0, μ_0, ϵ_0, q_0, m_e
    export GaussianPointSource, GaussianWavePointSource, GaussianWavePointSource2D, GaussianPlaneWaveSource2D_y, GaussianPlaneWaveSource3D_y,  GaussianWavePointSource3D, SinusoidalPointSource, RickerPointSource, RickerPointSource3D, sourceE!, sourceH!, intensity2amplitude, laserfluence
    export Fields1D, Fields2D, Fields3D, MaterialFields1D, LorentzFields1D, DrudeFields1D, TunnelFields1D, MaterialFields2D, TunnelFields2D, DrudeFields2D, LorentzFields2D, MaterialFields3D, LorentzFields3D, DrudeFields3D, TunnelFields3D
    export PointDetector, LineDetector, ZSliceDetector, YSliceDetector, XSliceDetector, BlockDetector, PlaneDetector, safeE!, safeP!, safeJ!, safePNl!, safeΓ_ADK!, safeJ_tunnel!, safeJ_bound!, safeJ_free!
    export Grid1D, Grid2D, Grid3D
    export CPML_Ψ_Fields_3D, CPML_Ψ_Fields_2D, CPML_Ψ_Fields_1D, CPML_Parameters_3D, CPML_Parameters_2D, CPML_Parameters_1D
    export GridCoefficients1D, GridCoefficients2D, GridCoefficients3D,  GridCoefficients1D_w_CPML, GridCoefficients2D_w_CPML, GridCoefficients3D_WIP, GridCoefficients3D_w_CPML, FieldIonizationCoefficients1D, FieldIonizationCoefficients2D, FieldIonizationCoefficients3D
    export updateE!, updateH!, updateEWIP!, updateHWIP!, updateJ!, updateJbound!, updatePbound!, updatePNl!, updatePlasma!, updatePlasmaTangent!, updateJfree!, updateJtunnel!, update_Ψ_E!, update_Ψ_H!, apply_Ψ_E!, apply_Ψ_H!
    export ABC!, LeftSideMurABC, RightSideMurABC, saveFields!, stepABC!, PeriodicBoundaryX, updateH!, updateE!
    export LorentzMedium1D, StaticMedium1D, TunnelMedium1D, DrudeMedium1D, StaticMedium2D, TunnelMedium2D, LorentzMedium2D, DrudeMedium2D, StaticMedium3D, TunnelMedium3D, LorentzMedium3D, DrudeMedium3D
    export mytukey, analytic_signal, isolate_around_freq
    export Γ_ADK, Γ_Tangent, ω_plasma, linear_predictor, multinomial_degen, effective_nonlinearity_m, binomial_degen
    export record_arr_over_time, slide_arr_over_time, plot_amplitude_spectrum, plot_log10_amplitude_spectrum, plot_log10_power_spectrum, plot_Reflection_spectrum, permutation_plot, timeseries_plot, plot_log10_power_spectrum_current_and_E
    export E_reflection_for_χ, χ_brunel_stat, χ_injection_stat, r_p, r_s, t_s, t_p, r_s_nl, n_bound, epsilon_complex, polarization_rotation_matrix_e, polarization_rotation_matrix_h
end