struct parameters
#====basic parameters====#
    debug_prior     ::Int64     # 0: normal inversion; 1: run without input data.
    plot_voronoi    ::Int64     # plot each voronoi diagram in the generated models or not.
    add_yVec        ::Int64     # 0: only work in 2D: x- and z- axis; 1: 3D tomography.

#====Voronoi diagram parameters====#
    sig             ::Int64     # in percent of scales for all parameters.
    zeta_scale      ::Int64     # std of the prior for Normal. +/- bounds for uniform. 
    max_cells       ::Int64     # don't go over 100
    min_cells       ::Int64     # minimum for a box shaped anomaly
    max_sig         ::Float64   # on t*
    interp_style    ::Int64     # nearest:1, Inverse Distance Weighting:2. Similar results for each
    enforce_discon  ::Int64     # load a discontinuity, and use it in the inversion
    prior           ::Int64     # Uniform:1, normal:2, or exponential:3. exponential not recommended
# check whether they are used
    event_statics   ::Int64     # much better if you use this if using relative measurements. 
    demean          ::Int64     # for plotting only.  
###
#=====Monte Carlo parameters=====#
    n_chains        ::Int64     # Per script
    n_iter          ::Float64   
    burn_in         ::Float64
    keep_each       ::Float64
    print_each      ::Float64

#=====map parameters=====#
    max_depth       ::Int64
    min_depth       ::Int64
#check whether it's used
    rotation        ::Int64
###
    ZnodeSpacing    ::Int64
    buffer          ::Int64
    XYnodeSpacing   ::Int64

#====Cross section parameters====#
    xyMap           ::Bool      # work on x-y cross section with a given z value
    zSlice          ::Int64     # the given z value (vertical axis)
    xzMap           ::Bool      # work on x-z cross section with a given y value
    ySlice          ::Int64     # the given y value (horizontal axis, E-W afer rotation of β)
    # yzMap           ::Bool      # work on y-z cross section with a given x value, currently not useful as it is along-strike cross section
    # xSlice          ::Int64     # the given x value (horizontal axis, N-S afer rotation of β)
end

function define_TDstructrure()

    TD_parameters = parameters(
        #====basic parameters====#
            0, 1, 1,
        #====Voronoi diagram parameters====#
            10, 50, 100, 5, 0.1, 1, 0, 1, 1, 1,
        #====Monte Carlo parameters====#
            2, 1e3, 5e2, 1e1, 1e2,
        #====map parameters====#
            660, 0, 20, 20, 100, 20,
        #====Cross section parameters====#
            true, 100, true, 150
        )
    
    return TD_parameters

end
