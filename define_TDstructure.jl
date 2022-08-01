function define_TDstructrure()

    TD_parameters = Dict(
        "debug_prior" => 1,     # decide whether the data will be used or not
        "xyMap" => false,       # work on x-y cross section with a given z value
        "zSlice" => [100],         # the given z value (vertical axis)
        "xzMap" => true,        # work on x-z cross section with a given y value
        "ySlice" => [150],         # the given y value (horizontal axis, E-W afer rotation of β)
        "yzMap" => false,       # work on y-z cross section with a given x value
        "xSlice" => [860],         # the given x value (horizontal axis, N-S afer rotation of β)
        "plot_voronoi" => 0,    # plot each voronoi diagram in the generated models or not
        "add_yVec" => 1,        # 0: only work in 2D: x- and z- axis; 1: 3D tomography
        "STdata3" => 0,         # use STdata3
        "sig" => 10,            # in percent of scales for all parameters.
        "zeta_scale" => 50,     # std of the prior for Normal. +/- bounds for uniform. 
         "max_cells" => 100,     # don't go over 100
         "min_cells" => 5,      # minimum for a box shaped anomaly
         "max_sig" => 0.1,     # on t*
         "interp_style" => 1,   # nearest:1, Inverse Distance Weighting:2. Similar results for each
         "enforce_discon" => 0, # load a discontinuity, and use it in the inversion
         "prior" => 1,          # Uniform:1, normal:2, or exponential:3. exponential not recommended
         "event_statics" =>1,   # much better if you use this if using relative measurements. 
         "demean" => 1,         # for plotting only.  

         "n_chains" => 2,       # Per script, not in total. See iter field in run_batch.pbs
         "n_iter" => 1e3,       # How many times you iterate each chain
         "burn_in" => 5e2,    # Don't save until you get to this iteration
        #  "keep_each" => 2.5e1,  # how often to save a model post burnin
         "keep_each" => 1e1,  # how often to save a model post burnin
         "print_each" => 1e2,   # how often to print to the screen if you are impatient and chekcing it (like I do)
         
        #  "n_chains" => 10,       # Per script, not in total. See iter field in run_batch.pbs
        #  "n_iter" => 1e4,       # How many times you iterate each chain
        #  "burn_in" => 5e3,    # Don't save until you get to this iteration
        #  "keep_each" => 1e2,  # how often to save a model post burnin
        # #  "keep_each" => 1e1,  # how often to save a model post burnin
        #  "print_each" => 1e3,   # how often to print to the screen if you are impatient and chekcing it (like I do)
 
    #     "n_chains" => 10,       # Per script, not in total. See iter field in run_batch.pbs
    #     "n_iter" => 1e5,       # How many times you iterate each chain
    #     "burn_in" => 5e4,    # Don't save until you get to this iteration
    #     "keep_each" => 1e2,  # how often to save a model post burnin
    #    #  "keep_each" => 1e1,  # how often to save a model post burnin
    #     "print_each" => 1e4,   # how often to print to the screen if you are impatient and chekcing it (like I do)
 

        #  "n_chains" => 10,       # Per script, not in total. See iter field in run_batch.pbs
        #  "n_iter" => 5e4,       # How many times you iterate each chain
        #  "burn_in" => 2.5e4,    # Don't save until you get to this iteration
        #  "keep_each" => 2.5e2,  # how often to save a model post burnin
        # #  "keep_each" => 1e1,  # how often to save a model post burnin
        #  "print_each" => 1e3,   # how often to print to the screen if you are impatient and chekcing it (like I do)
        
        # "n_chains" => 24,       # Per script, not in total. See iter field in run_batch.pbs
        # "n_iter" => 5e5,       # How many times you iterate each chain
        # "burn_in" => 2.5e5,    # Don't save until you get to this iteration
        # #  "keep_each" => 2.5e2,  # how often to save a model post burnin  
        # "keep_each" => 1e2,  # how often to save a model post burnin         
        # "print_each" => 1e4,   # how often to print to the screen if you are impatient and chekcing it (like I do)
        

         #### larger parameter
        #  "n_chains" => 10,       # Per script, not in total. See iter field in run_batch.pbs
        #  "n_iter" => 5e6,       # How many times you iterate each chain
        #  "burn_in" => 2.5e6,    # Don't save until you get to this iteration
        #  "keep_each" => 2.5e3,  # how often to save a model post burnin
        #  "print_each" => 1e5,   # how often to print to the screen if you are impatient and chekcing it (like I do)
        #####
         "max_depth" => 660,    # in km
         "min_depth" => 0,      # in km, should just be zero
         "rotation" => 20,      # used to rotate the grid. Helpful when using linear arrays so everthing is on one axis
         "ZnodeSpacing" => 20, #10   # in km, how often in the z-direction do you grid the ray paths
         "buffer" => 100,       # in km. If zero, edges of the modeling domain are [max_x min_x max_y min_y] station locations. Add buffer to broader the imaging domain. 
          ######map parameters
         "arc_inc" => 5,        # for Pn. Ignore it. 
         "XYnodeSpacing" => 20  #5  # not actually part of the inversion, but used for synthetic models.
    )
    
    return TD_parameters

end
