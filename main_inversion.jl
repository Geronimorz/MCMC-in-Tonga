using JLD, HDF5
using Distributed 
using Glob
@everywhere include("DefStruct.jl")
@everywhere include("define_TDstructure.jl")

@everywhere include("MCsub.jl")
@everywhere include("load_data_Tonga.jl")
@everywhere include("TD_inversion_function.jl")

@time TD_parameters = define_TDstructrure()
@time dataStruct = load_data_Tonga(TD_parameters)
println("--------data loaded-------")

@time models = pmap(x -> TD_inversion_function(TD_parameters, dataStruct, x), 1:TD_parameters.n_chains)

plot_model_hist(models, dataStruct, TD_parameters, 20.0)
@time save("model.jld","model",models)
println("--------finish-------")

model_checkpoint_lists = glob("chain*jld")
[rm(model_checkpoint_lists[i]) for i in 1:length(model_checkpoint_lists)]