using JLD, HDF5
# using SharedArrays
# @everywhere global modelnum = []
using Distributed
@everywhere include("DefStruct.jl")
@everywhere include("define_TDstructure.jl")
@everywhere include("MCsub.jl")
@everywhere include("load_data_Tonga.jl")
@everywhere include("TD_inversion_function.jl")

@time TD_parameters = define_TDstructrure()
println("*******loading Tonga data********")
@time dataStruct = load_data_Tonga(TD_parameters)

if TD_parameters.add_yVec == 0
    println("*******2d********")
else
    println("*******3d********")
end

@time models = pmap(x -> TD_inversion_function(TD_parameters, dataStruct, x), 1:TD_parameters.n_chains)

@time save("model.jld","model",models)
println("--------finish-------")

