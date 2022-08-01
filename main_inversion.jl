using JLD, HDF5
# using SharedArrays
# @everywhere global modelnum = []
using Distributed
@everywhere include("define_TDstructure.jl")
@everywhere include("MCsub.jl")
@everywhere include("load_data_Tonga.jl")
@everywhere include("load_data_STdata3.jl")
@everywhere include("TD_inversion_function.jl")

# models::Array{Tuple{Int64,Array{Float64,2}},1}

if isdir("voronoi")==false
    mkdir("voronoi")
end

if isdir("IndividualModels")==false
    mkdir("IndividualModels")
end

@time TD_parameters = define_TDstructrure()
if TD_parameters["STdata3"] == 1
    println("*******loading STdata3********")
    @time dataStruct = load_data_STdata3()
else
    println("*******loading Tonga data********")
    @time dataStruct = load_data_Tonga(TD_parameters)
end
if TD_parameters["add_yVec"] == 0
    println("*******2d********")
else
    println("*******3d********")
end

# m = zeros(2)
# num = SharedArray{Float64}(2)
# # zeta = []
# @time @distributed for nchain = 1:TD_parameters["n_chains"]
# # for nchain = 1:TD_parameters["n_chains"]
#     global (used_num, zeta) = TD_inversion_function(TD_parameters, dataStruct, nchain)
#     # println(used_num)
#     num[nchain] = used_num
#     # m[nchain] = zeta
#     # push!(m, zeta)
#     # println(num)
# end
# println(num)


# println(num)

# @time pmap(x -> TD_inversion_function(TD_parameters, dataStruct, x), 1:TD_parameters["n_chains"])
@time models = pmap(x -> TD_inversion_function(TD_parameters, dataStruct, x), 1:TD_parameters["n_chains"])

@time save("model.jld","model",models)
println("--------finish-------")
# println("in all",varinfo())
# used_num = pmap(x -> models[x][1],1:TD_parameters["n_chains"])
# averaged_model = pmap(x -> models[x][2],1:TD_parameters["n_chains"])

# model_mean = mean(averaged_model)
# # num_model_per_chain = Int((TD_parameters["n_iter"] - TD_parameters["burn_in"]) / TD_parameters["keep_each"])
# # bestmodels = models[1][1:num_model_per_chain]

# # @time for i = 1:length(models)
# #     num_model_per_chain = length(models[i])
# #     println(num_model_per_chain)
# #     for k = 1:num_model_per_chain
# #         push!(bestmodels, models[i][k])
# #     end
# # end

# # @time model_mean = MakeDistributions(bestmodels, dataStruct, TD_parameters)

# # if TD_parameters["STdata3"] == 1
# #     savename = "STdata3_inversion_n_iter_" * string(TD_parameters["n_iter"]) * "_try_1.jld"
# # else
# #     savename = "Tonga_inversion_n_iter_" * string(TD_parameters["n_iter"]) * "_try_1.jld"
# # end
# # # i = 1
# # # savename = "Tonga_inversion_n_iter_" * string(TD_parameters["n_iter"]) * "_try_" * string(i) * ".jld"
# # # while isfile(savename) == true
# # #     global i += 1
# # #     savename = "Tonga_inversion_n_iter_" * string(TD_parameters["n_iter"]) * "_try_" * string(i) * ".jld"
# # # end
# # @time save(savename, "bestmodels", bestmodels, "dataStruct", dataStruct, "TD_parameters", TD_parameters, "model_mean", model_mean)


# @time Plot_model_mean(dataStruct, model_mean, TD_parameters)

##########read inverted file############
# aaa = load("Tonga_inversion_n_iter_1000.0_try_1.jld")
# aaa = load("STdata3_inversion_n_iter_1000.0_try_1.jld")
# allsigs = pmap(x -> aaa["bestmodels"][x]["allSig"][1],1:length(aaa["bestmodels"]))
# phi = pmap(x -> aaa["bestmodels"][x]["phi"][1],1:length(aaa["bestmodels"]))
# action = pmap(x -> aaa["bestmodels"][x]["action"][1],1:length(aaa["bestmodels"]))




