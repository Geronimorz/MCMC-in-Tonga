using HDF5, JLD
using Plots, Distributed

@everywhere include("define_TDstructure.jl")
@everywhere include("MCsub.jl")
@everywhere include("load_data_Tonga.jl")
# @everywhere include("loadsub.jl")

ENV["GKSwstype"] = "nul"

if isdir("figures") == false
    mkdir("figures")
end
if isdir("./figures/xzMap") == false
    mkdir("./figures/xzMap")
end
if isdir("./figures/xzUncertainty") == false
    mkdir("./figures/xzUncertainty")
end
if isdir("./figures/xzMask") == false
    mkdir("./figures/xzMask")
end
if isdir("./figures/xzVoronoi") == false
    mkdir("./figures/xzVoronoi")
end
if isdir("./figures/xyMap") == false
    mkdir("./figures/xyMap")
end
if isdir("./figures/xyUncertainty") == false
    mkdir("./figures/xyUncertainty")
end
if isdir("./figures/xyMask") == false
    mkdir("./figures/xyMask")
end
if isdir("./figures/xyVoronoi") == false
    mkdir("./figures/xyVoronoi")
end

if isdir("./figures/xyAverage") == false
    mkdir("./figures/xyAverage")
end
if isdir("./figures/xzAverage") == false
    mkdir("./figures/xzAverage")
end


TD_parameters = define_TDstructrure()
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

num_models_per_chain = Int64((TD_parameters["n_iter"] - TD_parameters["burn_in"]) / TD_parameters["keep_each"])
models = load("model.jld")
println("####models loaded####")


plot_per_models = 50
cmax = 20

# models = pmap(x -> CrossSectionInterpolation(models["model"][x][1:plot_per_models:num_models_per_chain], dataStruct, TD_parameters, plot_per_models, x), 1:TD_parameters["n_chains"])
models = pmap(x -> CrossSectionInterpolation(models["model"][x], dataStruct, TD_parameters, x), 1:TD_parameters["n_chains"])
models = permutedims(hcat(models...))[:]
println("####models interpolated####")

n = length(models)

if TD_parameters["xzMap"] == 1
    # zeta_xz = pmap(x -> models[x]["zeta_xz"], 1:n)
    zeta_xz = []
    for i = 1:n
        append!(zeta_xz, [models[i]["zeta_xz"]])
    end
    model_mean_xz = mean(zeta_xz)
    poststd_xz = std(zeta_xz)

    mask_xz = ones(length(poststd_xz))
    for i = 1:length(poststd_xz)
        if poststd_xz[i] > 5
            mask_xz[i] = NaN
        end
    end
    mask_model_xz = mask_xz .* model_mean_xz

    # TD_parameters["xyMap"] = 0 
    cd("./figures/xzAverage")
    @time Plot_model_mean(dataStruct, model_mean_xz, TD_parameters, cmax)
    @time Plot_model_uncertainty(dataStruct, poststd_xz, TD_parameters)
    @time Plot_model_mask(dataStruct, mask_model_xz, TD_parameters, cmax)
    cd("../..")
end

# TD_parameters["xyMap"] = 1
# if TD_parameters["xyMap"] == 1
#     zeta_xy = pmap(x -> models[x]["zeta_xy"], 1:n)
#     model_mean_xy = mean(zeta_xy)
#     poststd_xy = std(zeta_xy)

#     mask_xy = ones(length(poststd_xy))
#     for i = 1:length(poststd_xy)
#         if poststd_xy[i] > 5
#             mask_xy[i] = NaN
#         end
#     end
#     mask_model_xy = mask_xy .* model_mean_xy
    
#     TD_parameters["xzMap"] = 0
#     cd("./figures/xyAverage")
#     @time Plot_model_mean(dataStruct, model_mean_xy, TD_parameters)
#     @time Plot_model_uncertainty(dataStruct, poststd_xy, TD_parameters)
#     @time Plot_model_mask(dataStruct, mask_model_xy, TD_parameters)
#     cd("../..")
# end

# TD_parameters["xzMap"] = 1

pmap(x -> PlotModelsOverIterations(models[Int64((x - 1) * num_models_per_chain + 1):5:Int64(x * num_models_per_chain)], dataStruct, TD_parameters, x,cmax), 1:TD_parameters["n_chains"])
# @time for x = 1:TD_parameters["n_chains"]
#     ind1 = Int64((x - 1) * num_models_per_chain + 1)
#     ind2 = Int64(x * num_models_per_chain)
#     PlotModelsOverIterations(models[ind1:ind2], dataStruct, TD_parameters, x, cmax)
# end




# global models = loadmodel["model"][1][1:num_models_per_chain]
# @time for i in 2:length(loadmodel["model"])
#     global models = vcat(models,loadmodel["model"][i])
#     # for k = 1:num_models_per_chain  
#     #     push!(models, loadmodel["model"][i][k])
#     # end
# end

# n = length(models)

# aaa = pmap(x -> CrossSectionInterpolation(models[Int64((x-1)*50+1):Int64(x*50)], dataStruct, TD_parameters), 1:TD_parameters["n_chains"])

# println(size(aaa))

