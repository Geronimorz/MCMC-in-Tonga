using HDF5, JLD
using Plots, Distributed

@everywhere include("DefStruct.jl")
@everywhere include("define_TDstructure.jl")
@everywhere include("MCsub.jl")
@everywhere include("load_data_Tonga.jl")

ENV["GKSwstype"] = "nul"

cmax = 20.0

FigDirLst = ["figures", "voronoi", 
        "./figures/xzMap", "./figures/xzUncertainty", "./figures/xzMask", "./figures/xzVoronoi", "./figures/xzMean",
        "./figures/xyMap", "./figures/xyUncertainty", "./figures/xyMask", "./figures/xyVoronoi", "./figures/xyMean"
            ]

for iDir in FigDirLst
    if isdir(iDir) == false
        mkdir(iDir)
    end
end

TD_parameters = define_TDstructrure()

println("*******loading Tonga data********")
@time dataStruct = load_data_Tonga(TD_parameters)
if TD_parameters.add_yVec == 0
    println("*******2d********")
else
    println("*******3d********")
end

num_models_per_chain = Int64((TD_parameters.n_iter - TD_parameters.burn_in) / TD_parameters.keep_each)
models = load("model.jld")
println("####models loaded####")

models = pmap(x -> CrossSectionInterpolation(models["model"][x], dataStruct, TD_parameters, x), 1:TD_parameters.n_chains)
models = permutedims(hcat(models...))[:]
println("####models interpolated####")

n = length(models)

# CrossSectionTypeList = ["xz", "xy"]
# for crosssectionType in CrossSectionTypeList
#     if (crosssectionType == "xz" && TD_parameters.xzMap == false) || 
#         (crosssectionType == "xy" && TD_parameters.xyMap == false)
#         continue
#     end
#     zeta_xz = []
#     for i = 1:n
#         ###############How to link a field name to a variable???#########
#         append!(zeta_xz, [models[i].zeta_xz]) 
#     end
#     model_mean_xz = mean(zeta_xz)
#     poststd_xz = std(zeta_xz)

#     mask_xz = ones(length(poststd_xz))
#     for i = 1:length(poststd_xz)
#         if poststd_xz[i] > 5
#             mask_xz[i] = NaN
#         end
#     end
#     mask_model_xz = mask_xz .* model_mean_xz

#     @time Plot_model(dataStruct, model_mean_xz, TD_parameters, cmax, "Mean", crosssectionType)
#     @time Plot_model(dataStruct, poststd_xz, TD_parameters, maximum(poststd_xz), "Uncertainty", crosssectionType)
#     @time Plot_model(dataStruct, mask_model_xz, TD_parameters, cmax, "Mask", crosssectionType)

# end

if TD_parameters.xzMap == true
    zeta_xz = []
    for i = 1:n
        append!(zeta_xz, [models[i].zeta_xz])
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

    @time Plot_model(dataStruct, model_mean_xz, TD_parameters, cmax, "Mean", "xz")
    @time Plot_model(dataStruct, poststd_xz, TD_parameters, maximum(poststd_xz), "Uncertainty", "xz")
    @time Plot_model(dataStruct, mask_model_xz, TD_parameters, cmax, "Mask", "xz")
end

if TD_parameters.xyMap == true
    zeta_xy = []
    for i = 1:n
        append!(zeta_xy, [models[i].zeta_xy])
    end
    model_mean_xy = mean(zeta_xy)
    poststd_xy = std(zeta_xy)

    mask_xy = ones(length(poststd_xy))
    for i = 1:length(poststd_xy)
        if poststd_xy[i] > 5
            mask_xy[i] = NaN
        end
    end
    mask_model_xy = mask_xy .* model_mean_xy

    @time Plot_model(dataStruct, model_mean_xy, TD_parameters, cmax, "Mean", "xy")
    @time Plot_model(dataStruct, poststd_xy, TD_parameters, maximum(poststd_xz), "Uncertainty", "xy")
    @time Plot_model(dataStruct, mask_model_xy, TD_parameters, cmax, "Mask", "xy")
end

if TD_parameters.plot_voronoi == 1
    pmap(x -> PlotModelsOverIterations(models[Int64((x - 1) * num_models_per_chain + 1):5:Int64(x * num_models_per_chain)], dataStruct, TD_parameters, x,cmax), 1:TD_parameters.n_chains)
end