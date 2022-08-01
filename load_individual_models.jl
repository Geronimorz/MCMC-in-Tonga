using HDF5, JLD
using Plots, Distributed

@everywhere include("define_TDstructure.jl")
@everywhere include("MCsub.jl")
@everywhere include("load_data_Tonga.jl")

ENV["GKSwstype"] = "nul"


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

models = []
# models::Array{Any,1}

# model1 = []

model1 = pmap(x -> load_model(x,num_models_per_chain),1:TD_parameters["n_chains"])
# models = model1[1]
@time for i = 1:length(model1)
    append!(models,model1[i])
end
# models = pmap(x -> load_model(x,num_models_per_chain),1:TD_parameters["n_chains"])
# @time for i in 1:TD_parameters["n_chains"]
#     # for j in 1:num_models_per_chain
#     for j in 1:3
#         println("loading chain#" * string(i) * ('_') * string(j) * (".jld"))
#         loadname = "./IndividualModels/chain#" * string(i) * ('_') * string(j) * (".jld")
#         model = load(loadname)
#         push!(models, model)
#     end
# end

n = length(models)

zeta = pmap(x -> models[x]["zeta"], 1:n)
poststd = std(zeta)

model_all = pmap(x -> models[x]["model"], 1:n)

nCells = pmap(x -> model_all[x]["nCells"][1], 1:n)
likelihood = pmap(x -> model_all[x]["likelihood"][1], 1:n)
phi = pmap(x -> model_all[x]["phi"][1], 1:n)
# allSig = pmap(x -> model_all[x]["allSig"][1], 1:n)

###Plot as we want~~~
###Plot nCells
ll = TD_parameters["n_chains"]
maxcells = TD_parameters["max_cells"]
gr()
p = scatter(1:length(nCells), nCells)
for i = 1:ll
    y = [0,maxcells]
    x = ones(length(y)) .* (i * num_models_per_chain)
    plot!(p, x, y, color=:pink)
end
savefig(p,"nCells")

@time Plot_model_uncertainty(dataStruct, poststd, TD_parameters)

mask = ones(length(poststd))

for i = 1:length(poststd)
    if poststd[i] > 5
        mask[i] = NaN
    end
end

mask_model = mask .* mean(zeta)

@time Plot_model_mask(dataStruct, mask_model, TD_parameters)





# if makeplot == 1 && (mod(used_num, 10) == 0 || used_num == 1)
#     ENV["GKSwstype"] = "nul"
    if TD_parameters["xzMap"] == 1
    @time p = contourf(vec(dataStruct["xVec"]), vec(dataStruct["zVec"]), vec(zeta), c=:jet, yflip=true)
    @time savefig(p, "./voronoi/chain#" * string(chain) * ('_') * string(used_num) * (".png"))
    end
#     if TD_parameters["xyMap"] == 1
#         p = contourf(vec(dataStruct["xVec"]), vec(dataStruct["yVec"]), vec(zeta), c=:jet)
#         savefig(p, "./voronoi/chain#" * string(chain) * ('_') * string(used_num) * (".png"))
#     end
# end



