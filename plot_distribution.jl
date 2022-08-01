using HDF5, JLD
using Plots, Distributed
pyplot()

@everywhere include("define_TDstructure.jl")
@everywhere include("MCsub.jl")
@everywhere include("load_data_Tonga.jl")
# @everywhere include("loadsub.jl")

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
all_model_num = Int64(num_models_per_chain * TD_parameters["n_chains"])
models = load("model.jld")
println("####models loaded####")



# julia> models[1]
# Dict{String,AbstractArray{T,1} where T} with 15 entries:
#   "xCell"      => [687.153, 983.883, 915.012, 751.487, -58.4778, 89.7843, 869.4…
#   "zCell"      => [648.463, 270.21, 28.8729, 620.605, 168.802, 598.539, 362.763…
#   "xVec"       => -79.47730270810919:20.0:1060.5226972918908
#   "action"     => [4]
#   "zVec"       => 0:20:660
#   "zeta_scale" => [9.09684, 45.3684, 23.2372, 31.2416, 0.67402, 16.4487, 8.0592…
#   "yVec"       => -164.40158664642206:20.0:495.59841335357794
#   "yCell"      => [-8.9557, 332.824, 445.33, 233.326, 215.729, 67.5656, -57.695…
#   "layer"      => [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
#   "id"         => 1.0:1.0:16.0
#   "phi"        => [1]
#   "zeta"       => [9.09684, 45.3684, 23.2372, 31.2416, 0.67402, 16.4487, 8.0592…
#   "likelihood" => [1]
#   "accept"     => [0]
#   "nCells"     => [9.0]

#plot nCells vs iteration for 1 chain
chain1 = models["model"][1]
nCells1 = pmap(x -> chain1[x]["nCells"][1], 1:num_models_per_chain) 
p1 = scatter(1:num_models_per_chain,nCells1)
title!(p1, "nCells vs iteration (chain1)")
savefig(p1, "chain1_nCells")

zetas1 = pmap(x -> chain1[x]["zeta"], 1:num_models_per_chain) 
zeta1 = zetas1[1]
for i = 2:num_models_per_chain
    zeta1 = vcat(zeta1,zetas1[i])
end
p3 = scatter(1:num_models_per_chain, zeta1)
title!(p3, "zeta vs iteration (chain1)")
savefig(p3, "chain1_zeta")

#plot nCells and zeta posterior distribution
models = permutedims(hcat(models["model"]...))[:]
nCells2 = pmap(x -> models[x]["nCells"][1], 1:all_model_num)
p2 = histogram(nCells2)
title!(p2, "nCells posterior distribution (in debug mode)")
savefig(p2, "nCells_distribution")

zetas2 = pmap(x -> models[x]["zeta"], 1:all_model_num)
zeta2 = zetas2[1]
for i = 2:all_model_num
    zeta2 = vcat(zeta2,zetas2[i])
end
p4 = histogram(zeta2)
title!(p4, "zeta posterior distribution (in debug mode)")
savefig(p4, "zeta_distribution")