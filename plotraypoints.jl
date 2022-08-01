using HDF5,JLD
using Plots

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

model = load("./IndividualModels/chain#7_18.jld")

p = contourf(vec(dataStruct["xVec"]), vec(dataStruct["zVec"]), vec(model["zeta"]), c=:jet, yflip=true,legend=false)
for i in 1:381
    ray = dataStruct["rayX"][:,i]
    ll = 131 - length(findall(x -> isnan(x), ray))
    if ll >1
        println(ll)
        scatter!(p, vec(dataStruct["rayX"][:,i]), vec(dataStruct["rayZ"][:,i]))
    end
end
# scatter!(p,vec(dataStruct["rayX"][:,379]),vec(dataStruct["rayZ"][:,379]))
# scatter!(p,vec(dataStruct["rayX"][:,18]),vec(dataStruct["rayZ"][:,18]))
# scatter!(p,vec(dataStruct["rayX"][:,35]),vec(dataStruct["rayZ"][:,35]))
# scatter!(p,vec(dataStruct["rayX"][:,99]),vec(dataStruct["rayZ"][:,99]))
# scatter!(p,vec(dataStruct["rayX"][:,129]),vec(dataStruct["rayZ"][:,129]))
# scatter!(p,vec(dataStruct["rayX"][:,78]),vec(dataStruct["rayZ"][:,78]))
savefig(p,"voronoinraypoints.png")

