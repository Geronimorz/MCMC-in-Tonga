using Random, Distributions, Distributed
using Plots
using HDF5,JLD

function lonlat2xy(
    lon0::Float64, 
    lat0::Float64, 
    beta::Float64, 
    lon1, 
    lat1)
    # function lonlat2xy(lon0::Float64, lat0::Float64, beta::Float64, lon1::Array{Float64,2}, lat1::Array{Float64,2})
    # Convert lon, lat to x, y
    # lon0, lat0: reference location
    # beta: angle of x axis from east (counterclockwise)
    re = 6371
    # pi = 4.*atan(1.)
    r2d = 180.0 / pi
    xx = (lon1 .- lon0) .* re ./ r2d
    yy = (lat1 .- lat0) .* re ./ r2d
    x1 = (xx .- yy .* tan.(beta)) .* cos.(beta)
    y1 = x1 .* tan.(beta) .+ yy ./ cos.(beta)

    ###################CANNOT WORK!!!! HOW TO DEAL WITH IT?!!##############
    # if abs.(x1) < 0.01
    #     x1 = 0
    # end
    # if abs.(y1) < 0.01
    #     y1 = 0
    # end
    return x1, y1
end

function xy2lonlat(
    lon0::Float64, 
    lat0::Float64, 
    beta::Float64, 
    x2::Array{Any,2}, 
    y2::Array{Any,2}
    )
    # Convert x, y to lon, lat 
    # lon0, lat0: reference location
    # beta: angle of x axis from east (counterclockwise)
    re = 6371
    # pi = 4.*atan(1.)
    r2d = 180.0 / pi
    yy = (y2 .- x2 .* tan.(beta)) .* cos.(beta)
    xx = yy .* tan.(beta) .+ x2 ./ cos.(beta)
    lon2 = xx .* r2d ./ re .+ lon0
    lat2 = yy .* r2d ./ re .+ lat0
    
    return lon2, lat2
end

function interp1(
    x::Array{Float64,1}, 
    y::Array{Float64,1}, 
    xx::Array{Float64,1}
    )

    xx_length = length(xx)
    x_length = length(x)
    # yy::Array{Float64,1} = zeros(xx_length)
    yy = NaN .* xx
    for i in 1:xx_length
        for j in 1:x_length
            if xx[i] >= x[j] && xx[i] < x[j + 1]
                # y1 = y[j] + (xx[i] - x[j]) / (x[j + 1] - x[j]) * (y[j + 1] - y[j])
                # println(xx[i], ' ', y1)
                yy[i] = y[j] + (xx[i] - x[j]) / (x[j + 1] - x[j]) * (y[j + 1] - y[j])
            end
        end
    end
    return yy
end

function build_starting(
    TD_parameters::parameters, 
    dataStruct::MutabledataStruct
    )

    # log uniform distribution, in Byrnes and Bezada, 2020, eq. 11

    xVec = dataStruct.xVec
    yVec = dataStruct.yVec
    zVec = dataStruct.zVec
    nCells = floor.(exp.(rand(1) * log(TD_parameters.max_cells / 
            TD_parameters.min_cells) .+ log(TD_parameters.min_cells)))
    # nCells = floor.(rand(1)*(TD_parameters.max_cells- TD_parameters.min_cells)) .+ TD_parameters.min_cells

    # ATTENTION! node_change hasn't been added

    xCell = min(xVec...) .+ (max(xVec...) - min(xVec...)) .* rand(Int(nCells[1]))
    yCell = min(yVec...) .+ (max(yVec...) - min(yVec...)) .* rand(Int(nCells[1]))
    zCell = min(zVec...) .+ (max(zVec...) - min(zVec...)) .* rand(Int(nCells[1]))


    if TD_parameters.prior == 1 
        # Uniform 
            # for absolute t*
        zeta = rand(Int(nCells[1])) .* TD_parameters.zeta_scale
            # for relative t*
        # zeta = rand(Int(nCells[1])) .* TD_parameters["zeta_scale"] * 2 .- TD_parameters["zeta_scale"]  
    elseif TD_parameters.prior == 2
        # Normal 
        zeta = rand(Normal(0, TD_parameters.zeta_scale), Int(nCells[1]))
    elseif TD_parameters.prior == 3
        # Exponential 
        zeta = -log.(rand(Int(nCells[1]))) .* TD_parameters.zeta_scale
    end

    # ATTENTION! If discontinuity is not a field:
    layer = ones(length(xCell))

    model = Model(nCells[1],
        xVec, yVec, zVec,
        xCell, yCell, zCell,
        zeta, zeta,
        layer, -1, -1, -1, -1, -1, -1, -1, -1
    )
    # setindex!(dataStruct, model["allSig"] .* ones(length(dataStruct["tS"])), "allSig")
    
    (model, dataStruct, valid) = evaluate(model, dataStruct, TD_parameters)
    valid = 1
    return model, dataStruct, valid
end

function evaluate(
    model_original::Model, 
    dataStruct::MutabledataStruct, 
    TD_parameters::parameters
    )
    valid            = 1
    likelihood = 1
    phi = 1
    model = deepcopy(model_original)
    model.phi = phi
    model.likelihood = likelihood

    if TD_parameters.debug_prior == 1
        return model, dataStruct, valid
    end

    (m, n) = size(dataStruct.rayX)

    ptS = zeros(n, 1)
    for i in 1:n
        zeta0 = Interpolation(TD_parameters, model, dataStruct.rayX[:,i], dataStruct.rayY[:,i], dataStruct.rayZ[:,i]) 
        rayzeta = []
 
        endpoint = length(zeta0)
        rayzeta = 0.5 .* (zeta0[1:endpoint - 1] + zeta0[2:endpoint])
    
        rayl = dataStruct.rayL[:,i]
        index = findall(x -> isnan(x), rayl)
        if isempty(index)
            ptS[i] = sum(rayl .* dataStruct.rayU[:,i] .* (rayzeta ./ 1000))
        else
            rayl = dataStruct.rayL[1:index[1] - 1,i]  
            rayu = dataStruct.rayU[1:index[1] - 1,i]
            ptS[i] = sum(rayl .* rayu .* (rayzeta ./ 1000))
        end
    end

    (m, n) = size(dataStruct.rayZ)

    tS = dataStruct.tS


    C = 0

    for k in 1:length(dataStruct.allSig)
        C += (ptS - tS)[k].^2 .* 1.0 / dataStruct.allSig[k][1]^2
        # println(dataStruct["allSig"][k][1])
    end
    model.phi = C
    model.ptS = ptS
    model.tS = tS

    likelihood = sum(-log.(dataStruct.allSig * sqrt(2 * pi)) * length(tS)) 
    - sum(0.5 * ((vec(ptS) .- vec(tS)) ./ vec(dataStruct.allSig)).^2)

    model.likelihood = likelihood
    
    return model, dataStruct, valid
end

function IDW(
    TD_parameters::parameters, 
    xCell::Array{Float64,2}, 
    yCell::Array{Float64,2}, 
    zCell::Array{Float64,2}, 
    zetaCell::Array{Float64,2}, 
    X, #::Array{Float64,1},#StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}, 
    Y, #::Array{Float64,1}, 
    Z, #::Array{Float64,1}, 
    nCells::Int64, 
    npoints::Int64)

    zeta = zeros(npoints, 1) 
    
    if nCells == 0
        return
    end

    X = X[1:npoints]

    xCell = xCell[1:nCells]
    yCell = yCell[1:nCells]
    zCell = zCell[1:nCells]
    zetaCell = zetaCell[1:nCells]

    if TD_parameters.add_yVec == 0
        for k in 1:npoints
            distance = (X[k] .- xCell).^2 + (Z[k] .- zCell).^2
            v_sum = sum(zetaCell ./ distance)
            inv_sum = sum(1 ./ distance)
            zeta[k] = v_sum / inv_sum
        end
    else
        for k in 1:npoints
            distance = sqrt.((X[k] .- xCell).^2 + (Y[k] .- yCell).^2 + (Z[k] .- zCell).^2)
            v_sum = sum(zetaCell ./ distance)
            inv_sum = sum(1 ./ distance)
            zeta[k] = v_sum / inv_sum
        end
    end
    return zeta
end

function NearestInterpolation(
    TD_parameters::parameters, 
    xCell::Array{Float64,2}, 
    yCell::Array{Float64,2}, 
    zCell::Array{Float64,2}, 
    zetaCell::Array{Float64,2}, 
    X, #::Array{Float64,1}, 
    Y, #::Array{Float64,1}, 
    Z, #::Array{Float64,1}, 
    nCells::Int64, 
    npoints::Int64)
    
    zeta = zeros(npoints, 1) 
    
    if nCells == 0
        return
    end

    X = X[1:npoints]

    xCell = xCell[1:nCells]
    yCell = yCell[1:nCells]
    zCell = zCell[1:nCells]
    zetaCell = zetaCell[1:nCells]


    # find the nearest cell to the points in each ray
    if TD_parameters.add_yVec == 0
        for k in 1:npoints
            tmp = (X[k] .- xCell).^2 + (Z[k] .- zCell).^2
            zeta[k] = zetaCell[findmin(tmp)[2]]
        end
    else
        for k in 1:npoints
            tmp = (X[k] .- xCell).^2 + (Y[k] .- yCell).^2 + (Z[k] .- zCell).^2
            zeta[k] = zetaCell[findmin(tmp)[2]]
        end
    end
    return zeta
end

function Interpolation(
    TD_parameters::parameters,
    model::Model, 
    X, Y, Z
    )
    
    xCell        = [ model.xCell;zeros(200 - Int(model.nCells), 1) ]
    yCell        = [ model.yCell;zeros(200 - Int(model.nCells), 1) ]
    zCell        = [ model.zCell;zeros(200 - Int(model.nCells), 1) ]
    # ATTENTION! IDK why doing this
    zeta    = [ model.zeta; zeros(200 - Int(model.nCells), 1) ]

    if .~isempty(findall(x -> isnan.(x), X))
        npoints = findall(x -> isnan.(x), X)[1] - 1
    else
        npoints = length(X)
    end
    if length(Y) == 1 # xzMap
        Y = Y .* ones(npoints)
    end
    if length(Z) == 1 # xyMap
        Z = Z .* ones(npoints)
    end
    # if length(X) == 1 # yzMap
    #     X = X .* ones(npoints)
    # end
    if TD_parameters.interp_style == 1  # Nearest Interpolation
        X[isnan.(X)] .= 99999
        Y[isnan.(Y)] .= 99999
        Z[isnan.(Z)] .= 99999
        zeta = NearestInterpolation(TD_parameters, xCell, yCell, zCell, zeta, X, Y, Z, Int(model.nCells), npoints)
    elseif TD_parameters.interp_style == 2  # Inverse Distance Weighting(IDW) 
        # X[isnan.(X)] .= 99999
        # Y[isnan.(Y)] .= 99999
        # Z[isnan.(Z)] .= 99999
        zeta = IDW(TD_parameters, xCell, yCell, zCell, zeta, X, Y, Z, Int(model.nCells), npoints)
    end

    return zeta
end

function CrossSectionInterpolation(
    models::Array{Any,1}, # models in an individual chain
    dataStruct::MutabledataStruct, 
    TD_parameters::parameters,
    chain::Int64
    )

    makeplot = TD_parameters.plot_voronoi
    ENV["GKSwstype"] = "nul"

    for k in 1:length(models)
        if TD_parameters.xzMap == 1
            zeta_xz = []
            for i in 1:length(dataStruct.zVec)
                tmp = Interpolation(TD_parameters, models[k], dataStruct.xVec, TD_parameters.ySlice, ones(length(dataStruct.xVec)) .* dataStruct.zVec[i])
                zeta_xz = vcat(zeta_xz, tmp)
            end
            models[k].zeta_xz = vec(zeta_xz)
            if makeplot == 1
                p = contourf(vec(dataStruct.xVec), vec(dataStruct.zVec), vec(zeta_xz), c=:jet,linewidth=0.001, clims=(0,20), yflip=true)
                scatter!(p, dataStruct.elonsX, vec(dataStruct.edep), marker=:o, color=:lightblue, label="events", markersize=4,legend=false)
                savefig(p, "./voronoi/chain#" * string(chain) * ('_') * string(k) * ("_xz_") * string(TD_parameters.ySlice) * (".png"))
            end
        end
        if TD_parameters.xyMap == 1
            zeta_xy = []
            for i in 1:length(dataStruct.yVec)
                tmp = Interpolation(TD_parameters, models[k], dataStruct.xVec, ones(length(dataStruct.xVec)) .* dataStruct.yVec[i], TD_parameters.zSlice)
                zeta_xy = vcat(zeta_xy, tmp)
            end
            models[k].zeta_xy = vec(zeta_xy)
            if makeplot == 1
                p = contourf(vec(dataStruct.xVec), vec(dataStruct.yVec), vec(zeta_xy), linewidth=0.001, clims=(0,20), c=:jet)
                savefig(p, "./voronoi/chain#" * string(chain) * ('_') * string(k) * ("_xy_") * string(TD_parameters.zSlice) * (".png"))
            end
        end
    end
    return models
end

function Plot_model(
    dataStruct::MutabledataStruct, 
    model_mean::Array{Float64}, 
    TD_parameters::parameters,
    cmax::Float64,
    Type::String,
    crosssectionType::String
    )

    ENV["GKSwstype"] = "nul"
    gr()
    closeenough = 2.0
    if Type == "Mask" 
        cmap = :jet
    elseif Type == "Uncertainty"
        cmap = :bone
    elseif Type == "Mean"
        cmap = :jet
    end

    dirname = "./figures/"*crosssectionType*Type

    if TD_parameters.add_yVec == 0 # Tonga 2D
        println("****2D Plotting****")
        tmpx = vcat(vec(dataStruct.dataX)', vec(dataStruct.elonsX)')
        tmpy = vcat(zeros(size(dataStruct.dataX))', vec(dataStruct.edep)') 
        tmpy = Array{Float64,2}(tmpy)
        p = contourf(vec(dataStruct.xVec), vec(dataStruct.zVec), vec(model_mean), xlabel="distance(km)", ylabel="depth(km)", yflip=true, c=cmap)

        title!("model_mean_2d_Tonga")
        savefig(p, "model_mean_2d_Tonga")

        title!("model_mean_2d_Tonga_0-300km")
        ylims!(p, (0, 300))
        savefig(p, "model_mean_2d_Tonga_0-300km") 

        ylims!(p, (0, 660))
        plot!(p, tmpx, tmpy, color=:white, legend=false)

        scatter!(p, dataStruct.dataX, zeros(size(dataStruct.dataX)), marker=:utri, color=:pink, label="station", markersize=6)
        scatter!(p, dataStruct.elonsX, vec(dataStruct.edep), marker=:o, color=:lightblue, label="events", markersize=4)
        title!("model_mean_2d_Tonga_ray")       
        savefig(p, dirname*"/model_mean_2d_Tonga_ray")   

        title!("model_mean_2d_Tonga__0-300km_ray")
        ylims!(p, (0, 300))
        savefig(p, dirname*"/model_mean_2d_Tonga__0-300km_ray") 
    else # 3d
        if crosssectionType == "xz"
            (m, n) = size(dataStruct.rayY)
            y0 = TD_parameters.ySlice
            nearrays = zeros(n)
            for i in 1:n
                yvec = dataStruct.rayY[:,i] .- y0    
                if sum(abs.(yvec)) - abs(sum(yvec)) > 1e-7
                    nearrays[i] = 1
                end  
                if min((abs.(yvec) .- closeenough)...) < 1e-7
                    nearrays[i] = 1
                end
            end
            nearrays = Array{Bool,1}(nearrays)
            tmpx = vcat(vec(dataStruct.dataX)', vec(dataStruct.elonsX)')
            tmpy = vcat(zeros(size(dataStruct.dataX))', vec(dataStruct.edep)') 
            tmpy = Array{Float64,2}(tmpy)
            nearraysX = tmpx[:,nearrays]
            nearraysY = tmpy[:,nearrays]

            p = contourf(vec(dataStruct.xVec), vec(dataStruct.zVec), vec(model_mean), linewidth=0.001, xlabel="distance(km)", ylabel="depth(km)", yflip=true, clims=(0,cmax), c=cmap)
           
            title!("Model_"*Type*"_xzMap" * string(TD_parameters.ySlice) * "km")
            savefig(p, dirname*"/Model_"*Type*"_xzMap" * string(TD_parameters.ySlice) * "km")

            title!("Model_"*Type*"_xzMap_0-300km" * string(TD_parameters.ySlice) * "km")
            # ylims!(p, (0, 300))
            # savefig(p, "model_mean_xzMap_0-300km" * string(TD_parameters["ySlice"][1]) * "km") 

            ylims!(p, (0, 660))
            scatter!(p, dataStruct.elonsX, vec(dataStruct.edep), marker=:o, color=:lightblue, label="events", markersize=4,legend=false)
            savefig(p, dirname*"/Model_"*Type*"_xzMap_events" * string(TD_parameters.ySlice) * "km") 

            plot!(p, tmpx, tmpy, color=:white, legend=false)
            plot!(p, nearraysX, nearraysY, color=:forestgreen)
            scatter!(p, dataStruct.dataX, zeros(size(dataStruct.dataX)), marker=:utri, color=:pink, label="station", markersize=6)
            
            title!("Model_"*Type*"_xzMap_ray" * string(TD_parameters.ySlice) * "km")       
            savefig(p, dirname*"/Model_"*Type*"_xzMap_ray" * string(TD_parameters.ySlice) * "km")   
  
            # title!("model_mean_xzMap_0-300km_ray" * string(TD_parameters["ySlice"][1]) * "km")
            # ylims!(p, (0, 300))
            # savefig(p, "model_mean_xzMap_0-300km_ray" * string(TD_parameters["ySlice"][1]) * "km") 
        end

        if crosssectionType == "xy"
            (m, n) = size(dataStruct.rayZ)
            z0 = TD_parameters.zSlice
            nearrays = zeros(n)
            for i in 1:n
                zvec = dataStruct.rayZ[:,i] .- z0    
                if sum(abs.(zvec)) - abs(sum(zvec)) > 1e-7
                    nearrays[i] = 1
                end  
                if min((abs.(zvec) .- closeenough)...) < 1e-7
                    nearrays[i] = 1
                end
            end
            nearrays = Array{Bool,1}(nearrays)
            tmpx = vcat(vec(dataStruct.dataX)', vec(dataStruct.elonsX)')
            tmpy = vcat(vec(dataStruct.dataY)', vec(dataStruct.elatsY)') 
            # tmpy = Array{Float64,2}(tmpy)
            nearraysX = tmpx[:,nearrays]
            nearraysY = tmpy[:,nearrays]
        
            p = contourf(dataStruct.xVec, dataStruct.yVec, vec(model_mean), xlabel="X(km)", ylabel="Y(km)", clims=(0,30), c=cmap)
            title!("Model_"*Type*"_xyMap" * string(TD_parameters.zSlice) * "km")        
            savefig(p, dirname*"/Model_"*Type*"_xyMap" * string(TD_parameters.zSlice) * "km")

            plot!(p, tmpx, tmpy, color=:white, legend=false)
            plot!(p, nearraysX, nearraysY, color=:forestgreen)
            scatter!(p, dataStruct.dataX, dataStruct.dataY, shape=:utri, color=:pink, label="station", markersize=6)
            scatter!(p, dataStruct.elonsX, dataStruct.elatsY, shape=:o, color=:lightblue, label="events", markersize=4)
            title!("Model_"*Type*"_xyMap_ray" * string(TD_parameters.zSlice) * "km")
        
            savefig(p, dirname*"/Model_"*Type*"_xyMap_ray" * string(TD_parameters.zSlice) * "km")       
        end

    end


end

function PlotModelsOverIterations(
    models, # models in an individual chain
    dataStruct::MutabledataStruct, 
    TD_parameters::parameters,
    chain::Int64,
    cmax::Float64
    )
    ENV["GKSwstype"] = "nul"
    gr()
    max_std = 5 

    
    n = length(models)

    if TD_parameters.xzMap == true

        voronoi_xz = vec(models[1].zeta_xz)
        p = contourf(vec(dataStruct.xVec), vec(dataStruct.zVec), voronoi_xz, linewidth=0.001, c=:jet, clims=(0,cmax), yflip=true)
        ylims!(p, (0, 660))
        scatter!(p, dataStruct.elonsX, vec(dataStruct.edep), marker=:o, color=:lightblue, label="events", markersize=4,legend=false)
        savefig(p, "./figures/xzVoronoi/chain#" * string(chain) * ("_1") * (".png"))

        average_xz = vec(models[1].zeta_xz)
        std_xz = zeros(length(average_xz))
        mask = ones(length(std_xz))
        for i = 1:length(std_xz)
            if std_xz[i] > max_std
                mask[i] = NaN
            end
        end
        mask_xz = mask .* average_xz
        Plot_Contours(dataStruct, average_xz, TD_parameters, chain, 1, 1, 1, cmax)
        Plot_Contours(dataStruct, std_xz, TD_parameters, chain, 1, 1, 2, cmax)
        Plot_Contours(dataStruct, mask_xz, TD_parameters, chain, 1, 1, 3, cmax)
    end

    if TD_parameters.xyMap == true
        voronoi_xy = vec(models[1].zeta_xy)
        p = contourf(vec(dataStruct.xVec), vec(dataStruct.yVec), voronoi_xy, clims=(0,cmax), c=:jet)
        savefig(p, "./figures/xyVoronoi/chain#" * string(chain) * ("_1") * (".png"))

        average_xy = vec(models[1].zeta_xy)
        std_xy = zeros(length(average_xy))
        mask = ones(length(std_xy))
        for i = 1:length(std_xy)
            if std_xy[i] > max_std
                mask[i] = NaN
            end
        end
        # mask_xy = mask .* mean(average_xy)
        mask_xy = mask .* average_xy
        Plot_Contours(dataStruct, average_xy, TD_parameters, chain, 1, 2, 1, cmax)
        Plot_Contours(dataStruct, std_xy, TD_parameters, chain, 1, 2, 2, cmax)
        Plot_Contours(dataStruct, mask_xy, TD_parameters, chain, 1, 2, 3, cmax)
    end

    for i in 2:n
        if TD_parameters.xzMap == true
            voronoi_xz = vec(models[i].zeta_xz)
            p = contourf(vec(dataStruct.xVec), vec(dataStruct.zVec), voronoi_xz,linewidth=0.001, c=:jet,clims=(0,cmax),  yflip=true)
            ylims!(p, (0, 660))
            scatter!(p, dataStruct.elonsX, vec(dataStruct.edep), marker=:o, color=:lightblue, label="events", markersize=4,legend=false)    
            savefig(p, "./figures/xzVoronoi/chain#" * string(chain) * ('_') * string(i) * (".png"))

            all_xz = pmap(x -> models[x].zeta_xz, 1:i)
            average_xz = mean(all_xz)
            std_xz = std(all_xz)
            mask = ones(length(std_xz))
            for j = 1:length(std_xz)
                if std_xz[j] > max_std
                    mask[j] = NaN
                end
            end
            mask_xz = mask .* average_xz
            Plot_Contours(dataStruct, average_xz, TD_parameters, chain, i, 1, 1, cmax)
            Plot_Contours(dataStruct, std_xz, TD_parameters, chain, i, 1, 2, cmax)
            Plot_Contours(dataStruct, mask_xz, TD_parameters, chain, i, 1, 3, cmax)    
        end
        if TD_parameters.xyMap == true
            voronoi_xy = vec(models[i].zeta_xy)
            p = contourf(vec(dataStruct.xVec), vec(dataStruct.yVec), voronoi_xy, clims=(0,cmax), c=:jet)
            savefig(p, "./figures/xyVoronoi/chain#" * string(chain) * ('_') * string(i) * (".png"))

            all_xy = pmap(x -> models[x].zeta_xy, 1:i)
            average_xy = mean(all_xy)
            std_xy = std(all_xy)
            mask = ones(length(std_xy))
            for j = 1:length(std_xy)
                if std_xy[j] > max_std
                    mask[j] = NaN
                end
            end
            mask_xy = mask .* average_xy
            Plot_Contours(dataStruct, average_xy, TD_parameters, chain, i, 2, 1, cmax)
            Plot_Contours(dataStruct, std_xy, TD_parameters, chain, i, 2, 2, cmax)
            Plot_Contours(dataStruct, mask_xy, TD_parameters, chain, i, 2, 3, cmax)    
        end
    end



end

function Plot_Contours(
    dataStruct::MutabledataStruct, 
    plot_model, # ::Array{Float64,2}, 
    TD_parameters::parameters,
    chain::Int64,
    k::Int64,
    crosssectionType::Int64,
    type::Int64, # 1: model, 2: model uncertainty, 3: masked model
    cmax::Float64
    )

    ENV["GKSwstype"] = "nul"
    gr()

    closeenough = 2.0

    if TD_parameters.add_yVec == 0 # Tonga 2D
        println("2dTonga")
        tmpx = vcat(vec(dataStruct.dataX)', vec(dataStruct.elonsX)')
        tmpy = vcat(zeros(size(dataStruct.dataX))', vec(dataStruct.edep)') 
        tmpy = Array{Float64,2}(tmpy)
        p = contourf(vec(dataStruct.xVec), vec(dataStruct.zVec), vec(plot_model), xlabel="distance(km)", ylabel="depth(km)", yflip=true, clims=(0,cmax), c=:jet)

        title!("model_mean_2d_Tonga")
        savefig(p, "model_mean_2d_Tonga")

        title!("model_mean_2d_Tonga_0-300km")
        ylims!(p, (0, 300))
        savefig(p, "model_mean_2d_Tonga_0-300km") 

        ylims!(p, (0, 660))
        plot!(p, tmpx, tmpy, color=:white, legend=false)

        scatter!(p, dataStruct.dataX, zeros(size(dataStruct.dataX)), marker=:utri, color=:pink, label="station", markersize=6)
        scatter!(p, dataStruct.elonsX, vec(dataStruct.edep), marker=:o, color=:lightblue, label="events", markersize=4)
        title!("model_mean_2d_Tonga_ray")       
        savefig(p, "model_mean_2d_Tonga_ray")   

        title!("model_mean_2d_Tonga__0-300km_ray")
        ylims!(p, (0, 300))
        savefig(p, "model_mean_2d_Tonga__0-300km_ray") 

    else # 3d
        if crosssectionType == 1 # xz

            if type == 1
                appendname = "./figures/xzMap/model_mean_chain#" * string(chain) * "_" * string(k)
                titlename = "model_mean_chain#" * string(chain) * "_" * string(k)
            elseif type == 2
                appendname = "./figures/xzUncertainty/model_uncertainty_chain#" * string(chain) * "_" * string(k)
                titlename = "model_uncertainty_chain#" * string(chain) * "_" * string(k)
            elseif type == 3
                appendname = "./figures/xzMask/model_mask_chain#" * string(chain) * "_" * string(k)
                titlename = "model_mask_chain#" * string(chain) * "_" * string(k)
            end

            (m, n) = size(dataStruct.rayY)
            y0 = TD_parameters.ySlice
            nearrays = zeros(n)
            for i in 1:n
                yvec = dataStruct.rayY[:,i] .- y0    
                if sum(abs.(yvec)) - abs(sum(yvec)) > 1e-7
                    nearrays[i] = 1
                end  
                if min((abs.(yvec) .- closeenough)...) < 1e-7
                    nearrays[i] = 1
                end
            end
            if type == 2 # uncertainty map
                p = contourf(vec(dataStruct.xVec), vec(dataStruct.zVec), vec(plot_model), linewidth=0.001, xlabel="distance(km)", ylabel="depth(km)", yflip=true, c=:bone)
            else # model, masked model 
                p = contourf(vec(dataStruct.xVec), vec(dataStruct.zVec), vec(plot_model), linewidth=0.001, xlabel="distance(km)", ylabel="depth(km)", yflip=true, clims=(0,cmax), c=:jet)
            end 
            ylims!(p, (0, 660))
            scatter!(p, dataStruct.elonsX, vec(dataStruct.edep), marker=:o, color=:lightblue, label="events", markersize=4,legend=false)
            title!(titlename * "_xzMap")      
            savefig(p, appendname * "_xzMap")   
        end

        if crosssectionType == 2 # xy
            if type == 1
                appendname = "./figures/xyMap/model_mean_chain#" * string(chain) * "_" * string(k)
                titlename = "model_mean_chain#" * string(chain) * "_" * string(k)
            elseif type == 2
                appendname = "./figures/xyUncertainty/model_uncertainty_chain#" * string(chain) * "_" * string(k)
                titlename = "model_uncertainty_chain#" * string(chain) * "_" * string(k)
            elseif type == 3
                appendname = "./figures/xyMask/model_mask_chain#" * string(chain) * "_" * string(k)
                titlename = "model_mask_chain#" * string(chain) * "_" * string(k)
            end
            (m, n) = size(dataStruct.rayZ)
            z0 = TD_parameters.zSlice
            nearrays = zeros(n)
            for i in 1:n
                zvec = dataStruct.rayZ[:,i] .- z0    
                if sum(abs.(zvec)) - abs(sum(zvec)) > 1e-7
                    nearrays[i] = 1
                end  
                if min((abs.(zvec) .- closeenough)...) < 1e-7
                    nearrays[i] = 1
                end
            end
            nearrays = Array{Bool,1}(nearrays)
            tmpx = vcat(vec(dataStruct.dataX)', vec(dataStruct.elonsX)')
            tmpy = vcat(vec(dataStruct.dataY)', vec(dataStruct.elatsY)') 
            # tmpy = Array{Float64,2}(tmpy)
            nearraysX = tmpx[:,nearrays]
            nearraysY = tmpy[:,nearrays]
            
            if type == 2 #uncertainty map
                p = contourf(dataStruct.xVec, dataStruct.yVec, vec(plot_model), xlabel="X(km)", ylabel="Y(km)", c=:bone)
            else # model, masked model
                p = contourf(dataStruct.xVec, dataStruct.yVec, vec(plot_model), xlabel="X(km)", ylabel="Y(km)",clims=(0,cmax), c=:jet)
            end
            
            title!(titlename * "_xyMap" * string(TD_parameters.zSlice) * "km")        
            savefig(p, appendname * "_xyMap" * string(TD_parameters.zSlice) * "km")

            plot!(p, tmpx, tmpy, color=:white, legend=false)
            plot!(p, nearraysX, nearraysY, color=:forestgreen)
            scatter!(p, dataStruct.dataX, dataStruct.dataY, shape=:utri, color=:pink, label="station", markersize=6)
            scatter!(p, dataStruct.elonsX, dataStruct.elatsY, shape=:o, color=:lightblue, label="events", markersize=4)
            title!(titlename * "_xyMap_ray" * string(TD_parameters.zSlice) * "km")
        
            savefig(p, appendname * "_xyMap_ray" * string(TD_parameters.zSlice) * "km")       
        end

    end
end