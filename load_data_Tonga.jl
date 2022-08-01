using MAT,DelimitedFiles,Interpolations,Plots
include("MCsub.jl")
struct Ray
    x
    y
    z
    lat
    lon
    U
    zeta
    id
end
function load_data_Tonga(TD_parameters::Dict{String,Any})
    allLats   = []
    allLons   = []
    allTS     = []
    allSig    = []
    dataE     = []
    elats     = []
    elons     = []
    edep      = []
    allSta    = []
    ray       = []
    #############If the data are stored in a .mat file##############
    #############add sth later if I would like to read the data from 3D ray tracing result file########
    traces    = matread("./Data/aaa.mat")

    elats = traces["ts_run"]["EventLatitude"]
    elons = traces["ts_run"]["EventLongitude"]
    edep = traces["ts_run"]["EventDepth"]

    allSta = traces["ts_run"]["station"]
    allLats = traces["ts_run"]["latitude"]
    allLons = traces["ts_run"]["longitude"]
    allTS = traces["ts_run"]["tStar"]

    # allSig  = 0.2 * ones(size(allTS))
    allSig  = traces["ts_run"]["error"]


    ##############unique dataE hasn't been considered#############

    ##############define the origin of the map and convert data from latlon to xy##########
    lat0 = -18.1404
    lon0 = 177.0798
    beta = 0.463647609

    # lat0 = -23.1000
    # lon0 = 174.6000
    # beta = 0.463647609
    
    # scatter(vec(allLons),vec(allLats))
    # scatter!([lon0],[lat0])
    # title!("station distribution before rotation (in latlon coordinate)")
    # savefig("station.png")


    (dataX, dataY) = lonlat2xy(lon0, lat0, beta, allLons, allLats)
    (elonsX, elatsY) = lonlat2xy(lon0, lat0, beta, elons, elats)

    # scatter(vec(dataX),vec(dataY))
    # scatter!([0],[0])
    # title!("station distribution after rotation (in xy coordinate)")
    # savefig("station_after.png")

    #############Read coastline data###################
    coastlines = matread("./Data/coastlines.mat")
    
    coastlon = coastlines["coastlon"]
    coastlat = coastlines["coastlat"]
    (coastX, coastY) = lonlat2xy(lon0, lat0, beta, coastlon, coastlat)

    minX = min(dataX[:,1]...) - TD_parameters["buffer"]
    maxX = max(dataX[:,1]...) + TD_parameters["buffer"]
    minY = min(dataY[:,1]...) - TD_parameters["buffer"]
    maxY = max(dataY[:,1]...) + TD_parameters["buffer"]
    # scatter(vec(dataX),vec(dataY))
    # scatter!([0],[0])
    # scatter!(vec([minX minX maxX maxX]),vec([minY maxY minY maxY]),markersize = 6)
    # title!("station distribution after rotation (in xy coordinate) with buffer = 100 km")
    # savefig("station_after_buffer.png")
    xVec = minX:TD_parameters["XYnodeSpacing"]:maxX
    yVec = minY:TD_parameters["XYnodeSpacing"]:maxY
    zVec = TD_parameters["min_depth"]:TD_parameters["ZnodeSpacing"]:TD_parameters["max_depth"]

    
    dataStruct = Dict(
        "tS" => allTS,
        "allLats" => allLats,    #Latitude for the stations in each event-station pair
        "allLons" => allLons,    #Longitude for the stations in each event-station pair
        "allSig" => allSig,      #ATTENTION!!! remain unknown!
        "allSta" => allSta,
        "dataX" => dataX,        #Station position in each event-station pair
        "dataY" => dataY,
        # "dataE" => dataE,
        "xVec" => xVec,
        "yVec" => yVec,
        "zVec" => zVec,
        "elonsX" => elonsX,
        "elatsY" => elatsY,
        "elons" => elons,
        "elats" => elats,
        "edep" => edep,
        "coastX" => coastX,
        "coastY" => coastY
    )


    ak135mod = readdlm("./Data/ak135f.txt", ',', Float64)
    # vp=LinearInterpolation(vp[:,1],vp[:,2]; xVec)
    # scatter(vp[:,1],vp[:,2])
    # plot!(vp)
    # savefig()
    # vp = interpolate(vp[:,1:2],BSpline(Linear()))
    # vp(6200,2)
    # etp = LinearInterpolation(xVec, vp[:,2]; extrapolation_bc=Throw())
    # vp = interp1(ak135mod[:,1], ak135mod[:,2], zVec)

    LAB = matread("./Data/LAB_discontinuity.mat")

    GLON = LAB["GLON"]
    GLAT = LAB["GLAT"]
    GDEP = LAB["GDEP"]

    #############################
    # nray=0
    # npt=0
    # x=zeros(200,1000)
    # y=zeros(200,1000)
    # z=zeros(200,1000)
    # for line in readlines("./Data/raypaths.p")
    #     if split(line)[1] == "1234567"
    #         nray=nray+1
    #         npt=0
    #     else
    #         npt=npt+1
    #         x[nray,npt]=parse(Float64,split(line)[1])
    #         y[nray,npt]=parse(Float64,split(line)[2])
    #         z[nray,npt]=parse(Float64,split(line)[3])
    #     end
    # end
    raypath=matread("./Data/raypathsp_nan_New_origin.mat")
    x=raypath["x_n"]
    y=raypath["y_n"]
    z=raypath["z_n"]

    for i in 1:length(dataStruct["tS"])
        elat = elats[i]
        elon = elons[i] 

        # zf = -a *ln(r/a)
        eidep = edep[i]

        vp = interp1(ak135mod[:,1], ak135mod[:,2], z[i,:])

        # Earth flattening transformation, (4.49) Peter Shearer
        vp = 6371/(6371-eidep) .* vp

        iray=Ray(x[i,:],y[i,:],z[i,:],ones(1,67),ones(1,67),
        1.0 ./vp,zeros(size(vp)),zeros(size(vp)))
        ray=[ray;iray]
    end

    

    setindex!(dataStruct,ray,"ray")

    x = [];y = [];z = [];U = []
    for k in 1:length(dataStruct["ray"])
        x1 = dataStruct["ray"][k].x
        y1 = dataStruct["ray"][k].y
        z1 = dataStruct["ray"][k].z
        # ATTENTION! why don't we use rayz? : related to rayL, * rayU, rayU is defined in zVec
        # z1 = dataStruct["zVec"]
        # U1 = dataStruct["ray"][k].U
        vp1 = interp1(ak135mod[:,1], ak135mod[:,2], z1)
        U1 = 1.0./vp1
        # U1 = dataStruct["ray"][k].z
        if k == 1
            x = x1;y = y1;z = z1;U = U1
        else
            x = [x x1];y = [y y1];z = [z z1];U = [U U1]
        end
    end
    setindex!(dataStruct, x, "rayX")
    setindex!(dataStruct, y, "rayY")
    setindex!(dataStruct, z, "rayZ")
    setindex!(dataStruct, U, "U")

    ###########ATTENTION! discontinuity～～～～rayz raylat rayX

end