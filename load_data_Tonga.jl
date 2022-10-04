using MAT,DelimitedFiles,Interpolations,Plots,JLD
include("MCsub.jl")

function load_data_Tonga(TD_parameters::parameters)
    allLats, allLons = [], []
    allTS, allSig = [], []
    dataE     = []
    elats, elons, edep = [], [], []
    allSta    = []
    ray       = []

    #############add sth later if I would like to read the data from 3D ray tracing result file########
    traces    = load("./Data/381traces.jld")

    elats = traces["EventLatitude"]
    elons = traces["EventLongitude"]
    edep = traces["EventDepth"]

    allSta = traces["station"]
    allLats = traces["latitude"]
    allLons = traces["longitude"]
    allTS = traces["tStar"]

    # allSig  = 0.2 * ones(size(allTS))
    allSig  = traces["error"]


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
    coastlines = load("./Data/coastlines.jld")

    coastlon = coastlines["coastlon"]
    coastlat = coastlines["coastlat"]
    (coastX, coastY) = lonlat2xy(lon0, lat0, beta, coastlon, coastlat)

    minX = min(dataX[:,1]...) - TD_parameters.buffer
    maxX = max(dataX[:,1]...) + TD_parameters.buffer
    minY = min(dataY[:,1]...) - TD_parameters.buffer
    maxY = max(dataY[:,1]...) + TD_parameters.buffer
    # scatter(vec(dataX),vec(dataY))
    # scatter!([0],[0])
    # scatter!(vec([minX minX maxX maxX]),vec([minY maxY minY maxY]),markersize = 6)
    # title!("station distribution after rotation (in xy coordinate) with buffer = 100 km")
    # savefig("station_after_buffer.png")
    xVec = minX:TD_parameters.XYnodeSpacing:maxX
    yVec = minY:TD_parameters.XYnodeSpacing:maxY
    zVec = TD_parameters.min_depth:TD_parameters.ZnodeSpacing:TD_parameters.max_depth

    ak135mod = readdlm("./Data/ak135f.txt", ',', Float64)
    # vp=LinearInterpolation(vp[:,1],vp[:,2]; xVec)
    # scatter(vp[:,1],vp[:,2])
    # plot!(vp)
    # savefig()
    # vp = interpolate(vp[:,1:2],BSpline(Linear()))
    # vp(6200,2)
    # etp = LinearInterpolation(xVec, vp[:,2]; extrapolation_bc=Throw())
    # vp = interp1(ak135mod[:,1], ak135mod[:,2], zVec)

    LAB = load("./Data/LAB_discontinuity.jld")

    GLON = LAB["GLON"]
    GLAT = LAB["GLAT"]
    GDEP = LAB["GDEP"]

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
    raypath=load("./Data/381raypaths.jld")
    x=raypath["x_n"]
    y=raypath["y_n"]
    z=raypath["z_n"]

    for i in 1:length(allTS)
        elat = elats[i]
        elon = elons[i] 

        # zf = -a *ln(r/a)
        eidep = edep[i]
        vp = interp1(ak135mod[:,1], ak135mod[:,2], z[i,:])

        # Earth flattening transformation, (4.49) Peter Shearer
        vp = 6371/(6371-eidep) .* vp

        # iray=Ray(x[i,:],y[i,:],z[i,:],ones(1,67),ones(1,67),
        # 1.0 ./vp,zeros(size(vp)),zeros(size(vp)))
        iray=Ray(x[i,:],y[i,:],z[i,:])
        ray=[ray;iray]
    end

    
    x = [];y = [];z = [];U = []
    for k in 1:length(ray)
        x1, y1, z1 = ray[k].x, ray[k].y, ray[k].z
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
   ###########ATTENTION! 最后的discontinuity～～～～rayz raylat rayX啥的没加

    dataStruct = MutabledataStruct(
        allTS, allLats, allLons, allSig, allSta,
        dataX, dataY,
        xVec, yVec, zVec,
        elonsX, elatsY,
        elons, elats, edep,
        coastX, coastY,
        # ray, 
        x, y, z, 
        0, 0, U
    )
    return dataStruct
end