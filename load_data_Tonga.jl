using DelimitedFiles,Interpolations,JLD
include("MCsub.jl")

function load_data_Tonga(TD_parameters::parameters)
    allLats, allLons    = [], [] #station latitudes and longitudes
    allTS, allSig, allaveatten  = [], [], []
    elats, elons, edep          = [], [], [] #event latitudes and longitudes
    
    traces    = load("./Data/traces.jld")
    elats       = traces["EventLatitude"]
    elons       = traces["EventLongitude"]
    edep        = traces["EventDepth"]
    allLats     = traces["latitude"]
    allLons     = traces["longitude"]
    allTS       = traces["tStar"]
    allaveatten = traces["aveatten"]
    allSig      = traces["error"]

    #drop to singleton dimensions
    allTS       = Vector{Float64}(vec(allTS))
    allLats     = Vector{Float64}(vec(allLats))
    allLons     = Vector{Float64}(vec(allLons))
    allaveatten = Vector{Float64}(vec(allaveatten))
    allSig      = Vector{Float64}(vec(allSig))

    lat0 = -23.1000
    lon0 = 174.6000
    beta = 0.463647609

    #stations and events coordicates in cartesian system
    (dataX, dataY)      = lonlat2xy(lon0, lat0, beta, allLons, allLats)
    (elonsX, elatsY)    = lonlat2xy(lon0, lat0, beta, elons, elats)

    #load coastline data
    #????not used yet
    coastlines = load("./Data/coastlines.jld")
    coastlon = coastlines["coastlon"]
    coastlat = coastlines["coastlat"]
    (coastX, coastY) = lonlat2xy(lon0, lat0, beta, coastlon, coastlat)

    #study area
    minX = minimum(dataX) - TD_parameters.buffer
    maxX = maximum(dataX) + TD_parameters.buffer
    minY = minimum(dataY) - TD_parameters.buffer
    maxY = maximum(dataY) + TD_parameters.buffer

    xVec = minX:TD_parameters.XYnodeSpacing:maxX
    yVec = minY:TD_parameters.XYnodeSpacing:maxY
    zVec = TD_parameters.min_depth:TD_parameters.ZnodeSpacing:TD_parameters.max_depth

    #load LAB discontinuity
    #??? not used yet
    LAB = load("./Data/LAB_discontinuity.jld")
    GLON = LAB["GLON"]
    GLAT = LAB["GLAT"]
    GDEP = LAB["GDEP"]

    #load raypaths
    raypath=load("./Data/raypaths.jld")
    x = raypath["x"]
    y = raypath["y"]
    z = raypath["z"]
    U = raypath["u"]

    #raylength and slowness for each segment
    rayl = sqrt.((x[1:end - 1,:] - x[2:end,:]).^2 +
    (y[1:end - 1,:] - y[2:end,:]).^2 +
    (z[1:end - 1,:] - z[2:end,:]).^2)
    rayu = 0.5 .* (U[1:end - 1,:] + U[2:end,:])

    dataStruct = DataStruct(
        allTS, allaveatten, allLats, allLons, allSig, 
        dataX, dataY,
        xVec, yVec, zVec,
        elonsX, elatsY,
        elons, elats, edep,
        coastX, coastY,
        # ray, 
        x, y, z, 
        rayl, rayu, U
    )

    return dataStruct
end

function load_synthetic_data_Tonga(TD_parameters::parameters)
    allLats, allLons    = [], [] #station latitudes and longitudes
    allTS, allSig, allaveatten  = [], [], []
    elats, elons, edep          = [], [], [] #event latitudes and longitudes
    
    traces    = load("./Data/synthetic_traces.jld")
    elats       = traces["EventLatitude"]
    elons       = traces["EventLongitude"]
    edep        = traces["EventDepth"]
    allLats     = traces["latitude"]
    allLons     = traces["longitude"]
    allTS       = traces["tStar"]
    allaveatten = traces["aveatten"]
    allSig      = traces["error"]

    #drop to singleton dimensions
    allTS       = Vector{Float64}(vec(allTS))
    allLats     = Vector{Float64}(vec(allLats))
    allLons     = Vector{Float64}(vec(allLons))
    allaveatten = Vector{Float64}(vec(allaveatten))
    allSig      = Vector{Float64}(vec(allSig))

    lat0 = -23.1000
    lon0 = 174.6000
    beta = 0.463647609

    #stations and events coordicates in cartesian system
    (dataX, dataY)      = lonlat2xy(lon0, lat0, beta, allLons, allLats)
    (elonsX, elatsY)    = lonlat2xy(lon0, lat0, beta, elons, elats)

    #load coastline data
    #????not used yet
    coastlines = load("./Data/coastlines.jld")
    coastlon = coastlines["coastlon"]
    coastlat = coastlines["coastlat"]
    (coastX, coastY) = lonlat2xy(lon0, lat0, beta, coastlon, coastlat)

    #study area
    minX = minimum(dataX) - TD_parameters.buffer
    maxX = maximum(dataX) + TD_parameters.buffer
    minY = minimum(dataY) - TD_parameters.buffer
    maxY = maximum(dataY) + TD_parameters.buffer

    xVec = minX:TD_parameters.XYnodeSpacing:maxX
    yVec = minY:TD_parameters.XYnodeSpacing:maxY
    zVec = TD_parameters.min_depth:TD_parameters.ZnodeSpacing:TD_parameters.max_depth

    #load LAB discontinuity
    #??? not used yet
    LAB = load("./Data/LAB_discontinuity.jld")
    GLON = LAB["GLON"]
    GLAT = LAB["GLAT"]
    GDEP = LAB["GDEP"]

    #load raypaths
    raypath=load("./Data/synthetic_raypaths.jld")
    x = raypath["x"]
    y = raypath["y"]
    z = raypath["z"]
    U = raypath["u"]

    #raylength and slowness for each segment
    rayl = sqrt.((x[1:end - 1,:] - x[2:end,:]).^2 +
    (y[1:end - 1,:] - y[2:end,:]).^2 +
    (z[1:end - 1,:] - z[2:end,:]).^2)
    rayu = 0.5 .* (U[1:end - 1,:] + U[2:end,:])

    dataStruct = DataStruct(
        allTS, allaveatten, allLats, allLons, allSig, 
        dataX, dataY,
        xVec, yVec, zVec,
        elonsX, elatsY,
        elons, elats, edep,
        coastX, coastY,
        # ray, 
        x, y, z, 
        rayl, rayu, U
    )

    return dataStruct
end