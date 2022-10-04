struct Ray
    x; y; z
end

mutable struct MutabledataStruct
    tS          ::Array{Any,2}
    allLats     ::Array{Any,2}      #Latitude for the stations in each event-station pair
    allLons     ::Array{Any,2}      #Longitude for the stations in each event-station pair
    allSig      ::Array{Any,2}      #error terms from t* inversion
    allSta      ::Array{Any,2}
    dataX       ::Array{Float64,2}  #Station position in each event-station pair
    dataY       ::Array{Float64,2}
    xVec        ::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    yVec        ::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    zVec        ::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    elonsX      ::Array{Float64,2}
    elatsY      ::Array{Float64,2}
    elons       ::Array{Any,2}
    elats       ::Array{Any,2}
    edep        ::Array{Any,2}
    coastX      ::Array{Float64,2}
    coastY      ::Array{Float64,2}
    # ray         ::Array{Any,1}
    rayX        ::Array{Float64,2}
    rayY        ::Array{Float64,2}
    rayZ        ::Array{Float64,2}
    rayL
    rayU
    U           ::Array{Float64,2}
end

mutable struct Model
    nCells
    xVec
    yVec
    zVec
    xCell
    yCell
    zCell
    zeta
    zeta_scale
    layer
    # allSig
    phi
    ptS
    tS
    likelihood
    action          ::Int64
    accept          ::Int64
    zeta_xz
    zeta_xy
end