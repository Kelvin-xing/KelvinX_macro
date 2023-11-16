# Initial BGP calibration: a struct for targeted moments
@with_kw struct BgpTargets
    growth::Float64       # Growth rate
    entry::Float64        # Firm entry rate
    rnd::Float64          # R&D to GDP ratio
    profit::Float64       # Profit share
    markup::Float64       # Average markup
    gd::Float64           # Dispersion of firm growth rates
    ec::Float64           # Net entry contribution to growth
end

# Transition calibration: a struct for targeted moments
@with_kw struct TransitionTargets
    entry::Float64              # Firm entry rate in 2015
    gd::Float64                 # Dispersion of firm growth rates in 2015
    markup::Float64             # Average markup in 2015
    entrydeclinein10::Float64   # Decline in firm entry rate in first 10 years
    entrydeclinein20::Float64   # Decline in firm entry rate in first 20 years
end

# It performs a transformation on parameters (For calibration purposes)
function TransformParameters(untransformedguess::Vector{Float64}, ranges::Vector{Tuple{Float64,Float64}})

    transformedguess = similar(untransformedguess)

    @inbounds for i = eachindex(untransformedguess)

        untransformedvariable = untransformedguess[i]
        range = ranges[i]
        low = range[1]
        high = range[2]

        transformedguess[i] = log((untransformedvariable - low) / (high - untransformedvariable))
    end

    return transformedguess
end

# It untransforms already transformed parameters (For calibration purposes)
function UntransformParameters(transformedguess::Vector{Float64}, ranges::Vector{Tuple{Float64,Float64}})

    untransformedguess = similar(transformedguess)

    @inbounds for i = eachindex(transformedguess)

        transformedvariable = transformedguess[i]
        range = ranges[i]
        low = range[1]
        high = range[2]

        untransformedguess[i] = (exp(transformedvariable) * high + low) / (1 + exp(transformedvariable))
    end

    return untransformedguess
end

# This function calculates the difference between model generated moments and targets for initial bgp calibration
function DiffCalc(moments::BgpTargets, targets::BgpTargets, weights::BgpTargets, momentnames::NTuple{7, Symbol})
    diff = 0.0

    for i = eachindex(momentnames)
        moment = momentnames[i]

        diff += getproperty(weights, moment) * (100 * (getproperty(moments, moment) - getproperty(targets, moment)) / ((getproperty(moments, moment) + getproperty(targets, moment)) / 2.0))^2
    end

    return diff
end

# This function calculates the difference between model generated moments and targets for transition calibration
function DiffCalc(moments::TransitionTargets, targets::TransitionTargets, weights::TransitionTargets, momentnames::NTuple{5, Symbol})

    diff = 0.0

    for i = eachindex(momentnames)

        moment = momentnames[i]

        diff += getproperty(weights, moment) * (100 * (getproperty(moments, moment) - getproperty(targets, moment)) / ((getproperty(moments, moment) + getproperty(targets, moment)) / 2.0))^2
    end

    return diff
end

# Initial BGP calibration. This function calculates the difference between model generated moments and targets.
# rangeall: bounds on parameter values. sectordata and panel are preallocated. momentnames specify what moments are targeted for the calibration. 
function BgpCalibration!(transformedparameters::Vector{Float64}, rangeall::Vector{Tuple{Float64,Float64}}, numericp::NumParam, targets::BgpTargets, sectordata::Matrix{Sector}, panel::Matrix{Float64}, momentnames::NTuple{7, Symbol}, momentweights::BgpTargets, τ::Float64, s::Float64)

    # Untransform parameter guesses
    untransformedparameters = UntransformParameters(transformedparameters, rangeall)

    # Create bgp structural parameters
    p = StrParam(τ = τ, s = s, α = untransformedparameters[1], αTilde = untransformedparameters[2], λ = untransformedparameters[3], δ = untransformedparameters[4], ϕ = untransformedparameters[5], β = untransformedparameters[6], ψ = untransformedparameters[7])

    # Solve for bgp
    initialbgp = BgpSolver(p, numericp)

    # Simulate bgp
    InitialBgpSimulation!(panel, sectordata, initialbgp)

    # Calculate moments for calibration
    moments = BgpCalibrationMoments(initialbgp, panel)

    return DiffCalc(moments, targets, momentweights, momentnames)
end

# Initial BGP calibration: This function minimizes BgpCalibration! function defined above in order to find initial bgp parameters
function BgpCalibration(numericp::NumParam, targets::BgpTargets, guess::Vector{Float64}, iterations::Int)

    panel, sectordata = SimPrep(numericp)

    αrange = (0.0001, 5.0000)
    αTilderange = (0.0001, 5.0000)
    λrange = (1.0001, 1.1000)
    δrange = (0.0001, 1.0000)
    ϕrange = (0.0001, 1.0000)
    βrange = (0.9000, 0.9999)
    ψrange = (0.0001, 1.0000)

    rangeall = [αrange, αTilderange, λrange, δrange, ϕrange, βrange, ψrange]

    momentnames = propertynames(targets)

    momentweights0 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    momentweights0 = momentweights0 ./ sum(momentweights0)
    momentweights = BgpTargets(growth = momentweights0[1], entry = momentweights0[2], rnd = momentweights0[3], profit = momentweights0[4], markup = momentweights0[5], gd = momentweights0[6], ec = momentweights0[7])

    transformedguess = TransformParameters(guess, rangeall)

    calibration = optimize(transformedparameters -> BgpCalibration!(transformedparameters, rangeall, numericp, targets, sectordata, panel, momentnames, momentweights, 0.30, 0.05), transformedguess, NelderMead(), Optim.Options(iterations = iterations, show_trace = true, store_trace = false))

    calibrated_param = UntransformParameters(calibration.minimizer, rangeall)

    return calibrated_param, calibration
end

# Transition calibration. This function calculates the difference between model generated moments and targets. For full calibration, we need to minimize this function.
# rangeall: bounds on parameter values. initialbgp: already calibrated initial bgp. entry1980: model generated (from initialbgp) firm entry rate in 1980. τpol and spol are 2015 values of corporate taxes and subsidy rates. panel and sectordata are preallocated. panelinitialbgplastyear, sectordatainitialbgplastperiod, and entrantid come from initial bgp calibration. 
function TransitionCalibration!(transformedparameters::Vector{Float64}, rangeall::Vector{Tuple{Float64,Float64}}, numericp::NumParam, targets::TransitionTargets, initialbgp::BGP, entry1980::Float64, momentnames::NTuple{5, Symbol}, momentweights::TransitionTargets, τpol::Float64, spol::Float64, panel::Matrix{Float64}, sectordata::Matrix{Sector}, entrantid::Int, panelinitialbgplastyear::Matrix{Float64}, sectordatainitialbgplastperiod::Vector{Sector})

    # Untransform parameter guesses
    untransformedparameters = UntransformParameters(transformedparameters, rangeall)

    # Create transition parameters
    transitionp = TransParam(τ = τpol/initialbgp.p.τ, s = spol/initialbgp.p.s, ντ = 1.0, νs = 1.0, αTilde = untransformedparameters[1], δ = untransformedparameters[2], ναTilde = untransformedparameters[3], νδ = untransformedparameters[4])

    # Solve for transition
    transition, _ = TransitionSolver(initialbgp, transitionp, numericp)

    # Simulate transition
    TransitionSimulation!(panel, sectordata, transition, entrantid, panelinitialbgplastyear, sectordatainitialbgplastperiod)
    
    # Calculate moments for calibration
    moments = TransitionCalibrationMoments(transition, numericp, entry1980, panel)

    return DiffCalc(moments, targets, momentweights, momentnames)
end

# Transition calibration: This function minimizes TransitionCalibration! function defined above in order to find transition parameters
# We need to start this function with an initialbgp object. targets are moment targets, and guess is the set of parameters that initiates minimization routine. entry1980 is the calibrated value of entry rate in initial bgp. We need this to calculate declines in entry rate over time. The remaining arguments entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod are used to sustain consistency between initial bgp simulations and transition simulations. 
function TransitionCalibration(initialbgp::BGP, targets::TransitionTargets, guess::Vector{Float64}, iterations::Int, entry1980::Float64, entrantid_initialbgp::Int, panelinitialbgplastyear::Matrix{Float64}, sectordatainitialbgplastperiod::Vector{Sector})
    
    @unpack numericp = initialbgp

    panel, sectordata = SimPrep(numericp)

    αTilderange = (1.0, 6.0)
    δrange = (0.01, 1.0)
    ναTilderange = (0.01, 10.0)
    νδrange = (0.01, 10.0)

    rangeall = [αTilderange, δrange, ναTilderange, νδrange]

    momentnames = propertynames(targets)

    momentweights0 = [1.0, 1.0, 1.0, 1.0, 1.0]
    momentweights0 = momentweights0 ./ sum(momentweights0)
    momentweights = TransitionTargets(entry = momentweights0[1], gd = momentweights0[2], markup = momentweights0[3], entrydeclinein10 = momentweights0[4], entrydeclinein20 = momentweights0[5])

    transformedguess = TransformParameters(guess, rangeall)

    calibration = optimize(transformedparameters -> TransitionCalibration!(transformedparameters, rangeall, numericp, targets, initialbgp, entry1980, momentnames, momentweights, 0.20, 0.20, panel, sectordata, entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod), transformedguess, NelderMead(), Optim.Options(iterations = iterations, show_trace = true, store_trace = false))

    transitionp_calibrated = UntransformParameters(calibration.minimizer, rangeall)

    return transitionp_calibrated, calibration
end