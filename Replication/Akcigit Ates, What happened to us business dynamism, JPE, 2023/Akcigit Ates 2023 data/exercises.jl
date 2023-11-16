# Section 5 - Understanding mechanisms in isolation: Calibrate αTilde to match the entry profile
function TransitionCalibration_alphaTilde!(transformedparameters::Vector{Float64}, rangeall::Vector{Tuple{Float64,Float64}}, numericp::NumParam, initialbgp::BGP, entry1980::Float64, entry2015_target::Float64, entrydeclinein10_target::Float64, entrydeclinein20_target::Float64)

    # Untransform x
    untransformedparameters = UntransformParameters(transformedparameters, rangeall)

    # Create transition parameters
    transitionp = TransParam(τ = 1.0, s = 1.0, αTilde = untransformedparameters[1], ναTilde = untransformedparameters[2], δ = 1.0, νδ = 1.0)

    # Solve for transition
    transition, _ = TransitionSolver(initialbgp, transitionp, numericp)

    # Calculate moments for calibration
    period2015 = numericp.tperiodpol
    period1990 = (1990 - 1980) * numericp.periodsperyear
    period2000 = (2000 - 1980) * numericp.periodsperyear

    entry2015 = EntryRate(transition.xe[:, period2015], transition.xenn[period2015], transition.μ[:,period2015], transition.μnn[period2015])
    entry1990 = EntryRate(transition.xe[:, period1990], transition.xenn[period1990], transition.μ[:,period1990], transition.μnn[period1990])
    entry2000 = EntryRate(transition.xe[:, period2000], transition.xenn[period2000], transition.μ[:,period2000], transition.μnn[period2000])
    
    entrydeclinein10 = 100.0 * (1.0 - entry1990/entry1980)
    entrydeclinein20 = 100.0 * (1.0 - entry2000/entry1980)

    diff = (100 * (entry2015 - entry2015_target) / ((entry2015 + entry2015_target) / 2.0))^2 + (100 * (entrydeclinein10 - entrydeclinein10_target) / ((entrydeclinein10 + entrydeclinein10_target) / 2.0))^2 + (100 * (entrydeclinein20 - entrydeclinein20_target) / ((entrydeclinein20 + entrydeclinein20_target) / 2.0))^2

    return diff
end

# Section 5 - Understanding mechanisms in isolation: Calibrate δ to match the entry profile
function TransitionCalibration_delta!(transformedparameters::Vector{Float64}, rangeall::Vector{Tuple{Float64,Float64}}, numericp::NumParam, initialbgp::BGP, entry1980::Float64, entry2015_target::Float64, entrydeclinein10_target::Float64, entrydeclinein20_target::Float64)

    # Untransform x
    untransformedparameters = UntransformParameters(transformedparameters, rangeall)

    # Create transition parameters
    transitionp = TransParam(τ = 1.0, s = 1.0, αTilde = 1.0, ναTilde = 1.0, δ = untransformedparameters[1], νδ = untransformedparameters[2])

    # Solve for transition
    transition, _ = TransitionSolver(initialbgp, transitionp, numericp)

    # Calculate moments for calibration
    period2015 = numericp.tperiodpol
    period1990 = (1990 - 1980) * numericp.periodsperyear
    period2000 = (2000 - 1980) * numericp.periodsperyear

    entry2015 = EntryRate(transition.xe[:, period2015], transition.xenn[period2015], transition.μ[:,period2015], transition.μnn[period2015])
    entry1990 = EntryRate(transition.xe[:, period1990], transition.xenn[period1990], transition.μ[:,period1990], transition.μnn[period1990])
    entry2000 = EntryRate(transition.xe[:, period2000], transition.xenn[period2000], transition.μ[:,period2000], transition.μnn[period2000])
    
    entrydeclinein10 = 100.0 * (1.0 - entry1990/entry1980)
    entrydeclinein20 = 100.0 * (1.0 - entry2000/entry1980)

    diff = (100 * (entry2015 - entry2015_target) / ((entry2015 + entry2015_target) / 2.0))^2 + (100 * (entrydeclinein10 - entrydeclinein10_target) / ((entrydeclinein10 + entrydeclinein10_target) / 2.0))^2 + (100 * (entrydeclinein20 - entrydeclinein20_target) / ((entrydeclinein20 + entrydeclinein20_target) / 2.0))^2

    return diff
end

# Section 5 - Understanding mechanisms in isolation: Following function first calibrates entry cost and diffusion channels to match declining entry profile. Then for each channel, transitions and corresponding moments are calculated. Finally, Figure 3, Figure 4, and Table 3 are returned. In matrix table3 numbers -1, 0, and 1 represent down, horizontal, and up arrows, respectively.
function Section5WithCalibration(initialbgp::BGP, initialbgp_moments::BgpMoments, numericp::NumParam, entrantid_initialbgp::Int, panelinitialbgplastyear::Matrix{Float64}, sectordatainitialbgplastperiod::Vector{Sector}, guess_αTilde::Vector{Float64}, guess_δ::Vector{Float64}, iterations::Int)
    
    @unpack_NumParam numericp
    entry1980 = initialbgp_moments.entry

    # Calibration of αTilde and ναTilde to match the declining profile of entry over the transition
    αTilderange = (1.0, 6.0)
    ναTilderange = (0.01, 10.0)

    rangeall = [αTilderange, ναTilderange]

    transformedguess = TransformParameters(guess_αTilde, rangeall)

    αTilde_calibration = optimize(transformedparameters -> TransitionCalibration_alphaTilde!(transformedparameters, rangeall, numericp, initialbgp, entry1980, 8.0, 15.0, 25.0), transformedguess, NelderMead(), Optim.Options(iterations = iterations, show_trace = true, store_trace = false))

    αTilde_calibrated, ναTilde_calibrated = UntransformParameters(αTilde_calibration.minimizer, rangeall)

    # Calibration of δ and νδ to match the declining profile of entry over the transition
    δrange = (0.01, 1.0)
    νδrange = (0.01, 10.0)

    rangeall = [δrange, νδrange]

    transformedguess = TransformParameters(guess_δ, rangeall)

    δ_calibration = optimize(transformedparameters -> TransitionCalibration_delta!(transformedparameters, rangeall, numericp, initialbgp, entry1980, 8.0, 15.0, 25.0), transformedguess, NelderMead(), Optim.Options(iterations = iterations, show_trace = true, store_trace = false))

    δ_calibrated, νδ_calibrated = UntransformParameters(δ_calibration.minimizer, rangeall)

    # Create figures and the table
    figure3, table3, figure4 = Section5(initialbgp, initialbgp_moments, numericp, entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod, αTilde_calibrated, ναTilde_calibrated, δ_calibrated, νδ_calibrated)

    return figure3, table3, figure4
end

# Section 5 - Understanding mechanisms in isolation: Contrary to above function, following function does not calibrate entry cost and diffusion channels. For given calibrated values, it calculates transitions and corresponding moments for each channel. Finally, Figure 3, Figure 4, and Table 3 are returned. In matrix table3 numbers -1, 0, and 1 represent down, horizontal, and up arrows, respectively.
function Section5(initialbgp::BGP, initialbgp_moments::BgpMoments, numericp::NumParam, entrantid_initialbgp::Int, panelinitialbgplastyear::Matrix{Float64}, sectordatainitialbgplastperiod::Vector{Sector}, αTilde_calibrated::Float64, ναTilde_calibrated::Float64, δ_calibrated::Float64, νδ_calibrated::Float64)
    
    @unpack_NumParam numericp

    # Solve four different transitions corresponding to each channel: corporate tax, R&D subsidy, entry cost, diffusion
    panel, sectordata = SimPrep(numericp)

    # Corporate tax channel: τ is reduced to 0% as explained in the paper
    transitionp_τ = TransParam(τ = 0, s = 1, αTilde = 1, δ = 1, ντ = 1, νs = 1, ναTilde = 1, νδ = 1)
    transition_τ, _ = TransitionSolver(initialbgp, transitionp_τ, numericp)
    TransitionSimulation!(panel, sectordata, transition_τ, entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod)
    transition_moments_τ = Moments(transition_τ, panel, sectordata)

    # R&D subsidy channel: s is increased to 50% as explained in the paper
    transitionp_s = TransParam(τ = 1, s = 0.50/0.05, αTilde = 1, δ = 1, ντ = 1, νs = 1, ναTilde = 1, νδ = 1)
    transition_s, _ = TransitionSolver(initialbgp, transitionp_s, numericp)
    TransitionSimulation!(panel, sectordata, transition_s, entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod)
    transition_moments_s = Moments(transition_s, panel, sectordata)

    # Entry cost channel: αTilde and ναTilde were calibrated to match entry rate decline. We use these calibrated values
    transitionp_αTilde = TransParam(τ = 1, s = 1, αTilde = αTilde_calibrated, δ = 1, ντ = 1, νs = 1, ναTilde = ναTilde_calibrated, νδ = 1)
    transition_αTilde, _ = TransitionSolver(initialbgp, transitionp_αTilde, numericp)
    TransitionSimulation!(panel, sectordata, transition_αTilde, entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod)
    transition_moments_αTilde = Moments(transition_αTilde, panel, sectordata)

    # Diffusion channel: δ and νδ were calibrated to match entry rate decline. We use these calibrated values
    transitionp_δ = TransParam(τ = 1, s = 1, αTilde = 1, δ = δ_calibrated, ντ = 1, νs = 1, ναTilde = 1, νδ = νδ_calibrated)
    transition_δ, _ = TransitionSolver(initialbgp, transitionp_δ, numericp)
    TransitionSimulation!(panel, sectordata, transition_δ, entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod)
    transition_moments_δ = Moments(transition_δ, panel, sectordata)

    # Table 3 in the paper
    momentnames = [:conc, :markup, :profit, :labor, :entry, :yes, :prodgap, :jr, :gd]
    channelnames = [:τ, :s, :αTilde, :δ]
    numericalchangeindata = Dict(
        "conc" => 5.0, 
        "markup" => 22.0 - 10.0, 
        "profit" => 14.0 - 6.0,
        "labor" => (100 - 14.0) - (100 - 6.0),
        "entry" => 8.0 - 12.0,
        "yes" =>  10.5 - 17.0,
        "prodgap" => 25.0,
        "jr" => -7.0,
        "gd" => -8.0)
    changeindata = [1; 1; 1; -1; -1; -1; 1; -1; -1] # Directional changes in data

    table3 = Matrix{Int}(undef, length(momentnames), length(channelnames) + 1)
    table3[:, 1] = changeindata

    for j = eachindex(channelnames)

        channel = channelnames[j]

        if channel == :τ
            moments = transition_moments_τ
        elseif channel == :s
            moments = transition_moments_s
        elseif channel == :αTilde
            moments = transition_moments_αTilde
        elseif channel == :δ
            moments = transition_moments_δ
        end

        for i = eachindex(momentnames)

            moment = momentnames[i]

            if  moment != :labor
                x2015 = getproperty(moments, moment)[tyearpol]
                x1980 = getproperty(initialbgp_moments, moment)
            else
                x2015 = 100 - getproperty(moments, :profit)[tyearpol]
                x1980 = 100 - getproperty(initialbgp_moments, :profit)
            end

            change = x2015 - x1980

            if abs(change) < 0.2 * abs(numericalchangeindata[String(moment)])
                table3[i, j+1] = 0
            elseif change > 0
                table3[i, j+1] = 1
            elseif change < 0 
                table3[i, j+1] = -1
            end

        end
    end

    # Figure 4 in the paper
    years = collect(1981:2015)

    pentry = plot(title = "Entry", legend = :bottomleft, grid = false, framestyle = :box, ylabel = "Entry", xlabel = "Years", ylims = (7, 14), yticks = 7:1:14)
    xticks!(1980:5:2015)
    plot!(pentry, years, transition_moments_δ.entry[1:tyearpol], lw = 3.5, label = "Diffusion")
    plot!(pentry, years, transition_moments_αTilde.entry[1:tyearpol], lw = 3.5, label = "Entry cost", ls = :dash, color = :red)
    plot!(pentry, years, transition_moments_s.entry[1:tyearpol], lw = 3.5, label = "Subsidy", ls = :dot, color = :orange)
    plot!(pentry, years, transition_moments_τ.entry[1:tyearpol], lw = 3.5, label = "Corp. tax", marker = :x, markerstrokewidth = 1, color = :purple)

    pyes = plot(title = "Young emp. share", legend = :bottomleft, grid = false, framestyle = :box, ylabel = "Employment share of firms age <= 5", xlabel = "Years", ylims = (8, 28), yticks = 8:2:28)
    xticks!(1980:5:2015)
    plot!(pyes, years, transition_moments_δ.yes[1:tyearpol], lw = 3.5, label = "Diffusion")
    plot!(pyes, years, transition_moments_αTilde.yes[1:tyearpol], lw = 3.5, label = "Entry cost", ls = :dash, color = :red)
    plot!(pyes, years, transition_moments_s.yes[1:tyearpol], lw = 3.5, label = "Subsidy", ls = :dot, color = :orange)
    plot!(pyes, years, transition_moments_τ.yes[1:tyearpol], lw = 3.5, label = "Corp. tax", marker = :x, markerstrokewidth = 1, color = :purple)

    pgd = plot(title = "Growth dispersion", legend = :bottomleft, grid = false, framestyle = :box, ylabel = "Firm growth rate standard deviation", xlabel = "Years", ylims = (0.3, 0.6), yticks = 0.3:0.05:0.6)
    xticks!(1980:5:2015)
    plot!(pgd, years, transition_moments_δ.gd[1:tyearpol] / 100, lw = 3.5, label = "Diffusion")
    plot!(pgd, years, transition_moments_αTilde.gd[1:tyearpol] / 100, lw = 3.5, label = "Entry cost", ls = :dash, color = :red)
    plot!(pgd, years, transition_moments_s.gd[1:tyearpol] / 100, lw = 3.5, label = "Subsidy", ls = :dot, color = :orange)
    plot!(pgd, years, transition_moments_τ.gd[1:tyearpol] / 100, lw = 3.5, label = "Corp. tax", marker = :x, markerstrokewidth = 1, color = :purple)

    pmarkup = plot(title = "Markup", legend = :topleft, grid = false, framestyle = :box, ylabel = "Markup", xlabel = "Years", ylims = (10, 30), yticks = 10:2:30)
    xticks!(1980:5:2015)
    plot!(pmarkup, years, transition_moments_δ.markup[1:tyearpol], lw = 3.5, label = "Diffusion")
    plot!(pmarkup, years, transition_moments_αTilde.markup[1:tyearpol], lw = 3.5, label = "Entry cost", ls = :dash, color = :red)
    plot!(pmarkup, years, transition_moments_s.markup[1:tyearpol], lw = 3.5, label = "Subsidy", ls = :dot, color = :orange)
    plot!(pmarkup, years, transition_moments_τ.markup[1:tyearpol], lw = 3.5, label = "Corp. tax", marker = :x, markerstrokewidth = 1, color = :purple)

    figure4 = plot(pentry, pyes, pgd, pmarkup, layout = (2, 2), size = (1000, 750), margin = 5mm)

    # Figure 3 in the paper
    figure3 = plot(title = "Path of entry rate, model vs data", ylims = (7, 14), ylabel = "Entry", xlabel = "Years", grid = false, framestyle = :box, legend = :bottomleft)
    xticks!(1980:5:2015)
    plot!(figure3, years, transition_moments_δ.entry[1:tyearpol], label = "Diffusion", lw = 3.5)
    plot!(figure3, years, transition_moments_αTilde.entry[1:tyearpol], label = "Entry cost", lw = 3.5, ls = :dash)

    return figure3, table3, figure4
end

# Section 6 - Investigation of Joint Forces. Given calibrated transition of the model, this function shuts off each channel one-by-one, and calculates corresponding transitions and moments. Finally, it produces Table 6 of the paper.
function Section6(numericp::NumParam, initialbgp::BGP, initialmoments::BgpMoments, panelinitialbgplastyear::Matrix{Float64}, sectordatainitialbgplastperiod::Vector{Sector}, entrantid_initialbgp::Int, transitionp::TransParam, transitionmoments::TransitionMoments)

    @unpack p = initialbgp
    initialp = p

    panel, sectordata = SimPrep(numericp)

    # Turn off τ
    transitionp_τ = TransParam(transitionp, τ = 1.0)
    transition_τ, terminal_τ = TransitionSolver(initialbgp, transitionp_τ, numericp)
    TransitionSimulation!(panel, sectordata, transition_τ, entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod)
    moments_τ = Moments(transition_τ, panel, sectordata)

    # Turn off s
    transitionp_s = TransParam(transitionp, s = 1.0)
    transition_s, terminal_s = TransitionSolver(initialbgp, transitionp_s, numericp)
    TransitionSimulation!(panel, sectordata, transition_s, entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod)
    moments_s = Moments(transition_s, panel, sectordata)

    # Turn off αTilde
    transitionp_αTilde = TransParam(transitionp, αTilde = 1.0)
    transition_αTilde, terminal_αTilde = TransitionSolver(initialbgp, transitionp_αTilde, numericp)
    TransitionSimulation!(panel, sectordata, transition_αTilde, entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod)
    moments_αTilde = Moments(transition_αTilde, panel, sectordata)

    # Turn off δ
    transitionp_δ = TransParam(transitionp, δ = 1.0)
    transition_δ, terminal_δ = TransitionSolver(initialbgp, transitionp_δ, numericp)
    TransitionSimulation!(panel, sectordata, transition_δ, entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod)
    moments_δ = Moments(transition_δ, panel, sectordata)

    # Generate Table 6
    momentnames = [:entry, :labor, :markup, :profit, :conc, :yes, :prodgap, :jr, :gd]
    channelnames = [:τ, :s, :αTilde, :δ]

    table6 = zeros(length(momentnames), length(channelnames))

    for j = eachindex(channelnames)
        channel = channelnames[j]
        
        if channel == :τ
            momentstemp = moments_τ
        elseif channel == :s
            momentstemp = moments_s
        elseif channel == :αTilde
            momentstemp = moments_αTilde
        elseif channel == :δ
            momentstemp = moments_δ
        end

        for i = eachindex(momentnames)
            moment = momentnames[i]

            if moment != :labor
                contr = round(100 * (getproperty(transitionmoments, moment)[numericp.tyearpol] - getproperty(momentstemp, moment)[numericp.tyearpol]) / (getproperty(transitionmoments, moment)[numericp.tyearpol] - getproperty(initialmoments, moment)), digits = 1)
            else
                contr = round(100 * ((100 - getproperty(transitionmoments, :profit)[numericp.tyearpol]) - (100 - getproperty(momentstemp, :profit)[numericp.tyearpol])) / ((100 - getproperty(transitionmoments, :profit)[numericp.tyearpol]) - (100 - getproperty(initialmoments, :profit))), digits = 1)
            end

            table6[i,j] = contr
        end        
    end

    return table6
end

# Section 7 - Market Power and Welfare. Function Section7_Figure7a and Section7_Figure7b below draws figures 7a and 7b, respectively.
@inline WelfareBgp(initialY::Float64, ρ::Float64, g::Float64) = log(initialY)/ρ + g/(ρ^2)

function Section7_Figure7a(initialbgp::BGP)
     
    # δrange = collect(range(0.3, 3, step = 0.3))
    δrange = collect(range(0.1, 3, step = 0.1))
    δrange = [0.03; δrange]

    ## Welfare across BGPs for different δs
    welfare_bgp = zeros(length(δrange))

    Threads.@threads for i = eachindex(δrange)

        δ = δrange[i]

        p′ = StrParam(initialbgp.p, δ = δ)
        bgp′ = BgpSolver(p′, initialbgp.numericp)

        # Initial period output: Q = 1 across BGPs. Equation 19
        first_term = Step(0, bgp′.p.ψ)*bgp′.μnn
        for m = 1:bgp′.numericp.mbar
            first_term += Step(m, bgp′.p.ψ)*bgp′.μ[m]
        end
        first_term = bgp′.p.λ^(-(first_term))/bgp′.ω
    
        second_term = bgp′.μnn * log(((bgp′.p.λ^(Step(0, bgp′.p.ψ)) * (bgp′.p.β*bgp′.snn*(1-bgp′.snn))/(1 - bgp′.p.β*bgp′.snn))^bgp′.p.β + ((bgp′.p.β*bgp′.snn*(1-bgp′.snn))/(1 - bgp′.p.β*(1-bgp′.snn)))^bgp′.p.β)^(1/bgp′.p.β))
        for m = 1:bgp′.numericp.mbar
            second_term += bgp′.μ[m] * log(((bgp′.p.λ^(Step(m, bgp′.p.ψ)) * (bgp′.p.β*bgp′.sl[m]*(1-bgp′.sl[m]))/(1 - bgp′.p.β*bgp′.sl[m]))^bgp′.p.β + ((bgp′.p.β*bgp′.sl[m]*(1-bgp′.sl[m]))/(1 - bgp′.p.β*(1-bgp′.sl[m])))^bgp′.p.β)^(1/bgp′.p.β))
        end
        second_term = exp(second_term)
    
        initialY = first_term * second_term    
        
        # Growth rate
        growth_rate = 0.0
        for m = 1:initialbgp.numericp.mbar
            growth_rate += bgp′.xl[m]*bgp′.μ[m]*(Step(m+1, bgp′.p.ψ) - Step(m, bgp′.p.ψ))
        end
        growth_rate = log(bgp′.p.λ)*((bgp′.xnn + bgp′.xnn + bgp′.xenn)*bgp′.μnn*(Step(1, bgp′.p.ψ) - Step(0, bgp′.p.ψ)) + growth_rate)*bgp′.numericp.dt

        growth_rate = (1+growth_rate)^(1/bgp′.numericp.dt)-1.0

        welfare_bgp[i] = WelfareBgp(initialY, bgp′.p.ρ, growth_rate)
    end

    ## Plot figure 7a
    figure7a = plot(title = "Aggregate welfare across BGPs", framestyle = :box, grid = false, xlabel = "Knowledge diffusion δ", ylabel = "Aggregate welfare", legend = false)
    plot!(figure7a, δrange, welfare_bgp, lw = 3.5, xticks = 0:0.3:3)

    return figure7a
end

# ATTENTION: Drawing figure 7b requires simulating a large number of sectors. Therefore, it can take a lot of time.  
function Section7_Figure7b(initialbgp::BGP)
        
    ## Number of sectors is increased to 100_000 from the benchmark value 10_000 to reduce statistical errors
    numericp′ = NumParam(initialbgp.numericp, numofsectors = 100_000)
    initialbgp′ = BgpSolver(initialbgp.p, numericp′)
    δrange2 = [0.01, 0.02, 0.04, 0.06, initialbgp.p.δ, 0.10, 0.12, 0.14, 0.16]

    ## Simulate initial bgp
    panel, sectordata = SimPrep(numericp′)
    entrantid_initialbgp = InitialBgpSimulation!(panel, sectordata, initialbgp′)
    panelinitialbgplastyear = panel[(initialbgp′.numericp.tyears - 1)*initialbgp′.numericp.numoffirms + 1:end, :]
    sectordatainitialbgplastperiod = sectordata[:, end]

    ## Simulate aggregate output: rows = periods, columns = different δ counterfactuals
    Ymatrix = zeros(numericp′.tperiods, length(δrange2))

    for (i, δ) = enumerate(δrange2)

        # Consider a shock to knowledge diffusion beginning from initial bgp
        transitionp = TransParam(
            τ = 1.0,
            s = 1.0,
            αTilde = 1.0,
            δ = δ/initialbgp.p.δ,
            ντ = 1.0,
            νs = 1.0,
            ναTilde = 1.0,
            νδ = 1.0
        )
    
        transition, terminalbgp = TransitionSolver(initialbgp′, transitionp, numericp′)
        TransitionSimulation!(panel, sectordata, transition, entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod)
        _, _, Ymatrix[:, i] = GrowthRateAndAggOutput(sectordata, transition.ω, transition.numericp, transition.sequencep, transition.sl, transition.sf, transition.snn, transition.zl, transition.zf, transition.znn)
    end

    ## Calculate welfares for different δ counterfactuals given simulated output streams
    e = Base.MathConstants.e # Euler's number e
    
    welfare_dict = Dict{Float64, Float64}()

    for (i, δ) = enumerate(δrange2)

        welfare_value = 0.0

        for t = 1:numericp′.tperiodpol
            welfare_value += e^(-initialbgp.p.ρ*t*numericp′.dt) * log(Ymatrix[t, i]) * numericp′.dt
        end

        welfare_dict[δ] = welfare_value
    end

    # η is the consumption equivalent welfare change for a particular δ counterfactual. That is, (1 + η) is the factor which must multiply the baseline consumption stream in order to get the same welfare as in a particular δ counterfactual (only until 2015)
    η_vec = zeros(length(δrange2))

    for (i,δ) = enumerate(δrange2)
        η_vec[i] = exp((welfare_dict[δ] - welfare_dict[initialbgp.p.δ])/((1 - e^(-initialbgp.p.ρ*numericp′.tyearpol))/initialbgp.p.ρ)) - 1
    end

    ## Plot figure 7b
    figure7b = plot(title = "Consumption-equivalent welfare across transitions", framestyle = :box, grid = false, xlabel = "Knowledge diffusion δ", ylabel = "Consumption-equivalent welfare change", legend = false)
    plot!(figure7b, δrange2, η_vec, lw = 2.5, xticks = 0:0.02:0.20, ylims = (-0.8, 0.4), yticks = (-0.8:0.2:0.4))
    vline!(figure7b, [initialbgp.p.δ], ls = :dash, lw = 1.75, color = :black)
    vline!(figure7b, [transitionp.δ*initialbgp.p.δ], ls = :dash, lw = 1.75, color = :black)
    # quiver!(figure7b, [0.075], [0.22], quiver=([-0.037], [0.0]), lw = 2, color = :black)
    # annotate!(figure7b, [(0.087, 0.3, ("Calibrated value", 9, :left)), (0.087, 0.24, ("of diffusion", 9, :left)), (0.057, 0.18, ("Estimated", 9)), (0.057, 0.12, ("decline in", 9)), (0.057, 0.06, ("diffusion", 9))])

    return figure7b
end

# Appendix D - Additional Quantitative Results. Given an initial BGP calibration, this function Sensitivity below calculates elasiticity of each moment with respect to parameters. In particular, it calculates the percentage change in each target in response to a 1 percent increase in each parameter, keeping the other parameters constant at their calibrated values.
function Sensitivity(initialbgp::BGP, initialbgp_moments::BgpMoments)
    
    # Calibrated parameter values
    initialp = initialbgp.p

    # Parameters of interest
    parameters = [:α, :αTilde, :λ, :δ, :ϕ, :β, :ψ]
    sensitivity_table = Dict() #Empty dictionary to store the results
    panel, sectordata = SimPrep(initialbgp.numericp)
    
    for param in parameters

        local paramvalue::Float64

        if param == :λ
            stepsize = initialp.λ - 1
            paramvalue = 1 + 1.01 * stepsize
        else
            paramvalue = 1.01 * getproperty(initialp, param)
        end

        local p::StrParam

        if param == :ϕ
            p = StrParam(initialp, Dict(:ϕ => paramvalue, :ϕTilde => paramvalue)) #Drastic innovation rate increases for both followers and entrants
        else
            p = StrParam(initialp, Dict(param => paramvalue))
        end

        bgp = BgpSolver(p, initialbgp.numericp)
        InitialBgpSimulation!(panel, sectordata, bgp)
        moments = Moments(bgp, panel)

        sensitivity_table[param] = BgpTargets(
            growth = 100 * (moments.growth - initialbgp_moments.growth)/initialbgp_moments.growth,
            entry = 100 * (moments.entry - initialbgp_moments.entry)/initialbgp_moments.entry,
            rnd = 100 * (moments.rnd - initialbgp_moments.rnd)/initialbgp_moments.rnd,
            profit = 100 * (moments.profit - initialbgp_moments.profit)/initialbgp_moments.profit,
            markup = 100 * (moments.markup - initialbgp_moments.markup)/initialbgp_moments.markup,
            gd = 100 * (moments.gd - initialbgp_moments.gd)/initialbgp_moments.gd,
            ec = 100 * (moments.netentrycont - initialbgp_moments.netentrycont)/initialbgp_moments.netentrycont
        )
    end

    return sensitivity_table
end

# Section 8 - Alternative mechanisms: Declining interest rate - The following function performs the anaylsis in which time discount rate ρ declines from 0.05 to 0.01 between 1980 and 2015. Other parameters are constant at their initial BGP values. We produce column 6 of Table 7. ρpol corresponds to a value of ρ in 2015 which should be specified. To replicate column 6, we take ρpol = 0.01.
function DecliningInterestRate(ρpol::Float64, initialbgp::BGP, initialbgp_moments::BgpMoments, entrantid_initialbgp::Int, panelinitialbgplastyear::Matrix{Float64}, sectordatainitialbgplastperiod::Vector{Sector})
    
    transition, terminal = TransitionSolver_DecliningInterestRate(initialbgp, ρpol)

    # Simulation
    panel, sectordata = SimPrep(initialbgp.numericp)
    TransitionSimulation!(panel, sectordata, transition, entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod)
    transition_moments = Moments(transition, panel, sectordata)

    # Generate 6th column of Table 7
    momentnames = [:conc, :markup, :profit, :labor, :entry, :yes, :prodgap, :jr, :gd]
    numericalchangeindata = Dict(
        "conc" => 5.0, 
        "markup" => 22.0 - 10.0, 
        "profit" => 14.0 - 6.0,
        "labor" => (100 - 14.0) - (100 - 6.0),
        "entry" => 8.0 - 12.0,
        "yes" =>  10.5 - 17.0,
        "prodgap" => 25.0,
        "jr" => -7.0,
        "gd" => -8.0)

    column6 = Vector{Int}(undef, length(momentnames))

    for i = eachindex(momentnames)
        moment = momentnames[i]

        if  moment != :labor
            x2015 = getproperty(transition_moments, moment)[initialbgp.numericp.tyearpol]
            x1980 = getproperty(initialbgp_moments, moment)
        else
            x2015 = 100 - getproperty(transition_moments, :profit)[initialbgp.numericp.tyearpol]
            x1980 = 100 - getproperty(initialbgp_moments, :profit)
        end

        change = x2015 - x1980

        if abs(change) < 0.2 * abs(numericalchangeindata[String(moment)])
            column6[i] = 0
        elseif change > 0
            column6[i] = 1
        elseif change < 0 
            column6[i] = -1
        end
    end

    return column6
end

# Section 8 - Alternative mechanisms: Ideas getting harder to find - The following function performs the analysis in which R&D cost parameters α and αTilde increase 10 times from 1980 to 2015. Other parameters are constant at their initial BGP values. We produce column 7 of Table 7. αpol and αTildepol corresponds to values of α and αTilde in 2015, respectively, which should be specified. To replicate column 7, we take αpol and αTildepol being 10 times the values of α and αTilde in initial BGP, respectively. 
function IdeasGettingHarder(αpol::Float64, αTildepol::Float64, initialbgp::BGP, initialbgp_moments::BgpMoments, entrantid_initialbgp::Int, panelinitialbgplastyear::Matrix{Float64}, sectordatainitialbgplastperiod::Vector{Sector})
    
    transition, terminal = TransitionSolver_IdeasGettingHarder(initialbgp, αpol, αTildepol)
    
    # Simulation
    panel, sectordata = SimPrep(initialbgp.numericp)
    TransitionSimulation!(panel, sectordata, transition, entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod)
    transition_moments = Moments(transition, panel, sectordata)

    # Generate 7th column of Table 7
    momentnames = [:conc, :markup, :profit, :labor, :entry, :yes, :prodgap, :jr, :gd]
    numericalchangeindata = Dict(
        "conc" => 5.0, 
        "markup" => 22.0 - 10.0, 
        "profit" => 14.0 - 6.0,
        "labor" => (100 - 14.0) - (100 - 6.0),
        "entry" => 8.0 - 12.0,
        "yes" =>  10.5 - 17.0,
        "prodgap" => 25.0,
        "jr" => -7.0,
        "gd" => -8.0)

    column7 = Vector{Int}(undef, length(momentnames))

    for i = eachindex(momentnames)
        moment = momentnames[i]

        if  moment != :labor
            x2015 = getproperty(transition_moments, moment)[initialbgp.numericp.tyearpol]
            x1980 = getproperty(initialbgp_moments, moment)
        else
            x2015 = 100 - getproperty(transition_moments, :profit)[initialbgp.numericp.tyearpol]
            x1980 = 100 - getproperty(initialbgp_moments, :profit)
        end

        change = x2015 - x1980

        if abs(change) < 0.2 * abs(numericalchangeindata[String(moment)])
            column7[i] = 0
        elseif change > 0
            column7[i] = 1
        elseif change < 0 
            column7[i] = -1
        end
    end

    return column7
end

# Section8 - Alternative mechanisms: Weaker market power of labor - We shock innovation step size, which is defined as λ - 1, so as to match the declining profile of labor share observed in the data. In particular, we calibrated this parameter that it increased by approximately 2.17 times between 1980 and 2015. In order to replicate column 8, we use the following exact value λpol for λ in 2015, λpol = 1 + 2.1669517156287204*(λ_1980 - 1). 
function WeakerLaborPower(λpol::Float64, initialbgp::BGP, initialbgp_moments::BgpMoments, entrantid_initialbgp::Int, panelinitialbgplastyear::Matrix{Float64}, sectordatainitialbgplastperiod::Vector{Sector})
    
    transition, terminal = TransitionSolver_WeakerLaborPower(initialbgp, λpol)
    
    # Simulation
    panel, sectordata = SimPrep(initialbgp.numericp)
    TransitionSimulation!(panel, sectordata, transition, entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod)
    transition_moments = Moments(transition, panel, sectordata)

    # Generate 8th column of Table 7
    momentnames = [:conc, :markup, :profit, :labor, :entry, :yes, :prodgap, :jr, :gd]
    numericalchangeindata = Dict(
        "conc" => 5.0, 
        "markup" => 22.0 - 10.0, 
        "profit" => 14.0 - 6.0,
        "labor" => (100 - 14.0) - (100 - 6.0),
        "entry" => 8.0 - 12.0,
        "yes" =>  10.5 - 17.0,
        "prodgap" => 25.0,
        "jr" => -7.0,
        "gd" => -8.0)

    column8 = Vector{Int}(undef, length(momentnames))

    for i = eachindex(momentnames)
        moment = momentnames[i]

        if  moment != :labor
            x2015 = getproperty(transition_moments, moment)[initialbgp.numericp.tyearpol]
            x1980 = getproperty(initialbgp_moments, moment)
        else
            x2015 = 100 - getproperty(transition_moments, :profit)[initialbgp.numericp.tyearpol]
            x1980 = 100 - getproperty(initialbgp_moments, :profit)
        end

        change = x2015 - x1980

        if abs(change) < 0.2 * abs(numericalchangeindata[String(moment)])
            column8[i] = 0
        elseif change > 0
            column8[i] = 1
        elseif change < 0 
            column8[i] = -1
        end
    end

    return column8
end