## ANALYTIC MOMENTS

function GrowthRate(xL::Vector{Float64}, xNN::Float64, xENN::Float64, μ::Vector{Float64}, μNN::Float64, λ::Float64, ψ::Float64, dt::Float64)

    # Leader formulation (per period)
    sum = 0.0
    for m = 1:length(μ)
        sum = sum + xL[m]*μ[m]*(Step(m+1, ψ) - Step(m, ψ))
    end
    gperiod = log(λ)*((xNN + xNN + xENN)*μNN*(Step(1, ψ) - Step(0, ψ)) + sum)*dt

    # Annual growth rate
    gann = 100*((1+gperiod)^(1/dt)-1.0)

    return gann
end

function InterestRate(g::Float64, ρ::Float64, dt::Float64)
    r = g + ρ*dt # per period
    rAnn = 100*((1+r)^(1/dt)-1.0)

    return rAnn, r
end

function EntryRate(xE::Vector{Float64}, xENN::Float64, μ::Vector{Float64}, μNN::Float64)

    sum = 0.5*xENN*μNN
    for m = 1:length(μ)
        sum += 0.5*xE[m]*μ[m]
    end

    return 100*sum
end

function RndToGdp(xL::Vector{Float64}, xF::Vector{Float64}, xNN::Float64, μ::Vector{Float64}, μNN::Float64, ω::Float64, α::Float64, γ::Float64)

    sum = (RndLabor(xNN, α, γ)*ω + RndLabor(xNN, α, γ)*ω)*μNN
    for m = 1:length(μ)
        sum += (RndLabor(xL[m], α, γ)*ω + RndLabor(xF[m], α, γ)*ω)*μ[m]
    end

    return 100*sum
end

@inline ProfitShare(ω::Float64) = 100*(1.0 - ω)

function Markup(zL::Vector{Float64}, zF::Vector{Float64}, zNN::Float64, μ::Vector{Float64}, μNN::Float64, sl::Vector{Float64}, sf::Vector{Float64}, snn::Float64)

    sum = (snn*zNN + snn*zNN)*μNN
    for m = 1:length(μ)
        sum += (sl[m]*zL[m] + sf[m]*zF[m])*μ[m]
    end

    return 100*(sum-1)
end

function Concentration(sL::Vector{Float64}, sF::Vector{Float64}, sNN::Float64, μ::Vector{Float64}, μNN::Float64)

    sum = (sNN^2 + sNN^2)*μNN
    for m = 1:length(μ)
        sum += (sL[m]^2 + sF[m]^2)*μ[m]
    end

    return 100*sum
end

function ProdGapDisp(μ::Vector{Float64}, μNN::Float64, λ::Float64, ψ::Float64)

    sum = log(λ^(0^ψ))*μNN
    for m = 1:length(μ)
        sum += log(λ^(m^ψ))*μ[m]
    end

    return 100*sum
end

function ProdGapDispStd(zl::Vector{Float64}, zf::Vector{Float64}, znn::Float64, μ::Vector{Float64}, μnn::Float64)

    logzl = log.(zl)
    logzf = log.(zf)
    logznn = log(znn)

    prodgap = μnn*(100 * StatsBase.std([logznn; logznn], corrected = true))

    for m = eachindex(μ)
        prodgap += μ[m]*(100 * StatsBase.std([logzl[m]; logzf[m]], corrected = true))
    end

    return prodgap
end



## SIMULATION MOMENTS

# Growth rate and aggregate output (on period-by-period panel), and welfares
function GrowthRateAndAggOutput(sectordata::Matrix{Sector}, ω::Vector{Float64}, numericp::NumParam, sequencep::Vector{StrParam}, sl::Matrix{Float64}, sf::Matrix{Float64}, snn::Vector{Float64}, zl::Matrix{Float64}, zf::Matrix{Float64}, znn::Vector{Float64})

    @unpack_NumParam numericp

    yperiods = zeros(tperiods) # Aggregate output per period
    yannual = zeros(tyears) # Aggregate output per year

    # Calculate aggregate output per period
    @inbounds for t = 1:tperiods
        
        @unpack_StrParam sequencep[t]

        @inbounds for id = 1:numofsectors

            sector = sectordata[id ,t]
            m = sector.m
            ql = sector.leader.q
            qf = sector.follower.q

            if m > 0
                yl = FirmOutput(sl[m, t], zl[m, t], ql, ω[t])
                yf = FirmOutput(sf[m, t], zf[m, t], qf, ω[t])
            else
                yl = FirmOutput(snn[t], znn[t], ql, ω[t])
                yf = FirmOutput(snn[t], znn[t], qf, ω[t])
            end

            yperiods[t] += log(SectorOutput(yl, yf, β)) * (1/numofsectors)
        end

        yperiods[t] = exp(yperiods[t])
    end

    # Calculate aggregate output per year
    for tyear = 1:tyears
        yannual[tyear] = sum(yperiods[(tyear - 1)*periodsperyear + 1 : tyear*periodsperyear] * dt)
    end

    # Annual growth rates
    g = zeros(tyears-1)

    for t = 1:tyears-1
        g[t] = 100 * (yannual[t+1] - yannual[t]) / yannual[t]
    end

    return g, yannual, yperiods
end

# Young employment share (yes)
function YoungEmploymentShare(panel::Matrix{Float64}, numoffirms::Int, tyears::Int, laborcol::Int, agecol::Int)

    # Pre-allocate yes vector
    yes = zeros(tyears)

    @inbounds for tyear = 1:tyears
        subpanel = @view panel[(tyear-1)*numoffirms+1:tyear*numoffirms, :]
        totalemployment = sum(subpanel[:, laborcol])
        youngemployment = sum(subpanel[subpanel[:, agecol] .<= 5.0, laborcol])
        yes[tyear] = 100.0 * youngemployment/totalemployment
    end

    return yes
end

# Job reallocation and Growth dispersion: First function calculates jr (job reallocation) and gd (growth dispersion) between tyear and tyear + 1, where tyear >= 1 (temp is tempforjr)
function JRbetweenTwoYears!(tyear::Int, temp::Matrix{Float64}, idxentrants::BitVector, panel::Matrix{Float64}, numofsectors::Int, numoffirms::Int, firmidcol::Int, laborcol::Int, agecol::Int)

    # temp has 3 columns to hold for: 1=x, 2=g, 3=jr. temp has 2*numoffirms rows. First half is for entrants, second half is for incumbents
    xcol = 1
    gcol = 2
    jrcol = 3

    # Reset temp due to correct weighting with x
    temp[:, xcol] .= 0.0

    # Split temp
    tempentrant = @view temp[1:numoffirms, :]
    tempincumbent = @view temp[numoffirms + 1:end, :]

    # Split panel
    subpaneltoday = @view panel[(tyear-1)*numoffirms+1:tyear*numoffirms, :]
    subpaneltomorrow = @view panel[tyear*numoffirms+1:(tyear+1)*numoffirms, :]

    # First, handle with today's entrants. Note that today's entrants have age = 1 tomorrow
    @. idxentrants = subpaneltomorrow[:, agecol] == 1.0
    @. tempentrant[idxentrants, xcol] = 0.5*subpaneltomorrow[idxentrants, laborcol]
    @. tempentrant[idxentrants, gcol] = 2.0
    @. tempentrant[idxentrants, jrcol] = 2.0

    # Second, handle with today's incumbents. Note that some of incumbents are continuing tomorrow, the rest are exiters.
    @inbounds for j = 1:numofsectors

        leaderrownumber = (j-1)*2 + 1       # leader's row number
        followerrownumber = (j-1)*2 + 2     # follower's row number

        # We know that an incumbent today cannot leave its sector. Therefore, we will try to find it in the same sector tomorrow
        leaderidtoday = subpaneltoday[leaderrownumber, firmidcol]
        followeridtoday = subpaneltoday[followerrownumber, firmidcol]

        leaderidtomorrow = subpaneltomorrow[leaderrownumber, firmidcol]
        followeridtomorrow = subpaneltomorrow[followerrownumber, firmidcol]

        # Allocate labors today and tomorrow
        leaderlabortoday = subpaneltoday[leaderrownumber, laborcol]
        followerlabortoday = subpaneltoday[followerrownumber, laborcol]

        leaderlabortomorrow = subpaneltomorrow[leaderrownumber, laborcol]
        followerlabortomorrow = subpaneltomorrow[followerrownumber, laborcol]

        # Now try to match them. First today's leader
        if leaderidtoday == leaderidtomorrow # Today's leader became tomorrow's leader
            x = 0.5*(leaderlabortomorrow + leaderlabortoday)
            g = (leaderlabortomorrow - leaderlabortoday)/x

            tempincumbent[leaderrownumber, xcol] = x
            tempincumbent[leaderrownumber, gcol] = g
            tempincumbent[leaderrownumber, jrcol] = abs(g)

        elseif leaderidtoday == followeridtomorrow # Today's leader became tomorrow's follower
            x = 0.5*(followerlabortomorrow + leaderlabortoday)
            g = (followerlabortomorrow - leaderlabortoday)/x

            tempincumbent[leaderrownumber, xcol] = x
            tempincumbent[leaderrownumber, gcol] = g
            tempincumbent[leaderrownumber, jrcol] = abs(g)

        else # Today's leader exited
            tempincumbent[leaderrownumber, xcol] = 0.5*leaderlabortoday
            tempincumbent[leaderrownumber, gcol] = -2.0
            tempincumbent[leaderrownumber, jrcol] = 2.0

        end

        # Now today's follower's turn
        if followeridtoday == leaderidtomorrow # Today's follower became tomorrow's leader
            x = 0.5*(leaderlabortomorrow + followerlabortoday)
            g = (leaderlabortomorrow - followerlabortoday)/x

            tempincumbent[followerrownumber, xcol] = x
            tempincumbent[followerrownumber, gcol] = g
            tempincumbent[followerrownumber, jrcol] = abs(g)

        elseif followeridtoday == followeridtomorrow # Today's follower became tomorrow's follower
            x = 0.5*(followerlabortomorrow + followerlabortoday)
            g = (followerlabortomorrow - followerlabortoday)/x

            tempincumbent[followerrownumber, xcol] = x
            tempincumbent[followerrownumber, gcol] = g
            tempincumbent[followerrownumber, jrcol] = abs(g)

        else # Today's follower exited
            tempincumbent[followerrownumber, xcol] = 0.5*followerlabortoday
            tempincumbent[followerrownumber, gcol] = -2.0
            tempincumbent[followerrownumber, jrcol] = 2.0

        end

    end

    jr = 100*mean(temp[:, jrcol], weights(temp[:, xcol]))
    gd = 100*StatsBase.std(temp[:, gcol], weights(temp[:, xcol]), corrected=false)

    return jr, gd
end

function JobReallocation!(temp::Matrix{Float64}, idxentrants::BitVector, panel::Matrix{Float64}, numofsectors::Int, numoffirms::Int, firmidcol::Int, laborcol::Int, agecol::Int)

    totalyears = Int(size(panel, 1) / numoffirms)

    # Pre-allocate results. 
    jr = zeros(totalyears - 1)
    gd = zeros(totalyears - 1)

    @inbounds for tyear = 1:totalyears - 1
        jr[tyear], gd[tyear] = JRbetweenTwoYears!(tyear, temp, idxentrants, panel, numofsectors, numoffirms, firmidcol, laborcol, agecol)
    end

    return jr, gd
end

# Growth decomposition (Based on Foster, Haltiwanger, Krizan (NBER, 2001). They propose two methods. Method 1 is their preferred one). k refers to number of years over which entry contributions are calculated
function GrowthDecompMethod1(panel::Matrix{Float64}, numericp::NumParam, k::Int)

    @unpack_NumParam numericp

    # Preallocation
    Qgrowth = zeros(tyears - k)
    ΔQ = zeros(tyears - k)
    within = zeros(tyears - k)
    between = zeros(tyears - k)
    cross = zeros(tyears - k)
    entry = zeros(tyears - k)
    exit = zeros(tyears - k)

    # More preallocation
    sharetoday = zeros(numoffirms)
    sharetomorrow = zeros(numoffirms)
    entering = BitVector(zeros(numoffirms))
    exiting = BitVector(zeros(numoffirms))

    # Run the loop over all years except the last k years. Calculate decomposition between t and t+k and assign it to t.
    @inbounds for t = 1:tyears-k

        # Reset bit vectors
        @. entering = 0
        @. exiting = 0

        # Reset sums
        sumwithin = 0.0
        sumbetween = 0.0
        sumcross = 0.0

        # Construct years
        thisyear = t
        nextyear = t + k

        # Extract relevant parts of panel and productivity
        qtoday = @view panel[(thisyear - 1)*numoffirms+1:thisyear*numoffirms, prodcol]
        qtomorrow = @view panel[(nextyear - 1)*numoffirms+1:nextyear*numoffirms, prodcol]

        firmidtoday = @view panel[(thisyear - 1)*numoffirms+1:thisyear*numoffirms, firmidcol]
        firmidtomorrow = @view panel[(nextyear - 1)*numoffirms+1:nextyear*numoffirms, firmidcol]

        @. sharetoday = panel[(thisyear - 1)*numoffirms+1:thisyear*numoffirms, laborcol]
        @. sharetomorrow = panel[(nextyear - 1)*numoffirms+1:nextyear*numoffirms, laborcol]
        temp = sum(sharetoday)
        @. sharetoday = sharetoday / temp
        temp = sum(sharetomorrow)
        @. sharetomorrow = sharetomorrow / temp

        # Calculate average productivities
        Qtoday = sum(qtoday .* sharetoday)
        Qtomorrow = sum(qtomorrow .* sharetomorrow)
        ΔQ[t] = Qtomorrow - Qtoday
        Qgrowth[t] = 100 * (1/k) * (Qtomorrow - Qtoday) / Qtoday # Growth rate with level q's

        # Determine continuing, entering and exiting firms (sets C, N, X respectively)
        @inbounds for j = 1:numofsectors

            # Row index of leader and follower
            leaderrownumber = (j-1)*numoffirmsinasector + 1
            followerrownumber = (j-1)*numoffirmsinasector + 2

            # Leader and follower firm ids today and tomorrow
            leaderidtoday = firmidtoday[leaderrownumber]
            followeridtoday = firmidtoday[followerrownumber]

            leaderidtomorrow = firmidtomorrow[leaderrownumber]
            followeridtomorrow = firmidtomorrow[followerrownumber]

            # Leader and follower quantities (q and shares) today and tomorrow
            leaderqtoday = qtoday[leaderrownumber]
            followerqtoday = qtoday[followerrownumber]

            leaderqtomorrow = qtomorrow[leaderrownumber]
            followerqtomorrow = qtomorrow[followerrownumber]

            leadersharetoday = sharetoday[leaderrownumber]
            followersharetoday = sharetoday[followerrownumber]

            leadersharetomorrow = sharetomorrow[leaderrownumber]
            followersharetomorrow = sharetomorrow[followerrownumber]

            # Determine entering firms tomorrow
            if leaderidtomorrow != leaderidtoday && leaderidtomorrow != followeridtoday
                entering[leaderrownumber] = 1
            end

            if followeridtomorrow != leaderidtoday && followeridtomorrow != followeridtoday
                entering[followerrownumber] = 1
            end

            # Calculate continuing components: withing, between, cross; and determine exiting firms

            # Leader
            if leaderidtoday == leaderidtomorrow # Leader today is leader tomorrow. Continuing

                sumwithin += leadersharetoday * (leaderqtomorrow - leaderqtoday)
                sumbetween += (leaderqtoday - Qtoday) * (leadersharetomorrow - leadersharetoday)
                sumcross += (leaderqtomorrow - leaderqtoday) * (leadersharetomorrow - leadersharetoday)

            elseif leaderidtoday == followeridtomorrow # Leader today is follower tomorrow. Continuing

                sumwithin += leadersharetoday * (followerqtomorrow - leaderqtoday)
                sumbetween += (leaderqtoday - Qtoday) * (followersharetomorrow - leadersharetoday)
                sumcross += (followerqtomorrow - leaderqtoday) * (followersharetomorrow - leadersharetoday)

            else # Leader today exits

                exiting[leaderrownumber] = 1

            end

            # Follower
            if followeridtoday == leaderidtomorrow # Follower today is leader tomorrow. Continuing

                sumwithin += followersharetoday * (leaderqtomorrow - followerqtoday)
                sumbetween += (followerqtoday - Qtoday) * (leadersharetomorrow - followersharetoday)
                sumcross += (leaderqtomorrow - followerqtoday) * (leadersharetomorrow - followersharetoday)

            elseif followeridtoday == followeridtomorrow # Follower today is follower tomorrow. Continuing

                sumwithin += followersharetoday * (followerqtomorrow - followerqtoday)
                sumbetween += (followerqtoday - Qtoday) * (followersharetomorrow - followersharetoday)
                sumcross += (followerqtomorrow - followerqtoday) * (followersharetomorrow - followersharetoday)

            else # Follower today exits

                exiting[followerrownumber] = 1

            end
        end

        # Finalize continuing component for year t
        within[t] = sumwithin
        between[t] = sumbetween
        cross[t] = sumcross

        # Calculate entry and exit components
        entry[t] = sum(@. entering * sharetomorrow * (qtomorrow - Qtoday))
        exit[t] = (-1) * sum(@. exiting * sharetoday * (qtoday - Qtoday))
    end

    # Contribution shares (at each point in time t, they have to sum up to 100%)
    withincont = @. 100 * within / ΔQ
    betweencont = @. 100 * between / ΔQ
    crosscont = @. 100 * cross / ΔQ
    entrycont = @. 100 * entry / ΔQ
    exitcont = @. 100 * exit / ΔQ

    return Qgrowth, withincont, betweencont, crosscont, entrycont, exitcont
end

function GrowthDecompMethod2(panel::Matrix{Float64}, numericp::NumParam, k::Int)

    @unpack_NumParam numericp

    # Preallocation (there is no cross term different than method 1)
    Qgrowth = zeros(tyears - k)
    ΔQ = zeros(tyears - k)
    within = zeros(tyears - k)
    between = zeros(tyears - k)
    entry = zeros(tyears - k)
    exit = zeros(tyears - k)

    # More preallocation
    sharetoday = zeros(numoffirms)
    sharetomorrow = zeros(numoffirms)
    entering = BitVector(zeros(numoffirms))
    exiting = BitVector(zeros(numoffirms))

    # Run the loop over all years except the last k years. Calculate decomposition between t and t+k and assign it to t.
    @inbounds for t = 1:tyears-k

        # Reset bit vectors
        @. entering = 0
        @. exiting = 0

        # Reset sums
        sumwithin = 0.0
        sumbetween = 0.0

        # Construct years
        thisyear = t
        nextyear = t + k

        # Extract relevant parts of panel and productivity
        qtoday = @view panel[(thisyear - 1)*numoffirms+1:thisyear*numoffirms, prodcol]
        qtomorrow = @view panel[(nextyear - 1)*numoffirms+1:nextyear*numoffirms, prodcol]

        firmidtoday = @view panel[(thisyear - 1)*numoffirms+1:thisyear*numoffirms, firmidcol]
        firmidtomorrow = @view panel[(nextyear - 1)*numoffirms+1:nextyear*numoffirms, firmidcol]

        @. sharetoday = panel[(thisyear - 1)*numoffirms+1:thisyear*numoffirms, laborcol]
        @. sharetomorrow = panel[(nextyear - 1)*numoffirms+1:nextyear*numoffirms, laborcol]
        temp = sum(sharetoday)
        @. sharetoday = sharetoday / temp
        temp = sum(sharetomorrow)
        @. sharetomorrow = sharetomorrow / temp

        # Calculate average productivities and growths
        Qtoday = sum(qtoday .* sharetoday)
        Qtomorrow = sum(qtomorrow .* sharetomorrow)
        Qaverage = 0.5 * (Qtomorrow + Qtoday)
        ΔQ[t] = Qtomorrow - Qtoday
        Qgrowth[t] = 100 * (1/k) * (Qtomorrow - Qtoday) / Qtoday

        # Look into each sector. Determine continuing firms and continuing components within, between. Also determine entering firms
        @inbounds for j = 1:numofsectors

            # Row index of leader and follower
            leaderrownumber = (j-1)*numoffirmsinasector + 1
            followerrownumber = (j-1)*numoffirmsinasector + 2

            # Leader and follower firm ids today and tomorrow
            leaderidtoday = firmidtoday[leaderrownumber]
            followeridtoday = firmidtoday[followerrownumber]

            leaderidtomorrow = firmidtomorrow[leaderrownumber]
            followeridtomorrow = firmidtomorrow[followerrownumber]

            # Leader and follower quantities (q and shares) today and tomorrow
            leaderqtoday = qtoday[leaderrownumber]
            followerqtoday = qtoday[followerrownumber]

            leaderqtomorrow = qtomorrow[leaderrownumber]
            followerqtomorrow = qtomorrow[followerrownumber]

            leadersharetoday = sharetoday[leaderrownumber]
            followersharetoday = sharetoday[followerrownumber]

            leadersharetomorrow = sharetomorrow[leaderrownumber]
            followersharetomorrow = sharetomorrow[followerrownumber]

            # Determine entering firms tomorrow
            if leaderidtomorrow != leaderidtoday && leaderidtomorrow != followeridtoday
                entering[leaderrownumber] = 1
            end

            if followeridtomorrow != leaderidtoday && followeridtomorrow != followeridtoday
                entering[followerrownumber] = 1
            end

            # Calculate continuing components: withing, between, cross; and determine exiting firms

            # Leader
            if leaderidtoday == leaderidtomorrow # Leader today is leader tomorrow. Continuing

                sumwithin += 0.5 * (leadersharetomorrow + leadersharetoday) * (leaderqtomorrow - leaderqtoday)
                sumbetween += (0.5 * (leaderqtomorrow + leaderqtoday) - Qaverage) * (leadersharetomorrow - leadersharetoday)

            elseif leaderidtoday == followeridtomorrow # Leader today is follower tomorrow. Continuing

                sumwithin += 0.5 * (followersharetomorrow + leadersharetoday) * (followerqtomorrow - leaderqtoday)
                sumbetween += (0.5 * (followerqtomorrow + leaderqtoday) - Qaverage) * (followersharetomorrow - leadersharetoday)

            else # Leader today exits

                exiting[leaderrownumber] = 1

            end

            # Follower
            if followeridtoday == leaderidtomorrow # Follower today is leader tomorrow. Continuing

                sumwithin += 0.5 * (leadersharetomorrow + followersharetoday) * (leaderqtomorrow - followerqtoday)
                sumbetween += (0.5 * (leaderqtomorrow + followerqtoday) - Qaverage) * (leadersharetomorrow - followersharetoday)

            elseif followeridtoday == followeridtomorrow # Follower today is follower tomorrow. Continuing

                sumwithin += 0.5 * (followersharetomorrow + followersharetoday) * (followerqtomorrow - followerqtoday)
                sumbetween += (0.5 * (followerqtomorrow + followerqtoday) - Qaverage) * (followersharetomorrow - followersharetoday)

            else # Follower today exits

                exiting[followerrownumber] = 1

            end
        end

        # Finalize continuing component for year t
        within[t] = sumwithin
        between[t] = sumbetween

        # Calculate entry and exit components
        entry[t] = sum(@. entering * sharetomorrow * (qtomorrow - Qaverage))
        exit[t] = (-1) * sum(@. exiting * sharetoday * (qtoday - Qaverage))
    end

    # Contribution shares (at each point in time t, they have to sum up to 100%)
    withincont = @. 100 * within / ΔQ
    betweencont = @. 100 * between / ΔQ
    entrycont = @. 100 * entry / ΔQ
    exitcont = @. 100 * exit / ΔQ

    return Qgrowth, withincont, betweencont, entrycont, exitcont
end

## ALL MOMENTS TOGETHER

# BGP moments
@with_kw struct BgpMoments
    growth::Float64
    entry::Float64
    rnd::Float64
    profit::Float64
    markup::Float64
    yes::Float64
    conc::Float64
    prodgap::Float64
    prodgapstd::Float64
    jr::Float64
    gd::Float64
    within::Float64
    between::Float64
    cross::Float64
    entrycont::Float64
    exitcont::Float64
    netentrycont::Float64
end

function Moments(bgp::BGP, panel::Matrix{Float64})

    @unpack_BGP bgp
    @unpack_NumParam numericp
    @unpack_StrParam p

    # Analytical moments
    gAnn        = GrowthRate(xl, xnn, xenn, μ, μnn, λ, ψ, dt)
    entry       = EntryRate(xe, xenn, μ, μnn)
    rnd         = RndToGdp(xl, xf, xnn, μ, μnn, ω, α, γ)
    profit      = ProfitShare(ω)
    markup      = Markup(zl, zf, znn, μ, μnn, sl, sf, snn)
    conc        = Concentration(sl, sf, snn, μ, μnn)
    prodGap     = ProdGapDisp(μ, μnn, λ, ψ)
    prodGapStd  = ProdGapDispStd(zl, zf, znn, μ, μnn)

    # Simulation moments
    yes = YoungEmploymentShare(panel, numoffirms, tyears, laborcol, agecol)
    tempforjr = zeros(2*numoffirms, 3)
    idxentrants = BitVector(zeros(numoffirms))
    jr, gd = JobReallocation!(tempforjr, idxentrants, panel, numofsectors, numoffirms, firmidcol, laborcol, agecol)
    k = 5
    Qgrowth, within, between, cross, entrycont, exitcont = GrowthDecompMethod1(panel, numericp, k)
    netentrycont = entrycont + exitcont # Net entry contribution which is the sum of two terms: entry contribution and exit contribution

    yesmean = mean(yes[tyearsburnin : end])
    jrmean = mean(jr[tyearsburnin : end])
    gdmean = mean(gd[tyearsburnin : end])
    withinmean = mean(within[numericp.tyearsburnin : end])
    betweenmean = mean(between[tyearsburnin : end])
    crossmean = mean(cross[tyearsburnin : end])
    entrycontmean = mean(entrycont[tyearsburnin : end])
    exitcontmean = mean(exitcont[tyearsburnin : end])
    netentrycontmean = mean(netentrycont[tyearsburnin : end])

    println()
    println("BGP MOMENTS")
    println("-----------")
    println("TFP growth                 = $(round(gAnn, digits=2))%")
    println("R&D to GDP                 = $(round(rnd, digits=2))%")
    println("Firm entry                 = $(round(entry, digits=1))%")
    println("Markup                     = $(round(markup, digits=1))%")
    println("Profit share               = $(round(profit, digits=2))%")
    println("Net entry contribution     = $(round(netentrycontmean, digits=1))%")
    println("Firm growth dispersion     = $(round(gdmean, digits=1))%")

    return BgpMoments(growth = gAnn, entry = entry, rnd = rnd, profit = profit, markup = markup, yes = yesmean, conc = conc, prodgap = prodGap, prodgapstd = prodGapStd, jr = jrmean, gd = gdmean, within = withinmean, between = betweenmean, cross = crossmean, entrycont = entrycontmean, exitcont = exitcontmean, netentrycont = netentrycontmean)
end

function BgpCalibrationMoments(bgp::BGP, panel::Matrix{Float64})

    @unpack_BGP bgp
    @unpack_NumParam numericp
    @unpack_StrParam p

    # Calculate only moments needed for calibration
    growth = GrowthRate(xl, xnn, xenn, μ, μnn, λ, ψ, dt)
    entry = EntryRate(xe, xenn, μ, μnn)
    rnd = RndToGdp(xl, xf, xnn, μ, μnn, ω, α, γ)
    profit = ProfitShare(ω)
    markup = Markup(zl, zf, znn, μ, μnn, sl, sf, snn)

    tempforjr = zeros(2*numoffirms, 3)
    idxentrants = BitVector(zeros(numoffirms))

    _, gdvector = JobReallocation!(tempforjr, idxentrants, panel, numofsectors, numoffirms, firmidcol, laborcol, agecol)
    gd = mean(gdvector[tyearsburnin : end])

    _, _, _, _, entrycontrvector, exitcontrvector = GrowthDecompMethod1(panel, numericp, 5)
    entrycontr = mean(entrycontrvector[tyearsburnin : end])
    exitcontr = mean(exitcontrvector[tyearsburnin : end])
    ec = entrycontr + exitcontr

    return BgpTargets(growth = growth, entry = entry, rnd = rnd, profit = profit, markup = markup, gd = gd, ec = ec)
end


# Transition moments
@with_kw struct TransitionMoments
    growth::Vector{Float64}                                                     # tyears - 1
    entry::Vector{Float64}                                                      # tyears
    rnd::Vector{Float64}                                                        # tyears
    profit::Vector{Float64}                                                     # tyears
    markup::Vector{Float64}                                                     # tyears
    yes::Vector{Float64}                                                        # tyears
    conc::Vector{Float64}                                                       # tyears
    prodgap::Vector{Float64}                                                    # tyears
    prodgapstd::Vector{Float64}                                                 # tyears
    jr::Vector{Float64}                                                         # tyears - 1
    gd::Vector{Float64}                                                         # tyears - 1
    within::Vector{Float64}                                                     # tyears - k
    between::Vector{Float64}                                                    # tyears - k
    cross::Vector{Float64}                                                      # tyears - k
    entrycont::Vector{Float64}                                                  # tyears - k
    exitcont::Vector{Float64}                                                   # tyears - k
    netentrycont::Vector{Float64}                                               # tyears - k
end

function AnnualizePeriodicMomentByAverage(moment::Vector{Float64}, tyears::Int, periodsperyear::Int)
    result = zeros(tyears)

    for tyear = 1:tyears
        result[tyear] = mean(moment[(tyear-1)*periodsperyear+1:tyear*periodsperyear])
    end

    return result
end

function Moments(transition::Transition, panel::Matrix{Float64}, sectordata::Matrix{Sector})

    @unpack_Transition transition
    @unpack_NumParam numericp

    # Pre-allocate
    entry = zeros(tperiods)
    rnd = similar(entry)
    profit = similar(entry)
    markup = similar(entry)
    conc = similar(entry)
    prodGap = similar(entry)
    prodGapStd = similar(entry)

    # Calculate analytic moments
    @inbounds for t = 1:tperiods

        @unpack α, γ, λ, ψ = sequencep[t]

        entry[t]       = EntryRate(xe[:,t], xenn[t], μ[:,t], μnn[t])
        rnd[t]         = RndToGdp(xl[:,t], xf[:,t], xnn[t], μ[:,t], μnn[t], ω[t], α, γ)
        profit[t]      = ProfitShare(ω[t])
        markup[t]      = Markup(zl[:, t], zf[:, t], znn[t], μ[:,t], μnn[t], sl[:, t], sf[:, t], snn[t])
        conc[t]        = Concentration(sl[:, t], sf[:, t], snn[t], μ[:,t], μnn[t])
        prodGap[t]     = ProdGapDisp(μ[:,t], μnn[t], λ, ψ)
        prodGapStd[t]     = ProdGapDispStd(zl[:,t], zf[:,t], znn[t], μ[:,t], μnn[t])
    end

    # Annualize
    entryann = AnnualizePeriodicMomentByAverage(entry, tyears, periodsperyear)
    rndann = AnnualizePeriodicMomentByAverage(rnd, tyears, periodsperyear)
    profitann = AnnualizePeriodicMomentByAverage(profit, tyears, periodsperyear)
    markupann = AnnualizePeriodicMomentByAverage(markup, tyears, periodsperyear)
    concann = AnnualizePeriodicMomentByAverage(conc, tyears, periodsperyear)
    prodGapann = AnnualizePeriodicMomentByAverage(prodGap, tyears, periodsperyear)
    prodGapStdann = AnnualizePeriodicMomentByAverage(prodGapStd, tyears, periodsperyear)

    # Calculate simulation moments
    yes = YoungEmploymentShare(panel, numoffirms, tyears, 4, 5)
    tempforjr = zeros(2*numoffirms, 3)
    idxentrants = BitVector(zeros(numoffirms))
    jr, gd = JobReallocation!(tempforjr, idxentrants, panel, numofsectors, numoffirms, 2, 4, 5)
    k = 5
    Qgrowth, within, between, cross, entrycont, exitcont = GrowthDecompMethod1(panel, numericp, k)
    netentrycont = entrycont + exitcont
    gann, _, _ = GrowthRateAndAggOutput(sectordata, ω, numericp, sequencep, sl, sf, snn, zl, zf, znn)

    return TransitionMoments(growth = gann, entry = entryann, rnd = rndann, profit = profitann, markup = markupann, yes = yes, conc = concann, prodgap = prodGapann, prodgapstd = prodGapStdann, jr = jr, gd = gd, within = within, between = between, cross = cross, entrycont = entrycont, exitcont = exitcont, netentrycont = netentrycont)
end

function TransitionCalibrationMoments(transition::Transition, numericp::NumParam, entry1980::Float64, panel::Matrix{Float64})

    @unpack_Transition transition
    @unpack_NumParam numericp

    # Corresponding periods
    tperiod2000 = (2000 - 1980) * periodsperyear
    tperiod1990 = (1990 - 1980) * periodsperyear

    # Entry rates
    entry2015 = EntryRate(xe[:,tperiodpol], xenn[tperiodpol], μ[:,tperiodpol], μnn[tperiodpol])
    entry2000 = EntryRate(xe[:,tperiod2000], xenn[tperiod2000], μ[:,tperiod2000], μnn[tperiod2000])
    entry1990 = EntryRate(xe[:,tperiod1990], xenn[tperiod1990], μ[:,tperiod1990], μnn[tperiod1990])

    entrydeclinein10 = 100.0 * (1.0 - entry1990/entry1980)
    entrydeclinein20 = 100.0 * (1.0 - entry2000/entry1980)

    # Growth dispersion
    tempforjr = zeros(2*numoffirms, 3)
    idxentrants = BitVector(zeros(numoffirms))
    _, gdvec = JobReallocation!(tempforjr, idxentrants, panel, numofsectors, numoffirms, firmidcol, laborcol, agecol)

    # Markup
    markup2015 = Markup(zl[:, tperiodpol], zf[:, tperiodpol], znn[tperiodpol], μ[:,tperiodpol], μnn[tperiodpol], sl[:, tperiodpol], sf[:, tperiodpol], snn[tperiodpol])

    return TransitionTargets(entry = entry2015, gd = gdvec[tyearpol], markup = markup2015, entrydeclinein10 = entrydeclinein10, entrydeclinein20 = entrydeclinein20)
end

