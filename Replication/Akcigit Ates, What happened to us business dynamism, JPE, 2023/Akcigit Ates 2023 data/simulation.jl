#=
sectordata = Matrix with the size (numericp.numofsectors × numericp.tperiods): Each column is the snapshot of the economy at the end of the period. Take a column t, i.e. sectordata[:,t]. This column is a snapshot of the economy at the end of the period t. That is, all the firms in this column produced in period t. Some of them are replaced "at the end of the period t" so that these firms will not appear in next period's data sectordata[:,t+1]. Therefore, all the firms in sectordata[:,t] aged "dt" years, since they operated in a period of length "dt".

Another interpretation of this is that, if en entrant replaces an incumbent in period t, that entrant will always appear in t+1's data set sectordata[:,t+1], since this entrant will produce and then innovate. Therefore, if an entrant enters in t, it will always age dt. In sectordata[:,t+1], its age will be dt.
=#


## FIRM AND SECTOR TYPES
@with_kw struct Firm
    firmid::Int             # Unique firm id
    n::Int                  # Firm step number
    q::Float64              # Firm productivity
end

@with_kw struct Sector
    sectorid::Int                       # Unique sector id
    leader::Firm                        # firmid, n
    follower::Firm                      # firmid, n
    m::Int = leader.n - follower.n      # Productivity gap between leader and follower
end

import Base.==
import Base.in

in(firm::Firm, sector::Sector) = (firm == sector.leader || firm == sector.follower) ? true : false
in(firmid::Float64, sector::Sector) = (firmid == sector.leader.firmid || firmid == sector.follower.firmid) ? true : false
==(firm::Firm, sector::Sector) = throw(ArgumentError("Cannot compare a firm with a sector!"))
==(firm1::Firm, firm2::Firm) = firm1.firmid == firm2.firmid ? true : false

## COMMON FUNCTIONS
# Preallocate
function SimPrep(numericp::NumParam)

    @unpack tyears, numoffirms, panelnumofcolumns, numofsectors, tperiods = numericp

    panel = zeros(tyears * numoffirms, panelnumofcolumns)
    sectordata = Matrix{Sector}(undef, numofsectors, tperiods)

    return panel, sectordata
end

# Calculate probabilities (weights) of events for unleveled sectors
function WeightsLevConstructer!(weightslev::Vector{ProbabilityWeights{Float64,Float64,Array{Float64,1}}}, xl::Vector{Float64}, xf::Vector{Float64}, xe::Vector{Float64}, mbar::Int, dt::Float64, ϕ::Float64, δ::Float64, ϕTilde::Float64)

    @inbounds for m = 1:mbar
        xL = xl[m]
        xF = xf[m]
        xE = xe[m]

        # In a small period of time dt, all possible events in an unleveled sector are as follows
        # 1. Leader innovates -> Probability = xL*dt
        # 2. Follower innovates drastically or knowledge diffusion happens -> Probability = (ϕ*xF + δ)*dt
        # 3. Follower incrementally innovates -> Probability = (1-ϕ)*xF*dt
        # 4. Entrant innovates drastically -> Probability = ϕTilde*xE*dt
        # 5. Entrant innovates incrementally -> Probability = (1-ϕTilde)*xE*dt
        # 6. Nothing happens -> Probability = One minus sum of all probabilities above

        temp = [xL, ϕ*xF + δ, (1-ϕ)*xF, ϕTilde*xE, (1-ϕTilde)*xE]*dt
        push!(temp, 1-sum(temp))

        weightslev[m] = ProbabilityWeights(temp)
    end
end

# Calculate probabilities (weights) of events for neck-and-neck sectors
function WeightsNNConstructer(xnn::Float64, xenn::Float64, dt::Float64)

    # In a small period of time dt, all possible events in a neck-and-neck sector are as follows
    # 1. One of firms innovates -> Probability = xnn*dt
    # 2. Other firm innovates -> Probability = xnn*dt
    # 3. Entrant replaces one of the firms -> Probability = 0.5*xenn*dt
    # 4. Entrant replaces other firm -> Probability = 0.5*xenn*dt
    # 5. Nothing happens -> Probability = One minus sum of all probabilities above

    temp = [xnn; xnn; 0.5*xenn; 0.5*xenn]*dt
    push!(temp, 1-sum(temp))

    return ProbabilityWeights(temp)
end

# Simulate events for a large number of sectors from time period t to t+1. Entrantid is the latest id of a unique entrant to the economy
function SectorsTransition!(sectordata::Matrix{Sector}, t::Int, mbar::Int, numOfSectors::Int, numofeventslev::Int, numofeventsnn::Int, weightslev::Vector{ProbabilityWeights{Float64,Float64,Array{Float64,1}}}, weightsNN::ProbabilityWeights{Float64,Float64,Array{Float64,1}}, entrantid::Int, λ::Float64, ψ::Float64)

    @inbounds for sectorid = 1:numOfSectors

        sector = sectordata[sectorid,t]
        @unpack_Sector sector

        if m == 0 #Neck-and-neck sectors

            event = sample(1:numofeventsnn, weightsNN)

            if event == 1 #Leader innovation
                sectordata[sectorid, t+1] = Sector(sectorid = sectorid, leader = Firm(firmid = leader.firmid, n = leader.n + 1, q = leader.q * increase_in_productivity(m, λ, ψ)), follower = Firm(firmid = follower.firmid, n = follower.n, q = follower.q), m = 1)
            elseif event == 2 #Follower innovation
                sectordata[sectorid, t+1] = Sector(sectorid = sectorid, leader = Firm(firmid = follower.firmid, n = follower.n + 1, q = follower.q * increase_in_productivity(m, λ, ψ)), follower = Firm(firmid = leader.firmid, n = leader.n, q = leader.q), m = 1)
            elseif event == 3 #Entrant replaces leader
                sectordata[sectorid, t+1] = Sector(sectorid = sectorid, leader = Firm(firmid = entrantid, n = leader.n + 1, q = leader.q * increase_in_productivity(m, λ, ψ)), follower = Firm(firmid = follower.firmid, n = follower.n, q = follower.q), m = 1)
                entrantid += 1
            elseif event == 4 #Entrant replaces follower
                sectordata[sectorid, t+1] = Sector(sectorid = sectorid, leader = Firm(firmid = entrantid, n = follower.n + 1, q = follower.q * increase_in_productivity(m, λ, ψ)), follower = Firm(firmid = leader.firmid, n = leader.n, q = leader.q), m = 1)
                entrantid += 1
            else #Nothing happens
                sectordata[sectorid, t+1] = sector
            end

        else #Unleveled sectors

            # weightss = weightslev[m]
            event = sample(1:numofeventslev, weightslev[m])

            if event == 1 #Leader innovation
                if m < mbar
                    sectordata[sectorid, t+1] = Sector(sectorid = sectorid, leader = Firm(firmid = leader.firmid, n = leader.n + 1, q = leader.q * increase_in_productivity(m, λ, ψ)), follower = Firm(firmid = follower.firmid, n = follower.n, q = follower.q), m = m + 1)
                else
                    sectordata[sectorid, t+1] = Sector(sectorid = sectorid, leader = Firm(firmid = leader.firmid, n = leader.n, q = leader.q), follower = Firm(firmid = follower.firmid, n = follower.n, q = follower.q), m = m)
                end
            elseif event == 2 #Follower drastic innovation
                sectordata[sectorid, t+1] = Sector(sectorid = sectorid, leader = Firm(firmid = leader.firmid, n = leader.n, q = leader.q), follower = Firm(firmid = follower.firmid, n = leader.n, q = leader.q), m = 0)
            elseif event == 3 #Follower incremental innovation
                sectordata[sectorid, t+1] = Sector(sectorid = sectorid, leader = Firm(firmid = leader.firmid, n = leader.n, q = leader.q), follower = Firm(firmid = follower.firmid, n = follower.n + 1, q = follower.q * increase_in_productivity(m-1, λ, ψ)), m = m - 1)
            elseif event == 4 #Entrant replaces follower, drastic innovation
                sectordata[sectorid, t+1] = Sector(sectorid = sectorid, leader = Firm(firmid = leader.firmid, n = leader.n, q = leader.q), follower = Firm(firmid = entrantid, n = leader.n, q = leader.q), m = 0)
                entrantid += 1
            elseif event == 5 #Entrant replaces follower, incremental innovation
                sectordata[sectorid, t+1] = Sector(sectorid = sectorid, leader = Firm(firmid = leader.firmid, n = leader.n, q = leader.q), follower = Firm(firmid = entrantid, n = follower.n + 1, q = follower.q * increase_in_productivity(m-1, λ, ψ)), m = m - 1)
                entrantid += 1
            else #Nothing happens
                sectordata[sectorid, t+1] = sector
            end
        end
    end

    return entrantid
end

# Initiate simulation of sectors from a given set of sectors (sectordatainitial). Entrantid is the latest id of a unique entrant to the economy
function SectorsTransition!(sectordata::Matrix{Sector}, sectordatainitial::Vector{Sector}, mbar::Int, numOfSectors::Int, numofeventslev::Int, numofeventsnn::Int, weightslev::Vector{ProbabilityWeights{Float64,Float64,Array{Float64,1}}}, weightsNN::ProbabilityWeights{Float64,Float64,Array{Float64,1}}, entrantid::Int, λ::Float64, ψ::Float64)

    @inbounds for sectorid = 1:numOfSectors

        sector = sectordatainitial[sectorid]
        @unpack_Sector sector
        t = 0

        if m == 0 #Neck-and-neck sectors

            event = sample(1:numofeventsnn, weightsNN)

            if event == 1 #Leader innovation
                sectordata[sectorid, t+1] = Sector(sectorid = sectorid, leader = Firm(firmid = leader.firmid, n = leader.n + 1, q = leader.q * increase_in_productivity(m, λ, ψ)), follower = Firm(firmid = follower.firmid, n = follower.n, q = follower.q), m = 1)
            elseif event == 2 #Follower innovation
                sectordata[sectorid, t+1] = Sector(sectorid = sectorid, leader = Firm(firmid = follower.firmid, n = follower.n + 1, q = follower.q * increase_in_productivity(m, λ, ψ)), follower = Firm(firmid = leader.firmid, n = leader.n, q = leader.q), m = 1)
            elseif event == 3 #Entrant replaces leader
                sectordata[sectorid, t+1] = Sector(sectorid = sectorid, leader = Firm(firmid = entrantid, n = leader.n + 1, q = leader.q * increase_in_productivity(m, λ, ψ)), follower = Firm(firmid = follower.firmid, n = follower.n, q = follower.q), m = 1)
                entrantid += 1
            elseif event == 4 #Entrant replaces follower
                sectordata[sectorid, t+1] = Sector(sectorid = sectorid, leader = Firm(firmid = entrantid, n = follower.n + 1, q = follower.q * increase_in_productivity(m, λ, ψ)), follower = Firm(firmid = leader.firmid, n = leader.n, q = leader.q), m = 1)
                entrantid += 1
            else #Nothing happens
                sectordata[sectorid, t+1] = sector
            end

        else #Unleveled sectors

            # weightss = weightslev[m]
            event = sample(1:numofeventslev, weightslev[m])

            if event == 1 #Leader innovation
                if m < mbar
                    sectordata[sectorid, t+1] = Sector(sectorid = sectorid, leader = Firm(firmid = leader.firmid, n = leader.n + 1, q = leader.q * increase_in_productivity(m, λ, ψ)), follower = Firm(firmid = follower.firmid, n = follower.n, q = follower.q), m = m + 1)
                else
                    sectordata[sectorid, t+1] = Sector(sectorid = sectorid, leader = Firm(firmid = leader.firmid, n = leader.n, q = leader.q), follower = Firm(firmid = follower.firmid, n = follower.n, q = follower.q), m = m)
                end
            elseif event == 2 #Follower drastic innovation
                sectordata[sectorid, t+1] = Sector(sectorid = sectorid, leader = Firm(firmid = leader.firmid, n = leader.n, q = leader.q), follower = Firm(firmid = follower.firmid, n = leader.n, q = leader.q), m = 0)
            elseif event == 3 #Follower incremental innovation
                sectordata[sectorid, t+1] = Sector(sectorid = sectorid, leader = Firm(firmid = leader.firmid, n = leader.n, q = leader.q), follower = Firm(firmid = follower.firmid, n = follower.n + 1, q = follower.q * increase_in_productivity(m-1, λ, ψ)), m = m - 1)
            elseif event == 4 #Entrant replaces follower, drastic innovation
                sectordata[sectorid, t+1] = Sector(sectorid = sectorid, leader = Firm(firmid = leader.firmid, n = leader.n, q = leader.q), follower = Firm(firmid = entrantid, n = leader.n, q = leader.q), m = 0)
                entrantid += 1
            elseif event == 5 #Entrant replaces follower, incremental innovation
                sectordata[sectorid, t+1] = Sector(sectorid = sectorid, leader = Firm(firmid = leader.firmid, n = leader.n, q = leader.q), follower = Firm(firmid = entrantid, n = follower.n + 1, q = follower.q * increase_in_productivity(m-1, λ, ψ)), m = m - 1)
                entrantid += 1
            else #Nothing happens
                sectordata[sectorid, t+1] = sector
            end
        end
    end

    return entrantid

end

## INITIAL BGP FUNCTIONS
# Construct initial set of sectors from distribution of sectors
function InitialSectorsConstructer(μ::Vector{Float64}, μnn::Float64, mbar::Int, numofsectors::Int, λ::Float64, ψ::Float64)

    # Initilaize sectors
    initialsectors = Array{Sector, 1}(undef, numofsectors);

    # Productivity gaps that we sample from
    M = 0:mbar

    # Corresponding weights
    μWeights = ProbabilityWeights([μnn; μ])

    # Fill initial sectors one by one
    @inbounds for sectorid = 1:numofsectors

        # Draw gap
        m = sample(M, μWeights)

        # Calculate leader productivity until reaching the gap m
        leaderq = 1.0
        for gap = 0:(m-1)
            leaderq = leaderq * increase_in_productivity(gap, λ, ψ)
        end

        # Create firms
        leader = Firm(firmid = sectorid, n = 1 + m, q = leaderq)
        follower = Firm(firmid = numofsectors + sectorid, n = 1, q = 1.0)

        # Create sector
        initialsectors[sectorid] = Sector(sectorid = sectorid, leader = leader, follower = follower, m = m)
    end

    return initialsectors
end

# Given functions above, simulate sectors for a specified time period (numericp.tperiods)
function SimulateSectorDataInitialBgp!(sectordata::Matrix{Sector}, weightslev::Vector{ProbabilityWeights{Float64,Float64,Array{Float64,1}}}, initial::BGP, numericp::NumParam)

    @unpack_BGP initial
    @unpack_NumParam numericp

    ## Initialize entrand ids. entrantid lets us assign unique ids for each entrant
    entrantid = numoffirms + 1

    ## Create constant weights
    WeightsLevConstructer!(weightslev, xl, xf, xe, mbar, dt, p.ϕ, p.δ, p.ϕTilde)
    weightsNN = WeightsNNConstructer(xnn, xenn, dt)

    ## Initialize simulation
    sectordatainitial = InitialSectorsConstructer(μ, μnn, mbar, numofsectors, p.λ, p.ψ)
    entrantid = SectorsTransition!(sectordata, sectordatainitial, mbar, numofsectors, numofeventslev, numofeventsnn, weightslev, weightsNN, entrantid, p.λ, p.ψ)

    ## Continue simulation for t > 1. Every period, recreate weights. Note that SectorsTransition! simulates sectors in t+1 given sectors in t (t => t+1)
    for t = 1:tperiods-1

        # Simulate t => t+1 by updating entrantid
        entrantid = SectorsTransition!(sectordata, t, mbar, numofsectors, numofeventslev, numofeventsnn, weightslev, weightsNN, entrantid, p.λ, p.ψ)
    end

    return entrantid
end

# This function creates an annual panel dataset given simulated sectordata. This panel is essentially a matrix with the size ((numericp.tyears * numericp.numoffirms) × numericp.panelnumofcolumns). It contains all the firms and relevant data such as age. Each year corresponds to the snapshot of simulated sectors at the end of the year.
function AnnualPanelDataInitialBgp!(panel::Matrix{Float64}, sectordata::Matrix{Sector}, numericp::NumParam, ll::Vector{Float64}, lf::Vector{Float64}, lnn::Float64, hl::Vector{Float64}, hf::Vector{Float64}, hnn::Float64)

    # Panel columns are as follows: 1 = year, 2 = firm id, 3 = firm step, 4 = firm total labor, 5 = age, 6 = productivity

    @unpack_NumParam numericp

    obsid = 0

    # Visit all years, and the last period of that year
    for tyear = 1:tyears

        # Last period of the year
        tperiod = tyear * periodsperyear

        @inbounds for i = 1:numofsectors # i is also the sector id

            # Sector id = i in period tperiod
            sector = sectordata[i,tperiod]
            leader = sector.leader
            follower = sector.follower
            m = sector.m

            # Extract variables. Where labors are extracted depends on gap m

            # If unleveled
            if m > 0

                # Leader
                obsid += 1
                panel[obsid, yearcol] = tyear + 1830
                panel[obsid, firmidcol] = leader.firmid
                panel[obsid, stepcol] = leader.n
                panel[obsid, laborcol] = ll[m] + hl[m]
                panel[obsid, prodcol] = leader.q

                # Follower
                obsid += 1
                panel[obsid, yearcol] = tyear + 1830
                panel[obsid, firmidcol] = follower.firmid
                panel[obsid, stepcol] = follower.n
                panel[obsid, laborcol] = lf[m] + hf[m]
                panel[obsid, prodcol] = follower.q

            else # If neck-and-neck

                # Leader
                obsid += 1
                panel[obsid, yearcol] = tyear + 1830
                panel[obsid, firmidcol] = leader.firmid
                panel[obsid, stepcol] = leader.n
                panel[obsid, laborcol] = lnn + hnn
                panel[obsid, prodcol] = leader.q

                # Follower
                obsid += 1
                panel[obsid, yearcol] = tyear + 1830
                panel[obsid, firmidcol] = follower.firmid
                panel[obsid, stepcol] = follower.n
                panel[obsid, laborcol] = lnn + hnn
                panel[obsid, prodcol] = follower.q
            end
        end
    end
end

# Calculate firm ages given the annual panel dataset panel.
function CalculateAgesInitialBgp!(panel::Matrix{Float64}, numofsectors::Int, numoffirms::Int, tyears::Int, firmidcol::Int, agecol::Int)

    # All the firms in year 1 have age 1
    panel[1:numoffirms, agecol] .= 1.0

    # Years > 1
    @inbounds for tyear = 2:tyears

        # Extract part of the panel corresponding to this year tyear, and the last year tyear-1
        subpanelyesterday = @view panel[(tyear-2)*numoffirms+1:(tyear-1)*numoffirms, :]
        subpaneltoday = @view panel[(tyear-1)*numoffirms+1:tyear*numoffirms, :]

        # Run over sectors, since a firm can belong to only one sector throughout its life cycle.
        @inbounds for j = 1:numofsectors

            # In each sector, there are two kinds of firms. In the panel their position is known due to its construction
            leaderindex = (j-1)*2 + 1
            followerindex = (j-1)*2 + 2

            # Firm ids this year in the sector that we try to match with last year firm ids in the same sector
            leaderid = subpaneltoday[leaderindex, firmidcol]
            followerid = subpaneltoday[followerindex, firmidcol]

            # Firm ids of the same sector last year
            leaderidyesterday = subpanelyesterday[leaderindex, firmidcol]
            followeridyesterday = subpanelyesterday[followerindex, firmidcol]

            # Firm ages of the sector last year
            leaderageyesterday = subpanelyesterday[leaderindex, agecol]
            followerageyesterday = subpanelyesterday[followerindex, agecol]

            # Try to match firms
            if leaderid == leaderidyesterday # Leader this year was the leader of the sector last year
                subpaneltoday[leaderindex, agecol] = leaderageyesterday + 1.0
            elseif leaderid == followeridyesterday # Leader this year was the follower of the sector last year
                subpaneltoday[leaderindex, agecol] = followerageyesterday + 1.0
            else # Leader this year is a new firm
                subpaneltoday[leaderindex, agecol] = 1.0
            end

            if followerid == leaderidyesterday # Follower this year was the leader of the sector last year
                subpaneltoday[followerindex, agecol] = leaderageyesterday + 1.0
            elseif followerid == followeridyesterday # Follower this year was the follower of the sector last year
                subpaneltoday[followerindex, agecol] = followerageyesterday + 1.0
            else # Follower this year is a new firm
                subpaneltoday[followerindex, agecol] = 1.0
            end
        end

    end
end

# Final simulation function for an initial bgp
function InitialBgpSimulation!(panel::Matrix{Float64}, sectordata::Matrix{Sector}, initial::BGP)

    @unpack ω, xl, xf, xnn, numericp = initial
    @unpack α, γ = initial.p
    @unpack numofsectors, numoffirms, tyears, firmidcol, agecol, mbar, numofeventslev = numericp
    @unpack sl, sf, snn, zl, zf, znn = initial

    # Create an empty vector of probabilities (weights)
    weightslev = [ProbabilityWeights(zeros(numofeventslev)) for m = 1:mbar]

    # Labors
    ll = ProdLabor.(sl, zl, ω)
    lf = ProdLabor.(sf, zf, ω)
    lnn = ProdLabor(snn, znn, ω)
    hl = RndLabor.(xl, α, γ)
    hf = RndLabor.(xf, α, γ)
    hnn = RndLabor(xnn, α, γ)

    # Simulation
    entrantid = SimulateSectorDataInitialBgp!(sectordata, weightslev, initial, numericp)

    # Create annual panel data
    AnnualPanelDataInitialBgp!(panel, sectordata, numericp, ll, lf, lnn, hl, hf, hnn)

    # Calculate ages
    CalculateAgesInitialBgp!(panel, numofsectors, numoffirms, tyears, firmidcol, agecol)

    return entrantid
end



## TRANSITION FUNCTIONS
# Simulate sectors over the transition starting from previous simulation of the baseline initial bgp (sectordatainitialbgplastperiod)
function SimulateSectorDataTransition!(sectordata::Matrix{Sector}, weightslev::Vector{ProbabilityWeights{Float64,Float64,Array{Float64,1}}}, transition::Transition, initial::BGP, numericp::NumParam, sequencep::Vector{StrParam}, sectordatainitialbgplastperiod::Vector{Sector}, entrantid::Int)

    @unpack_Transition transition
    @unpack_NumParam numericp

    ## Create bgp weights
    WeightsLevConstructer!(weightslev, initial.xl, initial.xf, initial.xe, mbar, dt, initial.p.ϕ, initial.p.δ, initial.p.ϕTilde)
    weightsNN = WeightsNNConstructer(initial.xnn, initial.xenn, dt)

    ## Initialize sectordata from the last period of initial bgp's simulated sectordata (sectordatainitialbgplastperiod) in order to find the initial distribution of transition
    entrantid = SectorsTransition!(sectordata, sectordatainitialbgplastperiod, mbar, numofsectors, numofeventslev, numofeventsnn, weightslev, weightsNN, entrantid, initial.p.λ, initial.p.ψ)

    ## Continue simulation for t > 1. Every period, recreate weights. Note that SectorsTransition! simulates sectors in t+1 given sectors in t (t => t+1).
    for t = 1:tperiods-1

        # Get parameters in t
        p = sequencep[t]

        # Create weights
        WeightsLevConstructer!(weightslev, xl[:,t], xf[:,t], xe[:,t], mbar, dt, p.ϕ, p.δ, p.ϕTilde)
        weightsNN = WeightsNNConstructer(xnn[t], xenn[t], dt)

        # Simulate t => t+1 by updating entrantid
        entrantid = SectorsTransition!(sectordata, t, mbar, numofsectors, numofeventslev, numofeventsnn, weightslev, weightsNN, entrantid, p.λ, p.ψ)
    end

    return entrantid
end

# Similar to AnnualPanelDataInitialBgp!. Only calendar years are different
function AnnualPanelDataTransition!(panel::Matrix{Float64}, sectordata::Matrix{Sector}, numericp::NumParam, ll::Matrix{Float64}, lf::Matrix{Float64}, lnn::Vector{Float64}, hl::Matrix{Float64}, hf::Matrix{Float64}, hnn::Vector{Float64})

    # Panel columns are as follows: 1 = year, 2 = firm id, 3 = firm step, 4 = firm total labor, 5 = age

    @unpack_NumParam numericp

    # obsid = numoffirms # Initial bgp's simulated panel is appended to the beginning. So obsid starts accordingly. obsid = row number
    obsid = 0

    # Visit all years, and the last period of that year
    for tyear = 1:tyears

        # Last period of the year
        tperiod = tyear * periodsperyear

        @inbounds for i = 1:numofsectors # i is also the sector id

            # Sector id = i in period tperiod
            sector = sectordata[i,tperiod]
            leader = sector.leader
            follower = sector.follower
            m = sector.m

            # Extract variables. Where labors are extracted depends on gap m

            # If unleveled
            if m > 0

                # Leader
                obsid += 1
                panel[obsid, yearcol] = tyear + 1980
                panel[obsid, firmidcol] = leader.firmid
                panel[obsid, stepcol] = leader.n
                panel[obsid, laborcol] = ll[m, tperiod] + hl[m, tperiod]
                panel[obsid, prodcol] = leader.q

                # Follower
                obsid += 1
                panel[obsid, yearcol] = tyear + 1980
                panel[obsid, firmidcol] = follower.firmid
                panel[obsid, stepcol] = follower.n
                panel[obsid, laborcol] = lf[m, tperiod] + hf[m, tperiod]
                panel[obsid, prodcol] = follower.q

            else # If neck-and-neck

                # Leader
                obsid += 1
                panel[obsid, yearcol] = tyear + 1980
                panel[obsid, firmidcol] = leader.firmid
                panel[obsid, stepcol] = leader.n
                panel[obsid, laborcol] = lnn[tperiod] + hnn[tperiod]
                panel[obsid, prodcol] = leader.q

                # Follower
                obsid += 1
                panel[obsid, yearcol] = tyear + 1980
                panel[obsid, firmidcol] = follower.firmid
                panel[obsid, stepcol] = follower.n
                panel[obsid, laborcol] = lnn[tperiod] + hnn[tperiod]
                panel[obsid, prodcol] = follower.q
            end
        end
    end
end

# Continue calculating ages starting from previous simulation of the underlying initial bgp
function CalculateAgesTransition!(panel::Matrix{Float64}, numofsectors::Int, numoffirms::Int, tyears::Int, firmidcol::Int, agecol::Int, panelinitialbgplastyear::Matrix{Float64})

    # Calculate ages for the first year (1981) given panelinitialbgplastyear (ages in 1980)
    subpaneltoday = @view panel[1:numoffirms, :] #panel1981

    @inbounds for j = 1:numofsectors

        # In each sector, there are two kinds of firms. In the panel their position is known due to its construction
        leaderindex = (j-1)*2 + 1
        followerindex = (j-1)*2 + 2

        # Firm ids this year in the sector that we try to match with last year firm ids in the same sector
        leaderid = subpaneltoday[leaderindex, firmidcol]
        followerid = subpaneltoday[followerindex, firmidcol]

        # Firm ids of the same sector last year
        leaderidyesterday = panelinitialbgplastyear[leaderindex, firmidcol]
        followeridyesterday = panelinitialbgplastyear[followerindex, firmidcol]

        # Firm ages of the sector last year
        leaderageyesterday = panelinitialbgplastyear[leaderindex, agecol]
        followerageyesterday = panelinitialbgplastyear[followerindex, agecol]

        # Try to match firms
        if leaderid == leaderidyesterday # Leader this year was the leader of the sector last year
            subpaneltoday[leaderindex, agecol] = leaderageyesterday + 1.0
        elseif leaderid == followeridyesterday # Leader this year was the follower of the sector last year
            subpaneltoday[leaderindex, agecol] = followerageyesterday + 1.0
        else # Leader this year is a new firm
            subpaneltoday[leaderindex, agecol] = 1.0
        end

        if followerid == leaderidyesterday # Follower this year was the leader of the sector last year
            subpaneltoday[followerindex, agecol] = leaderageyesterday + 1.0
        elseif followerid == followeridyesterday # Follower this year was the follower of the sector last year
            subpaneltoday[followerindex, agecol] = followerageyesterday + 1.0
        else # Follower this year is a new firm
            subpaneltoday[followerindex, agecol] = 1.0
        end
    end

    # Years >= 2 (since 1981 until 2129 - if numericp.tyears=150)
    @inbounds for tyear = 2:tyears

        # Extract part of the panel corresponding to this year tyear, and the last year tyear-1
        subpanelyesterday = @view panel[(tyear-2)*numoffirms+1:(tyear-1)*numoffirms, :]
        subpaneltoday = @view panel[(tyear-1)*numoffirms+1:tyear*numoffirms, :]

        # Run over sectors, since a firm can belong to only one sector throught its life cycle.
        @inbounds for j = 1:numofsectors

            # In each sector, there are two kinds of firms. In the panel their position is known due to its construction
            leaderindex = (j-1)*2 + 1
            followerindex = (j-1)*2 + 2

            # Firm ids this year in the sector that we try to match with last year firm ids in the same sector
            leaderid = subpaneltoday[leaderindex, firmidcol]
            followerid = subpaneltoday[followerindex, firmidcol]

            # Firm ids of the same sector last year
            leaderidyesterday = subpanelyesterday[leaderindex, firmidcol]
            followeridyesterday = subpanelyesterday[followerindex, firmidcol]

            # Firm ages of the sector last year
            leaderageyesterday = subpanelyesterday[leaderindex, agecol]
            followerageyesterday = subpanelyesterday[followerindex, agecol]

            # Try to match firms
            if leaderid == leaderidyesterday # Leader this year was the leader of the sector last year
                subpaneltoday[leaderindex, agecol] = leaderageyesterday + 1.0
            elseif leaderid == followeridyesterday # Leader this year was the follower of the sector last year
                subpaneltoday[leaderindex, agecol] = followerageyesterday + 1.0
            else # Leader this year is a new firm
                subpaneltoday[leaderindex, agecol] = 1.0
            end

            if followerid == leaderidyesterday # Follower this year was the leader of the sector last year
                subpaneltoday[followerindex, agecol] = leaderageyesterday + 1.0
            elseif followerid == followeridyesterday # Follower this year was the follower of the sector last year
                subpaneltoday[followerindex, agecol] = followerageyesterday + 1.0
            else # Follower this year is a new firm
                subpaneltoday[followerindex, agecol] = 1.0
            end
        end

    end
end

# Final simulation function for the transition of the model.
function TransitionSimulation!(panel::Matrix{Float64}, sectordata::Matrix{Sector}, transition::Transition, initial::BGP, entrantid::Int)

    @unpack_Transition transition
    @unpack_StrParam initial.p # Constant parameters over transition
    @unpack_NumParam numericp
    @unpack sl, sf, snn, zl, zf, znn = initial # Revenue shares and markups are constant over time

    # Create an empty vector of probabilities (weights)
    weightslev = [ProbabilityWeights(zeros(numofeventslev)) for m = 1:mbar]

    # Back out latest period and latest year from previous sectordata and panel
    sectordatainitialbgplastperiod = sectordata[:, end]
    panelinitialbgplastyear = panel[end - (numoffirms-1) : end, :]

    # Fundamentals: Correct way to do is to use parameter sequence. But for now, we know that these parameters are constant
    ll = ProdLabor.(sl, zl, ω')
    lf = ProdLabor.(sf, zf, ω')
    lnn = ProdLabor.(snn, znn, ω)
    hl = RndLabor.(xl, α, γ)
    hf = RndLabor.(xf, α, γ)
    hnn = RndLabor.(xnn, α, γ)

    # Simulation
    entrantid = SimulateSectorDataTransition!(sectordata, weightslev, transition, initial, numericp, sequencep, sectordatainitialbgplastperiod, entrantid)

    # Create annual panel data
    AnnualPanelDataTransition!(panel, sectordata, numericp, ll, lf, lnn, hl, hf, hnn)

    # Calculate ages
    CalculateAgesTransition!(panel, numofsectors, numoffirms, tyears, firmidcol, agecol, panelinitialbgplastyear)

    return entrantid
end

# This is the simulation function for transition where continuity between transition and initial bgp simulations are achieved
function TransitionSimulation_old!(panel::Matrix{Float64}, sectordata::Matrix{Sector}, transition::Transition, initial::BGP, entrantid::Int, panelinitialbgplastyear::Matrix{Float64}, sectordatainitialbgplastperiod::Vector{Sector})

    @unpack_Transition transition
    @unpack_StrParam initial.p # Constant parameters over transition
    @unpack_NumParam numericp
    @unpack sl, sf, snn, zl, zf, znn = initial # Revenue shares and markups are constant over time

    # Create an empty vector of probabilities (weights)
    weightslev = [ProbabilityWeights(zeros(numofeventslev)) for m = 1:mbar]

    # Fundamentals: Correct way to do is to use parameter sequence. But for now, we know that these parameters are constant
    ll = ProdLabor.(sl, zl, ω')
    lf = ProdLabor.(sf, zf, ω')
    lnn = ProdLabor.(snn, znn, ω)
    hl = RndLabor.(xl, α, γ)
    hf = RndLabor.(xf, α, γ)
    hnn = RndLabor.(xnn, α, γ)

    # Simulation
    entrantid = SimulateSectorDataTransition!(sectordata, weightslev, transition, initial, numericp, sequencep, sectordatainitialbgplastperiod, entrantid)

    # Create annual panel data
    AnnualPanelDataTransition!(panel, sectordata, numericp, ll, lf, lnn, hl, hf, hnn)

    # Calculate ages
    CalculateAgesTransition!(panel, numofsectors, numoffirms, tyears, firmidcol, agecol, panelinitialbgplastyear)

    return entrantid
end

function TransitionSimulation!(panel::Matrix{Float64}, sectordata::Matrix{Sector}, transition::Transition, entrantid::Int, panelinitialbgplastyear::Matrix{Float64}, sectordatainitialbgplastperiod::Vector{Sector})

    @unpack_Transition transition
    @unpack_NumParam numericp

    # Create an empty vector of probabilities (weights)
    weightslev = [ProbabilityWeights(zeros(numofeventslev)) for m = 1:mbar]

    # Labors
    ll = Matrix{Float64}(undef, mbar, tperiods)
    lf = similar(ll)
    lnn = Vector{Float64}(undef, tperiods)
    hl = similar(ll)
    hf = similar(ll)
    hnn = similar(lnn)

    for t = 1:tperiods
        @unpack α, γ = sequencep[t]
        
        ll[:, t] = ProdLabor.(sl[:, t], zl[:, t], ω[t]) 
        lf[:, t] = ProdLabor.(sf[:, t], zf[:, t], ω[t]) 
        lnn[t] = ProdLabor(snn[t], znn[t], ω[t]) 
        hl[:, t] = RndLabor.(xl[:, t], α, γ) 
        hf[:, t] = RndLabor.(xf[:, t], α, γ) 
        hnn[t] = RndLabor(xnn[t], α, γ)
    end

    # Simulation
    entrantid = SimulateSectorDataTransition!(sectordata, weightslev, transition, initialbgp, numericp, sequencep, sectordatainitialbgplastperiod, entrantid)

    # Create annual panel data
    AnnualPanelDataTransition!(panel, sectordata, numericp, ll, lf, lnn, hl, hf, hnn)

    # Calculate ages
    CalculateAgesTransition!(panel, numofsectors, numoffirms, tyears, firmidcol, agecol, panelinitialbgplastyear)

    return entrantid
end



## TERMINAL BGP FUNCTIONS
function SimulateSectorDataTerminalBgp!(sectordata::Matrix{Sector}, weightslev::Vector{ProbabilityWeights{Float64,Float64,Array{Float64,1}}}, terminal::BGP, transition::Transition, numericp::NumParam, sequencep::Vector{StrParam}, sctordatatransitionlastperiod::Vector{Sector}, entrantid::Int)

    @unpack_BGP terminal
    @unpack_NumParam numericp

    ## Create weights for the transition from last period of 2129 to the first period of terminal bgp which is 2130
    WeightsLevConstructer!(weightslev, transition.xl[:,end], transition.xf[:,end], transition.xe[:,end], mbar, dt, sequencep[end].ϕ, sequencep[end].δ, sequencep[end].ϕTilde)
    weightsNN = WeightsNNConstructer(transition.xnn[end], transition.xenn[end], dt)

    ## Initialize simulation
    entrantid = SectorsTransition!(sectordata, sctordatatransitionlastperiod, mbar, numofsectors, numofeventslev, numofeventsnn, weightslev, weightsNN, entrantid, sequencep[end].λ, sequencep[end].ψ)

    ## Create constant weights for terminal bgp transition from 2130 and on
    WeightsLevConstructer!(weightslev, xl, xf, xe, mbar, dt, p.ϕ, p.δ, p.ϕTilde)
    weightsNN = WeightsNNConstructer(xnn, xenn, dt)

    ## Continue simulation for t > 1. Every period, recreate weights. Note that SectorsTransition! simulates sectors in t+1 given sectors in t (t => t+1)
    for t = 1:tperiods-1

        # Simulate t => t+1 by updating entrantid
        entrantid = SectorsTransition!(sectordata, t, mbar, numofsectors, numofeventslev, numofeventsnn, weightslev, weightsNN, entrantid, p.λ, p.ψ)
    end

    return entrantid
end

function AnnualPanelDataTerminalBgp!(panel::Matrix{Float64}, sectordata::Matrix{Sector}, numericp::NumParam, ll::Vector{Float64}, lf::Vector{Float64}, lnn::Float64, hl::Vector{Float64}, hf::Vector{Float64}, hnn::Float64)

    # Panel columns are as follows: 1 = year, 2 = firm id, 3 = firm step, 4 = firm total labor, 5 = age

    @unpack_NumParam numericp

    obsid = 0

    # Visit all years, and the last period of that year
    for tyear = 1:tyears

        # Last period of the year
        tperiod = tyear * periodsperyear

        @inbounds for i = 1:numofsectors # i is also the sector id

            # Sector id = i in period tperiod
            sector = sectordata[i,tperiod]
            leader = sector.leader
            follower = sector.follower
            m = sector.m

            # Extract variables. Where labors are extracted depends on gap m

            # If unleveled
            if m > 0

                # Leader
                obsid += 1
                panel[obsid, yearcol] = tyear + 1980 + tyears
                panel[obsid, firmidcol] = leader.firmid
                panel[obsid, stepcol] = leader.n
                panel[obsid, laborcol] = ll[m] + hl[m]
                panel[obsid, prodcol] = leader.q

                # Follower
                obsid += 1
                panel[obsid, yearcol] = tyear + 1980 + tyears
                panel[obsid, firmidcol] = follower.firmid
                panel[obsid, stepcol] = follower.n
                panel[obsid, laborcol] = lf[m] + hf[m]
                panel[obsid, prodcol] = follower.q

            else # If neck-and-neck

                # Leader
                obsid += 1
                panel[obsid, yearcol] = tyear + 1980 + tyears
                panel[obsid, firmidcol] = leader.firmid
                panel[obsid, stepcol] = leader.n
                panel[obsid, laborcol] = lnn + hnn
                panel[obsid, prodcol] = leader.q

                # Follower
                obsid += 1
                panel[obsid, yearcol] = tyear + 1980 + tyears
                panel[obsid, firmidcol] = follower.firmid
                panel[obsid, stepcol] = follower.n
                panel[obsid, laborcol] = lnn + hnn
                panel[obsid, prodcol] = follower.q
            end
        end
    end
end

function CalculateAgesTerminalBgp!(panel::Matrix{Float64}, numofsectors::Int, numoffirms::Int, tyears::Int, firmidcol::Int, agecol::Int, paneltransitionlastyear::Matrix{Float64})

    # Calculate ages for the first year (2131) given paneltransitionlastyear (ages in 2130)
    subpaneltoday = @view panel[1:numoffirms, :] #panel2131

    @inbounds for j = 1:numofsectors

        # In each sector, there are two kinds of firms. In the panel their position is known due to its construction
        leaderindex = (j-1)*2 + 1
        followerindex = (j-1)*2 + 2

        # Firm ids this year in the sector that we try to match with last year firm ids in the same sector
        leaderid = subpaneltoday[leaderindex, firmidcol]
        followerid = subpaneltoday[followerindex, firmidcol]

        # Firm ids of the same sector last year
        leaderidyesterday = paneltransitionlastyear[leaderindex, firmidcol]
        followeridyesterday = paneltransitionlastyear[followerindex, firmidcol]

        # Firm ages of the sector last year
        leaderageyesterday = paneltransitionlastyear[leaderindex, agecol]
        followerageyesterday = paneltransitionlastyear[followerindex, agecol]

        # Try to match firms
        if leaderid == leaderidyesterday # Leader this year was the leader of the sector last year
            subpaneltoday[leaderindex, agecol] = leaderageyesterday + 1.0
        elseif leaderid == followeridyesterday # Leader this year was the follower of the sector last year
            subpaneltoday[leaderindex, agecol] = followerageyesterday + 1.0
        else # Leader this year is a new firm
            subpaneltoday[leaderindex, agecol] = 1.0
        end

        if followerid == leaderidyesterday # Follower this year was the leader of the sector last year
            subpaneltoday[followerindex, agecol] = leaderageyesterday + 1.0
        elseif followerid == followeridyesterday # Follower this year was the follower of the sector last year
            subpaneltoday[followerindex, agecol] = followerageyesterday + 1.0
        else # Follower this year is a new firm
            subpaneltoday[followerindex, agecol] = 1.0
        end
    end

    # Years >= 2 (since 2132 until 2280 - if numericp.tyears = 150)
    @inbounds for tyear = 2:tyears

        # Extract part of the panel corresponding to this year tyear, and the last year tyear-1
        subpanelyesterday = @view panel[(tyear-2)*numoffirms+1:(tyear-1)*numoffirms, :]
        subpaneltoday = @view panel[(tyear-1)*numoffirms+1:tyear*numoffirms, :]

        # Run over sectors, since a firm can belong to only one sector throught its life cycle.
        @inbounds for j = 1:numofsectors

            # In each sector, there are two kinds of firms. In the panel their position is known due to its construction
            leaderindex = (j-1)*2 + 1
            followerindex = (j-1)*2 + 2

            # Firm ids this year in the sector that we try to match with last year firm ids in the same sector
            leaderid = subpaneltoday[leaderindex, firmidcol]
            followerid = subpaneltoday[followerindex, firmidcol]

            # Firm ids of the same sector last year
            leaderidyesterday = subpanelyesterday[leaderindex, firmidcol]
            followeridyesterday = subpanelyesterday[followerindex, firmidcol]

            # Firm ages of the sector last year
            leaderageyesterday = subpanelyesterday[leaderindex, agecol]
            followerageyesterday = subpanelyesterday[followerindex, agecol]

            # Try to match firms
            if leaderid == leaderidyesterday # Leader this year was the leader of the sector last year
                subpaneltoday[leaderindex, agecol] = leaderageyesterday + 1.0
            elseif leaderid == followeridyesterday # Leader this year was the follower of the sector last year
                subpaneltoday[leaderindex, agecol] = followerageyesterday + 1.0
            else # Leader this year is a new firm
                subpaneltoday[leaderindex, agecol] = 1.0
            end

            if followerid == leaderidyesterday # Follower this year was the leader of the sector last year
                subpaneltoday[followerindex, agecol] = leaderageyesterday + 1.0
            elseif followerid == followeridyesterday # Follower this year was the follower of the sector last year
                subpaneltoday[followerindex, agecol] = followerageyesterday + 1.0
            else # Follower this year is a new firm
                subpaneltoday[followerindex, agecol] = 1.0
            end
        end

    end
end

function TerminalBgpSimulation!(panel::Matrix{Float64}, sectordata::Matrix{Sector}, terminal::BGP, transition::Transition, entrantid::Int)

    @unpack ω, xl, xf, xnn, numericp = terminal
    @unpack α, γ = terminal.p
    @unpack numofsectors, numoffirms, tyears, firmidcol, agecol, mbar, numofeventslev = numericp
    @unpack sl, sf, snn, zl, zf, znn = terminal
    @unpack sequencep = transition

    # Create an empty vector of probabilities (weights)
    weightslev = [ProbabilityWeights(zeros(numofeventslev)) for m = 1:mbar]

    # Back out latest period and latest year from previous sectordata and panel
    sectordatatransitionlastperiod = sectordata[:, end]
    paneltransitionlastyear = panel[end-(numoffirms - 1) : end, :]

    # Labors
    ll = ProdLabor.(sl, zl, ω)
    lf = ProdLabor.(sf, zf, ω)
    lnn = ProdLabor(snn, znn, ω)
    hl = RndLabor.(xl, α, γ)
    hf = RndLabor.(xf, α, γ)
    hnn = RndLabor(xnn, α, γ)

    # Simulation
    entrantid = SimulateSectorDataTerminalBgp!(sectordata, weightslev, terminal, transition, numericp, sequencep, sectordatatransitionlastperiod, entrantid)

    # Create annual panel data
    AnnualPanelDataTerminalBgp!(panel, sectordata, numericp, ll, lf, lnn, hl, hf, hnn)

    # Calculate ages
    CalculateAgesTerminalBgp!(panel, numofsectors, numoffirms, tyears, firmidcol, agecol, paneltransitionlastyear)

    return entrantid
end