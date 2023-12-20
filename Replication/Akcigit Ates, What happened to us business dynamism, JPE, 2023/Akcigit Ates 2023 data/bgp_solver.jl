## BGP (Balanced Growth Path) struct that contains endogenous variables
@with_kw struct BGP
    p::StrParam                 # Parameters of the BGP
    numericp::NumParam          # Numeric parameters of the BGP
    
    ω::Float64                  # Normalized wage or the labor share

    μ::Vector{Float64}          # Sector distribution for the unleveled sectors
    μnn::Float64                # Density of leveled (neck-and-neck) sectors 

    xl::Vector{Float64}         # Innovation rates of leaders
    xf::Vector{Float64}         # Innovation rates of followers
    xnn::Float64                # Innovation rate of neck-and-neck firms
    xe::Vector{Float64}         # Innovation rates of potential entrants in unleveled sectors
    xenn::Float64               # Innovation rates of potential entrants in neck-and-neck sectors

    vl::Vector{Float64}         # Value function (vector) of leaders
    vf::Vector{Float64}         # Value function of followers
    vnn::Float64                # Value of being neck-and-neck

    sl::Vector{Float64}         # Revenue share of leaders
    sf::Vector{Float64}         # Revenue share of followers
    snn::Float64                # Revenue share of neck-and-neck firms
    zl::Vector{Float64}         # Markups of leaders
    zf::Vector{Float64}         # Markups of followers
    znn::Float64                # Markups of neck-and-neck firms
end

## INNOVATION FUNCTIONS: Following functions calculate innovation rates given value functions

# Leader innovation
function InnL!(xL::AbstractVector{Float64}, vL::Vector{Float64}, m::Int, ω::Float64, s::Float64, α::Float64, γ::Float64, mbar::Int)

    if m == mbar
        vGain = 0.0
    else
        temp = vL[m+1] - vL[m]
        vGain = temp >= 0.0 ? temp : 0.0
    end

    xL[m] = (vGain/((1-s)*α*ω))^(1/(γ-1))
end

# Follower innovation
function InnF!(xF::AbstractVector{Float64}, vF::Vector{Float64}, vNN::Float64, m::Int, ω::Float64, s::Float64, α::Float64, γ::Float64, ϕ::Float64)

    if m == 1
        temp = vNN - vF[m]
        vGain = temp >= 0.0 ? temp : 0.0
    else
        temp = ϕ*vNN + (1-ϕ)*vF[m-1] - vF[m]
        vGain = temp >= 0.0 ? temp : 0.0
    end

    xF[m] = (vGain/((1-s)*α*ω))^(1/(γ-1))
end

# Neck-and-neck innovation
function InnNN(vL::Float64, vNN::Float64, ω::Float64, s::Float64, α::Float64, γ::Float64) #vL = vL[1]

    temp = vL - vNN
    vGain = temp >= 0.0 ? temp : 0.0

    xNN = (vGain/((1-s)*α*ω))^(1/(γ-1))

    return xNN
end

# Entrant innovation in unleveled sectors
function InnE!(xE::AbstractVector{Float64}, vL::Vector{Float64}, vF::Vector{Float64}, vNN::Float64, m::Int, ω::Float64, αTilde::Float64, γTilde::Float64, ϕTilde::Float64, mbar::Int)

    if m != mbar && m > 1
        temp = ϕTilde*vNN + (1-ϕTilde)*vF[m-1]
        vGain = temp >= 0.0 ? temp : 0.0
    elseif m == mbar
        temp = ϕTilde*vNN + (1-ϕTilde)*vF[m-1]
        vGain = temp >= 0.0 ? temp : 0.0
    elseif m == 1
        temp = vNN
        vGain = temp >= 0.0 ? temp : 0.0
    end

    xE[m] = (vGain/(αTilde*ω))^(1/(γTilde-1))
end

# Entrant innovation in neck-and-neck sectors
function InnENN(vL::Float64, ω::Float64, αTilde::Float64, γTilde::Float64) #vL = vL[1]
    vGain = vL >= 0.0 ? vL : 0.0

    xENN = (vGain/(αTilde*ω))^(1/(γTilde-1))

    return xENN
end

## BELLMAN EQUATIONS: ρ*v′ - (v′-v)/dt = A′ => v = (A′ - (ρ-1/dt)*v′)*dt.
# Value functions: vL, vF, vNN (both Next (′) and Prev (without ′) of all of them).
# Innovations: xL′, xF′, xNN′, xE′, xENN′

function LeaderBellman!(vL::AbstractVector{Float64}, vL′::Vector{Float64}, vNN′::Float64, xL′::Vector{Float64}, xF′::Vector{Float64}, xE′::Vector{Float64}, ω::Float64, p::StrParam, numericp::NumParam, sL::Vector{Float64}, A′::Vector{Float64})

    @unpack_StrParam p
    @unpack_NumParam numericp

    @inbounds for m = eachindex(vL)
        value = vL′[m]
        valueUp = m < mbar ? vL′[m+1] : value
        valueDown = m > 1 ? vL′[m-1] : vNN′

        xL = xL′[m]
        xF = xF′[m]
        xE = xE′[m]

        rndCost = RndLabor(xL, α, γ)*ω

        A′[m] = (1-τ)*ProfitNorm(sL[m], β) - (1-s)*rndCost + xL*(valueUp - value) + (ϕ*xF + δ)*(vNN′ - value) + (1-ϕ)*xF*(valueDown - value) + xE*(ϕTilde*vNN′ + (1-ϕTilde)*valueDown - value)

        vL[m] = (A′[m] - (ρ - 1 / dt) * vL′[m]) * dt
    end
end

function FollowerBellman!(vF::AbstractVector{Float64}, vF′::Vector{Float64}, vNN′::Float64, xL′::Vector{Float64}, xF′::Vector{Float64}, xE′::Vector{Float64}, ω::Float64, p::StrParam, numericp::NumParam, sF::Vector{Float64}, A′::Vector{Float64})

    @unpack_StrParam p
    @unpack_NumParam numericp

    @inbounds for m = eachindex(vF)
        value = vF′[m]
        valueUp = m > 1 ? vF′[m-1] : vNN′
        valueDown = m < mbar ? vF′[m+1] : value

        xL = xL′[m]
        xF = xF′[m]
        xE = xE′[m]

        rndCost = RndLabor(xF, α, γ)*ω

        A′[m] = (1-τ)*ProfitNorm(sF[m], β) - (1-s)*rndCost + (ϕ*xF + δ)*(vNN′ - value) + (1-ϕ)*xF*(valueUp - value) + xL*(valueDown - value) + xE*(0.0 - value)

        vF[m] = (A′[m] - (ρ - 1 / dt) * vF′[m]) * dt
    end
end

function NNBellman(vL′::Float64, vF′::Float64, vNN′::Float64, xNN′::Float64, xENN′::Float64, ω::Float64, p::StrParam, numericp::NumParam, sNN::Float64) #vL′ = vL′[1], vF′ = vF′[1]

    @unpack_StrParam p
    @unpack_NumParam numericp

    rndCost = RndLabor(xNN′, α, γ)*ω

    temp = (1-τ)*ProfitNorm(sNN, β) - (1-s)*rndCost + xNN′*(vL′ - vNN′) + xNN′*(vF′ - vNN′) + xENN′*(0.5*(0.0 - vNN′) + 0.5*(vF′ - vNN′))

    vNN = (temp - (ρ-1/dt)*vNN′)*dt

    return vNN
end

## VALUE FUNCTION ITERATION
function ValueIteration!(vL::Vector{Float64}, vF::Vector{Float64}, vNN::Float64, xL′::Vector{Float64}, xF′::Vector{Float64}, xNN′::Float64, xE′::Vector{Float64}, xENN′::Float64, ω::Float64, p::StrParam, numericp::NumParam, sL::Vector{Float64}, sF::Vector{Float64}, sNN::Float64, vL′::Vector{Float64}, vF′::Vector{Float64}, vNN′::Float64, A′::Vector{Float64}, temp1::Vector{Float64}, temp2::Vector{Float64})

    @unpack_StrParam p
    @unpack_NumParam numericp

    M = 1:mbar

    @inbounds for i = 1:maxvalueiter
        #println("Iteration i = $i")

        # Find innovation rates and R&D labors wrt Next (′) values
        @inbounds for m = eachindex(M)
            InnL!(xL′, vL′, m, ω, s, α, γ, mbar)
            InnF!(xF′, vF′, vNN′, m, ω, s, α, γ, ϕ)
            InnE!(xE′, vL′, vF′, vNN′, m, ω, αTilde, γTilde, ϕTilde, mbar)
        end
        xNN′ = InnNN(vL′[1], vNN′, ω, s, α, γ)
        xENN′ = InnENN(vL′[1], ω, αTilde, γTilde)

        # Iterate Belmmans one time
        LeaderBellman!(vL, vL′, vNN′, xL′, xF′, xE′, ω, p, numericp, sL, A′)
        FollowerBellman!(vF, vF′, vNN′, xL′, xF′, xE′, ω, p, numericp, sF, A′)
        vNN = NNBellman(vL′[1], vF′[1], vNN′, xNN′, xENN′, ω, p, numericp, sNN)

        # Difference and update
        @. temp1 = abs(vL - vL′)
        @. temp2 = abs(vF - vF′)
        diff = max(maximum(temp1), maximum(temp2), abs(vNN - vNN′))

        if diff < valuetol
            #println("Value iteration ends at i = $(i)")
            break
        end

        @. vL′ = vL
        @. vF′ = vF
        vNN′ = vNN
    end

    # Finally update innovations wrt today's values
    @inbounds for m = eachindex(M)
        InnL!(xL′, vL, m, ω, s, α, γ, mbar)
        InnF!(xF′, vF, vNN, m, ω, s, α, γ, ϕ)
        InnE!(xE′, vL, vF, vNN, m, ω, αTilde, γTilde, ϕTilde, mbar)
    end
    xNN′ = InnNN(vL[1], vNN, ω, s, α, γ)
    xENN′ = InnENN(vL[1], ω, αTilde, γTilde)

    # Objects of main interest are updated values: vL, vF, vNN; and innovations xL′, xF′, xNN′, xE′, xENN′
    return vNN, xNN′, xENN′
end

## DISTRIBUTION ITERATION

# This function calculates μ′=μ(t+dt) given μ=μ(t) s.t. (μ(t+dt)-μ(t))/dt = A′(t). Note μ′NN = 1- sum(μ′). Note that A′ was preallocated
function KolmogorovForward!(μ′::AbstractVector{Float64}, μ::Vector{Float64}, μNN::Float64, xL::Vector{Float64}, xF::Vector{Float64}, xNN::Float64, xE::Vector{Float64}, xENN::Float64, p::StrParam, numericp::NumParam, A′::Vector{Float64})

    @unpack_StrParam p
    @unpack_NumParam numericp

    # m = 1
    A′[1] = μNN*(xNN + xNN + xENN) + μ[2]*((1-ϕ)*xF[2] + (1-ϕTilde)*xE[2]) - μ[1]*(xL[1] + xF[1] + δ + xE[1])
    μ′[1] = μ[1] + A′[1] * dt

    # m = 2,...,mbar-1
    @inbounds for m = 2:mbar-1
        A′[m] = μ[m-1]*xL[m-1] + μ[m+1]*((1-ϕ)*xF[m+1] + (1-ϕTilde)*xE[m+1]) - μ[m]*(xL[m] + xF[m] + δ + xE[m])
        μ′[m] = μ[m] + A′[m] * dt
    end

    # m = mbar
    A′[mbar] = μ[mbar-1]*xL[mbar-1] + μ[mbar]*xL[mbar] - μ[mbar]*(xF[mbar] + δ + xE[mbar])
    μ′[mbar] = μ[mbar] + A′[mbar] * dt
end

function DistributionIteration!(μ′::Vector{Float64}, μ::Vector{Float64}, μNN::Float64, xL::Vector{Float64}, xF::Vector{Float64}, xNN::Float64, xE::Vector{Float64}, xENN::Float64, p::StrParam, numericp::NumParam, A′::Vector{Float64}, temp3::Vector{Float64})

    @unpack valuetol, maxvalueiter = numericp

    @inbounds for i = 1:maxvalueiter
        # println("Iteration i = $i")

        # Iterate KF one time to replace μ′
        KolmogorovForward!(μ′, μ, μNN, xL, xF, xNN, xE, xENN, p, numericp, A′)
        μNN′ = 1-sum(μ′)

        # Difference and update
        @. temp3 = abs(μ′ - μ)
        diff = max(maximum(temp3), abs(μNN′ - μNN))

        if diff < valuetol
            #println("Distribution iteration ends at i = $(i)")
            break
        end

        @. μ = μ′
        μNN = μNN′
    end
end

## BGP SOLVER
# Computes implied wage rate ω′ that clears the labor market
function ImpliedWage!(μ::Vector{Float64}, μnn::Float64, hl::Vector{Float64}, hf::Vector{Float64}, hnn::Float64, he::Vector{Float64}, henn::Float64, sl::Vector{Float64}, sf::Vector{Float64}, snn::Float64, zl::Vector{Float64}, zf::Vector{Float64}, znn::Float64)

    sumNumerator = μnn*(snn/znn + snn/znn)
    sumDenominator = μnn*(hnn + hnn + henn)

    @inbounds for m = eachindex(μ)
        sumNumerator = sumNumerator + μ[m]*(sl[m]/zl[m] + sf[m]/zf[m])
        sumDenominator = sumDenominator + μ[m]*(hl[m] + hf[m] + he[m])
    end

    ω′ = sumNumerator/(1-sumDenominator)

    return ω′
end

# Solves the bgp of the model
function BgpSolver(p::StrParam, numericp::NumParam)

    @unpack_StrParam p
    @unpack_NumParam numericp

    # Shares and markups
    sl, sf, snn = Shares(mbar, λ, β, ψ)
    zl, zf, znn = Markups(sl, sf, snn, β)

    # Initial guess for values and distribution: vl′, vf′, vnn′, and μ, μnn
    vl′ = zeros(mbar)
    vf′ = zeros(mbar)
    vnn′ = 0.0

    μ = fill(1.0/(mbar+1), mbar)
    μnn = 1.0/(mbar+1)

    # Pre-allocation for vl, vf, vnn; xl, xf, xnn, xe, xenn; hl, hf, hnn, he, henn; μ′, μnn′
    vl = similar(vl′)
    vf = similar(vf′)
    vnn = 0.0

    xl = similar(vl′)
    xf = similar(vf′)
    xnn = 0.0
    xe = similar(vf′)
    xenn = 0.0

    hl = similar(xl)
    hf = similar(xf)
    hnn = 0.0
    he = similar(xe)
    henn = 0.0

    μ′ = similar(μ)
    μnn′ = 0.0

    # Computational pre-allocation: A′, temp1, temp2, temp3
    A′ = similar(vl′)
    temp1 = similar(vl′)
    temp2 = similar(vl′)
    temp3 = similar(vl′)

    # Wage guess
    ω = 0.9
    local ω′

    # Wage iteration
    @inbounds for i = 1:maxwageiter
        # println("Iteration = $i")

        # Given wage, solve for values vl, vf, vnn; and innovations xl, xf, xnn, xe, xenn
        vnn, xnn, xenn = ValueIteration!(vl, vf, vnn, xl, xf, xnn, xe, xenn, ω, p, numericp, sl, sf, snn, vl′, vf′, vnn′, A′, temp1, temp2)

        # Given policies, solve for distribution μ′, μnn′
        DistributionIteration!(μ′, μ, μnn, xl, xf, xnn, xe, xenn, p, numericp, A′, temp3)
        μnn′ = 1.0 - sum(μ′)

        # Calculate R&D labors
        @. hl = RndLabor(xl, α, γ)
        @. hf = RndLabor(xf, α, γ)
        hnn = RndLabor(xnn, α, γ)
        @. he = RndLabor(xe, αTilde, γTilde)
        henn = RndLabor(xenn, αTilde, γTilde)

        # Given distribution, clear labor market and find implied wage
        ω′ = ImpliedWage!(μ′, μnn′, hl, hf, hnn, he, henn, sl, sf, snn, zl, zf, znn)

        # Wage difference and update
        diff = abs(ω′ - ω)

        if diff < wagetol
            # println("Wage iteration ends at i = $(i)")
            break
        end

        # Update wage
        ω = wageweight*ω + (1-wageweight)*ω′
    end

    return BGP(p = p, numericp = numericp, ω = ω′, μ = μ′, μnn = μnn′, xl = xl, xf = xf, xnn = xnn, xe = xe, xenn = xenn, vl = vl, vf = vf, vnn = vnn, sl = sl, sf = sf, snn = snn, zl = zl, zf = zf, znn = znn)
end
