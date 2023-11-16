## FUNCTIONS RELATED TO CURVATURE IN STEP SIZES
@inline Step(m::Int, ψ::Float64) = m ^ ψ
@inline increase_in_productivity(m::Int, λ::Float64, ψ::Float64) = λ^((m+1)^ψ - m^ψ)

# If a firm with productivity q innovates, new productivity q′ is q′ = increase_in_productivity(m, λ, ψ) * q
# Above structure always implies the following: In a sector with gap m, we have q_leader = q_follower * λ^(m^ψ)

## REVENUE SHARES
function LeaderShareFunction!(F::Vector{Float64}, sl::Vector{Float64}, λ::Float64, β::Float64, ψ::Float64)

    for m = eachindex(F)
        F[m] = sl[m]*(λ^(-Step(m, ψ)) * sl[m]/(1-sl[m]) * (1-β*sl[m])/(1-β*(1-sl[m])))^(β/(1-β)) + sl[m] - 1.0
    end
end

function Shares(mbar::Int, λ::Float64, β::Float64, ψ::Float64)

    slinit = fill(0.80, mbar)
    sl = nlsolve((F, sl) -> LeaderShareFunction!(F, sl, λ, β, ψ), slinit).zero
    sf = 1.0 .- sl
    snn = 0.5

    return sl, sf, snn
end

## MARKUPS
@inline Markup(s, β) = (1 - β*s) / (β*(1 - s))

function Markups(sl::Vector{Float64}, sf::Vector{Float64}, snn::Float64, β::Float64)

    zl = Markup.(sl, β)
    zf = Markup.(sf, β)
    znn = Markup(snn, β)

    return zl, zf, znn
end

## TIME VARYING SHARES AND MARKUPS
function TimeVaryingSharesMarkups(sequencep::Vector{StrParam}, numericp::NumParam)
    @unpack_NumParam numericp

    sl = zeros(mbar, tperiods)

    sl[:,1] .= nlsolve((F, sl) -> LeaderShareFunction!(F, sl, sequencep[1].λ, sequencep[1].β, sequencep[1].ψ), fill(0.8, mbar)).zero

    @inbounds for t = 2:tperiods
        sl[:,t] .= nlsolve((F, sl) -> LeaderShareFunction!(F, sl, sequencep[t].λ, sequencep[t].β, sequencep[t].ψ), sl[:,t-1]).zero
    end

    sf = 1 .- sl
    snn = fill(0.5, tperiods)

    zl = Matrix{Float64}(undef, mbar, tperiods)
    zf = similar(zl)
    znn = Vector{Float64}(undef, tperiods)

    @inbounds for t = 1:tperiods
        zl[:,t], zf[:,t], znn[t] = Markups(sl[:,t], sf[:,t], snn[t], sequencep[t].β)
    end

    return sl, sf, snn, zl, zf, znn
end

## OUTPUTS
@inline FirmOutput(s, z, q, ω) = s / z * q / ω
@inline FirmOutput(q, l) = q * l
@inline SectorOutput(yl, yf, β) = (yl ^ β + yf ^ β) ^ (1 / β)
@inline AggregateOutput(ysector::Vector{Float64}, ysectornn::Float64, μ::Vector{Float64}, μnn::Float64) = exp(sum(log.(ysector) .* μ) + log(ysectornn)*μnn)

## LABORS
@inline ProdLabor(s, z, ω) = s / z / ω
@inline RndLabor(x, α, γ) = α * x ^ γ / γ

## PROFITS
@inline Profit(s, Y, β) = ((1 - β) * s) / (1 - β * s) * Y
@inline ProfitNorm(s, β) = ((1 - β) * s) / (1 - β * s)
