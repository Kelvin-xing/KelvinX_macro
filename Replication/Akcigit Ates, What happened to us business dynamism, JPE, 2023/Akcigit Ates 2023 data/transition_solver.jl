## Transition struct that contains endogenous variables over time. Columns represent each time period
@with_kw struct Transition
    initialbgp::BGP
    sequencep::Vector{StrParam}
    numericp::NumParam
    
    ω::Vector{Float64}

    μ::Matrix{Float64}
    μnn::Vector{Float64}

    xl::Matrix{Float64}
    xf::Matrix{Float64}
    xnn::Vector{Float64}
    xe::Matrix{Float64}
    xenn::Vector{Float64}

    vl::Matrix{Float64}
    vf::Matrix{Float64}
    vnn::Vector{Float64}

    sl::Matrix{Float64}
    sf::Matrix{Float64}
    snn::Vector{Float64}
    zl::Matrix{Float64}
    zf::Matrix{Float64}
    znn::Vector{Float64}
end

# To solve for the transition of the model starting from a BGP, we need a parameter sequence. The function below creates the parameter sequence given transitionp (the set of parameters that govern the change in parameters over time)
function CreateParameters(initialp::StrParam, transitionp::TransParam, numericp::NumParam)

    @unpack tperiodpol, tperiods = numericp
    @unpack_StrParam initialp

    # Preallocate
    sequencep = Vector{StrParam}(undef, tperiods)

    # Policy values of variables
    τpol = τ * transitionp.τ
    spol = s * transitionp.s
    αTildepol = αTilde * transitionp.αTilde
    δpol = δ * transitionp.δ

    local τtemp, stemp

    # Create parameters until tperiodpol
    @inbounds for t = 1:tperiodpol

        αTildetemp = αTilde + (αTildepol - αTilde) * ( 1.0 - ( ( exp( (1.0 - t / tperiodpol) * transitionp.ναTilde ) - 1.0 ) / ( exp( transitionp.ναTilde ) - 1.0 ) ) )
        
        δtemp = δ + (δpol - δ) * ( 1.0 - ( ( exp( (1.0 - t / tperiodpol) * transitionp.νδ ) - 1.0 ) / ( exp( transitionp.νδ ) - 1.0 ) ) )

        τtemp = τ + (τpol - τ) * ( 1.0 - ( ( exp( (1.0 - t / tperiodpol) * transitionp.ντ ) - 1.0 ) / ( exp( transitionp.ντ ) - 1.0 ) ) )

        stemp = s + (spol - s) * ( 1.0 - ( ( exp( (1.0 - t / tperiodpol) * transitionp.νs ) - 1.0 ) / ( exp( transitionp.νs ) - 1.0 ) ) )

        sequencep[t] = StrParam(initialp, τ = τtemp, s = stemp, αTilde = αTildetemp, δ = δtemp)
    end

    # Create parameters after tperiodpol
    @inbounds for t = tperiodpol+1:tperiods
        sequencep[t] = sequencep[t-1]
    end

    return sequencep
end

# Solver function for transition of the model (Shares and markups are constant during transition. This is because parameters that determine revenue shares such as λ, β are not time varying during transition in the benchmark case)
function TransitionSolver(initial::BGP, transitionp::TransParam, numericp::NumParam)

    @unpack_NumParam numericp
    @unpack sl, sf, snn, zl, zf, znn = initial

    ### Create sequence of parameters
    sequencep = CreateParameters(initial.p, transitionp, numericp)

    ### Solve for terminal bgp
    terminal = BgpSolver(sequencep[end], numericp)

    ### Pre-allocate
    vl = Matrix{Float64}(undef, mbar, tperiods)
    vf = similar(vl)
    vnn = Vector{Float64}(undef, tperiods)
    xl = similar(vl)
    xf = similar(vl)
    xnn = similar(vnn)
    xe = similar(vl)
    xenn = similar(vnn)
    μ = similar(vl)
    μnn = similar(vnn)
    A′ = similar(initial.vl)
    hl = similar(initial.vl)
    hf = similar(initial.vl)
    hnn = 0.0
    he = similar(initial.vl)
    henn = 0.0

    ### Guess wage
    ω = fill(initial.ω, tperiods)
    ω′ = similar(ω)
    temp = similar(ω)

    ### Wage iteration
    @inbounds for i = 1:maxwageiter

        # println("Wage iteration = $i")

        ## VALUE ITERATION

        # For t = tperiods, we need to calculate values from terminal bgp
        vL = @view vl[:,tperiods]
        vF = @view vf[:,tperiods]

        LeaderBellman!(vL, terminal.vl, terminal.vnn, terminal.xl, terminal.xf, terminal.xe, terminal.ω, terminal.p, numericp, sl, A′)
        FollowerBellman!(vF, terminal.vf, terminal.vnn, terminal.xl, terminal.xf, terminal.xe, terminal.ω, terminal.p, numericp, sf, A′)
        vnn[tperiods] = NNBellman(terminal.vl[1], terminal.vf[1], terminal.vnn, terminal.xnn, terminal.xenn, terminal.ω, terminal.p, numericp, snn)

        # For t < tperiods
        @inbounds for t = tperiods-1:-1:1

            # Parameters
            p = sequencep[t+1]
            @unpack_StrParam p

            # Views
            vL = @view vl[:,t]
            vF = @view vf[:,t]
            xL = @view xl[:,t+1]
            xF = @view xf[:,t+1]
            xE = @view xe[:,t+1]

            # Innovations in t+1
            @inbounds for m = 1:mbar
                InnL!(xL, vl[:,t+1], m, ω[t+1], s, α, γ, mbar)
                InnF!(xF, vf[:,t+1], vnn[t+1], m, ω[t+1], s, α, γ, ϕ)
                InnE!(xE, vl[:,t+1], vf[:,t+1], vnn[t+1], m, ω[t+1], αTilde, γTilde, ϕTilde, mbar)
            end
            xnn[t+1] = InnNN(vl[1,t+1], vnn[t+1], ω[t+1], s, α, γ)
            xenn[t+1] = InnENN(vl[1,t+1], ω[t+1], αTilde, γTilde)

            # Iterate Bellmans to find t values
            LeaderBellman!(vL, vl[:,t+1], vnn[t+1], xl[:,t+1], xf[:,t+1], xe[:,t+1], ω[t+1], p, numericp, sl, A′)
            FollowerBellman!(vF, vf[:,t+1], vnn[t+1], xl[:,t+1], xf[:,t+1], xe[:,t+1], ω[t+1], p, numericp, sf, A′)
            vnn[t] = NNBellman(vl[1,t+1], vf[1,t+1], vnn[t+1], xnn[t+1], xenn[t+1], ω[t+1], p, numericp, snn)
        end

        # For t = 1, we need to calculate innovations
        p = sequencep[1]
        @unpack_StrParam p

        xL = @view xl[:,1]
        xF = @view xf[:,1]
        xE = @view xe[:,1]

        @inbounds for m = 1:mbar
            InnL!(xL, vl[:,1], m, ω[1], s, α, γ, mbar)
            InnF!(xF, vf[:,1], vnn[1], m, ω[1], s, α, γ, ϕ)
            InnE!(xE, vl[:,1], vf[:,1], vnn[1], m, ω[1], αTilde, γTilde, ϕTilde, mbar)
        end
        xnn[1] = InnNN(vl[1,1], vnn[1], ω[1], s, α, γ)
        xenn[1] = InnENN(vl[1,1], ω[1], αTilde, γTilde)


        ## DISTRIBUTION ITERATION

        # For t = 1, we need to calculate distribution from initial bgp
        μview = @view μ[:,1]

        KolmogorovForward!(μview, initial.μ, initial.μnn, initial.xl, initial.xf, initial.xnn, initial.xe, initial.xenn, initial.p, numericp, A′)
        μnn[1] = 1.0 - sum(μ[:,1])

        # For t > 2
        @inbounds for t = 1:tperiods-1

            # Views
            μview = @view μ[:,t+1]

            # Iterate Kolmogorov Forward to find distribution in t+1
            KolmogorovForward!(μview, μ[:,t], μnn[t], xl[:,t], xf[:,t], xnn[t], xe[:,t], xenn[t], sequencep[t], numericp, A′)
            μnn[t+1] = 1.0 - sum(μ[:,t+1])
        end

        ## IMPLIED WAGES
        @inbounds for t = 1:tperiods

            @unpack α, γ, αTilde, γTilde = sequencep[t]

            # Calculate R&D labors
            @. hl = RndLabor(xl[:,t], α, γ)
            @. hf = RndLabor(xf[:,t], α, γ)
            hnn = RndLabor(xnn[t], α, γ)
            @. he = RndLabor(xe[:,t], αTilde, γTilde)
            henn = RndLabor(xenn[t], αTilde, γTilde)

            ω′[t] = ImpliedWage!(μ[:,t], μnn[t], hl, hf, hnn, he, henn, sl, sf, snn, zl, zf, znn)
        end

        ## DIFFERENCE AND UPDATE
        @. temp = abs(ω′ - ω)
        diff = maximum(temp)

        if diff < wagetol
            # println("Wage iteration ends at i = $(i)")
            break
        end

        @. ω = wageweight * ω + (1.0 - wageweight) * ω′
    end

    transition = Transition(initialbgp = initial, sequencep = sequencep, numericp = numericp, ω = ω′, μ = μ, μnn = μnn, xl = xl, xf = xf, xnn = xnn, xe = xe, xenn = xenn, vl = vl, vf = vf, vnn = vnn, sl = repeat(sl, 1, tperiods), sf = repeat(sf, 1, tperiods), snn = fill(snn, tperiods), zl = repeat(zl, 1, tperiods), zf = repeat(zf, 1, tperiods), znn = fill(znn, tperiods))

    return transition, terminal
end

## Other solver functions for the transition of the model

# 1. Section 8 - Alternative mechanisms: Declining interest rate
function TransitionSolver_DecliningInterestRate(initial::BGP, ρpol::Float64)
    
    @unpack_NumParam initial.numericp
    @unpack sl, sf, snn, zl, zf, znn = initial
    
    ### Create sequence of parameters for declining ρ
    sequencep = Vector{StrParam}(undef, tperiods)
    νρ = 1.0

    # Create parameters until tperiodpol
    @inbounds for t = 1:tperiodpol

        ρtemp = initial.p.ρ + (ρpol - initial.p.ρ) * ( 1.0 - ( ( exp( (1.0 - t / tperiodpol) * νρ ) - 1.0 ) / ( exp( νρ ) - 1.0 ) ) )
        
        sequencep[t] = StrParam(initial.p, ρ = ρtemp)
    end

    # Create parameters after tperiodpol
    @inbounds for t = tperiodpol+1:tperiods
        sequencep[t] = sequencep[t-1]
    end

    ### Solve for terminal bgp
    terminal = BgpSolver(sequencep[end], numericp)

    ### Pre-allocate
    vl = Matrix{Float64}(undef, mbar, tperiods)
    vf = similar(vl)
    vnn = Vector{Float64}(undef, tperiods)
    xl = similar(vl)
    xf = similar(vl)
    xnn = similar(vnn)
    xe = similar(vl)
    xenn = similar(vnn)
    μ = similar(vl)
    μnn = similar(vnn)
    A′ = similar(initial.vl)
    hl = similar(initial.vl)
    hf = similar(initial.vl)
    hnn = 0.0
    he = similar(initial.vl)
    henn = 0.0

    ### Guess wage
    ω = fill(initial.ω, tperiods)
    ω′ = similar(ω)
    temp = similar(ω)

    ### Wage iteration
    @inbounds for i = 1:maxwageiter

        # println("Wage iteration = $i")

        ## VALUE ITERATION

        # For t = tperiods, we need to calculate values from terminal bgp
        vL = @view vl[:,tperiods]
        vF = @view vf[:,tperiods]

        LeaderBellman!(vL, terminal.vl, terminal.vnn, terminal.xl, terminal.xf, terminal.xe, terminal.ω, terminal.p, numericp, sl, A′)
        FollowerBellman!(vF, terminal.vf, terminal.vnn, terminal.xl, terminal.xf, terminal.xe, terminal.ω, terminal.p, numericp, sf, A′)
        vnn[tperiods] = NNBellman(terminal.vl[1], terminal.vf[1], terminal.vnn, terminal.xnn, terminal.xenn, terminal.ω, terminal.p, numericp, snn)

        # For t < tperiods
        @inbounds for t = tperiods-1:-1:1

            # Parameters
            p = sequencep[t+1]
            @unpack_StrParam p

            # Views
            vL = @view vl[:,t]
            vF = @view vf[:,t]
            xL = @view xl[:,t+1]
            xF = @view xf[:,t+1]
            xE = @view xe[:,t+1]

            # Innovations in t+1
            @inbounds for m = 1:mbar
                InnL!(xL, vl[:,t+1], m, ω[t+1], s, α, γ, mbar)
                InnF!(xF, vf[:,t+1], vnn[t+1], m, ω[t+1], s, α, γ, ϕ)
                InnE!(xE, vl[:,t+1], vf[:,t+1], vnn[t+1], m, ω[t+1], αTilde, γTilde, ϕTilde, mbar)
            end
            xnn[t+1] = InnNN(vl[1,t+1], vnn[t+1], ω[t+1], s, α, γ)
            xenn[t+1] = InnENN(vl[1,t+1], ω[t+1], αTilde, γTilde)

            # Iterate Bellmans to find t values
            LeaderBellman!(vL, vl[:,t+1], vnn[t+1], xl[:,t+1], xf[:,t+1], xe[:,t+1], ω[t+1], p, numericp, sl, A′)
            FollowerBellman!(vF, vf[:,t+1], vnn[t+1], xl[:,t+1], xf[:,t+1], xe[:,t+1], ω[t+1], p, numericp, sf, A′)
            vnn[t] = NNBellman(vl[1,t+1], vf[1,t+1], vnn[t+1], xnn[t+1], xenn[t+1], ω[t+1], p, numericp, snn)
        end

        # For t = 1, we need to calculate innovations
        p = sequencep[1]
        @unpack_StrParam p

        xL = @view xl[:,1]
        xF = @view xf[:,1]
        xE = @view xe[:,1]

        @inbounds for m = 1:mbar
            InnL!(xL, vl[:,1], m, ω[1], s, α, γ, mbar)
            InnF!(xF, vf[:,1], vnn[1], m, ω[1], s, α, γ, ϕ)
            InnE!(xE, vl[:,1], vf[:,1], vnn[1], m, ω[1], αTilde, γTilde, ϕTilde, mbar)
        end
        xnn[1] = InnNN(vl[1,1], vnn[1], ω[1], s, α, γ)
        xenn[1] = InnENN(vl[1,1], ω[1], αTilde, γTilde)


        ## DISTRIBUTION ITERATION

        # For t = 1, we need to calculate distribution from initial bgp
        μview = @view μ[:,1]

        KolmogorovForward!(μview, initial.μ, initial.μnn, initial.xl, initial.xf, initial.xnn, initial.xe, initial.xenn, initial.p, numericp, A′)
        μnn[1] = 1.0 - sum(μ[:,1])

        # For t > 2
        @inbounds for t = 1:tperiods-1

            # Views
            μview = @view μ[:,t+1]

            # Iterate Kolmogorov Forward to find distribution in t+1
            KolmogorovForward!(μview, μ[:,t], μnn[t], xl[:,t], xf[:,t], xnn[t], xe[:,t], xenn[t], sequencep[t], numericp, A′)
            μnn[t+1] = 1.0 - sum(μ[:,t+1])
        end

        ## IMPLIED WAGES
        @inbounds for t = 1:tperiods

            @unpack α, γ, αTilde, γTilde = sequencep[t]

            # Calculate R&D labors
            @. hl = RndLabor(xl[:,t], α, γ)
            @. hf = RndLabor(xf[:,t], α, γ)
            hnn = RndLabor(xnn[t], α, γ)
            @. he = RndLabor(xe[:,t], αTilde, γTilde)
            henn = RndLabor(xenn[t], αTilde, γTilde)

            ω′[t] = ImpliedWage!(μ[:,t], μnn[t], hl, hf, hnn, he, henn, sl, sf, snn, zl, zf, znn)
        end

        ## DIFFERENCE AND UPDATE
        @. temp = abs(ω′ - ω)
        diff = maximum(temp)

        if diff < wagetol
            # println("Wage iteration ends at i = $(i)")
            break
        end

        @. ω = wageweight * ω + (1.0 - wageweight) * ω′
    end

    transition = Transition(initialbgp = initial, sequencep = sequencep, numericp = numericp, ω = ω′, μ = μ, μnn = μnn, xl = xl, xf = xf, xnn = xnn, xe = xe, xenn = xenn, vl = vl, vf = vf, vnn = vnn, sl = repeat(sl, 1, tperiods), sf = repeat(sf, 1, tperiods), snn = fill(snn, tperiods), zl = repeat(zl, 1, tperiods), zf = repeat(zf, 1, tperiods), znn = fill(znn, tperiods))

    return transition, terminal
end

# 2. Section 8 - Alternative mechanisms: Ideas getting harder to find
function TransitionSolver_IdeasGettingHarder(initial::BGP, αpol::Float64, αTildepol::Float64)
    
    @unpack_NumParam initial.numericp
    @unpack sl, sf, snn, zl, zf, znn = initial
    
    ### Create sequence of parameters for declining ρ
    sequencep = Vector{StrParam}(undef, tperiods)
    να = 1.0
    ναTilde = 1.0

    # Create parameters until tperiodpol
    @inbounds for t = 1:tperiodpol

        αtemp = initial.p.α + (αpol - initial.p.α) * ( 1.0 - ( ( exp( (1.0 - t / tperiodpol) * να ) - 1.0 ) / ( exp( να ) - 1.0 ) ) )

        αTildetemp = initial.p.αTilde + (αTildepol - initial.p.αTilde) * ( 1.0 - ( ( exp( (1.0 - t / tperiodpol) * ναTilde ) - 1.0 ) / ( exp( ναTilde ) - 1.0 ) ) )
        
        sequencep[t] = StrParam(initial.p, α = αtemp, αTilde = αTildetemp)
    end

    # Create parameters after tperiodpol
    @inbounds for t = tperiodpol+1:tperiods
        sequencep[t] = sequencep[t-1]
    end

    ### Solve for terminal bgp
    terminal = BgpSolver(sequencep[end], numericp)

    ### Pre-allocate
    vl = Matrix{Float64}(undef, mbar, tperiods)
    vf = similar(vl)
    vnn = Vector{Float64}(undef, tperiods)
    xl = similar(vl)
    xf = similar(vl)
    xnn = similar(vnn)
    xe = similar(vl)
    xenn = similar(vnn)
    μ = similar(vl)
    μnn = similar(vnn)
    A′ = similar(initial.vl)
    hl = similar(initial.vl)
    hf = similar(initial.vl)
    hnn = 0.0
    he = similar(initial.vl)
    henn = 0.0

    ### Guess wage
    ω = fill(initial.ω, tperiods)
    ω′ = similar(ω)
    temp = similar(ω)

    ### Wage iteration
    @inbounds for i = 1:maxwageiter

        # println("Wage iteration = $i")

        ## VALUE ITERATION

        # For t = tperiods, we need to calculate values from terminal bgp
        vL = @view vl[:,tperiods]
        vF = @view vf[:,tperiods]

        LeaderBellman!(vL, terminal.vl, terminal.vnn, terminal.xl, terminal.xf, terminal.xe, terminal.ω, terminal.p, numericp, sl, A′)
        FollowerBellman!(vF, terminal.vf, terminal.vnn, terminal.xl, terminal.xf, terminal.xe, terminal.ω, terminal.p, numericp, sf, A′)
        vnn[tperiods] = NNBellman(terminal.vl[1], terminal.vf[1], terminal.vnn, terminal.xnn, terminal.xenn, terminal.ω, terminal.p, numericp, snn)

        # For t < tperiods
        @inbounds for t = tperiods-1:-1:1

            # Parameters
            p = sequencep[t+1]
            @unpack_StrParam p

            # Views
            vL = @view vl[:,t]
            vF = @view vf[:,t]
            xL = @view xl[:,t+1]
            xF = @view xf[:,t+1]
            xE = @view xe[:,t+1]

            # Innovations in t+1
            @inbounds for m = 1:mbar
                InnL!(xL, vl[:,t+1], m, ω[t+1], s, α, γ, mbar)
                InnF!(xF, vf[:,t+1], vnn[t+1], m, ω[t+1], s, α, γ, ϕ)
                InnE!(xE, vl[:,t+1], vf[:,t+1], vnn[t+1], m, ω[t+1], αTilde, γTilde, ϕTilde, mbar)
            end
            xnn[t+1] = InnNN(vl[1,t+1], vnn[t+1], ω[t+1], s, α, γ)
            xenn[t+1] = InnENN(vl[1,t+1], ω[t+1], αTilde, γTilde)

            # Iterate Bellmans to find t values
            LeaderBellman!(vL, vl[:,t+1], vnn[t+1], xl[:,t+1], xf[:,t+1], xe[:,t+1], ω[t+1], p, numericp, sl, A′)
            FollowerBellman!(vF, vf[:,t+1], vnn[t+1], xl[:,t+1], xf[:,t+1], xe[:,t+1], ω[t+1], p, numericp, sf, A′)
            vnn[t] = NNBellman(vl[1,t+1], vf[1,t+1], vnn[t+1], xnn[t+1], xenn[t+1], ω[t+1], p, numericp, snn)
        end

        # For t = 1, we need to calculate innovations
        p = sequencep[1]
        @unpack_StrParam p

        xL = @view xl[:,1]
        xF = @view xf[:,1]
        xE = @view xe[:,1]

        @inbounds for m = 1:mbar
            InnL!(xL, vl[:,1], m, ω[1], s, α, γ, mbar)
            InnF!(xF, vf[:,1], vnn[1], m, ω[1], s, α, γ, ϕ)
            InnE!(xE, vl[:,1], vf[:,1], vnn[1], m, ω[1], αTilde, γTilde, ϕTilde, mbar)
        end
        xnn[1] = InnNN(vl[1,1], vnn[1], ω[1], s, α, γ)
        xenn[1] = InnENN(vl[1,1], ω[1], αTilde, γTilde)


        ## DISTRIBUTION ITERATION

        # For t = 1, we need to calculate distribution from initial bgp
        μview = @view μ[:,1]

        KolmogorovForward!(μview, initial.μ, initial.μnn, initial.xl, initial.xf, initial.xnn, initial.xe, initial.xenn, initial.p, numericp, A′)
        μnn[1] = 1.0 - sum(μ[:,1])

        # For t > 2
        @inbounds for t = 1:tperiods-1

            # Views
            μview = @view μ[:,t+1]

            # Iterate Kolmogorov Forward to find distribution in t+1
            KolmogorovForward!(μview, μ[:,t], μnn[t], xl[:,t], xf[:,t], xnn[t], xe[:,t], xenn[t], sequencep[t], numericp, A′)
            μnn[t+1] = 1.0 - sum(μ[:,t+1])
        end

        ## IMPLIED WAGES
        @inbounds for t = 1:tperiods

            @unpack α, γ, αTilde, γTilde = sequencep[t]

            # Calculate R&D labors
            @. hl = RndLabor(xl[:,t], α, γ)
            @. hf = RndLabor(xf[:,t], α, γ)
            hnn = RndLabor(xnn[t], α, γ)
            @. he = RndLabor(xe[:,t], αTilde, γTilde)
            henn = RndLabor(xenn[t], αTilde, γTilde)

            ω′[t] = ImpliedWage!(μ[:,t], μnn[t], hl, hf, hnn, he, henn, sl, sf, snn, zl, zf, znn)
        end

        ## DIFFERENCE AND UPDATE
        @. temp = abs(ω′ - ω)
        diff = maximum(temp)

        if diff < wagetol
            # println("Wage iteration ends at i = $(i)")
            break
        end

        @. ω = wageweight * ω + (1.0 - wageweight) * ω′
    end

    transition = Transition(initialbgp = initial, sequencep = sequencep, numericp = numericp, ω = ω′, μ = μ, μnn = μnn, xl = xl, xf = xf, xnn = xnn, xe = xe, xenn = xenn, vl = vl, vf = vf, vnn = vnn, sl = repeat(sl, 1, tperiods), sf = repeat(sf, 1, tperiods), snn = fill(snn, tperiods), zl = repeat(zl, 1, tperiods), zf = repeat(zf, 1, tperiods), znn = fill(znn, tperiods))

    return transition, terminal
end

# 3. Section 8 - Alternative mechanisms: Weaker market power of labor - The important distinction of the following transition solver from the others is the fact that when λ changes during the transition, revenue shares and markups are no longer constant over time.
function TransitionSolver_WeakerLaborPower(initial::BGP, λpol::Float64)
    
    @unpack_NumParam initial.numericp
    
    ### Create sequence of parameters for declining ρ
    sequencep = Vector{StrParam}(undef, tperiods)
    stepsize = initial.p.λ - 1
    stepsizepol = λpol - 1
    νstepsize = 1.0

    # Create parameters until tperiodpol
    @inbounds for t = 1:tperiodpol

        stepsizetemp = stepsize  + (stepsizepol - stepsize) * ( 1.0 - ( ( exp( (1.0 - t / tperiodpol) * νstepsize ) - 1.0 ) / ( exp( νstepsize ) - 1.0 ) ) )
        
        sequencep[t] = StrParam(initial.p, λ = 1 + stepsizetemp)
    end

    # Create parameters after tperiodpol
    @inbounds for t = tperiodpol+1:tperiods
        sequencep[t] = sequencep[t-1]
    end

    ### Solve for terminal bgp
    terminal = BgpSolver(sequencep[end], numericp)

    ### First calculate time varying revenue shares and markups as matrices sl, sf, snn, zl, zf, znn. As before, each row corresponds to productivity gap m, while each column corresponds a time period during transition.
    sl, sf, snn, zl, zf, znn = TimeVaryingSharesMarkups(sequencep, initial.numericp)

    ### Pre-allocate
    vl = Matrix{Float64}(undef, mbar, tperiods)
    vf = similar(vl)
    vnn = Vector{Float64}(undef, tperiods)
    xl = similar(vl)
    xf = similar(vl)
    xnn = similar(vnn)
    xe = similar(vl)
    xenn = similar(vnn)
    μ = similar(vl)
    μnn = similar(vnn)
    A′ = similar(initial.vl)
    hl = similar(initial.vl)
    hf = similar(initial.vl)
    hnn = 0.0
    he = similar(initial.vl)
    henn = 0.0

    ### Guess wage
    ω = fill(initial.ω, tperiods)
    ω′ = similar(ω)
    temp = similar(ω)

    ### Wage iteration
    @inbounds for i = 1:maxwageiter

        # println("Wage iteration = $i")

        ## VALUE ITERATION

        # For t = tperiods, we need to calculate values from terminal bgp. Note that we make use of shares and markups at terminal BGP in this step
        vL = @view vl[:,tperiods]
        vF = @view vf[:,tperiods]

        LeaderBellman!(vL, terminal.vl, terminal.vnn, terminal.xl, terminal.xf, terminal.xe, terminal.ω, terminal.p, numericp, terminal.sl, A′)
        FollowerBellman!(vF, terminal.vf, terminal.vnn, terminal.xl, terminal.xf, terminal.xe, terminal.ω, terminal.p, numericp, terminal.sf, A′)
        vnn[tperiods] = NNBellman(terminal.vl[1], terminal.vf[1], terminal.vnn, terminal.xnn, terminal.xenn, terminal.ω, terminal.p, numericp, terminal.snn)

        # For t < tperiods
        @inbounds for t = tperiods-1:-1:1

            # Parameters
            p = sequencep[t+1]
            @unpack_StrParam p

            # Views
            vL = @view vl[:,t]
            vF = @view vf[:,t]
            xL = @view xl[:,t+1]
            xF = @view xf[:,t+1]
            xE = @view xe[:,t+1]

            # Innovations in t+1
            @inbounds for m = 1:mbar
                InnL!(xL, vl[:,t+1], m, ω[t+1], s, α, γ, mbar)
                InnF!(xF, vf[:,t+1], vnn[t+1], m, ω[t+1], s, α, γ, ϕ)
                InnE!(xE, vl[:,t+1], vf[:,t+1], vnn[t+1], m, ω[t+1], αTilde, γTilde, ϕTilde, mbar)
            end
            xnn[t+1] = InnNN(vl[1,t+1], vnn[t+1], ω[t+1], s, α, γ)
            xenn[t+1] = InnENN(vl[1,t+1], ω[t+1], αTilde, γTilde)

            # Iterate Bellmans to find t values
            LeaderBellman!(vL, vl[:,t+1], vnn[t+1], xl[:,t+1], xf[:,t+1], xe[:,t+1], ω[t+1], p, numericp, sl[:, t+1], A′)
            FollowerBellman!(vF, vf[:,t+1], vnn[t+1], xl[:,t+1], xf[:,t+1], xe[:,t+1], ω[t+1], p, numericp, sf[:, t+1], A′)
            vnn[t] = NNBellman(vl[1,t+1], vf[1,t+1], vnn[t+1], xnn[t+1], xenn[t+1], ω[t+1], p, numericp, snn[t+1])
        end

        # For t = 1, we need to calculate innovations
        p = sequencep[1]
        @unpack_StrParam p

        xL = @view xl[:,1]
        xF = @view xf[:,1]
        xE = @view xe[:,1]

        @inbounds for m = 1:mbar
            InnL!(xL, vl[:,1], m, ω[1], s, α, γ, mbar)
            InnF!(xF, vf[:,1], vnn[1], m, ω[1], s, α, γ, ϕ)
            InnE!(xE, vl[:,1], vf[:,1], vnn[1], m, ω[1], αTilde, γTilde, ϕTilde, mbar)
        end
        xnn[1] = InnNN(vl[1,1], vnn[1], ω[1], s, α, γ)
        xenn[1] = InnENN(vl[1,1], ω[1], αTilde, γTilde)


        ## DISTRIBUTION ITERATION

        # For t = 1, we need to calculate distribution from initial bgp
        μview = @view μ[:,1]

        KolmogorovForward!(μview, initial.μ, initial.μnn, initial.xl, initial.xf, initial.xnn, initial.xe, initial.xenn, initial.p, numericp, A′)
        μnn[1] = 1.0 - sum(μ[:,1])

        # For t > 2
        @inbounds for t = 1:tperiods-1

            # Views
            μview = @view μ[:,t+1]

            # Iterate Kolmogorov Forward to find distribution in t+1
            KolmogorovForward!(μview, μ[:,t], μnn[t], xl[:,t], xf[:,t], xnn[t], xe[:,t], xenn[t], sequencep[t], numericp, A′)
            μnn[t+1] = 1.0 - sum(μ[:,t+1])
        end

        ## IMPLIED WAGES
        @inbounds for t = 1:tperiods

            @unpack α, γ, αTilde, γTilde = sequencep[t]

            # Calculate R&D labors
            @. hl = RndLabor(xl[:,t], α, γ)
            @. hf = RndLabor(xf[:,t], α, γ)
            hnn = RndLabor(xnn[t], α, γ)
            @. he = RndLabor(xe[:,t], αTilde, γTilde)
            henn = RndLabor(xenn[t], αTilde, γTilde)

            ω′[t] = ImpliedWage!(μ[:,t], μnn[t], hl, hf, hnn, he, henn, sl[:, t], sf[:, t], snn[t], zl[:, t], zf[:, t], znn[t])
        end

        ## DIFFERENCE AND UPDATE
        @. temp = abs(ω′ - ω)
        diff = maximum(temp)

        if diff < wagetol
            # println("Wage iteration ends at i = $(i)")
            break
        end

        @. ω = wageweight * ω + (1.0 - wageweight) * ω′
    end

    transition = Transition(initialbgp = initial, sequencep = sequencep, numericp = numericp, ω = ω′, μ = μ, μnn = μnn, xl = xl, xf = xf, xnn = xnn, xe = xe, xenn = xenn, vl = vl, vf = vf, vnn = vnn, sl = sl, sf = sf, snn = snn, zl = zl, zf = zf, znn = znn)

    return transition, terminal
end