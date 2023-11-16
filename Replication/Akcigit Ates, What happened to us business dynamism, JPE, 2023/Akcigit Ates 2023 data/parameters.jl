# Structural parameters of the model
@with_kw struct StrParam

    # Externally calibrated parameters
    ρ::Float64       = 0.05
    τ::Float64
    s::Float64
    γ::Float64       = 1/0.35
    γTilde::Float64  = 1/0.35

    # Internally calibrated parameters
    α::Float64
    αTilde::Float64
    λ::Float64
    δ::Float64
    ϕ::Float64
    ϕTilde::Float64  = ϕ
    β::Float64
    ψ::Float64
end

# These parameters give proportional changes in structural parameters and their speeds (ν's) over the transition
@with_kw struct TransParam 

    # Proportional change
    τ::Float64 = 0.20 / 0.30
    s::Float64 = 0.20 / 0.05
    αTilde::Float64
    δ::Float64

    # Speed parameters
    ντ::Float64 = 1.0
    νs::Float64 = 1.0
    ναTilde::Float64
    νδ::Float64
end

# Other numeric parameters
@with_kw struct NumParam

    mbar::Int                   = 100
    dt::Float64                 = 1/50

    numoffirmsinasector::Int    = 2                                             # Number of firms in a sector
    numofsectors::Int           = 10000                                         # Number of sectors in economy
    numoffirms::Int             = numoffirmsinasector*numofsectors              # Number of firms in economy

    periodsperyear              = Int(1/dt)              
    tyears::Int                 = 150                                           # Transition years (simulation length as well)
    tperiods::Int               = periodsperyear*tyears                         # Transition periods
    tyearsburnin::Int           = 30                                            # Years to burn-in
    tperiodsburnin::Int         = periodsperyear*tyearsburnin                   # Periods to burn-in
    tyearpol::Int               = 35                                            # Year after when parameter transitions stopped (2015)
    tperiodpol::Int             = periodsperyear*tyearpol                       # Last period of 2015

    maxvalueiter::Int           = 10000                                         # Maximum number of iterations for value function and distribution iterations
    maxwageiter::Int            = 60                                            # Maximum number of iterations for wage iteration
    valuetol::Float64           = 1e-8                                          # Tolerance for value and distribution convergence
    wagetol::Float64            = 1e-4                                          # Tolerance for wage convergence
    wageweight::Float64         = 0.75                                          # Weight for wage update

    panelnumofcolumns           = 6                                             # Number of columns in simulations
    numofeventslev              = 6                                             # Number of events in unleveled sectors
    numofeventsnn               = 5                                             # Number of events in neck-and-neck sectors

    yearcol::Int                = 1
    firmidcol::Int              = 2
    stepcol::Int                = 3
    laborcol::Int               = 4
    agecol::Int                 = 5
    prodcol::Int                = 6
end
