using Parameters, NLsolve, StatsBase, Plots, Optim

# Following codes define parameters and fundamental functions used throughout the solution
include("parameters.jl")
include("fundamentals.jl")

# Following codes define solver functions
include("bgp_solver.jl")
include("transition_solver.jl")

# Following codes define simulation functions and other functions that help compute several moments from the model
include("simulation.jl")
include("moments.jl")

## Define Parameters
numericp = NumParam(
    mbar = 100,
    dt = 1/50,
    tyears = 150,
    tyearsburnin = 30,
    tyearpol = 35
)

initialp = StrParam(
    ρ = 0.05,
    τ = 0.3,
    s = 0.05,
    γ = 1/0.35,
    γTilde = 1/0.35,
    α = 0.0065197905757763555,
    αTilde = 0.5654045445345105,
    λ = 1.0093991854956592,
    δ = 0.08378196909781864,
    ϕ = 0.00923568449429258,
    ϕTilde = 0.00923568449429258,
    β = 0.9781494840234898,
    ψ = 0.8649852401597248
)

transitionp = TransParam(
    τ = 0.2/0.3,
    s = 0.2/0.05,
    αTilde = 2.262104086275848,
    δ = 0.2936794269870019,
    ντ = 1.0,
    νs = 1.0,
    ναTilde = 0.03244688686336543,
    νδ = 7.1814073609357205
)

## Solve Initial BGP
initialbgp = BgpSolver(initialp, numericp)

## Simulate Initial BGP
panel, sectordata = SimPrep(numericp)
entrantid_initialbgp = InitialBgpSimulation!(panel, sectordata, initialbgp)

## Calculate Initial BGP Moments
initialbgp_moments = Moments(initialbgp, panel)

## Solve Transition
transition, terminalbgp = TransitionSolver(initialbgp, transitionp, numericp)

## Simulate Transition
TransitionSimulation!(panel, sectordata, transition, initialbgp, entrantid_initialbgp)

## Calculate Transition Moments
transition_moments = Moments(transition, panel, sectordata)

