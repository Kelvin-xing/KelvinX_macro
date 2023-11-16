# We first call some necessary packages
using Parameters, NLsolve, StatsBase, Plots, Plots.PlotMeasures, Optim

# Following codes define parameters and fundamental functions used throughout the solution
include("parameters.jl")
include("fundamentals.jl")

# Following codes define solver functions
include("bgp_solver.jl")
include("transition_solver.jl")

# Following codes define simulation functions and other functions that help us compute several moments from the model
include("simulation.jl")
include("moments.jl")

# Following codes define calibration functions and minimization procedure
include("calibration.jl")

# Following codes include several exercises of section 5, 6, and 7. 
include("exercises.jl")

# Define numerical parameters packed in numericp object
numericp = NumParam(
    mbar = 100,
    dt = 1/50,
    tyears = 150,
    tyearsburnin = 30,
    tyearpol = 35
)

#=
====================
PART 1 - CALIBRATION
====================
=#

## INITIAL BGP CALIBRATION
initialbgp_targets = BgpTargets(growth = 1.37, entry = 12.0, rnd = 2.4, profit = 6.0, markup = 10.0, gd = 54.0, ec = 30.0)
initialp_guess = [0.006220049811151615, 0.5783604690108269, 1.0096842992448012, 0.09379516372885797, 0.0200, 0.9784870390300153, 0.8513377282566288]
iterations = 500

initialp_calibrated_vector, initialbgp_calibration = BgpCalibration(numericp, initialbgp_targets, initialp_guess, iterations)

initialp_calibrated = StrParam(
    ρ = 0.05,
    τ = 0.30,
    s = 0.05,
    γ = 1/0.35,
    γTilde = 1/0.35,
    α = initialp_calibrated_vector[1],
    αTilde = initialp_calibrated_vector[2],
    λ = initialp_calibrated_vector[3],
    δ = initialp_calibrated_vector[4],
    ϕ = initialp_calibrated_vector[5],
    ϕTilde = initialp_calibrated_vector[5],
    β = initialp_calibrated_vector[6],
    ψ = initialp_calibrated_vector[7]
)

# Solve calibrated initial bgp
initialbgp_calibrated = BgpSolver(initialp_calibrated, numericp)

# Initial bgp simulation of sectors and firms
panel, sectordata = SimPrep(numericp)
entrantid_initialbgp_calibrated = InitialBgpSimulation!(panel, sectordata, initialbgp_calibrated)

# Calculate calibrated initial bgp moments
initialbgp_moments_calibrated = Moments(initialbgp_calibrated, panel)

## TRANSITION CALIBRATION
transition_targets = TransitionTargets(entry = 8.0, gd = 42.0, markup = 22.0, entrydeclinein10 = 15.0, entrydeclinein20 = 25.0)
transitionp_guess = [2.262104086275848, 0.2936794269870019, 0.03244688686336543, 7.1814073609357205]
iterations = 300
entry1980 = initialbgp_moments_calibrated.entry
panelinitialbgplastyear_calibrated = panel[(numericp.tyears - 1)*numericp.numoffirms + 1:end, :]
sectordatainitialbgplastperiod_calibrated = sectordata[:, end]

transitionp_calibrated_vector, transition_calibration = TransitionCalibration(initialbgp_calibrated, transition_targets, transitionp_guess, iterations, entry1980, entrantid_initialbgp_calibrated, panelinitialbgplastyear_calibrated, sectordatainitialbgplastperiod_calibrated)

transitionp_calibrated = TransParam(
    τ = 0.20/0.30,
    s = 0.20/0.05,
    αTilde = transitionp_calibrated_vector[1],
    δ = transitionp_calibrated_vector[2],
    ντ = 1.0,
    νs = 1.0,
    ναTilde = transitionp_calibrated_vector[3],
    νδ = transitionp_calibrated_vector[4]
)

# Solve calibrated transition
transition_calibrated, _ = TransitionSolver(initialbgp_calibrated, transitionp_calibrated, numericp)

# Transition simulation of sectors and firms
TransitionSimulation!(panel, sectordata, transition_calibrated, entrantid_initialbgp_calibrated, panelinitialbgplastyear_calibrated, sectordatainitialbgplastperiod_calibrated)

# Calculate calibrated transition moments as time series
transition_moments_calibrated = Moments(transition_calibrated, panel, sectordata)

#=
====================
PART 2 - REPLICATION
====================
=#

## CALIBRATED PARAMETERS USED IN THE PAPER - REPLICATION OF RESULTS
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
    τ = 0.20/0.30,
    s = 0.20/0.05,
    αTilde = 2.442926303158954,
    δ = 0.3565012860205644,
    ντ = 1.0,
    νs = 1.0,
    ναTilde = 0.019361861209644962,
    νδ = 7.360306879059708
)

# Solve initial bgp
initialbgp = BgpSolver(initialp, numericp)

# Table 1
begin
    println("TABLE 1: LIST OF PARAMETER VALUES (BGP)")
    println("---------------------------------------")
    println("Panel A: Externally calibrated parameters")
    println("ρ           = $(round(100*initialbgp.p.ρ, digits = 3))%")
    println("γ, γTilde   = $(round(initialbgp.p.γ, digits = 3))")    
    println("τ           = $(round(100*initialbgp.p.τ, digits = 3))%")    
    println("s           = $(round(100*initialbgp.p.s, digits = 3))%")
    println()
    println("Panel B: Internally calibrated parameters")
    println("β           = $(round(initialbgp.p.β, digits = 3))")    
    println("λ           = $(round(initialbgp.p.λ, digits = 3))")    
    println("ψ           = $(round(initialbgp.p.ψ, digits = 3))")    
    println("α           = $(round(initialbgp.p.α, digits = 3))")    
    println("αTilde      = $(round(initialbgp.p.αTilde, digits = 3))")    
    println("δ           = $(round(100*initialbgp.p.δ, digits = 2))%")    
    println("ϕ = ϕTilde  = $(round(100*initialbgp.p.ϕ, digits = 2))%")    
end

# Simulate initial bgp sectors and firms
panel, sectordata = SimPrep(numericp)
entrantid_initialbgp = InitialBgpSimulation!(panel, sectordata, initialbgp)
panelinitialbgplastyear = panel[(numericp.tyears - 1)*numericp.numoffirms + 1:end, :]
sectordatainitialbgplastperiod = sectordata[:, end]

# Calculate initial bgp moments: Note that net entry contribution and firm growth dispersion are calculated by simulation. Therefore, they may different than reported value in the paper
initialbgp_moments = Moments(initialbgp, panel)

# Solve transition
transition, terminalbgp = TransitionSolver(initialbgp, transitionp, numericp)

# Table 4
begin
    println("TABLE 4: LIST OF PARAMETER VALUES (TRANSITION)")
    println("----------------------------------------------")
    println("Panel A: Externally calibrated parameters")
    println("τ_T         = $(round(100*transition.sequencep[transition.numericp.tperiodpol].τ, digits = 3))%")    
    println("s_T         = $(round(100*transition.sequencep[transition.numericp.tperiodpol].s, digits = 3))%")
    println("ν_τ         = $(round(transitionp.ντ, digits = 3))")
    println("ν_s         = $(round(transitionp.νs, digits = 3))")    
    println()
    println("Panel B: Internally calibrated parameters")
    println("αTilde      = $(round(transition.sequencep[transition.numericp.tperiodpol].αTilde, digits = 3))")    
    println("δ           = $(round(100*transition.sequencep[transition.numericp.tperiodpol].δ, digits = 2))%")
    println("ν_αTilde    = $(round(transitionp.ναTilde, digits = 3))")
    println("ν_δ         = $(round(transitionp.νδ, digits = 3))")
end

# Simulate transition: Connecting this simulation to initial bgp simulation via giving panelinitialbgplastyear and sectordatainitialbgplastperiod as arguments of the simulation function is suggested for consistency of entry and firm ages after 1981, which is the first year of transition.
TransitionSimulation!(panel, sectordata, transition, entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod)

# Calculate transition moments
transition_moments = Moments(transition, panel, sectordata)

# Table 5: Level changes in untargeted moments
begin
    println("TABLE 5: LEVEL CHANGES IN UNTARGETED MOMENTS")
    println("--------------------------------------------")
    println("Δ Concentration            = $(round(transition_moments.conc[numericp.tyearpol] - initialbgp_moments.conc, digits = 1))%")
    println("Δ Profit share             = $(round(transition_moments.profit[numericp.tyearpol] - initialbgp_moments.profit, digits = 1))%")
    println("Δ Labor share              = $(round((1-transition_moments.profit[numericp.tyearpol]) - (1-initialbgp_moments.profit), digits = 1))%")
    println("Δ Young employment share   = $(round(transition_moments.yes[numericp.tyearpol] - initialbgp_moments.yes, digits = 1))%")
    println("Δ Job reallocation         = $(round(transition_moments.jr[numericp.tyearpol] - initialbgp_moments.jr, digits = 1))%")
    println("Δ Productivity gap         = $(round(transition_moments.prodgap[numericp.tyearpol] - initialbgp_moments.prodgap, digits = 1))%")
    println("Δ TFPR stdev.              = $(round(transition_moments.prodgapstd[numericp.tyearpol] - initialbgp_moments.prodgapstd, digits = 1))%")
end

# Replication of Figure 5
begin
    plot(title = "Transition paths of αTilde and δ", grid = false, xlabel = "Year", right_margin = 15mm, left_margin = 5mm, framestyle = :box)
    plot!(getproperty.(transition.sequencep, :δ)[1:numericp.tperiodpol], lw = 3.5, label = "Diffusion", legend = :topleft, ylabel = "Knowledge Diffusion Parameter Path", xticks = (0:5*numericp.periodsperyear:numericp.tyearpol*numericp.periodsperyear, 1980:5:2015), ylim = (0.02, 0.09), yticks = 0.02:0.01:0.09)
    plot!(twinx(), getproperty.(transition.sequencep, :αTilde)[1:numericp.tperiodpol], lw = 3.5, label = "Entry cost", legend = :topright, ls = :dash, color = :red, ylabel = "Entry Cost Parameter Path", xticks = (0:5*numericp.periodsperyear:numericp.tyearpol*numericp.periodsperyear, 1980:5:2015), ylims = (0.5, 1.4), yticks = 0.5:0.1:1.4)
end

# Replication of Figure D.4 in the Appendix
begin
    plot(title = "Transition paths of s and τ", grid = false, xlabel = "Year", right_margin = 15mm, left_margin = 5mm, framestyle = :box)
    plot!(getproperty.(transition.sequencep, :s)[1:numericp.tperiodpol], lw = 3.5, label = "Subsidy", legend = :topleft, ylabel = "Subsidy Rate Parameter Path", xticks = (0:5*numericp.periodsperyear:numericp.tyearpol*numericp.periodsperyear, 1980:5:2015), ylims = (0.05, 0.2), yticks = 0.05:0.05:0.2)
    plot!(twinx(), getproperty.(transition.sequencep, :τ)[1:numericp.tperiodpol], lw = 3.5, label = "Corporate", legend = :topright, ls = :dash, color = :red, ylabel = "Tax Rate Parameter Path", xticks = (0:5*numericp.periodsperyear:numericp.tyearpol*numericp.periodsperyear, 1980:5:2015), ylims = (0.2, 0.32), yticks = 0.2:0.02:0.32)
end

# Replication of Figure 6 (without data)
begin
    ll = @layout [a b; [_ c{0.5w} _]]
    
    plt_entry = plot(transition_moments.entry[1:numericp.tyearpol], lw = 3.5, legend = false, color = :black, xticks = (0:5:numericp.tyearpol, 1980:5:2015), frame = :box, grid = false, title = "Entry", titlefontsize = 12, margin = 5mm)
    plt_markup = plot(transition_moments.markup[1:numericp.tyearpol], lw = 3.5, legend = false, color = :black, xticks = (0:5:numericp.tyearpol, 1980:5:2015), frame = :box, grid = false, title = "Markup", titlefontsize = 12, margin = 5mm)
    plt_gd = plot(transition_moments.gd[1:numericp.tyearpol], lw = 3.5, legend = false, color = :black, xticks = (0:5:numericp.tyearpol, 1980:5:2015), frame = :box, grid = false, title = "Dispersion of firm growth", titlefontsize = 12, margin = 5mm)

    plot(plt_entry, plt_markup, plt_gd, layout = ll, size = (1000, 750), margin = 5mm)
end

# SECTION 5 - UNDERSTANDING THE MECHANISMS IN ISOLATION: In this section, we can replicate Figure 3, Table 3, and Figure 4 in the paper. In matrix table3, rows correspond to moments with the same order as in the paper. Similarly, first column corresponds to data column, and remaining columns correspond to channels with the same order as in the paper. 
# Note that last two columns, gross job reallocation and dispersion of firm growth, are calculated by simulation. Therefore, for these rows, there might be occasional differences between the paper and the replication here. 
# Note that the function Section5WithCalibration below calibrates the parameters from the beginning. Therefore, the results obtained might be different than the results in the paper. For exact replication, see below.
guess_αTilde = [transitionp.αTilde, transitionp.ναTilde]
guess_δ = [transitionp.δ, transitionp.νδ]
iterations = 100

figure3, table3, figure4 = Section5WithCalibration(initialbgp, initialbgp_moments, numericp, entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod, guess_αTilde, guess_δ, iterations)

# Below are calibrated parameters for the analysis in Section 5 of the paper. After defining these, we can simply call function Section5 to replicate Figure 3, Table 3, and Figure 4 in the paper. These objects are figure3_paper, table3_paper, and figure4_paper. In table3_paper, we have 1 = up arrow, 0 = horizontal arrow, -1 = down arrow.
section5_αTilde_calibrated = 2.9989384492968845
section5_ναTilde_calibrated = 0.2865843289484703
section5_δ_calibrated = 0.1012622345133649
section5_νδ_calibrated = 6.191182509094215

figure3_paper, table3_paper, figure4_paper = Section5(initialbgp, initialbgp_moments, numericp, entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod, section5_αTilde_calibrated, section5_ναTilde_calibrated, section5_δ_calibrated, section5_νδ_calibrated)

## SECTION 6 - INVESTIGATION OF JOINT FORCES: In matrix table6, rows correspond to moments with the same order as in the paper, while columns correspond to each channel. Similar to table3 object previously, last two rows are calculated by simulating the model. Therefore, numbers might differ slightly when replicating them.
table6 = Section6(numericp, initialbgp, initialbgp_moments, panelinitialbgplastyear, sectordatainitialbgplastperiod, entrantid_initialbgp, transitionp, transition_moments)

## SECTION 7 - MARKET POWER AND WELFARE: In this section, we replicate Figure 7 with its two panels a and b
figure7a = Section7_Figure7a(initialbgp)
# ATTENTION: Drawing figure 7b requires simulating a large number of sectors. Therefore, it can take a lot of time.  
figure7b = Section7_Figure7b(initialbgp)

## SECTION 8 - ALTERNATIVE MECHANISMS: In this section, we replicate Table 7 of the paper.
# Definitions and explanations for each exercise can be found in function definitions in exercise.jl
# As before, 1 means up arrow, 0 means horizontal arrow (no change), -1 means down error in the table. Order of rows are the same as in the paper.
column6_table7 = DecliningInterestRate(0.01, initialbgp, initialbgp_moments, entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod)
column7_table7 = IdeasGettingHarder(10*initialbgp.p.α, 10*initialbgp.p.αTilde, initialbgp, initialbgp_moments, entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod)
column8_table7 = WeakerLaborPower(1 + 2.1669517156287204*(initialbgp.p.λ - 1), initialbgp, initialbgp_moments, entrantid_initialbgp, panelinitialbgplastyear, sectordatainitialbgplastperiod)

## APPENDIX D - ADDITIONAL QUANTITATIVE RESULTS. In this section we replicate Table D.1. Note that elasticities for the moments that are obtained by simulating the model can be different than what is reported in the paper. These moments are entry contribution and growth dispersion.
table_d1 = Sensitivity(initialbgp, initialbgp_moments)