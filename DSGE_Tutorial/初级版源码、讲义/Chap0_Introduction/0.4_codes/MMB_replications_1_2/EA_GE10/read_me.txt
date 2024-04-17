EA_GE10 Replication

+ last change: 2011-08-24

+ replication: IRF to a monetary policy shock 
               IRF to a technology shock
               IRF to a investment specific shock
               IRF to a labor supply shock
	
Note: the IRFs in the paper are generated with an older version of the code (as being notified by the author), 
    	hence the replicated impulse responses do not match the IRFs in the paper precisely, 
	but they do match qualitatively

Note2: the published version of the model is implemented into the Macro Model Data Base 

        
+ replicated IRFs: 	EA_GE10_irf_monetary.pdf   
			EA_GE10_irf_investment.pdf
			EA_GE10_irf_labor_supply.pdf	
			EA_GE10_irf_technology.pdf 

+ software requirements: DYNARE version 4 

+ file to produce replicated IRFs: run.m (which calls EA_GE10_rep.mod in the folder EA_GE10_rep)

+ original IRFs: Figures 5-8 on pages 64-67 in Gelain (2010)

+ literature:
  - Gelain (2010): "The external finance premium in the Euro area:
	A dynamic stochastic general equilibrium analysis", North american Journal of Economics and Finance 21, 49-71.