% Compute the RHS for realization shock
function xxx=RHSrealization_sunspot2(epsi,par)

state1_next  = [ 1 log(par.K_t/par.K_ss) log(par.C_t/par.C_ss) ];
state2_next  = [ epsi*pi/2 ];

mu_next      = exp(state1_next*par.eta1)+par.eta_dzeta*sin(state2_next*par.eta2);
%C_next       = mu_next^(-1/par.nu);
L_next       = (par.b*mu_next*par.labda1*par.K_next^par.alpha/par.labda2)^(1/(par.chi+1-par.beta));
xxx          = par.dfactor*mu_next*(par.a*par.labda1*par.K_next^(par.alpha-1)*L_next^par.beta+1-par.delta);
%disp([C_next mu_next par.K_next L_next])
%disp([par.C_ss par.C_ss^(-par.nu) par.K_ss par.L_ss])
%pause


% **********************************************************************
% **********************************************************************

