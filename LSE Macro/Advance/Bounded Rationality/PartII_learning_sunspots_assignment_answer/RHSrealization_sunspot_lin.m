% Compute the RHS for realization shock
function xxx=RHSrealization_sunspot_lin(epsi,par)

c_tilde_next = par.JJ(2,1)*log(par.K_t/par.K_ss)   +par.JJ(2,2)*log(par.C_t/par.C_ss)+par.lin_dzeta*epsi;
C_next       = exp(c_tilde_next + log(par.C_ss));
mu_next      = C_next^(-par.nu);

L_next       = (par.b*mu_next*par.labda1*par.K_next^par.alpha/par.labda2)^(1/(par.chi+1-par.beta));
xxx          = par.dfactor*mu_next*(par.a*par.labda1*par.K_next^(par.alpha-1)*L_next^par.beta+1-par.delta);
%disp([C_next mu_next par.K_next L_next])
%disp([par.C_ss par.C_ss^(-par.nu) par.K_ss par.L_ss])
%pause


% **********************************************************************
% **********************************************************************

