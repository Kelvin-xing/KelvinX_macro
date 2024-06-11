%% Alvarez Lucas 2007 JME

close all;
warning off

% Endogenous
var

end;

% Exogenous
varexo

end;


% Parameters
parameters

end;

% Model
model;

end;


% steady state
steady_state_model;

end;

steady;
check;

% Shocks
shocks;
var z; stderr 0.01; 

end;

% IRFs
stoch_simul(irf=40,order=1, hp_filter=1600,pruning) w;



