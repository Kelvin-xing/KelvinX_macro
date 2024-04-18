clear;
clc;

format short;

global theta xai

theta = 0.25;
xai    = 0.1008;

W     = 0.9725;
r     = 0.0192;

lb = [0.95;  0.015];          
ub = [0.99;  0.025];

% guess for calibrated parameters 

xstart = [W; r];

objective_equilibrium(xstart)
break

%gaoptions = gaoptimset(gaoptions,'UseParallel','always','Display','iter','InitialPopulation',xstart);
%gaoptions = gaoptimset('UseParallel', 'always', 'Display','iter','InitialPopulation',xstart');
%x = ga(@objective_equilibrium,size(xstart,1),[],[],[],[],lb,ub,[],gaoptions); 

x = neldmead_bounds('objective_equilibrium',xstart,lb,ub);

save x x
