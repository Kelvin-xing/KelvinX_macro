clear;
clc;

format short;

W   = 1.5721;
r   = 0.04019;

lb  = [1.55;  0.036];          
ub  = [1.59;  0.044];

% guess for calibrated parameters 

xstart = [W; r];

objective_equilibrium(xstart)
break

%gaoptions = gaoptimset(gaoptions,'UseParallel','always','Display','iter','InitialPopulation',xstart);
gaoptions = gaoptimset('UseParallel', 'always', 'Display','iter','InitialPopulation',xstart');
x = ga(@objective_equilibrium,size(xstart,1),[],[],[],[],lb,ub,[],gaoptions); 

%x = neldmead_bounds('objective_equilibrium',xstart,lb,ub);

save x x
