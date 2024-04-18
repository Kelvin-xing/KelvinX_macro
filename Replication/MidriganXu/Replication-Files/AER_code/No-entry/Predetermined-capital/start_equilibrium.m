clear;
clc;

format short;

W      = 0.838838;
r      = 0.030724;

lb = [0.836;   0.030];          
ub = [0.840;   0.032];

% guess for calibrated parameters 

xstart = [W; r];

objective_equilibrium(xstart)
break


%gaoptions = gaoptimset(gaoptions,'UseParallel','always','Display','iter','InitialPopulation',xstart);
%gaoptions = gaoptimset('UseParallel', 'always', 'Display','iter','InitialPopulation',xstart');
%x = ga(@objective,size(xstart,1),[],[],[],[],lb,ub,[],gaoptions); 


x = neldmead_bounds('objective_equilibrium',xstart,lb,ub);
