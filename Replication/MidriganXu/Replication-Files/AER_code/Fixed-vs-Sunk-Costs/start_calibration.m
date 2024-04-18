clear;
clc;

format short;

kappau = 1.1930;              % fixed cost of joining u

theta  = 0.8648; 
xai    = 0.1008; 

rho    = 0.2457;             % persistence transitory productivity
se     = 0.4958;

W      = 1.1031149;
r      = 0.0467542;

lb = [1.15;  0.81;  0.09;  0.15;   0.48;  1.09;  0.043];          
ub = [1.25;  0.89;  0.11;  0.25;   0.52;  1.11;  0.049];

% guess for calibrated parameters 

xstart = [kappau; theta; xai; rho; se; W; r];

objective(xstart)
break

%gaoptions = gaoptimset(gaoptions,'UseParallel','always','Display','iter','InitialPopulation',xstart);
%gaoptions = gaoptimset('UseParallel', 'always', 'Display','iter','InitialPopulation',xstart');
%x = ga(@objective,size(xstart,1),[],[],[],[],lb,ub,[],gaoptions); 

x = neldmead_bounds('objective',xstart,lb,ub);

save x x
