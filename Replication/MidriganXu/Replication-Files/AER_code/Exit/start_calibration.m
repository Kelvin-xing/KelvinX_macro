clear;
clc;

format short;
% 17-inch: start

% kappau = 2.5;                % fixed cost of joining u
% theta  = 0.65; 
% xai    = 0.22; 
% rho    = 0.40;               % persistence transitory productivity
% se     = 0.90;
% F      = 0.3;
% W      = 1.01;

% 13 inch start

% kappau = 2.6;                % fixed cost of joining u
% theta  = 0.68; 
% xai    = 0.20; 
% rho    = 0.40;               % persistence transitory productivity
% se     = 0.95;
% F      = 0.3;
% W      = 1.017;

% 13 inch end: wins

% kappau = 2.6610;
% theta  = 0.6781; 
% xai    = 0.2190; 
% rho    = 0.4040;               % persistence transitory productivity
% se     = 0.9573;
% F      = 0.2671;
% W      = 1.0214;

% 17 inch end: 

kappau = 2.5215;
theta  = 0.6624; 
xai    = 0.2254; 
rho    = 0.4082;               % persistence transitory productivity
se     = 0.9269;
F      = 0.3058;
W      = 1.0176;



lb = [1.5;  0.5;  0.1;  0.2;   0.6;  0.05;  0.93];          
ub = [3.5;  0.8;  0.3;  0.5;   1.2;  0.45;  1.07];

% guess for calibrated parameters 

xstart = [kappau; theta; xai; rho; se; F; W];

 objective(xstart)
 break

%gaoptions = gaoptimset(gaoptions,'UseParallel','always','Display','iter','InitialPopulation',xstart);
gaoptions = gaoptimset('UseParallel', 'always', 'Display','iter','InitialPopulation',xstart');
x = ga(@objective,size(xstart,1),[],[],[],[],lb,ub,[],gaoptions); 

%x = neldmead_bounds('objective',xstart,lb,ub);

%save x x
