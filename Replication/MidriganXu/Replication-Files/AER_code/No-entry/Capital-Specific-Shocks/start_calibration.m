clear;
clc;

format short;


theta  = 0.5685; 

rho    = 0.3000;             % persistence transitory productivity
se     = 0.8310;

W      = 0.8381;
r      = 0.0311;

lb = [0.55;  0.3;  0.80;  0.82;   0.03];          
ub = [0.63;  0.4;  0.86;  0.85;   0.036];

% guess for calibrated parameters 

xstart = [theta; rho; se; W; r];

objective(xstart)
break


% Asset Demand vs. Supply       =     3.554      3.571 
% Labor Demand vs. Supply       =     1.001      1.000 
% 
% 
% 
% 
% Moments
% 
%                    1st col: model,      2nd col: data
% 
% s.d.  dlog(Y)                 =      0.58       0.59 
% s.d.   log(Y)                 =      1.31       1.31 
% 
% cor y y1                      =      0.90       0.90 
% cor y y3                      =      0.86       0.87 
% cor y y5                      =      0.86       0.85 
% 
% aggregate Debt to GDP         =      1.22       1.20 
%     0.0180
% 
% 
% ans =
% 
%     0.0180



%gaoptions = gaoptimset(gaoptions,'UseParallel','always','Display','iter','InitialPopulation',xstart);
%gaoptions = gaoptimset('UseParallel', 'always', 'Display','iter','InitialPopulation',xstart');
%x = ga(@objective,size(xstart,1),[],[],[],[],lb,ub,[],gaoptions); 



x = neldmead_bounds('objective',xstart,lb,ub);

save x x
