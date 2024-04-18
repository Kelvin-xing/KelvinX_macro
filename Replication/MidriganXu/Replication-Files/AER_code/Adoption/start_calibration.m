clear;
clc;

format short;

kappau  = 0.30272992403060;              % fixed cost of joining u

theta   = 0.78219600032876; 
xai     = 0.08091351214415; 

rho     = 0.11026196075805;              % persistence transitory productivity
se      = 0.5047725543333;              % stand. dev. transitory shocks

W       = 1.57427065589513;
r       = 0.04001885482174;

lb = [0.28;  0.76;  0.07;  0.1;   0.40;  1.55;  0.037];          
ub = [0.32;  0.80;  0.09;  0.4;   0.60;  1.59;  0.043];

% guess for calibrated parameters 

xstart = [kappau; theta; xai; rho; se; W; r];

 objective(xstart)
 break


% Equilibrium Conditions
% 
% Asset Demand vs. Supply       =     7.319      7.299 
% Labor Demand vs. Supply       =     0.997      1.000 
% 
% 
% Moments
% 
%                    1st col: model,      2nd col: data
% 
% s.d.  dlog(Y)                 =      0.56       0.59 
% s.d.   log(Y)                 =      1.37       1.31 
% 
% cor y y1                      =      0.90       0.90 
% cor y y3                      =      0.83       0.87 
% cor y y5                      =      0.80       0.85 
% 
% aggregate Debt to GDP         =      1.16       1.20 
% Market capitaliz to GDP       =      0.30       0.30 
% 
% Fixed Inv  to total Y, perc   =       4.6        4.6 
% 
% Variance exog component       =      1.43 
%     0.0267


if 1 

%gaoptions = gaoptimset(gaoptions,'UseParallel','always','Display','iter','InitialPopulation',xstart);
gaoptions = gaoptimset('UseParallel', 'always', 'Display','iter','InitialPopulation',xstart');
x = ga(@objective,size(xstart,1),[],[],[],[],lb,ub,[],gaoptions); 

else

x = neldmead_bounds('objective',xstart,lb,ub);

end

save x x
