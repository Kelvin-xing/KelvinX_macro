clear all;
clc;
tic
global rho alpha eta beta delta W r smin smax fspace cr egrid P theta lbar Pw  eta1 eta2 eta3 Perg 

beta      = 0.92;            

alpha     = 2/3;
%eta       = 0.85;               % decreasing returns
delta     = 0.06; 
Lbar      = 1;

spliorder = 1; 

lambdae  = 0.785;              % probability worker stays in e
lambdau  = 0.5;                % 1

%Data Targets

sdYdata   = 0.59; 
sYdata    = 1.31;
ac1Ydata  = 0.90;
ac3Ydata  = 0.87;
ac5Ydata  = 0.85;

DtoY      = 1.2;
EtoY      = 0.3;

eta1 = 0.55;
eta2 = 0.85;
eta3 = 0.95;

%calibrated parameters

theta  = 0.59; 
rho    = 0.25;               % persistence transitory productivity
se     = 0.72;

W      = 0.9284;
r      = 0.047;

% theta  = 0;
% W      = 0.857;

solve_workers;

disp('case I - low eta')
eta = eta1; 
solve_entrepreneurs;
ergodic;                       % compute distribution for these guys

Y1 = Ya;
K1 = Ka;
L1 = La; 
A1 = Aa;
C1 = n'*C;
n1 = n;
Debt1 = (K - state(:,1));

LK1 = L./K;
YLK1 = Y./(W*L + r*K); 

Y1i = Y; 
K1i = K;
L1i = L;

TFP1 = TFP;
TFP1best = TFPbest;
%save all relevant simulated micro-data for moments
simulate1
Ys1=Y;
Ls1=L;
Ks1=K;
es1=e;
rs1=rr;

disp('caseII - medium eta')
eta = eta2; 

solve_entrepreneurs;
ergodic;                       % compute distribution for these guys

Y2 = Ya;
K2 = Ka;
L2 = La; 
A2 = Aa;
C2 = n'*C;
n2 = n;
Debt2 = (K - state(:,1));

LK2 = L./K;
YLK2 = Y./(W*L + r*K); 

Y2i = Y; 
K2i = K;
L2i = L;

TFP2 = TFP;
TFP2best = TFPbest;
simulate1;
%save all relevant simulated micro-data for moments
Ys2=Y;
Ls2=L;
Ks2=K;
es2=e;
rs2=rr;

disp('caseIII - high eta')
eta = eta3;

solve_entrepreneurs;
ergodic;                       % compute distribution for these guys

Y3 = Ya;
K3 = Ka;
L3 = La; 
A3 = Aa;
C3 = n'*C;
n3 = n;
Debt3 = (K - state(:,1));

LK3 = L./K;
YLK3 = Y./(W*L + r*K); 

Y3i = Y; 
K3i = K;
L3i = L;

TFP3 = TFP;
TFP3best = TFPbest;

Yi = [Y1i; Y2i; Y3i];
Ki = [K1i; K2i; K3i];
Li = [L1i; L2i; L3i];
ni = [n1/3; n2/3; n3/3];
ei = [egrid(state(:,2)); egrid(state(:,2)); egrid(state(:,2))];

% Naive misallocation loss

etam = (eta1 + eta2 + eta3)/3;

TFPnaive = (ni'*(exp(ei).*(Yi./Ki).^(-(1-alpha)*etam/(1-etam)))).^(1-alpha*etam)/...
    (ni'*(exp(ei).*(Yi./Ki).^((alpha*etam-1)/(1-etam)))).^((1-alpha)*etam)/(ni'*(exp(ei))).^(1-etam);



simulate1;
%save all relevant simulated micro-data for moments
Ys3=Y;
Ls3=L;
Ks3=K;
es3=e;
rs3=rr;

Ya = 1/3*(Y1 + Y2 + Y3); 
Ka = 1/3*(K1 + K2 + K3); 
La = 1/3*(L1 + L2 + L3); 
DtoYm = 1/3*(n1'*Debt1 + n2'*Debt2 + n3'*Debt3)/Ya;

disp('Equilibrium Conditions')

fprintf('\n');

fprintf('Asset Demand vs. Supply       = %9.3f  %9.3f \n',   [Ka,  Aa + Aw]);
fprintf('Labor Demand vs. Supply       = %9.3f  %9.3f \n',   [La,  Lbar]);
fprintf('Debt to GDP                   = %9.3f  %9.3f \n',   [DtoYm,  DtoY]);

fprintf('\n');

sim_moments;

fprintf('\n');
disp('Model Implications')

fprintf('TFP  1                        = %9.3f \n',          TFP1);
fprintf('Misallocation Loss   1        = %9.3f \n',          log(TFP1best/TFP1));
fprintf('TFP  2                        = %9.3f \n',          TFP2);
fprintf('Misallocation Loss   2        = %9.3f \n',          log(TFP2best/TFP2));
fprintf('TFP  3                        = %9.3f \n',          TFP3);
fprintf('Misallocation Loss   3        = %9.3f \n',          log(TFP3best/TFP3));


% Compute total losses from misallocation

Kbar = Ka; 
Lbar = La;
Ybar = Ya;

optset('broyden', 'showiters', 1);
x = broyden('equilibrium', [W; r+delta], Kbar, Lbar);


W = x(1);
R = x(2); 

% Compute total losses from misallocation

eta = eta1; 

L1   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*R.^(-(1-alpha)*eta/(1-eta)).*exp(egrid);
K1   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*R.^((alpha*eta-1)/(1-eta)).*exp(egrid);
Y1   = exp(egrid).^(1-eta).*(L1.^alpha.*K1.^(1-alpha)).^eta;

eta = eta2; 

L2   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*R.^(-(1-alpha)*eta/(1-eta)).*exp(egrid);
K2   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*R.^((alpha*eta-1)/(1-eta)).*exp(egrid);
Y2   = exp(egrid).^(1-eta).*(L2.^alpha.*K2.^(1-alpha)).^eta;

eta = eta3; 

L3   = eta^(1/(1-eta))*alpha^((1-(1-alpha)*eta)/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*W^(((1-alpha)*eta-1)/(1-eta)).*R.^(-(1-alpha)*eta/(1-eta)).*exp(egrid);
K3   = eta^(1/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha*eta)/(1-eta))*W^(-alpha*eta/(1-eta)).*R.^((alpha*eta-1)/(1-eta)).*exp(egrid);
Y3   = exp(egrid).^(1-eta).*(L3.^alpha.*K3.^(1-alpha)).^eta;


fprintf('Planner L        = %9.3f %9.3f \n',     [1/3*Perg'*(L1 + L2 + L3), Lbar ]);
fprintf('Planner K        = %9.3f %9.3f \n',     [1/3*Perg'*(K1 + K2 + K3), Kbar ]);
fprintf('Planner Y        = %9.3f %9.3f \n',     [1/3*Perg'*(Y1 + Y2 + Y3), Ybar ]);
fprintf('Total TFP loss   = %9.3f  \n',          log(1/3*Perg'*(Y1 + Y2 + Y3)/ Ybar)*100);
fprintf('Naive Misallocation Loss Modern     = %9.3f \n',          -log(TFPnaive)*100);


