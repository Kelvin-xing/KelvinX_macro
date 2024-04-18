clear all;
clc;

global rhom rhot alpha eta beta delta W r smin smax fspace cr znode wnode  kappa pait ...
       theta mu phi

mu        = 1.08;              % growth rate
beta      = 0.92*mu;            

alpha     = 2/3;
eta       = 0.85;              % decreasing returns
delta     = 0.06; 
Lbar      = 1;

spliorder = 1; 

% Targets

sdYdata   = 0.59; 
sYdata    = 1.31;
ac1Ydata  = 0.90;
ac3Ydata  = 0.87;
ac5Ydata  = 0.85;

DtoY      =  1.2;
DtoE      = (30.2 - 3.96 - 3.03 - 2.15)/7;      % from  Outline of FSA Korea
DtoA      = (46.6 - 5.42 - 5.08 - 4.28)/7/10; 
EtoY      =  0.3;

% Version 1:

phi     = 1/(1-eta)*0.20;   % productivity gap modern
kappa   = 0.10;             % fixed cost of joining modern sector
rhom    = 0.9999;           % probability to keep z of modern guy
rhot    = 0.50;             % probability to keep z of trad. guys

sz      = 2.5;              % std. dev. of permanent component               

theta   = 0.9999; 
W       = 1.523;
r       = 0.0395;

% theta   = 0; 
% W       = 1.33;
% r       = -0.0599;


k       = 51; 

[znode, wnode]  = qnwnorm(k, 0, sz^2); 

znode = znode(wnode>1e-3); 
wnode = wnode(wnode>1e-3); 
wnode = wnode/sum(wnode);

k = size(wnode,1);


pait    = (1-eta)*(eta/W)^(eta/(1-eta));
paim    = (1-eta)*eta^(eta/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*...
           W^(-alpha*eta/(1-eta))*(r+delta)^(-(1-alpha)*eta/(1-eta));
             
%%%%%%%%%%%%%%%%%%%% START Entrepreneurs %%%%%%%%%%%%%%%%%%%%%

fprintf('Solve Entrepreneurs Problem...');

fprintf('\n');
fprintf('\n');

amin =  exp(-7)+.05;                      % tricky right at 0, so move a bit to be able to solve it        
amax =  3*exp(znode(end))*exp(phi);         % bounds for assets

zmin = 1;
zmax = k;

smin = [amin, zmin];
smax = [amax, zmax];

n = [101, k];

curv   = .25;                       % the close to 0, the closer to log-scale
agrid  =  nodeunif(n(1), 0, (amax - amin).^curv).^(1/curv) + amin;             % make grid have more nodes close to lower bound, then upper bound
agrid  =  exp(nodeunif(n(1), log(amin), log(amax)));

fspace = fundef({'spli', agrid, 0, spliorder},...
                {'spli', znode, 0, 1});
           
grid = funnode(fspace);            % 4d grid of collocation nodes 
s    = gridmake(grid);             % collection of  states
ns   = length(s);

solvem_static;

c  =  funfitxy(fspace, s, zeros(ns,4));      % guess for coefficients

load c

for i = 1 : 1000                               % do a couple to give better starting values for newton

cnew = c;

[v1, v2, v3, v4]     = saveBelmax(cnew, fspace, s);

c      = funfitxy(fspace, s, [v1, v2, v3, v4]); 

fprintf('%4i %6.2e \n',[i, norm((c-cnew))/norm(c)]);    

if norm((c-cnew))/norm(c)< 1e-5 , break, end

end

c = vec(c);

for i = 1 : 5*0       %newton iterations
    
cnew = c;   

[bel, beljac] = solvebel(cnew, fspace ,s);        % RHS - LHS + derivatives

c  = cnew - (beljac\bel);                         % Newton step

fprintf('%4i %6.2e \n',[i, norm((c-cnew))/norm(c)]);    
if norm((c-cnew))/norm(c) < 1e-7 , break, end

end

c = [c(1:ns,:), c(ns+1:2*ns,:), c(2*ns+1:3*ns,:), c(3*ns+1:4*ns,:)];


[~, ~, ~, ~,  xm, xt] = saveBelmax(c, fspace, s); 
cx = funfitxy(fspace, s, [xm, xt]);   

[~, vtm, vtt] = valfunc2t(c, fspace, s);
cst = funfitxy(fspace, s, [vtm, vtt]); 
    
save c c

clear bel beljac

clc
ergodic
simulate






