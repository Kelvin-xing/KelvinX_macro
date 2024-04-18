clear all;
clc;

global rho alpha eta beta delta W r smin smax fspacep fspaceu ...
       fspace cp cu cr egrid P kappap kappau pait sminp sminu ...
       theta mu lbar Pw phip phiu xai cinj fspacei

mu        = 1.08;               % growth rate
beta      = 0.92*mu;            

nuk       = 1;                  % fraction of kappa that is pledgeable

alpha     = 2/3;
eta       = 0.85;               % decreasing returns
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
DtoE      = (30.2 - 3.96 - 3.03 - 2.15)/7; % from  Outline of FSA Korea
DtoA      = (46.6 - 5.42 - 5.08 - 4.28)/7/10; 
EtoY      = 0.3;

lambdae  = 0.785;            % probability worker stays in e
lambdau  = 0.5;              % 1

% Version 1: 

phiu    = 1/(1-eta)*0.20;    % productivity gap unproductive
phip    = phiu;              % productivity gap productive

kappau = 1.1930;             % fixed cost of joining u

theta  = 0.8648; 
xai    = 0.1008; 

rho    = 0.2457;             % persistence transitory productivity
se     = 0.4958;

W      = 1.103114901451660;
r      = 0.046754247588991;

kappap  = kappau + 1;        % fixed cost of joining p

findvarz = 0;
bunp     = 0;    % need to raise slightly in some experiments the lower bound on assets to ensure convergence
curvs    = [];

% % Experiment without transitory shocks
% 
% se = 1e-7;
% r  = 0.058;
% W  = 1.035;
% 
% theta = 0.50;
% xai   = 0;
% r = 0.012;
% W = 0.935;

% CLOSED ECONOMY, xai = original

% theta = 0.75;
% W     = 1.084960504273230;
% r     = 0.042286399236787;
 
% theta = 0.50;
% W     = 1.023325366622690;
% r     = 0.032559950879239;

% theta = 0.25;      
% W     = 0.9725;
% r     = 0.0192;
% bunp  = 0.01;     % raise unproductive lower bound 0.01 (for computations to converge)

% theta = 0.00;
% W     = 0.8787;
% r     = -0.0599;
% bunp  = 0.01;

% OPEN ECONOMY, xai = original 

% theta = 0.75; 
% W     = 1.064602035768120; 

% theta = 0.50; 
% W     = 0.964761924297550; 
 
% theta = 0.25; 
% W     = 0.918225560351186; 
 
% theta = 0;
% W     = 0.8962;

% CLOSED ECONOMY, xai = 0

% theta = 0.75;
% xai   = 0;
% W     = 1.009979949035960;
% r     = 0.031028750964899;

% theta = 0.50;
% xai   = 0;
% W     = 0.957318767387261;
% r     = 0.010324465349904;

% theta = 0.25;
% xai   = 0;
% W     = 0.915587528197036;
% r     = -0.0599;
% bunp  = 0.05;
% curvs = 0.5; 

% OPEN ECONOMY, xai = 0

% theta = 0.75;
% xai   = 0;
% W     = 0.985;

% theta = 0.50;
% xai   = 0;
% W     = 0.9374;

% theta = 0.25;
% xai   = 0;
% W     = 0.9117;

% theta = 0;
% xai   = 0; 
% W     = 0.8962;

%%%%%% START WORKERS  %%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Solve Worker Problem...');

fprintf('\n');
fprintf('\n');

Pw    = [lambdau,     1-lambdau; 
         1-lambdae,   lambdae];
  
     
Pwerg = Pw^5000; Pwerg = Pwerg(1,:)';
lbar  = Lbar/Pwerg(end);                                        % make sure total sums up to 1

smin = 1e-5;
smax = 50;
 
n = [101, 2];

curv  = .25;                                                    % the closer to 0, the closer to log-scale
agrid =  nodeunif(n(1), smin.^curv, smax.^curv).^(1/curv);      % make grid have more nodes close to lower bound, then upper bound

fspace = fundef({'spli',  agrid, 0, spliorder},...
                {'spli', (1:2)', 0, 1});
           
grid = funnode(fspace);                 % 4d grid of collocation nodes 
s    = gridmake(grid);                  % collection of  states

c = zeros(size(s,1),1);

%load c

for i = 1 : 2                           % do a couple to give better starting values for newton

cnew = c;

v = saveBelmaxw(cnew,fspace,s);
c = funfitxy(fspace,s, v);              % update coefficients by projecting values on function space

fprintf('%4i %6.2e \n',[i,norm((c-cnew))/norm(c)]);    

if norm((c-cnew))/norm(c)< 1e-5 , break, end

end

c = vec(c);

for i = 1 : 50       % newton iterations
    
cnew = c;   

[bel, beljac] = solvebelw(cnew,fspace,s);        % RHS - LHS + derivatives

c  = cnew - (beljac\bel);                        % Newton step

fprintf('%4i %6.2e \n',[i, norm((c-cnew))/norm(c)]);    
if norm((c-cnew))/norm(c) < 1e-7 , break, end

end

[~, x] = saveBelmaxw(c,fspace,s); 
cx = funfitxy(fspace,s,x);   

fspacew = fspace;  % different state-space
cxw     = cx;

ergodicw; 

save c c
%simulatew;

%%%%%%%%%%%%%%%%%%%% START Productive Modern %%%%%%%%%%%%%%%%%%%%%
fprintf('\n');

fprintf('Solve Productive Sector Problem...');

fprintf('\n');
fprintf('\n');

k = 9;   % number of nodes for e

q  = (rho+1)/2;
nu = sqrt((k-1)/(1-rho^2))*se;

P = [q 1-q; 1-q q];

for i = 2:k-1
    
   P =     q*[P zeros(i,1); zeros(1,i+1)] + (1-q)*[zeros(i,1) P; zeros(1,i+1)] + ...
       (1-q)*[zeros(1,i+1); P zeros(i,1)] +     q*[zeros(1,i+1); zeros(i,1) P];
   
   P(2:i,:) = P(2:i,:)/2;
   
end

Perg  = P^5000; Perg = Perg(1,:)';
egrid = linspace(-nu, +nu ,k)';

pait  = (1-eta)*(eta/W)^(eta/(1-eta));
paim = (1-eta)*eta^(eta/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*...
        W^(-alpha*eta/(1-eta))*(r+delta)^(-(1-alpha)*eta/(1-eta));

fprintf('\n');
disp('Fixed cost to join unproductive sector relative to traditional sector profits')
disp(kappau/pait);
fprintf('\n')
disp('Fixed cost to join productive sector relative to unproductive sector profits')
disp(kappap/paim/exp(phiu));
fprintf('\n')
fprintf('\n')

solvep_static_min

fprintf('Lower bound on net worth in Productive Sector = %10.2f \n', Apmin);

amin =  Apmin + exp(-7);   % tricky right at 0, so move a bit to be able to solve it        
amax =  3*exp(nu)*exp(phip);      % bounds for assets (logs)

emin = 1;
emax = k;

smin = [amin, emin];
smax = [amax, emax];

sminp = smin;
smaxp = smax;
 
n = [101, k];

curv  = .1;                        % the close to 0, the closer to log-scale
agrid =  nodeunif(n(1), (amin-Apmin).^curv, (amax-Apmin).^curv).^(1/curv) + Apmin;             % make grid have more nodes close to lower bound, then upper bound

fspace = fundef({'spli', agrid, 0, spliorder},...
                {'spli', (1:1:k)', 0, 1});
           
grid = funnode(fspace);                  % 4d grid of collocation nodes 
s    = gridmake(grid);                   % collection of  states
ns   = length(s);

solvep_static;

paip = (1-eta)*eta^(eta/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*...
        W^(-alpha*eta/(1-eta))*(r+delta)^(-(1-alpha)*eta/(1-eta))*exp(phip);


v  =  zeros(ns,1);
c  =  funfitxy(fspace,s,[v,v,v]);         % guess for coefficients: almost never works, use homotopy (start with flex price,

%load cp; c = cp;

% Step 1: Solve productive modern producer problem

for i = 1 : 2                             % do a couple to give better starting values for newton

cnew = c;

[v1, ~, ~, x]     = saveBelmaxp(cnew, fspace, s);
c(:,1) = funfitxy(fspace, s, v1);         % update coefficients by projecting values on function space

v2     = valfunc2p(c, fspace, s);
c(:,2) = funfitxy(fspace, s, v2);

v3     = valfunc3p(c, fspace, s, x);      % price of claim to entire future stream of dividends
c(:,3) = funfitxy(fspace, s, v3);

fprintf('%4i %6.2e \n',[i, norm((c-cnew))/norm(c)]);    

if norm((c-cnew))/norm(c)< 1e-5 , break, end

end

c = vec(c);

for i = 1 : 50       %newton iterations
    
cnew = c;   

[bel, beljac] = solvebelp(cnew,fspace,s);        % RHS - LHS + derivatives

c  = cnew - (beljac\bel);                        % Newton step

fprintf('%4i %6.2e \n',[i, norm((c-cnew))/norm(c)]);    
if norm((c-cnew))/norm(c) < 1e-7 , break, end

end

c = [c(1:ns,:), c(ns+1:2*ns,:), c(2*ns+1:3*ns,:)];

[~, ~, ~, x] = saveBelmaxp(c,fspace,s); 
cx = funfitxy(fspace,s,x);   

fspacep = fspace;  % different state-space
cp      = c;       % value functions and the third is the price of a claim to dividends
cxp     = cx;      % savings choice
crp     = cr;      % multipliers

save cp cp

clear bel beljac c cx v1 v2 x

%%%%%%%%%%%%%%%%%%%% START Unproductive Modern %%%%%%%%%%%%%%%%%%%%%

fprintf('Solve Unproductive Sector Problem...');

fprintf('\n');
fprintf('\n');

solveu_static_min

fprintf('Lower bound on net worth in Unproductive Sector = %10.2f \n', Aumin);

amin =  Aumin + exp(-7) + bunp;     % tricky right at 0, so move a bit to be able to solve it        
amax =  3*exp(nu)*exp(phiu);      % bounds for assets (logs)

emin = 1;
emax = k;

smin = [amin, emin];
smax = [amax, emax];

sminu = smin;
smaxu = smax;
 
n = [101, k];

curv  = .5;    % the closet to 0, the closer to log-scale
if ~isempty(curvs)
    curv = curvs; 
end

agrid =  nodeunif(n(1), (amin-Aumin).^curv, (amax-Aumin).^curv).^(1/curv) + Aumin;             % make grid have more nodes close to lower bound, then upper bound

fspace = fundef({'spli', agrid, 0, spliorder},...
                {'spli', (1:1:k)', 0, 1});
           
grid = funnode(fspace);            % 4d grid of collocation nodes 
s    = gridmake(grid);             % collection of  states
ns   = length(s);

solveu_static;

v  =  zeros(ns,1);
c  =  funfitxy(fspace,s,[v,v,v]);         % guess for coefficients: almost never works, use homotopy (start with flex price,

%load cu; c = cu;

% Step 1: Solve unproductive modern producer problem

for i = 1 : 2                           % do a couple to give better starting values for newton

cnew = c;

[v1, ~, ~, x]     = saveBelmaxu(cnew, fspace, s);
c(:,1) = funfitxy(fspace, s, v1);         % update coefficients by projecting values on function space

v2     = valfunc2u(c, fspace, s);
c(:,2) = funfitxy(fspace, s, v2);

v3     = valfunc3u(c, fspace, s, x);
c(:,3) = funfitxy(fspace, s, v3);

fprintf('%4i %6.2e \n',[i, norm((c-cnew))/norm(c)]);    

if norm((c-cnew))/norm(c)< 1e-5 , break, end

end

c = vec(c);

for i = 1 : 50       %newton iterations
    
cnew = c;   

[bel, beljac] = solvebelu(cnew,fspace,s);        % RHS - LHS + derivatives

c  = cnew - (beljac\bel);                        % Newton step

fprintf('%4i %6.2e \n',[i, norm((c-cnew))/norm(c)]);    
if norm((c-cnew))/norm(c) < 1e-7 , break, end

end

c = [c(1:ns,:), c(ns+1:2*ns,:), c(2*ns+1:3*ns,:)];

[~, ~, ~, x] = saveBelmaxu(c,fspace,s); 
cx = funfitxy(fspace,s,x);   

fspaceu = fspace;  % different state-space
cu      = c;       % value functions
cxu     = cx;      % savings choice
cru     = cr;      % multipliers

[~, vup, vuu] = valfunc2u(c,fspace,s);

csu     = funfitxy(fspace,s,[vup,vuu]);      % expected continuation value for workers of switching

save cu cu

clear bel beljac c cx v1 v2 x


%%%%%%%%%%%%%%%%%%%% START Traditional %%%%%%%%%%%%%%%%%%%%%

fprintf('\n');

fprintf('Solve Traditional Sector Problem...');

fprintf('\n');
 
Amin = 0;

amin =  Amin + exp(-15) + 0.00;
amax =  exp(nu);          % bounds for assets (logs)

emin = 1;
emax = k;

smin  = [amin, emin];
smax  = [amax, emax];
smint = smin;
smaxt = smax;

n = [101, k];

curv  = .1;           % the close to 0, the closer to log-scale
agrid =  nodeunif(n(1), (amin - Amin).^curv, (amax-Amin).^curv).^(1/curv) + Amin;             % make grid have more nodes close to lower bound, then upper bound

fspace = fundef({'spli', agrid, 0, spliorder},...
                {'spli', (1:1:k)', 0, 1});
            
grid = funnode(fspace);                 % 4d grid of collocation nodes 
s    = gridmake(grid);                  % collection of  states
ns   = length(s);

% Solve for producer's net worth after joining the modern sector:
% a' = x - kappa_u + theta*xai*pu*(a',e), where pu is coded in cu(:,3)

aigrid = nodeunif(1500, (amin - Amin).^curv, (amax-Amin).^curv).^(1/curv) + Amin; 
fspacei = fundef({'spli', aigrid, 0, 1},...
                 {'spli', (1:1:k)', 0, 1});
            
gridi = funnode(fspacei);                 % 4d grid of collocation nodes 
si    = gridmake(gridi);                  % collection of  states
solvenetworth;

cinj = funfitxy(fspacei,si,inject); 


v  =  zeros(ns,1);
c  =  funfitxy(fspace,s,[v,v]);         % guess for coefficients: almost never works, use homotopy (start with flex price,

%load ct; c = ct;

for i = 1 : 2                           % do a couple to give better starting values for newton

cnew = c;

v1     = saveBelmaxt(cnew,fspace,s);
c(:,1) = funfitxy(fspace,s,v1);         % update coefficients by projecting values on function space
v2     = valfunc2t(c,fspace,s);         % careful here: want to use worker's cont. value to update coefficient
c(:,2) = funfitxy(fspace,s,v2);

fprintf('%4i %6.2e \n',[i,norm((c-cnew)./c)]);    

end


c = vec(c);

for i = 1 : 50       %newton iterations
    
cnew = c;   

[bel, beljac] = solvebelt(cnew,fspace,s);        % RHS - LHS + derivatives

c  = cnew - (beljac\bel);                        % Newton step

fprintf('%4i %6.2e \n',[i,norm((c-cnew)./c)]);    
if norm((c-cnew)./c) < 1e-7 , break, end

end

c = [c(1:ns,:), c(ns+1:2*ns,:)];
clear bel beljac

[~, ~, x] = saveBelmaxt(c,fspace,s); 
cx    = funfitxy(fspace,s,x);   


fspacet = fspace;  % different state-space
ct      = c;
cxt     = cx;

[~, vtu, vtt] = valfunc2t(ct,fspacet,s);

cst     = funfitxy(fspacet,s,[vtu,vtt]);      % expected continuation value for workers of switching
clear fspace c

save ct ct

fprintf('\n');

disp('Done with Value Functions.')

clear e w ez wz s n bel beljac

clc
ergodic
simulate

break






