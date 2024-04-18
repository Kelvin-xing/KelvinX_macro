clear all;
clc;

global rho alpha eta beta delta W r smin smax fspaceu ...
       fspace cu cr egrid P kappau pait sminu ...
       theta mu lbar phiu

mu        = 1.08;              % growth rate
beta      = 0.92*mu;            

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

% Version 1: 

phiu   = 1/(1-eta)*0.20;    % productivity gap unproductive

kappau = 1.1930;            % fixed cost of joining u

theta  = 0.8648; 

rho    = 0.2457;            % persistence transitory productivity
se     = 0.4958;

W      = 1.042530270930890;
r      = 0.038232184088732;

% After Crisis

theta = 0.59;
W     = 0.9592;

findvarz = 1;

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

%%%%%%%%%%%%%%%%%%%% START Unproductive Modern %%%%%%%%%%%%%%%%%%%%%

fprintf('Solve Unproductive Sector Problem...');

fprintf('\n');
fprintf('\n');

solveu_static_min

fprintf('Lower bound on net worth in Unproductive Sector = %10.2f \n', Aumin);

amin =  Aumin + exp(-7);          % tricky right at 0, so move a bit to be able to solve it        
amax =  3*exp(nu)*exp(phiu);      % bounds for assets (logs)

emin = 1;
emax = k;

smin = [amin, emin];
smax = [amax, emax];

sminu = smin;
smaxu = smax;
 
n = [101, k];

curv   = .1;                       % the close to 0, the closer to log-scale
agrid  =  nodeunif(n(1), (amin-Aumin).^curv, (amax-Aumin).^curv).^(1/curv) + Aumin;             % make grid have more nodes close to lower bound, then upper bound

fspace = fundef({'spli', agrid, 0, spliorder},...
                {'spli', (1:1:k)', 0, 1});
           
grid = funnode(fspace);            % 4d grid of collocation nodes 
s    = gridmake(grid);             % collection of  states
ns   = length(s);

solveu_static;

v  =  zeros(ns,1);
c  =  funfitxy(fspace,s,[v,v]);         % guess for coefficients: almost never works, use homotopy (start with flex price,

%load cu; c = cu;

% Step 1: Solve unproductive modern producer problem

for i = 1 : 2                           % do a couple to give better starting values for newton

cnew = c;

[v1, ~, x]     = saveBelmaxu(cnew, fspace, s);
c(:,1) = funfitxy(fspace, s, v1);         % update coefficients by projecting values on function space

v2     = valfunc2u(c, fspace, s);
c(:,2) = funfitxy(fspace, s, v2);

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

c = [c(1:ns,:), c(ns+1:2*ns,:)];

[~, ~, x] = saveBelmaxu(c,fspace,s); 
cx = funfitxy(fspace,s,x);   

fspaceu = fspace;  % different state-space
cu      = c;       % value functions
cxu     = cx;      % savings choice
cru     = cr;      % multipliers
su      = s;

save cu cu


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRANSITIONS MODERN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T     = 35;
Wsave = W*ones(T,1);

%load Wsave           % first and last periods are end-points

cu_irf  = zeros(size(cu,1),  size(cu,2),  T);    cu_irf(:,:,T)  = cu;
cxu_irf = zeros(size(cxu,1), size(cxu,2), T);    cxu_irf(:,:,T) = cxu;
cru_irf = zeros(size(cru,1), size(cru,2), T);    cru_irf(:,:,T) = cru;


for t = T-1: -1: 1
    
disp(t)
    
W = Wsave(t,:);

% Modern sector

fspace  = fspaceu;
s       = su;
smin    = sminu; 
smax    = smaxu;

solveu_static;        % get new multipliers on BC

cru_irf(:,:,t) = cr;  % save multipliers, notice new cr is the global
 
cu = cu_irf(:,:,t+1); 
 
[v1, ~, x] = saveBelmaxu(cu, fspaceu, su); 

cu(:,1) =  funfitxy(fspaceu, su, v1);            % update to compute v2

v2  = valfunc2u(cu,fspaceu,su);

cu_irf(:,:,t) = funfitxy(fspaceu, su, [v1,v2]);  
cxu_irf(:,:,t) = funfitxy(fspaceu, su, x);      


end


%%%%%%%%%%%%%%%%%%%% START Traditional %%%%%%%%%%%%%%%%%%%%%

fprintf('\n');

fprintf('Solve Traditional Sector Problem...');

fprintf('\n');
 
Amin = 0;

amin =  Amin + exp(-15);
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
st      = s;

[~, vtu, vtt] = valfunc2t(ct,fspacet,s);

cst     = funfitxy(fspacet,s,[vtu,vtt]);      % expected continuation value for workers of switching
clear fspace c

save ct ct

fprintf('\n');

disp('Done with Value Functions.')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRANSITIONS TRADITIONAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T     = 35;
Wsave = W*ones(T,1);

load Wsave           % first and last periods are end-points


ct_irf  = zeros(size(ct,1),  size(ct,2),  T);    ct_irf(:,:,T)  = ct;
cxt_irf = zeros(size(cxt,1), size(cxt,2), T);    cxt_irf(:,:,T) = cxt;
cst_irf = zeros(size(cst,1), size(cst,2), T);    cst_irf(:,:,T) = cst;

for t = T-1: -1: 1
      
disp(t)
    
W = Wsave(t,:);

% Traditional sector 

%(notice this needs cu(t+1) to compute continuation value to switch so update

cu = cu_irf(:,:,t+1); 
ct = ct_irf(:,:,t+1);

fspace = fspacet;
smin   = smint;
smax   = smaxt;

[v1, ~, x] = saveBelmaxt(ct,fspacet,st); 
cxt_irf(:,:,t)    = funfitxy(fspacet,st,x);  

ct(:,1) = funfitxy(fspacet, st, v1);

[v2, vtu, vtt] = valfunc2t(ct,fspacet,s);

ct_irf(:,:,t) = funfitxy(fspacet, st, [v1,v2]);  
cst_irf(:,:,t) = funfitxy(fspacet, st, [vtu,vtt]);

end


ergodic_transitions

simulate_transitions


