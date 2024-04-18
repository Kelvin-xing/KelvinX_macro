function y = objective_equilibrium(x)

global rho alpha eta beta delta W r smin smax fspace cr egrid P theta lbar Pw


beta      = 0.92;            

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

DtoY      = 1.2;
EtoY      = 0.3;

lambdae  = 0.785;            % probability worker stays in e
lambdau  = 0.5;              % 1

% Version 1: 

theta  = 0.5685; 
rho    = 0.3000;               % persistence transitory productivity
se     = 0.8310;

W      = x(1);
r      = x(2);

findvarz = 1;

% WORKERS

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


% ENTREPRENEURS: 

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

pai = (1-eta)*eta^(eta/(1-eta))*alpha^(alpha*eta/(1-eta))*(1-alpha)^((1-alpha)*eta/(1-eta))*...
        W^(-alpha*eta/(1-eta))*(r+delta)^(-(1-alpha)*eta/(1-eta));

amin =  exp(-7);        % tricky right at 0, so move a bit to be able to solve it        
amax =  2*exp(nu);      % bounds for assets (logs)

emin = 1;
emax = k;

smin = [amin, emin];
smax = [amax, emax];

n = [101, k];

curv  = .1;                        % the close to 0, the closer to log-scale
agrid =  nodeunif(n(1), amin.^curv, amax.^curv).^(1/curv);             % make grid have more nodes close to lower bound, then upper bound

fspace = fundef({'spli', agrid, 0, spliorder},...
                {'spli', (1:1:k)', 0, 1});
           
grid = funnode(fspace);                  % 4d grid of collocation nodes 
s    = gridmake(grid);                   % collection of  states
ns   = length(s);

solve_static;

v  =  zeros(ns,1);
c  =  funfitxy(fspace,s,[v,v]);         % guess for coefficients: almost never works, use homotopy (start with flex price,

%load c c

for i = 1 : 2                             % do a couple to give better starting values for newton

cnew = c;

[v1, ~, x]     = saveBelmax(cnew, fspace, s);
c(:,1) = funfitxy(fspace, s, v1);         % update coefficients by projecting values on function space

v2     = valfunc2(c, fspace, s);
c(:,2) = funfitxy(fspace, s, v2);

fprintf('%4i %6.2e \n',[i, norm((c-cnew))/norm(c)]);    

if norm((c-cnew))/norm(c)< 1e-5 , break, end

end

c = vec(c);

for i = 1 : 50       %newton iterations
    
cnew = c;   

[bel, beljac] = solvebel(cnew,fspace,s);        % RHS - LHS + derivatives

c  = cnew - (beljac\bel);                        % Newton step

fprintf('%4i %6.2e \n',[i, norm((c-cnew))/norm(c)]);    
if norm((c-cnew))/norm(c) < 1e-7 , break, end

end

c = [c(1:ns,:), c(ns+1:2*ns,:)];

[~, ~, x] = saveBelmax(c,fspace,s); 
cx = funfitxy(fspace,s,x);   

save c c

clear bel beljac v1 v2 x

clc

y = zeros(2,1); 

ergodic1;

y(1) = (Ka -  Aa - Aw)./Ka;
y(2) = (La - Lbar)./Lbar;

weight = ones(2,1); 

y = sqrt(sum(y.^2.*weight)/sum(weight));

disp(y); 