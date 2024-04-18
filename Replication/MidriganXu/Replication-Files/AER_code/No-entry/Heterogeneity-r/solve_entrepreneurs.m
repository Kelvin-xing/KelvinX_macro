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
amax =  2*exp(nu);      % bounds for assets

emin = 1;
emax = k;

smin = [amin, emin];
smax = [amax, emax];

n = [101, k];

curv  = .1;                        % the close to 0, the closer to log-scale
agrid =  nodeunif(n(1), amin.^curv, amax.^curv).^(1/curv);             % make grid have more nodes close to lower bound, then upper bound

fspace = fundef({'spli', agrid, 0, spliorder},...
                {'spli', (1:1:k)', 0, 1});
           
grid = funnode(fspace);                   % 4d grid of collocation nodes 
s    = gridmake(grid);                    % collection of  states
ns   = length(s);

solve_static;

v  =  zeros(ns,1);
c  =  funfitxy(fspace,s,[v,v]);           % guess for coefficients: almost never works, use homotopy (start with flex price,

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
