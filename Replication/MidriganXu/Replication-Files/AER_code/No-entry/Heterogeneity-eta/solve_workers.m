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