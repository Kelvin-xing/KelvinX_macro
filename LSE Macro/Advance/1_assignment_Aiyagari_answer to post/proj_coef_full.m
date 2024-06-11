function coef = proj_coef_full(r,par)

% Parameters
%--------------------------------------------------------------------------

par.r     = r;
par.kd    = (par.r/par.alpha)^(1/(par.alpha-1));
par.wss   = (1-par.alpha)*par.kd^par.alpha;
par.ass   = par.kd*(1+par.r-par.delta) + par.wss;

% Grid
%--------------------------------------------------------------------------

a_low  = 0.75*(par.ass);  
a_high = 1.25*(par.ass); 
siz.a  = 300;
% a_step = (a_high-a_low)/(siz.a-1);   % equidistant nodes
% grid.a = (a_low:a_step:a_high)';
cna = chebnode(siz.a);                  % Chebyshev nodes
grid.a = (a_high-a_low)/2*(cna+1)+a_low; 

% GH nodes
%--------------------------------------------------------------------------

siz.e = 10; 
[gh.e,gh.w] = hernodes(siz.e);

% Projection
%--------------------------------------------------------------------------

% initial values for coefficients of polynomial
% init = zeros(order+1,1); 

if par.order == 1
    
    init = [1.581429 - 0.0228279912193402*par.ass; 0.0228279912193402];
else
    init = [1.581429 - 0.0228279912193402*par.ass; 0.0228279912193402; ...
        zeros(par.order-1,1)];

end

%Minimization routine

options = optimset('Display','Iter','MaxFunEvals',1E5,'MaxIter',...
          1E4,'TolFun',1E-10,'TolX',1E-10);%,'UseParallel','Always');
coef    = fminsearch(@(coef) projerrs_full(coef,par,grid,gh,siz),init,options);


