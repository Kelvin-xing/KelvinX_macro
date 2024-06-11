function coef = proj_coef(r,par)

% Parameters
%--------------------------------------------------------------------------

par.r     = r;
par.kd    = (par.r/par.alpha)^(1/(par.alpha-1));
par.wss   = (1-par.alpha)*par.kd^par.alpha;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Below specify the steady state value for wealth
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
par.ass   =;

% Grid
%--------------------------------------------------------------------------

%!!! note that there is a one dimensional grid for wealth instead of one
%for income and one for savings. This is ok because one parameter value
%takes on a particular choice. Which parameter is that? Make sure you
%understand why this is ok.
%
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
coef    = fminsearch(@(coef) projerrs(coef,par,grid,gh,siz),init,options);