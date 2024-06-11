%**************************************************************************
%                   LSE Macroeconomics Summer Program
%                   Part II: Heterogeneous Agents
%                   Instructor: Wouter J. Den Haan
%
%                   use of this program in any fee-based program requires
%                   explicit permission (wjdenhaan@gmail.com)
%**************************************************************************
%
%==========================================================================
% This program takes the values from perturbation and solves the partial
% equilibrium problem using projection methods
%--------------------------------------------------------------------------

clc
clear

% loading parameters and results from perturbation
load par
load r_ge
load ks_ge

ks_partial = zeros(1,2);

par.order = 4;
sigs = [0.001;0.3];
par.T = 100000;

%--------------------------------------------------------------------------
% 2. generate innovations
%--------------------------------------------------------------------------

randn('seed',20100807)
innovations = randn(par.T,1);

r = r_ge(1);

for i = 1:2

% i = 1: we use r=r_ge(1), i.e., the low volatility equilibrium interest
%        rate, so we simply resolve at equilibrium value
% i = 2: we increase uncertainty but do not adjust r, i.e. partial eq.

    par.sigshock = sigs(i);
    shocks   = par.sigshock*innovations;

    coef = proj_coef_full(r,par);
    k = zeros(par.T,1);
    a = zeros(par.T,1);
    c = zeros(par.T,1);
    k(1) = (r/par.alpha)^(1/(par.alpha-1));
    w = (1-par.alpha)*k(1)^par.alpha;
    a(1) = k(1)*(1+r-par.delta) + w;

    for t = 2:par.T
        c(t-1) = polyn(a(t-1),coef,par.order);
        k(t) = a(t-1) - c(t-1);
        a(t) = (1+r-par.delta)*k(t) + w*(1+shocks(t));
    end

    ks_partial(i) = mean(k(par.T0:end));          % average capital supply (from HH problem)

end

disp(' ')
disp('Partial equilibrium capital supply')
disp([ks_partial(1), 100*(ks_partial(2)-ks_partial(1))/ks_partial(1)])
disp('General equilibrium capital supply')
disp([ks_ge(1), 100*(ks_ge(2)-ks_ge(1))/ks_ge(1)])