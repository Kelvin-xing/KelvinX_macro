%**************************************************************************
%                   LSE Macroeconomics Summer Program
%                   Part II: Heterogeneous Agents
%                   Instructor: Wouter J. Den Haan
%
%                   use of this program in any fee-based program requires
%                   explicit permission (wjdenhaan@gmail.com)
%**************************************************************************
%
%--------------------------------------------------------------------------
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

par.order = 2;
sigs = [0.001;0.3];

%--------------------------------------------------------------------------
% 2. generate innovations
%--------------------------------------------------------------------------

rng(20100807,'philox') % set the seed of the random number generator
innovations = randn(par.T,1);

r = r_ge(1);

for i = 1:2

% i = 1: we use r=r_ge(1), i.e., the low volatility equilibrium interest
%        rate, so we simply resolve at equilibrium value
% i = 2: we increase uncertainty but do not adjust r, i.e. partial eq.

    par.sigshock = sigs(i);
    shocks   = par.sigshock*innovations;

    coef = projmain(r,par);
    k = zeros(par.T,1);
    a = zeros(par.T,1);
    c = zeros(par.T,1);
    k(1) = (r/par.alpha)^(1/(par.alpha-1));
    w = (1-par.alpha)*k(1)^par.alpha;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Below specify the first element of a to be the steady state value of
% wealth and simulate the economy given your solution from the projection.
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    a(1) = ;

    for t = 2:par.T

    end

    ks_partial(i) = mean(k(par.T0:end));          % average capital supply (from HH problem)

end

disp(' ')
disp('Partial equilibrium capital supply')
disp([ks_partial(1), 100*(ks_partial(2)-ks_partial(1))/ks_partial(1)])
disp('General equilibrium capital supply')
disp([ks_ge(1), 100*(ks_ge(2)-ks_ge(1))/ks_ge(1)])