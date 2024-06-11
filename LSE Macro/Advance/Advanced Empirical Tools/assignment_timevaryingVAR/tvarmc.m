%%                      Time-varying VAR
%
%                     Monte Carlo exercise
%
%         Structure of this script: just a for loop around tvar.m
%
%       Interesting to play around with constant-of-proportionality
%
%                         Toy example:
%
% y(t) = const(t) + rho1(t) * y(t-1) + rho2(t) * y(t-2) + ... + epsilon(t)
%
%   - no stochastic volatility
%   - parameters assumed to follow random walk
%   - general number of lags
%
%       Joris de Wind (June 2014)
%--------------------------------------------------------------------------

%Cleaning
clear all; close all; clc

%Random numbers
rng(25052014,'twister')   % seed is the date Joris' son Magnus was born s

%% Monte Carlo exercise
%--------------------------------------------------------------------------

%Number of Monte Carlo repititions
M = 50;

%Allocating memory
mcTH = zeros(2,161,M);          %Dimensions hard-coded ...

%Main Monte Carlo for loop
progressbar(0);
for mc = 1:M,

    %Progressbar
    progressbar(mc/M);

%% Genearating some data
%--------------------------------------------------------------------------

% - data generating process used below is different from time-varying VAR
%   model specified above
% - feel free to try out other data generating processes or real data

%Number of time periods
T = 200;

%Trajectory AR(1) coefficient
rho1 = 0.9*ones(T,1);
rho1(floor(0.4*T):T,1) = 0.5;
rho1(floor(0.4*T)+(0:50),1) = 0.9 + (0.5 - 0.9) * 1./(1+exp(-0.5*(-25:25)));    %gradual change according to logistic function

%Trajectory AR(2) coefficient
% rho2 = 0.1*ones(T,1);
% rho2(floor(0.7*T):T,1) = 0.05;
rho2 = zeros(T,1);

%Trajectory constant
% cons = 0.3*ones(T,1);
% cons(floor(0.5*T):T,1) = 0.5;
% cons(floor(0.5*T)+(0:50),1) = 0.3 + (0.5 - 0.3) * 1./(1+exp(-0.5*(-25:25)));  %gradual change according to logistic function
cons = 1-rho1;	%so as to keep the unconditional mean constant

%Variance of shocks
sigma2 = 0.0005;

%Data generating process
data = zeros(T,1);
data(1,1) = cons(1)/(1-rho1(1)-rho2(1));        %start at local approximation of unconditional mean
data(2,1) = cons(1)/(1-rho1(1)-rho2(1));
for t = 3:T,
    data(t,1) = cons(t) + rho1(t) * data(t-1,1) + rho2(t) * data(t-2,1) + ...
        sqrt(sigma2)*randn;
end

%Scaling data
scaling = 1;          %scaling matters for the relative amount of time variation
                      %in the intercept versus the autoregressive coefficients
data = data/scaling;
cons = cons/scaling;

%% Preparing
%--------------------------------------------------------------------------

%Number of lags
p = 1;                                 %feel free to try out other settings

%Getting the right dimensions
T = length(data)-p;
k = 1+p;                               %number of time-varying parameters

%Notation as in Cogley and Sargent (2001) and De Wind and Gambetti (2014):
%                   y(t) = X(t)' * theta(t) + u(t)
y = data(1+p:end,:);                                    % T x 1 vector
X = zeros(k,T);                                         % k x T matrix
for t = 1:T,
    X(:,t) = [1,data(t+p-1:-1:t,:)'];
end

%Splitting the sample: training versus main
np = floor(0.2*T);   T = T - np;
yp = y(1:np,:);      y(1:np,:) = [];
Xp = X(:,1:np);      X(:,1:np) = [];

%% Priors
%--------------------------------------------------------------------------

% - checking robustness to different constant-of-proportionality crucial!!

%OLS on training sample
thetahat    = Xp'\yp;
err         = yp-Xp'*thetahat;
Sigma       = err'*err/np;
varthetahat = Sigma*inv(Xp*Xp');                                 %#ok<MINV>

%Prior TH
THbar      = thetahat;          %mean normal prior TH(0|0)
PeK        = zeros(k,k,T+1);    %allocating memory
PeK(:,:,1) = 4*varthetahat;     %variance normal prior TH(0|0)

%Prior Q
gam2 = 0.01;                    %constant-of-proportionality, crucial to check for robustness here ...
Qbar = gam2*varthetahat;        %scale matrix inverse-Wishart prior Q
T0Q  = k+1;                     %degrees-of-freedom inverse-Wishart prior Q

%Prior R
Rbar = Sigma;                   %scale matrix inverse-Wishart prior R
T0R  = 2;                       %degrees-of-freedom inverse-Wishart prior R

%% Bayesian estimation
%--------------------------------------------------------------------------

%Settings
g = 250;                                        %Length of Markov chain
b = 0.5;                                        %Burn-in period fraction

%Allocating memory and initialization
TH = zeros(k,T+1,g);
Q1 = zeros(k,k,g+1);        Q1(:,:,1) = Qbar;   %Actually better to randomize
R1 = zeros(1,1,g+1);        R1(:,:,1) = Rbar;   %Actually better to randomize

%Gibbs sampler
%progressbar(0);
for i = 1:g,
    
    %Progressbar
    %progressbar(i/g);
    
    %Sample TH
    TH(:,1,i) = THbar;
    TH(:,:,i) = sampleTH(y,X,R1(:,:,i),Q1(:,:,i),...
        T,PeK,zeros(k,k,T),TH(:,:,i),randn(k,T+1));
    
    %Sample Q1
    v = TH(:,2:end,i)-TH(:,1:end-1,i);
    Qtemp = v*v'+T0Q*Qbar; Qtemp = (Qtemp+Qtemp')/2;
    Q1(:,:,i+1) = iwishrnd(Qtemp,T+T0Q);
    
    %Sample R1
    u  = createyhat(y,X,TH(:,2:end,i),T);
    Rtemp = u'*u+T0R*Rbar; Rtemp = (Rtemp+Rtemp')/2;
    R1(:,:,i+1) = iwishrnd(Rtemp,T+T0R);
    
end

%Median
results.TH = median(TH(:,:,b*g+1:g),3);
results.Q1 = median(Q1(:,:,b*g+1:g),3);
results.R1 = median(R1(:,:,b*g+1:g),3);

    %Monte Carlo results
    mcTH(:,:,mc) = results.TH;

end

%Monte Carlo median
results.mcTH = median(mcTH,3);

%% Plotting results
%--------------------------------------------------------------------------

%Median trajectories of time-varying coefficients
figure('name','trajectories - Monte Carlo exericse')
subplot(2,1,1)
plot(cons,'LineWidth',1.5), hold all, plot(p+np+(0:T),results.mcTH(1,:),'LineWidth',1.5)
title('\bf{constant}'), l = legend('true','estimated'); set(l,'box','off','location','northwest')
subplot(2,1,2)
plot(rho1,'LineWidth',1.5), hold all, plot(p+np+(0:T),results.mcTH(2,:),'LineWidth',1.5)
title('\bf{AR(1) coeff}'), l = legend('true','estimated'); set(l,'box','off','location','northwest')
maximize