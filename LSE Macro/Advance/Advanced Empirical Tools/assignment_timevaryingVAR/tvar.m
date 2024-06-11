%%                      Time-varying VAR
%
%                         univariate example:
%
% y(t) = const(t) + rho1(t) * y(t-1) + rho2(t) * y(t-2) + ... + epsilon(t)
%
%   - no stochastic volatility
%   - parameters assumed to follow random walk
%   - general number of lags (although figures that compare true
%     time-varying parameters with estimates only make sense if p <= 2)
%
%   Joris de Wind (June 2014) plus some comments and modifications by Wouter
%--------------------------------------------------------------------------

%Cleaning
clear all; close all; clc

%Random numbers

rng(25052014,'twister')   % seed is the date Joris' son Magnus was born s
load shocks

%If rng(number) gives an error message, then you are using an older Matlab
%in this case you can use randn('state',25052014). 

%The results are quite senstive to the actual time series used, i.e. the vector data generated below. 
%The results are also somewhat sensitive to the shocks used in the Gibbs sampler, but less so.
%
%To make sure everybody gets somewhat similar results, we all use the same
%shocks to generate the actual observations. These are stored in shocks.mat
%


%Number of time periods
T = 200;

%% Section 1: Genearating some data 
%--------------------------------------------------------------------------

% This section specifies how the coefficients of the dgp, const(t), rho1(t), and rho2(t)
% vary over time. Do not read this section until you got the program running. Just treat
% it as something unknown, which is the way it would be in an actual empirical application.

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
scaling = 1;          %scaling matters for the relative amount of time variation in the 
                      %intercept versus time variation in the autoregressive coefficients
                      %with scaling = 1 there is obviously no scaling
data = data/scaling;
cons = cons/scaling;

%% Preparing
%--------------------------------------------------------------------------

%Number of lags
p = 2;                                 %graphs below rely on p<=2 but algorithm allows for higher values for p

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

%% Section 2: Priors
%--------------------------------------------------------------------------

%Splitting the sample: training (to form prior) versus main (to do the analysis)
np = floor(0.2*T);   T = T - np;
yp = y(1:np,:);      y(1:np,:) = [];
Xp = X(:,1:np);      X(:,1:np) = [];


%  results often depend on how you set priors, so it is important that you check
%  whether your results are robust to reasonable changes in the prior
%  here we use OLS estimates obtained with the training sample to form prior

thetahat    = XXX;              %calculate OLS estimate                 
err         = XXX;              %calculate residual          
Sigma       = XXX;              %calculate SSR
varthetahat = XXX;              %calculate variance-covariance matrix of estimate

%Prior TH
THbar      = thetahat;          %mean of the prior for theta, i.e. TH(0|0)
PeK        = zeros(k,k,T+1);    %allocating memory for the variance of theta
PeK(:,:,1) = 4*varthetahat;     %variance of the prior for theta TH(0|0)
                                %the number 4 is pretty arbitrary

%Prior Q
gam2 = 0.01;                    %constant-of-proportionality, crucial to check for robustness here ...
Qbar = gam2*varthetahat;        %scale matrix inverse-Wishart prior Q
T0Q  = k+1;                     %degrees-of-freedom inverse-Wishart prior Q
%(relative to the original version of this program, I have changed the
%notation such that it is consistent with the slides, Qbar is now like Vbar
%in the slides)

% About TOQ: The DGFs are set equal to “dimension + 1” which is the minimal degrees-of-freedom to have a finite first moment.
% So with this minimal prior information is used. For bigger models you may
% want to use a less diffuse prior
% About gam2: varthetahat says something about sampling uncertainty. This
% would be present even if theta was not time-varying. How estimated sampling uncertainty 
% (with a regression model assuming constant parameters) is a tough one and
% some judgement is required here.

%Prior R
T0R  = 2;                       %degrees-of-freedom inverse-Wishart prior R (See comments on prior above)
Rbar = TOR*Sigma;               %scale matrix inverse-Wishart prior R

%% Section 3: Bayesian estimation
%--------------------------------------------------------------------------

%Settings
g = 2000;                                       %Length of Markov chain
b = 0.5;                                        %Burn-in period fraction

%Allocating memory and initialization
TH = zeros(k,T+1,g);
Q1 = zeros(k,k,g+1);        Q1(:,:,1) = gam2*varthetahat;   %Actually better to randomize
R1 = zeros(1,1,g+1);        R1(:,:,1) = Sigma;   %Actually better to randomize
%Q1 & R1 are NOT priors. They are simply initial conditions of the Gibbs sampler
%Note that Q1 & R1 are averages of sum of squares

%Gibbs sampler
progressbar(0);   %progressbar is an included program that gives a time estimate of the remaining time in loop
for i = 1:g
    
    %Progressbar
    progressbar(i/g);
    
    %Sample TH step I of the Gibbs sampler   
    TH(:,1,i) = THbar;
    TH(:,:,i) = sampleTH(y,X,R1(:,:,i),Q1(:,:,i),...
        T,PeK,zeros(k,k,T),TH(:,:,i),randn(k,T+1));
    %note that TH(:,t,i) goes from t=1 to t=T+1
    %TH(:,t,i) is the estimate of theta_(t-1) conditional on info (more detailed info on this in function sampleTH)
    %TH(:,1,i) is theta_0 (here set equal to the mean of prior, but you can
    %also take a draw from prior). Either way we as econometricians do not
    %know this value
    %TH(:,T+1,i) corresponds to the last theta_t about which we can say something, i.e., theta_T

    %Step II of the Gibbs sampler
    %combine prior with information from the data (note that prior remains the same for different i
    
    %Sample Q1
    v = TH(:,2:end,i)-TH(:,1:end-1,i); %calculates prediction error for thetat
    Qtemp = XXX; %calculate scalings matrix posterior
    Q1(:,:,i+1) = iwishrnd(XXX,XXX); %iwishrnd(A,I) draws from inverted Wishart with scalings matrix A and I dgf
    % whereas Qtemp is like a sum of squares, Q1 is like an average sum of squares 
    
    %Sample R1
    u  = createyhat(y,X,TH(:,2:end,i),T); % calculate prediction error for data
    Rtemp = XXX;  %calculate scalings matrix posterior
    R1(:,:,i+1) = iwishrnd(XXX,XXX); 
    % whereas Rtemp is like a sum of squares, R1 is like an average sum of squares 
    
end

%Median
results.TH = median(TH(:,:,b*g+1:g),3);
results.Q1 = median(Q1(:,:,b*g+1:g),3);
results.R1 = median(R1(:,:,b*g+1:g),3);


%% Plotting results
%--------------------------------------------------------------------------

%Median trajectories of time-varying coefficients
figure('name','trajectories')
subplot(2,2,1)
plot(cons,'LineWidth',1.5), hold all, plot(p+np+(0:T),results.TH(1,:),'LineWidth',1.5)
title('\bf{constant}'), l = legend('true','estimated','Location','Northwest'); set(l,'box','off')
subplot(2,2,2)
plot(rho1,'LineWidth',1.5), hold all, plot(p+np+(0:T),results.TH(2,:),'LineWidth',1.5)
title('\bf{AR(1) coeff}'), l = legend('true','estimated','Location','Northeast'); set(l,'box','off')
if p == 2,
    subplot(2,2,3)
    plot(rho2,'LineWidth',1.5), hold all, plot(p+np+(0:T),results.TH(3,:),'LineWidth',1.5)
    title('\bf{AR(2) coeff}'), l = legend('true','estimated','Location','Northwest'); set(l,'box','off')
end
subplot(2,2,4)
pl = plot(data,'LineWidth',1.5,'Color','red');
title('\bf{data}')
maximize

%Posterior distribution of trace(Q)
Qtrace_prior     = trace(gam2*varthetahat);
Qtrace_posterior = zeros((1-b)*g,1);
for i = 1:k,
    Qtrace_posterior = Qtrace_posterior + squeeze(Q1(i,i,b*g+1:g));  
    %squeeze removes singletons, e.g., it would turn a N1x1xN3 object into a N1xN3 object
    %in this case it turns this 3 dimensional object into a one-dimensional object
end
figure('name','Qtrace')
hist(Qtrace_posterior,50), title('\bf{Posterior trace Q}')
patch('XData',Qtrace_prior+diff(get(gca,'XLim'))*[-0.5 0.5 0.5 -0.5]/100,'YData',...
reshape(repmat(get(gca,'YLim').*[1 0.98],2,1),1,4),'FaceColor','r'), box off
maximize



%% Plot percentiles

All.TH1     = squeeze(TH(1,:,b*g+1:g));     %squeeze gets rid of dimension with singleton
All.TH2     = squeeze(TH(2,:,b*g+1:g));
All.TH3     = squeeze(TH(3,:,b*g+1:g));
All.TH1     = sort(All.TH1,2);
All.TH2     = sort(All.TH2,2);
All.TH3     = sort(All.TH3,2);
[S1 S2]     = size(All.TH1);
p5          = floor(0.05*S2);
p10         = floor(0.10*S2);
p90         = floor(0.90*S2);
p95         = floor(0.95*S2);

%Percentiles
figure('name','percentiles')
subplot(2,2,1)
plot(cons,'LineWidth',1.5),hold all, plot(p+np+(0:T),All.TH1(:,p5),'LineWidth',1.5)
hold all, plot(p+np+(0:T),All.TH1(:,p10),'LineWidth',1.5)
hold all, plot(p+np+(0:T),All.TH1(:,p90),'LineWidth',1.5)
hold all, plot(p+np+(0:T),All.TH1(:,p95),'LineWidth',1.5)
title('\bf{constant}'), l = legend('true','5% percentile','10% percentile','90% percentile','95% percentile','Location','Best'); set(l,'box','off')
subplot(2,2,2)
plot(rho1,'LineWidth',1.5), hold all, plot(p+np+(0:T),All.TH2(:,5),'LineWidth',1.5)
hold all, plot(p+np+(0:T),All.TH2(:,p10),'LineWidth',1.5)
hold all, plot(p+np+(0:T),All.TH2(:,p90),'LineWidth',1.5)
hold all, plot(p+np+(0:T),All.TH2(:,p95),'LineWidth',1.5)
title('\bf{AR(1) coeff}'), l = legend('true','5% percentile','10% percentile','90% percentile','95% percentile','Location','Best'); set(l,'box','off')
if p == 2,
    subplot(2,2,3)
    plot(rho2,'LineWidth',1.5), hold all, plot(p+np+(0:T),All.TH3(:,5),'LineWidth',1.5)
hold all, plot(p+np+(0:T),All.TH3(:,p10),'LineWidth',1.5)
hold all, plot(p+np+(0:T),All.TH3(:,p90),'LineWidth',1.5)
hold all, plot(p+np+(0:T),All.TH3(:,p95),'LineWidth',1.5)
title('\bf{AR(2) coeff}'), l = legend('true','5% percentile','10% percentile','90% percentile','95% percentile','Location','Best'); set(l,'box','off')
end
subplot(2,2,4)
pl = plot(data,'LineWidth',1.5,'Color','red');
title('\bf{data}')
maximize
