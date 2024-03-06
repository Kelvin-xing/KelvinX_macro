%**************************************************************************
%                   LSE Macroeconomics Summer Program
%                   Part I: The Essentials
%                   Instructor: Wouter J. Den Haan
%
%                   use of this program in any fee-based program requires
%                   explicit permission (wjdenhaan@gmail.com)
%**************************************************************************


%==========================================================================
%           Simple matching model and business cycle statistics
%==========================================================================


%--------------------------------------------------------------------------
% This m-file sets values for entrepreneur share and then runs a loop and
% does the following:
%
%   1. Runs Dynare to obtain policy functions
%   2. Simulates the economy given the policy functions
%   3. HP-filters the data and computes some labour market statistics
%--------------------------------------------------------------------------

close all
clear all
clc

%--------------------------------------------------------------------------
% 0. Setting values for the entrepreneur share 
%--------------------------------------------------------------------------

ent_share = [0.25,0.1,0.05,0.025,0.01];
I = 5;

% allocating memory
s_ny    = zeros(1,I);

%--------------------------------------------------------------------------
% 0.1. Starting the loop over entrepreneur share values
%--------------------------------------------------------------------------

for i = 1:I

    % specifying entrepreneur share value for loop i
    om_e = ent_share(i);
    
    %----------------------------------------------------------------------
    % 1. Running dynare for a given value of entrepreneur share
    %----------------------------------------------------------------------
    
    save ent_share_value om_e       % saving value of entrepreneur share
    dynare SimpleMatching_answer.mod noclearall    % running dynare

    pert_order= 1;
    VarsToUse = {'c';'n';'v';'y';'m';'w';'Q';'g';'pf';'z'};
 
    X = get_policy_rule_coefs_fcn(pert_order,VarsToUse,1,M_,oo_);
    
    %----------------------------------------------------------------------
    % 2. Simulating the economy given policy functions
    %----------------------------------------------------------------------
    
    s = 0.007;           %Standard deviation of the shock
    T = 10000;           %Length of the series
    D = 500;             %Discarded periods in the beginning
    
    %Shock series:
    randn('seed',666);   % fixing random seed
    e = s*randn(T,1);

    % defining variables
    % Reserving space:
    n = zeros(T,1);
    z = zeros(T,1);
    y = zeros(T,1);
    
    %Initial values:
    y(1,1) = X(1,4);
    n(1,1) = X(1,2);
    z(1,1) = X(1,10);

    %Recursion that computes simulated data using the shock series and the
    %policy functions:

    for t = 2:T          
        y(t,1) = X(1,4) + X(2,4)*(n(t-1)-n(1)) + X(3,4)*z(t-1) + X(4,4)*e(t);
        n(t,1) = X(1,2) + X(2,2)*(n(t-1)-n(1)) + X(3,2)*z(t-1) + X(4,2)*e(t);
        z(t,1) = X(1,10) + X(2,10)*(n(t-1)-n(1)) + X(3,10)*z(t-1) + X(4,10)*e(t);

    end
    %Computing logs, ratios and HP-filtering the artificial data:
    lny  = log(y);
    lnn  = log(n);
    
    %----------------------------------------------------------------------
    % 3. HP filtering simulated data
    %----------------------------------------------------------------------

    hp_lny  = lny(D:end,:)  - hpfilter2(lny(D:end,:),1600);
    hp_lnn  = lnn(D:end,:)  - hpfilter2(lnn(D:end,:),1600);

    %Computing and displaying business cycle statistics
    
    %volatility of log employment relative to volatility of log output
    %Note that in the US data it is 0.466 

    s_ny(1,i)    = std(hp_lnn)/std(hp_lny);
    
    subplot(3,2,i)
    plot(1:101,hp_lnn(end-100:end),'k',1:101,hp_lny(end-100:end),':r','LineWidth',2)
    axis([1 101 -0.025 0.04])

    if i == 1
        legend('employment','output')
        set(legend,'Position',[0.5729 0.2649 0.2028 0.05018]);
%        annotation('textbox','string',{'Cyclical '},'FontWeight','bold', ...
%        'FontSize',12,'FitBoxToText','off','Position', ...
%        [0.4757 0.9903 0 0]);
    end
    
    xx1=sprintf('The share of revenues going to the entrepreneur is %8.3f',om_e);
    xx2=sprintf('The volatility of employment relative to output is %8.3f',s_ny(1,i));
    disp(' ')
    disp(xx1)
    disp(xx2)
    disp(' ')
    disp('Look at the graph and hit return to continue')
    disp(' ')
    pause
end

figure
plot(ent_share,s_ny,'--r',ent_share,ones(I)*0.466,'k','LineWidth',2)
legend('model','US data')
title('std(log(n))/std(log(y))','FontWeight','bold','FontSize',12)

clear n, clear y, clear z

%--------------------------------------------------------------------------
% 4. Computing impulse response functions
%--------------------------------------------------------------------------

s  = 0.007;           %Standard deviation of the shock
Ti = 24;              %Length of the series

%Shock series:
ei    = zeros(Ti,1);
ei(2) = s;

% defining variables
% Reserving space:
v = zeros(Ti,1);
n = zeros(Ti,1);
z = zeros(Ti,1);
y = zeros(Ti,1);

%Initial values:
v(1,1) = X(1,3);
y(1,1) = X(1,4);
n(1,1) = X(1,2);
z(1,1) = X(1,10);

%Recursion that computes simulated data using the shock series and the
%policy functions:

for t = 2:Ti

    v(t,1) = X(1,3) + X(2,3)*(n(t-1)-n(1)) + X(3,3)*z(t-1) + X(4,3)*ei(t);
    y(t,1) = X(1,4) + X(2,4)*(n(t-1)-n(1)) + X(3,4)*z(t-1) + X(4,4)*ei(t);
    n(t,1) = X(1,2) + X(2,2)*(n(t-1)-n(1)) + X(3,2)*z(t-1) + X(4,2)*ei(t);
    z(t,1) = X(1,10) + X(2,10)*(n(t-1)-n(1)) + X(3,10)*z(t-1) + X(4,10)*ei(t);

end

u = ones(Ti,1) - n;

figure
title('XXX')
subplot(2,2,1)
plot(1:Ti-1,z(2:end))
title('aggr. prod.')
subplot(2,2,2)
plot(1:Ti-1,(v(2:end) - v(1))/v(1))
title('vacancies')
subplot(2,2,3)
plot(1:Ti-1,(u(2:end) - u(1))/u(1))
title('unemployment')
subplot(2,2,4)
plot(1:Ti-1,(y(2:end) - y(1))/y(1))
title('output')
%annotation('textbox','string',{'IRFs'},'FontWeight','bold', ...
%    'FontSize',12,'FitBoxToText','on','Position', ...
%    [0.45 0.9903 0 0]);




