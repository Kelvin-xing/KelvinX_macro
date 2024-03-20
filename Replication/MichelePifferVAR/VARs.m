%% Example 1
clc
clear
close all

% Parameter matrics
A = [1, -5.12; 2.19, 1];
A_minus = [-1.01, -2.02; 1.52, -0.3];
B = A^(-1);
B_minus = B * A_minus;

% stationary check
[V, C] = eig(B_minus);
eigenvalues = max(C);
dominant_root = max(max(abs(C)));
assert(abs(dominant_root) < 1)

% pseudo data
T = 50; % periods
k = 2; % # of variables
randn('seed', 0)
S = mvnrnd(zeros(k,1),eye(k),T)'; % structural shocks with 0 mean and 1 sigma
R = B * S; % reduced form shocks

Y = zeros(k,T); % initial value
Y(:,1) = R(:,1);
for t = 2:T
    Y(:,t) = B_minus * Y(:,t-1) + R(:,t);
end

figure(1)
    subplot(3,2,1)
plot(S(1,:)), hold on
plot(zeros(T,1),'k')
title('Structural shocks to 1'), axis([1 T -Inf Inf]), set (gca, 'box','off')

    subplot(3,2,2)
plot(S(2,:)), hold on
plot(zeros(T,1),'k')
title('Structural shocks to 2'), axis([1 T -Inf Inf]), set (gca, 'box','off')

    subplot(3,2,3)
plot(R(1,:)), hold on
plot(zeros(T,1),'k')
title('Reduced shocks to 1'), axis([1 T -Inf Inf]), set (gca, 'box','off')

    subplot(3,2,4)
plot(R(2,:)), hold on
plot(zeros(T,1),'k')
title('Reduced shocks to 2'), axis([1 T -Inf Inf]), set (gca, 'box','off')

    subplot(3,2,5)
plot(Y(1,:)), hold on
plot(zeros(T,1),'k')
title('Variable 1'), axis([1 T -Inf Inf]), set (gca, 'box','off')

    subplot(3,2,6)
plot(Y(2,:)), hold on
plot(zeros(T,1),'k')
title('Variable 2'), axis([1 T -Inf Inf]), set (gca, 'box','off')

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[0 0 1000*.7 1200*.7])

% Shock correlation check
% within period
corr(S');
corr(R');

% across periods
corr([S(1,2:T)', S(1,1:T-1)']);
corr([S(1,3:T)', S(1,1:T-2)']);

corr([S(2,2:T)', S(2,1:T-1)']);
corr([S(2,3:T)', S(2,1:T-2)']);

corr([R(1,2:T)', R(1,1:T-1)']);
corr([R(1,3:T)', R(1,1:T-2)']);

corr([R(2,2:T)', R(2,1:T-1)']);
corr([R(2,3:T)', R(2,1:T-2)']);

% OLS
% regressor
X = Y(:,1:T-1)';
y = Y(:,2:T)';

B_minus_hat = (X'*X)^(-1)*X'*y;

% rediced form
display(B_minus_hat');
display(B_minus);

% Structural
y = Y(1,2:T)';
X = [Y(2,2:T)', Y(1,1:T-1)', Y(2,1:T-1)'];
step1 = (X'*X)^(-1)*X'*y;
y = Y(2,2:T)';
X = [Y(1,2:T)', Y(1,1:T-1)', Y(2,1:T-1)'];
step2 = (X'*X)^(-1)*X'*y;

A_hat = [1, -step1(1); step2(1), 1]
A
A_minus_hat = [step1(2:3)'; step2(2:3)']
A_minus

%% Forecasting

T = 100;
k = 2;

randn('seed', 0)
S = mvnrnd(zeros(k,1),eye(k),T)'; % structural shocks with 0 mean and 1 sigma
R = B * S; % reduced form shocks

A_minus = [-1.01, -2.02; 1.52, -0.3];
B = A^(-1);
B_minus = B * A_minus;
% original value
Y = zeros(k,T); % initial value
Y(:,1) = R(:,1);
for t = 2:T
    Y(:,t) = B_minus * Y(:,t-1) + R(:,t);
end

% OLS
X = Y(:,1:T-1)';
y = Y(:,2:T)';
B_minus_hat = (X'*X)^(-1)*X'*y;

T_predict = 10;
% Forecasting
Y_pre_estimatedpara = ones(k,T_predict);
Y_pre_truepara = ones(k,T_predict);

Y_pre_estimatedpara(:,1) = B_minus_hat * Y(:,T);
Y_pre_truepara(:,1) = B_minus * Y(:,T);

for t = 2:T_predict
    Y_pre_estimatedpara(:,t) = B_minus_hat * Y_pre_estimatedpara(:,t-1);
    Y_pre_truepara(:,t) = B_minus * Y_pre_truepara(:,t-1);
end

%realization
randn('seed',0)
S_future = mvnrnd(zeros(k,1),eye(k),T_predict);
R_future = B * S_future;

Y_future = zeros(k,T_predict);
Y_future(:,1) = B_minus * Y(:,T) + R_future(:,1);
for t = 2:T_predict
    Y_future(:,t) = B_minus * Y_future(:,t-1) + R_future(:,t);
end

% Ploting
% first 100 periods
figure(1)
    subplot(4,1,1)
plot(Y(1,:)), hold on
plot(zeros(T,1),'k')
title('Variable 1'), axis([1 T -Inf Inf]), set (gca, 'box','off')

    subplot(4,1,2)
plot(Y(2,:)), hold on
plot(zeros(T,1),'k')
title('Variable 2'), axis([1 T -Inf Inf]), set (gca, 'box','off')

