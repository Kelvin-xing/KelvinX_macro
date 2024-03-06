clear;

%% Parameters

sigma = 1/1.16;         % Elasticity if intertemporal substitution.
beta = 0.997;           % Discount factor.
alpha = 0.7747;         % Price stickiness (measure of firms that cannot reset their prices in a given period). 
theta = 12.7721;        % Elasticity of substitution.
omega = 1.5692;         % Frisch elasticity of labour supply.
phi_pi = 1.5;           % Coefficient on pi in Taylor rule.
phi_y = 0.5/4;          % Coefficient on y in Taylor rule.
rs = -0.0104;           % Shock to preference that brings the economy to a liquidity trap.
mu = 0.9030;            % Probability of remaining in a liquidity trap.

% Reduced form parameters

kappa = (1-alpha)*(1-alpha*beta)/alpha*(1/sigma+omega)/(1+omega*theta);

psi = 1/(1/sigma+omega);

% Theoretical multiplier according to equation (30) in Eggertson

dydg = ((1-mu)*(1-beta*mu)-mu*kappa*psi)/((1-mu)*(1-beta*mu)-sigma*mu*kappa);

% Transition matrix.

T = [1, 0;1-mu,mu];

%% Steady state values.

y_ss = 0;
pi_ss = 0;
i_ss = -log(beta);
r_ss = i_ss;
g_ss = 0;

% Set up the two systems symbolically.

syms ym y ypt ypn piem pie piept piepn gm g gpt gpn im i ip rm r rp

no_trap = [-y+T(1,1)*ypn+T(1,2)*ypt-sigma*(i-(T(1,1)*piepn+T(1,2)*piept)-r)+g-(T(1,1)*gpn+T(1,2)*gpt);
    -pie+kappa*y+kappa*psi*(-1/sigma*g)+beta*(T(1,1)*piepn+T(1,2)*piept);
    -i+r+phi_pi*pie+phi_y*y;
    -r+r_ss;
    -g+g_ss;];

trap = [-y+T(2,1)*ypn+T(2,2)*ypt-sigma*(i-(T(2,1)*piepn+T(2,2)*piept)-r)+g-(T(2,1)*gpn+T(2,2)*gpt);
    -pie+kappa*y+kappa*psi*(-1/sigma*g)+beta*(T(2,1)*piepn+T(2,2)*piept);
    -i+0;
    -r+rs;
    -g+gm;];

% If you uncomment the below code you will instead find the multiplier when
% the economy is not in a liquidity trap.

% trap = [-y+T(2,1)*ypn+T(2,2)*ypt-sigma*(i-(T(2,1)*piepn+T(2,2)*pipt)-rm)+g-(T(2,1)*gpn+T(2,2)*gpt);
%     -pie+kappa*y+kappa*psi*(-1/sigma*g)+beta*(T(2,1)*piepn+T(2,2)*pipt);
%     -i+rm+phi_pi*pie+phi_y*y;
%     -r+r_ss;
%     -g+gm;];



%% Get the matrices. Tedious a lot of coding. But repetitive.

variables = [ym y ypt ypn piem pie pipt piepn gm g gpt gpn im i ip rm r rp];
steady_states = [y_ss y_ss y_ss y_ss pi_ss pi_ss pi_ss pi_ss g_ss g_ss g_ss g_ss i_ss i_ss i_ss r_ss r_ss r_ss];

Xm = [ym piem gm im rm];
X = [y pie g i r];
Xpn = [ypn piepn gpn ip rp];
Xpt = [ypt piept gpt ip rp];

A_n = jacobian(no_trap,Xm);
A_n = subs(A_n,variables,steady_states);

B_n = %%XX
B_n = subs(B_n, variables,steady_states);

C_nn = %%XX
C_nn = subs(C_nn,variables,steady_states);

C_nt = %%XX
C_nt = subs(C_nt,variables,steady_states);

A_t = %%XX
A_t = subs(A_t,variables,steady_states);

B_t = %%XX
B_t = subs(B_t,variables,steady_states);

C_tn = %%XX
C_tn = subs(C_tn,variables,steady_states);

C_tt = %%XX
C_tt = subs(C_tt,variables,steady_states);

D_n = %%XX

D_t = %%XX


% Make sure the matrices are only numbers, not symbolic
A_n  = double(A_n);
B_n  = double(B_n);
C_nn = double(C_nn);
C_nt = double(C_nt);
A_t  = double(A_t);
B_t  = double(B_t);
C_tn = double(C_tn);
C_tt = double(C_tt);
D_n  = double(D_n);
D_t  = double(D_t);




%% Initial guesses. Stupid as usual.

E_n = %%XX 
E_t = %%XX

F_n = %%XX
F_t = %%XX

metric = 1;

%% Let's solve the problem.

while metric>1e-13
    
    % First calculate new guesses given the previous:

    E_n1 = %%XX
    E_t1 = %%XX

    F_n1 = %%XX
    F_t1 = %%XX
    
    % Check for convergence

    metric = %%XX
    
    % Then update the guesses for the next iteration:

    E_n = %%XX
    E_t = %%XX
    F_n = %%XX
    F_t = %%XX

end

% Done problem solved


%% IRF
% Set the length of the impulse response.

T_imp = 50;

v = [0,1];

for i = 1:T_imp
    v(i+1,:) = v(i,:)*T;
end

Ptn = v(1:end-1,2)./v(2:end,1)*T(2,1);
Pnt = v(1:end-1,1)./v(2:end,2)*T(1,2);
Ptt = v(1:end-1,2)./v(2:end,2)*T(2,2);
Pnn = v(1:end-1,1)./v(2:end,1)*T(1,1);

% Initial value:

% uss = zeros(5,1);

ussg = uss;
ussg(3) = 0.05;

Ex_n = uss;
Ex_t = uss;
Exp = uss;

Exg_n = ussg;
Exg_t = ussg;
Exp_g = ussg;

for i = 1:T_imp-1

    Ex_t(:,i+1) = Ptt(i)*(E_t+F_t*Ex_t(:,i))+Pnt(i)*(E_t+F_t*Ex_n(:,i));
    Ex_n(:,i+1) = %%XX
    Exp(:,i+1) = v(i,1)*Ex_n(:,i+1)+v(i,2)*Ex_t(:,i+1);    
    Exg_t(:,i+1) = %%XX
    Exg_n(:,i+1) = %%XX   
    Exp_g(:,i+1) = %%XX
    
end

uss = zeros(5,1);

% Ok done. Let's calculate the impulse responses when the negative shock 
% lasts for 10, 20, 30, and 40 periods

X10 = zeros(5,1);

X10(:,2) = E_t+F_t*X10;

X20 = X10;
X30 = X10;
X40 = X10;

X10g = zeros(5,1);

X10g(3) = 0.05;

X10g(:,2) = E_t+F_t*X10g;

X20g = X10g;
X30g = X10g;
X40g = X10g;

for i = 2:T_imp
    
%% XX
    
end

Exp_g(3,1) = 0;

subplot(2,2,1);
plot(100*(Exp(1,:)),'LineWidth',1.6);
hold on
plot(100*(Exp_g(1,:)),'LineWidth',1.6);
plot(100*(X10(1,:)),'LineWidth',1.0,'color',[0,0,0]+0.5);
plot(100*(X20(1,:)),'LineWidth',1.0,'color',[0,0,0]+0.5);
plot(100*(X30(1,:)),'LineWidth',1.0,'color',[0,0,0]+0.5);
plot(100*(X40(1,:)),'LineWidth',1.0,'color',[0,0,0]+0.5);
ylabel('Output, Y_{t}','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')

subplot(2,2,2);
plot(100*(Exp(2,:)),'LineWidth',1.6);
hold on
plot(100*(Exp_g(2,:)),'LineWidth',1.6);
plot(100*(X10(2,:)),'LineWidth',1.0,'color',[0,0,0]+0.5);
plot(100*(X20(2,:)),'LineWidth',1.0,'color',[0,0,0]+0.5);
plot(100*(X30(2,:)),'LineWidth',1.0,'color',[0,0,0]+0.5);
plot(100*(X40(2,:)),'LineWidth',1.0,'color',[0,0,0]+0.5);
ylabel('Inflation, \pi_t','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')


subplot(2,2,3);
plot(100*(Exp_g(3,:)),'LineWidth',1.6);
hold on
plot(100*(X10g(3,:)),'LineWidth',1.0,'color',[0,0,0]+0.5);
plot(100*(X20g(3,:)),'LineWidth',1.0,'color',[0,0,0]+0.5);
plot(100*(X30g(3,:)),'LineWidth',1.0,'color',[0,0,0]+0.5);
plot(100*(X40g(3,:)),'LineWidth',1.0,'color',[0,0,0]+0.5);
ylabel('Government spending, G_{t}','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
xlabel('Time (quarters)','FontSize',12,'fontname','times')

subplot(2,2,4);
plot((Exp_g(1,:)-Exp(1,:))./Exp_g(3,:),'LineWidth',1.6);
axis([0 50 0 3])
ylabel('Fiscal multiplier','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
xlabel('Time (quarters)','FontSize',12,'fontname','times')
