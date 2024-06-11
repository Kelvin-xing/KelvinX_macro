clear;

% This problem solves the stochastic growth model with a labor leisure
% choice.

% Define the parameters

alpha = 0.33;
beta = 1/(1+0.03/4);
delta = 0.025;
gamma = 2;
% Frisch elasticity (or something: I never get what's what)
psi = 2;
% Persistence of shocks
rho = 0.95;
% Standard deviation of shocks
sigma = 0.007;

% Length of stochastic simulation
T = 1000;
% Length of impulse response
TT = 100;

% Define the symbolic variables. I use an "m" to denoted t-1, and a "p" to
% denote t+1. The state is therefore km and zm.

syms cm c cp km k kp lm l lp zm z zp

% Set up the system of equations.

system = [c^(-gamma)-beta*(1+alpha*z*k^(alpha-1)*lp^(1-alpha)-delta)*cp^(-gamma);
    k+c-(zm*km^alpha*l^(1-alpha)+(1-delta)*km);
    c^(-gamma)*(1-alpha)*zm*km^(alpha)*l^(-alpha)-l^psi
    z-(1-rho)-rho*zm];

% Set up the system of equations to solve for the steady state.

ss = [c^(-gamma)-beta*(1+alpha*z*k^(alpha-1)*l^(1-alpha)-delta)*c^(-gamma);
    k+c-(z*k^alpha*l^(1-alpha)+(1-delta)*k);
    c^(-gamma)*(1-alpha)*z*k^(alpha)*l^(-alpha)-l^psi
    z-(1-rho)-rho*z];

% Linearize. 

A = jacobian(system,[km,cm,zm,lm]);

B = jacobian(system,[k,c,z,l]);

C = jacobian(system,[kp,cp,zp,lp]);

% Now find the steady state values of c k l and z, so that we can
% substitute them in to the matrices A B and C.

num_system = matlabFunction(ss,'vars',{[c k l z]});

ss = fsolve(num_system,[1,1,1,1]);

css = ss(1); kss=ss(2); lss=ss(3); zss=ss(4);

Xss = [kss;css;zss;lss];

% Now substitute them in.

A = subs(A,{km k kp cm c cp lm l lp zm z zp},{kss kss kss css css css lss lss lss zss zss zss});

B = subs(B,{km k kp cm c cp lm l lp zm z zp},{kss kss kss css css css lss lss lss zss zss zss});

C = subs(C,{km k kp cm c cp lm l lp zm z zp},{kss kss kss css css css lss lss lss zss zss zss});

% The model is now given by A*u_{t-1}+B*u_t+C*u_{t+1}=0, with u_t = [k c z
% l]'.

% The policy function is called "x". I start with a silly
% initial guess: x=0. That is, the economy reverts to its steady state
% immediately.

x = 0;

metric = 1;

% Time iteration!

while metric>1e-13

    x = inv(B+C*x)*(-A);
    
    metric = max(max(abs(A+B*x+C*x*x)));
    
end


% Simulate model for T periods

X = zeros(4,T+1);

X(1,1) = kss;
X(3,1) = zss;

e = randn(T,1)*sigma;

for t = 1:T
    X(:,t+1) = Xss+x*(X(:,t)-Xss)+[0;0;e(t);0];
end

output = X(3,1:end-1).*X(1,1:end-1).^(alpha).*X(4,2:end).^(1-alpha);
consumption = X(2,2:end);
investment = X(1,2:end)-(1-delta)*X(1,1:end-1);
hours = X(4,2:end);

% Plot the results.

subplot(2,2,1);
plot(output);
title('Output')

subplot(2,2,2);
plot(consumption);
title('Consumption')

subplot(2,2,3);
plot(investment);
title('Investment')

subplot(2,2,4);
plot(hours);
title('Hours')

pause;

% Calculate an impulse reponse:

X = zeros(4,TT);

X(1,1) = kss;
X(3,1) = zss;

X(:,2) = Xss+x*(X(:,1)-Xss)+[0;0;-sigma;0];

for t = 2:TT
    X(:,t+1) = Xss+x*(X(:,t)-Xss);
end

output = X(3,1:end-1).*X(1,1:end-1).^(alpha).*X(4,2:end).^(1-alpha);
consumption = X(2,2:end);
investment = X(1,2:end)-(1-delta)*X(1,1:end-1);
hours = X(4,2:end);

% Plot the results.

subplot(2,3,1);
plot(output);
title('Output')

subplot(2,3,2);
plot(consumption);
title('Consumption')

subplot(2,3,3);
plot(investment);
title('Investment')

subplot(2,3,4);
plot(hours);
title('Hours')

subplot(2,3,5);
plot(X(3,1:end-1));
title('TFP')

% Over and out.


