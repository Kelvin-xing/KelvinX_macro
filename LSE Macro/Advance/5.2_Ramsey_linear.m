clear;

% A simple example solving the stochastic growth model with endogenous
% labor supply using linear time iteration. Pontus Rendahl, 20/09/2017.

% Declare parameters

alpha = 1/3;        % Capital share of output
beta = 1.03^(-1/4); % Discount factor.
gamma = 2;          % Coefficient of risk aversion
eta = 2;            % Frisch elasticity of labor supply
delta = 0.025;      % Depreciation rate of capital
rho = 0.9;          % Persistence of TFP process.

% First solve for the deterministic steady state: Ee is the Euler equation;
% Rc the resource constraint; and Ls is labour supply.

Ee = @(x) -x(1).^(-gamma)+beta*(1+alpha*x(2).^(alpha-1).*x(3).^(1-alpha)-delta).*x(1).^(-gamma);
Rc = @(x) -x(1)-x(2)+x(2).^(alpha).*x(3).^(1-alpha)+(1-delta)*x(2);
Ls = @(x) -x(1).^(-gamma).*(1-alpha).*x(2).^(alpha).*x(3).^(-alpha)+x(3).^(eta);

% Collect the equations as a system of equations.

ss = @(x) [Ee(x);Rc(x);Ls(x)];

% Let fsolve do the job (but you can change it).

options = optimoptions('fsolve','Display','off','tolfun',1e-10,'tolx',1e-10);
xss = fsolve(ss,[1 1 1],options);

css = xss(1);
kss = xss(2);
lss = xss(3);

% Let's include output and investment too just for illustration

yss = kss^(alpha)*lss^(1-alpha);
Iss = kss-(1-delta)*kss;

% So the vector of variables at the steady state is:

xss = [yss,Iss,xss,1];

% Set up the stochastic system symbolically

syms ym y yp cm c cp Im I Ip km k kp lm l lp zm z zp

system = [-y+z*km^(alpha)*l^(1-alpha);
    -I+k-(1-delta)*km;
    -c.^(-gamma)+beta*(1+zp*alpha*k.^(alpha-1).*lp.^(1-alpha)-delta).*cp.^(-gamma);
    c+k-(z.*km.^(alpha).*l.^(1-alpha)+(1-delta)*km);
    c.^(-gamma).*(1-alpha).*km.^(alpha).*l.^(-alpha).*z-l.^(eta);
    -z+zm*rho];

% Useful vectors not to confuse things.

Xm = [ym Im cm km lm zm];
X = [y I c k l z];
Xp = [yp Ip cp kp lp zp];
Xss = [yss yss yss Iss Iss Iss css css css kss kss kss lss lss lss 1 1 1];
Vars = [ym y yp Im I Ip cm c cp km k kp lm l lp zm z zp];

% Linearize system.

A = jacobian(system,Xm); A = double(subs(A,Vars,Xss));
B = jacobian(system,X);  B = double(subs(B,Vars,Xss));
C = jacobian(system,Xp); C = double(subs(C,Vars,Xss));

% Convert to log-linear system (doesn't matter as long as you interpret things correctly).

M = ones(6,1)*xss;
A = A.*M; B = B.*M; C = C.*M;

% Solve.

metric = 1;
F = 0;

while metric>1e-12
    F = -(B+C*F)\A;
    metric = max(max(abs((A+B*F+C*F*F))));
end
Q = -(B+C*F)\eye(6);

% The problem is solved. The below is all about illustrating the results.
% We will do this in two ways: 1, calculate impulse response functions; and
% 2, show a stochastic simulation.

% 1, calculate impulse response functions

T = 40;

u = zeros(6,1); u(end) = 1;
x(:,1) = Q*u;

for t = 1:T-1
    x(:,t+1) = F*x(:,t);
end

% Plot the results

black = [0 0 0];

add = 0.015;
add3 = 0.015;
add2 = add3/2;
add4 = add/4;

figure;
h = subplot(3,2,1);
ph = get(h, 'pos');
ph(4) = ph(4) + add;
ph(3) = ph(3) + add2;
ph(2) = ph(2) - add4;
set(h, 'pos', ph);
p1 = plot(x(1,:),'linewidth',1.6,'color',black);
set(gca,'FontSize',8,'fontname','times')
ylabel('Percent deviation','FontSize',10,'fontname','times')
title('Output','FontSize',10,'fontname','times','FontWeight','Normal')

h = subplot(3,2,2);
ph = get(h, 'pos');
ph(4) = ph(4) + add;
ph(3) = ph(3) + add2;
ph(1) = ph(1) - add3;
ph(2) = ph(2) - add4;
set(h, 'pos', ph);
p1 = plot(x(3,:),'linewidth',1.6,'color',black);
set(gca,'FontSize',8,'fontname','times')
title('Consumption','FontSize',10,'fontname','times','FontWeight','Normal')

h = subplot(3,2,3);
ph = get(h, 'pos');
ph(4) = ph(4) + add;
ph(3) = ph(3) + add2;
ph(2) = ph(2) - add4;
set(h, 'pos', ph);
p1 = plot(x(2,:),'linewidth',1.6,'color',black);
set(gca,'FontSize',8,'fontname','times')
ylabel('Percent deviation','FontSize',10,'fontname','times')
title('Investment','FontSize',10,'fontname','times','FontWeight','Normal')

h = subplot(3,2,4);
ph = get(h, 'pos');
ph(4) = ph(4) + add;
ph(3) = ph(3) + add2;
ph(1) = ph(1) - add3;
ph(2) = ph(2) - add4;
set(h, 'pos', ph);
p1 = plot([0,x(4,1:end-1)],'linewidth',1.6,'color',black);
set(gca,'FontSize',8,'fontname','times')
title('Capital','FontSize',10,'fontname','times','FontWeight','Normal')

h = subplot(3,2,5);
ph = get(h, 'pos');
ph(4) = ph(4) + add;
ph(3) = ph(3) + add2;
ph(2) = ph(2) - add4;
set(h, 'pos', ph);
p1 = plot(x(5,:),'linewidth',1.6,'color',black);
set(gca,'FontSize',8,'fontname','times')
ylabel('Percent deviation','FontSize',10,'fontname','times')
title('Hours','FontSize',10,'fontname','times','FontWeight','Normal')
xlabel('Time (quarters)','FontSize',10,'fontname','times')

h = subplot(3,2,6);
ph = get(h, 'pos');
ph(4) = ph(4) + add;
ph(3) = ph(3) + add2;
ph(1) = ph(1) - add3;
ph(2) = ph(2) - add4;
set(h, 'pos', ph);
p1 = plot(x(6,:),'linewidth',1.6,'color',black);
set(gca,'FontSize',8,'fontname','times')
title('Productivity','FontSize',10,'fontname','times','FontWeight','Normal')
xlabel('Time (quarters)','FontSize',10,'fontname','times')

xSize = 20.5; 
ySize = 18;

set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
     ,'PaperPosition',[0 0 xSize ySize],'PaperSize',[xSize-3.2 ySize-1.2],'PaperPositionMode','auto')

print -dpdf -painters IRF.pdf

% 2, calculate a stochastic simulation

T = 200;

u = zeros(6,T); u(end,:) = randn(1,T);
x(:,1) = Q*u(:,1);

for t = 1:T-1
    x(:,t+1) = F*x(:,t)+Q*u(:,t+1);
end

% Plot the results

figure;
h = subplot(3,2,1);
ph = get(h, 'pos');
ph(4) = ph(4) + add;
ph(3) = ph(3) + add2;
ph(2) = ph(2) - add4;
set(h, 'pos', ph);
p1 = plot(x(1,:),'linewidth',1.6,'color',black);
set(gca,'FontSize',8,'fontname','times')
ylabel('Percent deviation','FontSize',10,'fontname','times')
title('Output','FontSize',10,'fontname','times','FontWeight','Normal')

h = subplot(3,2,2);
ph = get(h, 'pos');
ph(4) = ph(4) + add;
ph(3) = ph(3) + add2;
ph(1) = ph(1) - add3;
ph(2) = ph(2) - add4;
set(h, 'pos', ph);
p1 = plot(x(3,:),'linewidth',1.6,'color',black);
set(gca,'FontSize',8,'fontname','times')
title('Consumption','FontSize',10,'fontname','times','FontWeight','Normal')

h = subplot(3,2,3);
ph = get(h, 'pos');
ph(4) = ph(4) + add;
ph(3) = ph(3) + add2;
ph(2) = ph(2) - add4;
set(h, 'pos', ph);
p1 = plot(x(2,:),'linewidth',1.6,'color',black);
set(gca,'FontSize',8,'fontname','times')
ylabel('Percent deviation','FontSize',10,'fontname','times')
title('Investment','FontSize',10,'fontname','times','FontWeight','Normal')

h = subplot(3,2,4);
ph = get(h, 'pos');
ph(4) = ph(4) + add;
ph(3) = ph(3) + add2;
ph(1) = ph(1) - add3;
ph(2) = ph(2) - add4;
set(h, 'pos', ph);
p1 = plot([0,x(4,1:end-1)],'linewidth',1.6,'color',black);
set(gca,'FontSize',8,'fontname','times')
title('Capital','FontSize',10,'fontname','times','FontWeight','Normal')

h = subplot(3,2,5);
ph = get(h, 'pos');
ph(4) = ph(4) + add;
ph(3) = ph(3) + add2;
ph(2) = ph(2) - add4;
set(h, 'pos', ph);
p1 = plot(x(5,:),'linewidth',1.6,'color',black);
set(gca,'FontSize',8,'fontname','times')
ylabel('Percent deviation','FontSize',10,'fontname','times')
title('Hours','FontSize',10,'fontname','times','FontWeight','Normal')
xlabel('Time (quarters)','FontSize',10,'fontname','times')

h = subplot(3,2,6);
ph = get(h, 'pos');
ph(4) = ph(4) + add;
ph(3) = ph(3) + add2;
ph(1) = ph(1) - add3;
ph(2) = ph(2) - add4;
set(h, 'pos', ph);
p1 = plot(x(6,:),'linewidth',1.6,'color',black);
set(gca,'FontSize',8,'fontname','times')
title('Productivity','FontSize',10,'fontname','times','FontWeight','Normal')
xlabel('Time (quarters)','FontSize',10,'fontname','times')

xSize = 20.5; 
ySize = 18;

set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
     ,'PaperPosition',[0 0 xSize ySize],'PaperSize',[xSize-3.2 ySize-1.2],'PaperPositionMode','auto')

print -dpdf -painters stoch_sim.pdf

