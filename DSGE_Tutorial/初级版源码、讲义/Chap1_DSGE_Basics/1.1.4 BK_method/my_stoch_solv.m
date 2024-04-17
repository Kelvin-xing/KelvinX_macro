%This code is written by Xiangyang Li, 2014-1@UND
clear;
%the parameters and steady states
beta = .99;
alpha = .36;
sigma = 2;
delta = .025;
rho = .9 ;
ks = (alpha/(1/beta - 1 + delta))^(1/(1-alpha)); %steady state of capital stock
cs = ks^alpha - delta*ks; %steady state of consumption;
ck= cs/ks; %consumption capital ratio
R = 1/beta - 1 + delta; %simplifying parameter
ys = ks^alpha; %steady state of proudction
is = ys - cs; %steady state of investment

%preparation for solution
n = 1; %number of jumper
m =2; %number of states
M=zeros(3,3); %the coefficient matrix
M(1,1)= 1- ck*(alpha-1)*beta*R/sigma;
M(1,2)= (alpha-1)*R/sigma;
M(1,3)= beta*R*(rho + (alpha-1)*R/alpha)/sigma;
M(2,1)=-ck;
M(2,2) = 1/beta;
M(2,3) =R/alpha;
M(3,3)=rho;

[vv,lamb] = eig(M); %find the eigenvalues of M;
[lamb_sorted, index] = sort(abs(diag(lamb))); 

%sort the eigenvector matirx
for ii = 1:m + n
    vv_sorted(:,ii) = vv(:,index(ii));
end

%find the index of the first eigenvalue whose value >=1
first_unstable_index = find(abs(lamb_sorted)>=1, 1 );

% num of stable eigs
Q = first_unstable_index - 1; 
% num of unstable eigs
B = m+n - Q; 

G = inv(vv_sorted);
G21 = G(Q+1:Q+B,1:n); % low left, G21,size = B*n;
G22 = G(Q+1:Q+B,n+1:n+m);% lower right, G22, size = B*m;

%the policy function coefficients
pol= -inv(G21)*G22;

%% verify the state space representation of RBC model
M11 = M(1:n,1:n);
M12 = M(1:n,n+1:n+m);
M21 = M(n+1:n+m,1:n);
M22 = M(n+1:n+m,n+1:n+m);
C = M11*pol +M12;
A = M21 * pol + M22;
% exactly the same;
C == pol*A
%% do the IRF
%number of simulation;
%if you use great number,like H=10000, you will get strange results, this
%about why?
%but for small number it works.
H = 10000;

%n+m rows correspondent to the number of variables
IRF = zeros(n+m,H); 
%setting up the first period 
IRF(n+2,1) = 1; %one unit shock to technology shock.
IRF(n+1,1) = 0; % for capital, starts from steady state 

% for the controls, actually only for consumption, n=1, the policy function
IRF(1:n,1) = pol*IRF(n+1:n+m,1); 

%
for ii=2:H
    IRF(:,ii) = M*IRF(:,ii-1);
end

%including another two variables
IRF2 = zeros(2,H);
for jj=1:H
    IRF2(1,jj) = IRF(3,jj) + alpha*IRF(2,jj);%for output
    IRF2(2,jj) = ys/is*IRF2(1,jj) - cs/is*IRF(1,jj);%for investment
end

%plot the IRF
figure(1)
title('IRF of Technology Shock')
subplot(3,2,1)
plot(IRF(1,:),'-b','Linewidth',2)
title('Consumption')

subplot(3,2,2)
plot(IRF2(1,:),'-b','Linewidth',2)
title('Output ')

subplot(3,2,3)
plot(IRF2(2,:),'-b','Linewidth',2)
title('Investment ')

subplot(3,2,4)
plot(IRF(2,:),'-b','Linewidth',2)
title('Capital ')

subplot(3,2,5)
plot(IRF(3,:),'-b','Linewidth',2)
title('Technology')

%% do the IRF by state space representation 
%no matter how large of H, the IRFs are correct!, by construction is
%stable.
H = 10000;

%n+m rows correspondent to the number of variables
IRF = zeros(n+m,H); 
%setting up the first period, the inital states
IRF(n+1,1) = 0; % for capital, starts from steady state 
IRF(n+2,1) = 1; %one unit shock to technology shock.

%the initial controls
% for the controls, actually only for consumption, n=1, the policy function
IRF(1:n,1) = pol*IRF(n+1:n+m,1); 

%state sapce, by construction is stable
M1 = M(1:n,:); %the first n rows
M2 = M(n+1:n+m,:); % the last m rows;
M21=M2(:,1:n);
M22=M2(:,n+1:n+m);
A = M21*pol + M22;
B = [0;1]; %by definition, we only have one exogenous shock, technology shock. 

for ii=2:H
    IRF(n+1:n+m,ii)=A*IRF(n+1:n+m,ii-1);%transition equation
    IRF(1:n,ii)=pol*A*IRF(n+1:n+m,ii-1); %measurement equation
end 

%including another two variables
IRF2 = zeros(2,H);
for jj=1:H
    IRF2(1,jj) = IRF(3,jj) + alpha*IRF(2,jj);%for output
    IRF2(2,jj) = ys/is*IRF2(1,jj) - cs/is*IRF(1,jj);%for investment
end

%plot the IRF
figure(1)
title('IRF of Technology Shock')
subplot(3,2,1)
plot(IRF(1,:),'-b','Linewidth',2)
title('Consumption')

subplot(3,2,2)
plot(IRF2(1,:),'-b','Linewidth',2)
title('Output ')

subplot(3,2,3)
plot(IRF2(2,:),'-b','Linewidth',2)
title('Investment ')

subplot(3,2,4)
plot(IRF(2,:),'-b','Linewidth',2)
title('Capital ')

subplot(3,2,5)
plot(IRF(3,:),'-b','Linewidth',2)
title('Technology')
%% do simulation of the model by state space representation

%ensure the replication of simulation results
randn('state',1234567);

%number of data to be simulated
T = 500;

%sampling from standard normal distribution;
e = 0.1*randn(1,T);

%state sapce, by construction is stable
M1 = M(1:n,:); %the first n rows
M2 = M(n+1:n+m,:); % the last m rows;
M21=M2(:,1:n);
M22=M2(:,n+1:n+m);
A = M21*pol + M22;
B = [0;1]; %by definition, we only have one exogenous shock, technology shock. 
C = pol;

% [consumption, output ,investment]

C = [C;alpha 1;(ys/is)*alpha - (cs/is)*pol(1,1) (ys/is)*1 - (cs/is)*pol(1,2)];

s = zeros(2,T); % states, [capital , technology]
c = zeros(3,T); % controls, [consumption, output ,investment]

%s(:,1) = B*e(1,1);

for j = 2:T
    s(:,j) = A*s(:,j-1) + B*e(1,j); %the state equations
end

for j = 1:T
    c(:,j) = C*s(:,j); %the measurement equations, states map to controls;
end
% coverting to levels from log-linearized form
Simss = repmat( [ks 1]',1,T);
s =Simss.* (repmat([1;1],1,T) + s);
% 
Sim2ss = repmat([cs ys is]',1,T);
c = Sim2ss.*(1+c);

% do hp filtering
lambda = 1600;

%*t, trend part; *c, cycle part;
% [yt,yc] = hp_filter(c(2,:)',lambda);
% [ct,cc] = hp_filter(c(1,:)',lambda);
% [it,ic] = hp_filter(c(3,:)',lambda);
% 
% find out the simulated moments:standard deviations
% from cyclic part of the data
% SDy = std(yc);
% SDc = std(cc);
% SDi = std(ic);
% 
% relative volatility
% Rc = SDc/SDy;
% Ri = SDi/SDy;

% the following simulation will explode up.
%add the shocks to the system
% Sigma = [0; 0 ;1];
% Sim = zeros(n+m,T);
% %Sim(:,1) = Sigma*e(1,1);
% for ii=2:T
%     Sim(:,ii) = M*Sim(:,ii-1) + Sigma*e(1,ii);
% end
% % Y and Investment
% Sim2 = zeros(2,T);
% for jj=1:T
%     Sim2(1,jj) = Sim(3,jj) + alpha*Sim(2,jj);
%     Sim2(2,jj) = ys/is*Sim2(1,jj) - cs/is*Sim(1,jj);
% end
%% plot out the simulated variable in levels
%plot the IRF
figure(1)
subplot(2,2,1)
plot(c(1,:),'-b','Linewidth',1)
title('Consumption ')

subplot(2,2,2)
plot(c(2,:),'-b','Linewidth',1)
title('Output ')

subplot(2,2,3)
plot(c(3,:),'-b','Linewidth',1)
title('Investment')

subplot(2,2,4)
plot(e,'-b','Linewidth',1)
title('Technology Shock')
