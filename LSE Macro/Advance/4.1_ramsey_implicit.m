clear;

% Parameters
rho = 0.01;
alpha = 0.33;
delta = 0.1;
gamma = 2;
N = 101;
% Smoothing parameter
lambda = 0.5;

% Utility function
if gamma == 1
    u = @(c) log(c);
else
    u = @(c) (c.^(1-gamma)-1)/(1-gamma);
end
mu = @(c) c.^(-gamma);
mu_inv = @(x) x.^(-1/gamma);

% Grid
kss = (alpha/(rho+delta))^(1/(1-alpha));
kgrd = linspace(kss*0.5,kss*1.5,N)';
ygrd =  kgrd.^(alpha)-delta*kgrd;

% Derivative matrix
dk = kgrd(2)-kgrd(1);
D = zeros(N,N);
D(1,1) = -1/dk; D(1,2) = 1/dk; D(end,end-1) = -1/dk; D(end,end) = 1/dk;
for i = 2:N-1
    D(i,i-1) = -1/(2*dk); 
    D(i,i+1) = 1/(2*dk);
end

% Initial guess
c0 = 0.5*(ygrd+delta*kgrd);
V0 = u(c0)/rho;
I = eye(N);

metric = 1;
iter = 0;
% Let's go. 
tic
while metric>1e-9
    iter = iter+1;
    DV0 = D*V0;
    c1 = mu_inv(DV0);
    
    V1 = (rho*I-diag(ygrd-c1)*D)\u(c1);
    metric = (-rho*V0+u(c1)+DV0.*(ygrd-c1))'*(-rho*V0+u(c1)+DV0.*(ygrd-c1))
    V0 = lambda*V1+(1-lambda)*V0;
end
time = toc;

% Display results
fprintf('Problem is solved using %d iterations, taking %.2f seconds.\n',iter,time);
% And plot
plot(kgrd,ygrd-c1,'LineWidth',1.5);
hold on
plot(kss,0,'.','MarkerSize',15);


