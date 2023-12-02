function xt = simulate_AR1(S,mu,phi,sigeps,x0);

randn('state',0);   % reset seed of random number generator

E = rand(S,1);    %% S-vector of draws from independent uniform [0,1]  

epsilons = sigeps*norminv(E);

%%%%% iteratively construct sample path

xt    = zeros(S,1);

for s=2:S,

xt(1) = x0;
xt(s) = (1-phi)*mu+phi*xt(s-1)+epsilons(s);

end

time = (1:1:S)';







