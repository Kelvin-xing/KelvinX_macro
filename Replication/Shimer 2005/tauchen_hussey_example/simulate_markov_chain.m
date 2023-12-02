function [chain,state] = simulate_markov_chain(T,x,P,s0);

%  x     = the quantity corresponding to each state, typical element x(i)
%  P     = Markov transition matrix, typical element p(i,j) i,j=1,...n
%  s0    = initial state
%  T     = number of periods to simulate
%  
%  chain = sequence of realizations from the simulation

randn('state',0);   % reset seed of random number generator

n = length(x); %% what is the size of the state vector?
E = rand(T,1); %% T-vector of draws from independent uniform [0,1]  

cumsumP = P*triu(ones(size(P)));
               %% creates a matrix whose rows are the cumulative sums of
               %% the rows of P               
               
%%%%% SET INITIAL STATE               

s = zeros(n,1);

[~,arg] = min(abs(x-s0));

s(arg) = 1;

%%%%% ITERATE ON THE CHAIN

for t=1:T,
    state(:,t) = s;
    ppi        = [0,s'*cumsumP];
    s          = ((E(t)<=ppi(2:n+1)).*(E(t)>ppi(1:n)))';
end

chain = x'*state;

chain = chain';
state = state';


