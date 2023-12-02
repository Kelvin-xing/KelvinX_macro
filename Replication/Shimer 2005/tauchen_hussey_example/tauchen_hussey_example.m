clear all;
close all;
clc;

%%%%% AR1 process

phi    = 0.90; %% AR1 coefficient
sigeps = 0.10; %% innovation std deviation

% long run moments
mu     = 0;    
sigma  = sigeps/sqrt(1-phi^2);

%%%% discrete-state approximation to AR1 

N      = 33;  %% number of nodes 
floden = 1;   %% indicator for Floden correction 

[nodes,weights,P] = get_tauchen_hussey(mu,sigeps,phi,N,floden); 

%%%% stationary distribution for P

[eig_vectors,eig_values] = eig(P'); 
[~,arg] = min(abs(diag(eig_values)-1)); 
unit_eig_vector = eig_vectors(:,arg); 

pibar = unit_eig_vector/sum(unit_eig_vector); 


%%%%% simulate sample-path from AR1

x0 = 0;
T  = 250;

xt1 = simulate_AR1(T,mu,phi,sigeps,x0);

[chain,state] = simulate_markov_chain(T,nodes,P,x0);

xt2 = chain;

figure(1)
plot([xt1,xt2])
legend('AR1','markov chain')
xlabel('time')
ylabel('x(t)')




return











