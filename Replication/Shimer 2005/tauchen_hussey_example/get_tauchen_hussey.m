function [nodes,weights,P] = get_tauchen_hussey1(mu,sigeps,phi,N,floden)

% where mu     =  unconditional mean of AR1
%       sigeps = std deviation of innovations
%       phi    = AR1 coefficient
%       N      = number of node points desired
%       floden = 1 to turn on Floden correction for persistence

%%%%% INDICATOR FOR FLODEN CORRECTION
if floden==1,
    
w      = 0.5 + phi/4;
sigx   = sigeps/sqrt(1-phi^2); %% std deviation of unconditional distribution

flodensigma  = w*sigeps + (1-w)*sigx;    
    
else 

flodensigma  = sigeps;

end

%%%%% LOOKUP QUADRATURE NODES AND WEIGHTS
[nodes,weights] = qnwnorm(N,mu,flodensigma^2); 


P = zeros(N,N);
F = zeros(N,N);
regularity_function = zeros(N,1);

%%%%% CONSTRUCT TRANSITION MATRIX
%p[ij] = f[ij]*quadrature_weight(j)/regularity_function(j) 

for i=1:N,
    for j=1:N,

        %%%%% conditional mean
        mean    = (1-phi)*mu + phi*nodes(i);

        %%%%% given we are at node(i), what is likelihood of node(j)?
        F(i,j)  = normpdf(nodes(j),mean,sigeps); 
        
        %%%%% multiply by quadrature weights 
        P(i,j)  = F(i,j)*weights(j);
        
        %%%%% divide by regularity_function
        regularity_function(j) = normpdf(nodes(j),mu,flodensigma);
        
        P(i,j) = P(i,j)/regularity_function(j);
        
    end
end

%%%%% normalize so rows sum to 1
for i=1:N,
    
    P(i,:) = P(i,:) / sum(P(i,:),2);

end





