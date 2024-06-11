function TH = sampleTH(y,X,R,Q,T,Pe,Po,TH,draw)

%Simulation smoother as developed by Carter and Kohn (1994), see Cogley and
%Sargent (2001) or De Wind and Gambetti (2014) for the formulae
%
%   This version: univariate case only
%
%       Joris de Wind (June 2014)

%Kalman filter
for t = 1:T,
    Xtemp       = X(:,t);
    Po(:,:,t)   = XXX;
    K           = XXX;
    Pe(:,:,t+1) = XXX;
    TH(:,t+1)   = XXX;
end

%note that TH(:,t,i) goes from t=1 to t=T+1
%
%In the for loop above, TH(:,t,i) is the estimate of theta_(t-1) conditional on info available up to period t-1
%
%TH(:,1,i) is theta_0 (here set equal to the mean of prior, but you can also take a draw from prior)
%TH(:,T+1,i) corresponds to the last theta_t about which we can say something, i.e., theta_T
%
%Pe(:,:,t+1) is P_(t|t)
%Po(:,:,t)   is P_(t|t-1)


%Draw from TH(T|T)
TH(:,T+1) = TH(:,T+1)+XXX*draw(:,T+1);    
%HINT #1: XXX should be the cholesky decomposition of something
%HINT #2: Think carefully about transpose (see below)

%when A is positive definite then chol(A) produces an upper triangular R so that R'*R = A
%draw(:,t) is a kx1 vector with uncorrelated N(0,1) shocks
%we want that E[X*draw(:,t)*draw(:,t)'*X']=A
%thus X should be the transpose of chol(Pe)


%Draw from TH(t-1|t)
%
%   In the for loop below t goes from high to low.
%
%   at a particular t:
%   1. TH(:,t+1) it is a random draw from a normal that has
%   already been determined (either in this loop or for T above)
%   2. TH(:,t) on the RHS of the mean equation is equal to theta_(t-1)|(t-1)
%   3. TH(:,t) what we end up with is a random draw for theta(t-1) conditional on knowning theta in the next period
%

for t = T:-1:1,
    mean     = XXX;
    variance = XXX;
    TH(:,t) = mean +...
              chol(variance)'*draw(:,t);
end
% in the for loop above, understanding timing is key. Consider the first
% iteration when tt=T. In terms of the model, this corresponds to t = T-1,
% so the terms on the RHSs of the expressions are theta_{T-1|T-1}, P_{T-1|T-1} & P_{T|T-1}.,
% In terms of the program vectors this corresponds to TH(T), Pe(T) & Po(T),
% that is, TH(tt), Pe(tt) & Po(tt)



end

