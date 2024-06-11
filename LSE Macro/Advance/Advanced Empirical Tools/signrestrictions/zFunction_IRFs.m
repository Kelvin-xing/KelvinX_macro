function IRF = zFunction_IRFs(Ahat, P, T_irf, shocked_variab, shock)

% -----------------------------------------------------%

% DESCRIPTION OF THE FUNCTION
% computes impulse response functions in a VAR 

% Generates a k x T_irf matrix containing the impulse response functions of
% the k variables to a shock of value "shock" to variable "shocked_variable". 
% Ahat contains the coefficients of the autoregressive VAR model, P is the
% matrix that maps structural shocks into reduced form shocks

% -----------------------------------------------------%


k = size(Ahat,1); % number of variables in the underlying VAR model
p = size(Ahat,2)/k; % number of lags in the underlying VAR model

Spseudo = zeros(k,T_irf);
Spseudo(shocked_variab,1) = shock;
Rpseudo = P*Spseudo;

IRFestim = zeros(k,T_irf); 
IRFestim(:,1) = Rpseudo(:,1); % the impulse vector
for j=2:T_irf
    step1 = [zeros(k,p-1), IRFestim];
    step2 = step1(:,j-1:j-1+p-1);
    step3 = fliplr(step2);
    step4 = reshape(step3,k*p,1);
    IRFestim(:,j) =  Ahat*step4 + Rpseudo(:,j);  
end

IRF = IRFestim;

end