function [Ahat,Consthat,Residuals,T_used,Sigmahat] = zFunction_RVAR(Y,p,constantD)

    % -----------------------------------------------------%
    
    % DESCRIPTION 
    % estimates the reduced-form VAR model for the variables
    % included in the matrix Y adding p lags. A constant is added if the 
    % dummy constantD equals 1. It generates:
    
    % Ahat = autoregression parameters 
    % Consthat = the estimates of the constants
    % Residuals = the reduced form shocks estimated
    % T_used = number of observations used in the estimation 
    % Sigmahat = the variance covariance matrix of the 
    
    % -----------------------------------------------------%
            
    % PARAMETERS USED
    k = size(Y,1);
    T = size(Y,2);

    % OPERATIONAL VARIABLES
    y = Y(:,p+1:T)';

    X = Y(:,p:T-1)';
    for g = 1:p-1
       X = [X, Y(:,p-g:T-1-g)'];   % could preallocate to reduce speed
    end

    % constant
    if constantD == 1
        X = [X, ones(T-p,1)];
    end

    df_lost = size(X,2);
    T_used = T-p;

    % ESTIMATE THE MODEL
    Ahat = (X'*X)^(-1)*X'*y;     

    Resid = y-X*Ahat;
    RSS = Resid'*Resid;                       
    Sigmahat = RSS/(T_used-df_lost);    % variance of the reduced-form shocks

    % REARRANGE ESTIMATES
    if constantD == 1
        Consthat = Ahat((k*p)+constantD,:)';
    else
        Consthat = NaN*ones(k,1);
    end
    
    Ahat = Ahat(1:(k*p),:)';
    Residuals = Resid';

end