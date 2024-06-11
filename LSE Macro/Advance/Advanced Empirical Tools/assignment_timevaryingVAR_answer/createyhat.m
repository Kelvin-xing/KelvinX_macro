function yhat = createyhat(y,X,TH,T)

yhat = zeros(T,1);
for t = 1:T
    yhat(t)=X(:,t)'*TH(:,t);
end
yhat = y-yhat;    
end


%=========================================================================
% Joris' original createyhat file was more complex but allowed for more
% general problems and in particular dealt efficiently with the case when
% the DGP was an actual VAR the idea is to stack the different equations in
% the VAR (typically we program such that we have the RHS variables of the
% VAR in columns, here everything is stacked in one column)
%
%=========================================================================
%
% In the main program you would have
%
%                         %Sample R1
%                         yhat  = reshape(createyhat(y,X,mTH(:,2:end),T,n),n,T);
%                         Rtemp = yhat*yhat'+T0R*Rbar; Rtemp = (Rtemp+Rtemp')/2;
%                         mR1   = iwishrnd(Rtemp,T+T0R);
%
%=========================================================================
%
% And the function would be as follows
%
%
%=========================================================================
%
% function yhat = createyhat(y,X,TH,T,n)
% 
% %Declare persistent memory
% persistent index;   % 'persistent' ensures that index is known whenever
%                     % Matlab gets into this function
% 
% %Create index if not already done
% if isempty(index),
%     for i = 1:n,
%         index(i,:) = i+n*(T+1)*(0:T-1);
%     end
%     index = index(:);
% end
% 
% %Create yhat
% temp = X'*TH;
% 
% 
% yhat = y-temp(index);
% 
% end
%
%
%=========================================================================
%
% INPUTS:
% 
% - k = n + p * n^2, i.e. the number of time-varying parameters; not used as input for createyhat but used in the discussion below
% 
% - Input y is one long vector of dimension nT x 1; first all variables for first period, then all variables for second period, etc.
% 
% - Input X is one big matrix of dimension k x nT; details are a bit more complicated to describe in this email, see in chapter 2 of my dissertation for details, on page 12 I explain how the k x n matrix X(t) looks like, then to get the k x nT matrix X from the k x n matrix X(t) it is just stacked similarly as y (but then over the second dimension instead of the first)
% 
% - Input mTH is one big matrix of dimension k x T+1; so the input mTH actually given to the routine is k x T since the index starts from the second column; note that in the snapshot given above there is no third index for mTH, as opposed to the example code I gave you where there is a third dimension with index i that runs over the various Gibbs steps, it is just that my full-blown code does the saving of Gibbs draws at the end (to account for thinning directly and avoid memory issues for really long chains)
% 
% OUTPUT:
% 
% The output is a long vector of dimension nT x 1 with the prediction errors. I reshape them to make them in the appropriate n x T format. Of course, with n = 1 as in the example code I gave you there is no need for a reshape.
% 
