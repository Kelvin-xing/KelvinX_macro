function [IRF_median, IRF_median_target, IRF_low, IRF_high, mt] = zFunction_MedianTarget(IRFs,interval)

% -----------------------------------------------------%

% DESCRIPTION OF THE FUNCTION
% this function computes the median target a' la Fry and Pagan and other
% statistics
%
% IRFs is a kxTxN matrix containing the N candidate IRFs. Interval is
% the interval of the computed bottom and top impulse response. 95% interval 
% has interval = 95, not 0.95. The code also computes the median IRF. The
% representation that corresponds to the median target is denoted mt

% -----------------------------------------------------%


k = size(IRFs,1);       % number of regressors
T_irf = size(IRFs,2);   % length of the IRF
N = size(IRFs,3);       % number of IRFs for each variable

IRF_median = zeros(k,T_irf);
IRF_standdev = zeros(k,T_irf);
IRF_low = zeros(k,T_irf);
IRF_high = zeros(k,T_irf);
for kkk = 1:k
    for ttt = 1:T_irf   
        IRF_median(kkk,ttt)   =  prctile(IRFs(kkk,ttt,:),50);
        IRF_standdev(kkk,ttt) =  std(IRFs(kkk,ttt,:));  
        IRF_low(kkk,ttt)   =  prctile(IRFs(kkk,ttt,:),(100-interval)/2);
        IRF_high(kkk,ttt)   =  prctile(IRFs(kkk,ttt,:),100-(100-interval)/2);   
    end
end

% standardize IRFs, in order to account for the fact that they have
% different deviations around the median
IRF_standardized = zeros(k,T_irf,N); 
for n = 1:N
    for kkk = 1:k
        for ttt = 1:T_irf 
            IRF_standardized(kkk,ttt,n) =  (IRFs(kkk,ttt,n)-IRF_median(kkk,ttt))/IRF_standdev(kkk,ttt);   
        end
    end
end

% combine so that find the median target across the different k dimensions
% considered
vec_IRF_standardized = zeros(N,k*T_irf);
for n = 1:N
    vec_IRF_standardized(n,1:T_irf) = IRF_standardized(1,:,n);
    
    if k > 1
        for kkk = 2:k 
            vec_IRF_standardized(n,T_irf*(kkk-1)+1:T_irf*(kkk-1)+1+T_irf-1) = IRF_standardized(kkk,:,n);
        end
    end
end

step = vec_IRF_standardized*vec_IRF_standardized';
step = diag(step);
mt = find(step == min(step));
IRF_median_target = IRFs(:,:,mt);

end

