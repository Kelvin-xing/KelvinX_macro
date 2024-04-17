function [Rk,omega,G,F,Gamma,Gam_muG,Fprime,k,n,c] = steadystate(R,sigma,mu,alpha,gamma,delta,z)

% using Rk as a starting point to find the steady states. Xiangyang
% li,2013-5-23
% When Rk is given, then efficiency condition could be
% solved to find omega. 
RRk=R+.001:.001:1.1;
xRk=RRk(1);
[ff,sp1,omega,G,F,Gamma,Gam_muG,Fprime,k,n]=get_bank_zero_profit_diff_cond_on_Rk(xRk,R,sigma,mu,alpha,gamma,delta,z);

%where ff stand for the difference of the two sides of the bank's zero
%profit conditions, this condition are used to pinpoint Rk.
%first, we need a initial value of the difference to find zeros
gg(1)=ff;
for ii = 2:length(RRk)
    xRk=RRk(ii);
    [ff,sp1,omega,G,F,Gamma,Gam_muG,Fprime,k,n]=get_bank_zero_profit_diff_cond_on_Rk(xRk,R,sigma,mu,alpha,gamma,delta,z);    
    gg(ii)=ff;
    if gg(ii)*gg(ii-1) < 0 %trick to find zero
        II=ii;
        break
    end
end
RRk_low=RRk(II-1);
RRk_high=RRk(II);

%zeros range is minimized to [RRk_low RRk_high].
%fzero care about the first return value by get_bank_zero_profit_diff_cond_on_Rk.m.
%values in the interval [RRk_low RRk_high] will be treated as the first
%argument of function get_bank_zero_profit_diff_cond_on_Rk.m.
% the return value x, the first argument of get_bank_zero_profit_diff_cond_on_Rk,  
%that solves get_bank_zero_profit_diff_cond_on_Rk == 0, i.e. the first return value ff==0.

[x,fval,exitflag]=fzero('get_bank_zero_profit_diff_cond_on_Rk',[RRk_low RRk_high],optimset('Display','iter'),R,sigma,mu,alpha,gamma,delta,z);
if exitflag ~= 1 || abs(fval) > .1e-9
    error('fatal (bgg_steadystate) failed to find steady state')
end

%we find the ss of  captial return that solve the bank zero profit condition given
%other parameters.
Rk=x;

%using the steady state of capital return to find other ss.
[ff,sp1,omega,G,F,Gamma,Gam_muG,Fprime,k,n]=get_bank_zero_profit_diff_cond_on_Rk(Rk,R,sigma,mu,alpha,gamma,delta,z); 
c=z*k^alpha-delta*k-mu*G*Rk*k;


