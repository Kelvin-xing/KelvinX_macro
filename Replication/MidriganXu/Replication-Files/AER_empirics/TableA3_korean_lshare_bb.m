% ROBUSTNESS
% Blundell and Bond (1998)
% August 6th, 2012

% This matfile requires data created on the stata do file: data4matlab.do
% It performs Dynamic Panel Estimation by 2-digit KSIC industry.


% Fminsearch
% clc
% clear all
% 
% i = [15 17:36];
% 
% for j=1:size(i,2)
%     k=i(1,j);
%     a=['data_korea_bb_' num2str(k) '.mat'];
%     load(a)
% 
% eta=0.85;
% alpha=0.66;
% rho=0.9;
% par = [eta, alpha, rho];
% 
% options=optimset('Display','iter','Diagnostics','off','TolX',1e-13,'TolFun',1e-13,'MaxIter',1e3,'MaxFunEvals',1e13);
% [est,fv]=fminsearch(@gmm_korea3,par, options, lnv, lnv_1, lnv_2, lnv_3, lnl, lnl_1, lnl_2, lnl_3, lnk, lnk_1, lnk_2, lnk_3)
% 
% eta_est_bb(j,1)=exp(est(1,1))/(1+exp(est(1,1)))
% alpha_est_bb(j,1)=exp(est(1,2))/(1+exp(est(1,2)))
% rho_est_bb(j,1)=exp(est(1,3))/(1+exp(est(1,3)))
% mark_up_bb(j,1)=1/eta_est_bb(j,1)
% 
% end

%% 
% GA
clc
clear all
cd ('C:\Users\yx42\Documents\Research Projects\Korea_inv\Revision_2012\Blundell Bond')
i = [15 17:36];
global lnv lnv_1 lnv_2 lnv_3 lnl lnl_1 lnl_2 lnl_3 lnk lnk_1 lnk_2 lnk_3

for j=1:size(i,2)
    k=i(1,j);
    a=['data_korea_bb_' num2str(k) '.mat'];
    load(a)

est=ga(@gmm_korea3ga,3);

eta_est_bb(j,1)=exp(est(1,1))/(1+exp(est(1,1)))
alpha_est_bb(j,1)=exp(est(1,2))/(1+exp(est(1,2)))
rho_est_bb(j,1)=exp(est(1,3))/(1+exp(est(1,3)))
mark_up_bb(j,1)=1/eta_est_bb(j,1)

end

%%