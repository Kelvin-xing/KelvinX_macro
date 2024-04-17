%% This is the main file that recursively invokes cv.mod file

clear all;
close all;
clc; 

%the high and low volatility; we only focus on this unique parameter here;
%you could add more desired parameters if you want to study the
%sensitivity for the cv welfare metric of these parameters;
%however, the cv will generally depends on the size of sigmae, the
%volatility of exogenous shock.
sigmae_arr = [0.01, 0.02];

% unconditional and conditional welfare metric
uncond_mean = zeros(3,length(sigmae_arr));
cond = zeros(3,length(sigmae_arr));

%running the mod file and save the results
for ii=1:length(sigmae_arr)         
    sigmae=sigmae_arr(ii);
    save parameterfile_cv sigmae;
    dynare cv noclearall;
    
    %what does this mean here, explain more;
    %variables in delcaration order in oo_.mean
    %we extract unconditional welfare metric here; 
    %The welfare value function w,w_c,w_l  are at No. 1, 2,3 in declaration order
    uncond_mean(1,ii) = oo_.mean(1); %value function w
    uncond_mean(2,ii) = oo_.mean(2);%consumption part w_c
    uncond_mean(3,ii) = oo_.mean(3); %labor part w_l
    
    %we extract conditional welfare metric here;
    %oo_.steady_state in declaration order
    %oo_.dr_ghs2 in DR(Decision Rule) order:w n wage y i r k a w_c w_l c Rk.;
    %The welfare value function w,w_c,w_l  are at No. 1, 9,10 in DR order
    %The welfare value function w,w_c,w_l  are at No. 1, 2,3 in declaration order
    cond(1,ii) = oo_.steady_state(1) + oo_.dr.ghs2(1); %value function
    cond(2,ii) = oo_.steady_state(2) + oo_.dr.ghs2(9); %consumption part
    cond(3,ii) = oo_.steady_state(3) + oo_.dr.ghs2(10);%labor part
end

% calculating the unconditional compensation variation welfare metric
nomin= (uncond_mean(1,1) - uncond_mean(3,2) + cab_lab);
denomin = uncond_mean(2,2) +cab_lab;

%see the note for details, lambda_u, unconditional CV
% in section 2.2.2 Additively Separable; lambda_u is usually small, 
%we multiply it by 100;
%under the current parameter setting, cv = 0.0106. This mean that low
%volatility regime will be preferred.
lambda_u = 100*((nomin/denomin)^(1/(1-sigma)) -1);

% conditional compensation variation welfare metric
nomin_c= (cond(1,1) - cond(3,2) + cab_lab);
denomin_c = cond(2,2) +cab_lab;
lambda = 100*((nomin_c/denomin_c)^(1/(1-sigma)) -1);

%after calculation, we display it.
disp('conditional   unconditional');
disp([lambda lambda_u]);