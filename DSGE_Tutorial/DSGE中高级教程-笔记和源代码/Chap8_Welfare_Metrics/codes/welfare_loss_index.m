%This file is the main file that recursively invokes welfare_loss.mod 
clear all;
close all;
clc;

%we iterate on two parameters as in textbook P177,Table7.2
%phi: the reciprocal of elasticity of Frisch labor supply
%epsilon: the elasticity of substitution bwteen varieties produced within
%any given country.
% mu = epsilon/(epsilon-1);
epsilon_arr=[6,11];
phi_arr = [3,10];

beta =0.995;
theta = 0.75;
lambda = (1-beta*theta)*(1-theta)/theta;

%welfare loss metric
welfare_loss_arr=zeros(length(epsilon_arr),length(phi_arr));

var_inflation = zeros(length(epsilon_arr),length(phi_arr));
var_outputgap=zeros(length(epsilon_arr),length(phi_arr));

%running the mode file and save the desired results
for ii=1:length(epsilon_arr)
    for jj=1:length(phi_arr)
        epsilon = epsilon_arr(ii);
        phi = phi_arr(jj);
        save parameterfile_welfare epsilon phi;
        dynare welfare_loss noclearall
        var_outputgap(ii,jj) = oo_.var(1,1);
        var_inflation(ii,jj) = oo_.var(2,2);
        welfare_loss_arr(ii,jj) = epsilon/lambda*var_inflation(ii,jj)+...
            (1+phi)*var_outputgap(ii,jj);
    end
end

disp('The welfare loss');
disp(100*welfare_loss_arr);