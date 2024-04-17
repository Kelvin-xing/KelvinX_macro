%initial setup
clear all;
close all;
clc;

%sensitivity setting, only one parameter is set for illustration
alpha_d=[0.3 0.35 0.4];

%where the results are saved, here we only interested in the IRF of output
%y with respect to exogenous shock e, y_e;
%the default periods for IRF is 40;
save_results =zeros(40,length(alpha_d));

%delete old mat file
if exist('save_results.mat','file')
    delete save_results.mat
end

for ii=1:length(alpha_d)
    
    if exist('parametersaved.mat','file') 
        delete parametersaved.mat
    end
    
    %alpha will be saved in a mat file named parametersaved.mat
    alpha =alpha_d(ii);
    
    save parametersaved alpha;
    
    %new parameters will be loaded in the Dynare mod file
    %see mod file for how we load new parameters into;
    dynare recursively_running_dynare noclearall;
    
    %directly invoke the complied m file, you need use dynare command in the first place to generate m file 
    %recursively_running_dynare;
    
    %IRF of output w.r.t. exogenous shock e;
    save_results(:,ii) = y_e;
    
end

save save_results save_results;

%% Plot the IRF from the results saved
load save_results;

close all;
plot(100*save_results);
legend('\alpha=0.3','\alpha=0.35','\alpha=0.4');
title('IRF of output to technology shock')