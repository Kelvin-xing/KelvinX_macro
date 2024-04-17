%This is the main file which invoke the mod file
%initial setup
clear all;
close all;
clc;

%sensitivity setting, only one parameter is set for illustration.
alpha_d = [0.3,0.35,0.4];

%where results are saved, here we are interested in the 
%IRF of output y with respect to exogenous shock e, y_e;
%the default periods is 40 for IRF.
saved_results=zeros(40,length(alpha_d));

%delete old results mat file.
if exist('saved_results.mat','file')
delete saved_results.mat
end

for ii=1:length(alpha_d)
    %delete the old mat file
     if exist('parametersaved.mat','file')
    delete parametersaved.mat
     end
    
    %new parameter value are assigned
    alpha= alpha_d(ii);
    
    % alpha will be saved in a mat file named as parametersaved.mat
    % if more parameters are stored you may use the command like 
    % save parametersaved alpha beta ...;
    save parametersaved alpha;
    
    %new parameters will be loaded in the Dynare mod file;
    %see mod file for how we load new parameters into mod file;
    dynare recursively_running_dynare noclearall
    
    %directly invoke the complied m file, you need 
    %use the dynare command at first place.
    %recursively_running_dynare
    
    %IRF of output y w.r.t exogenous shock;
    saved_results(:,ii) = y_e;
end

%% Plot the IRF from the results saved.
%we are not going to detail the plot.
 plot(saved_results*100);
 legend('\alpha =0.3','\alpha =0.35','\alpha =0.4');
