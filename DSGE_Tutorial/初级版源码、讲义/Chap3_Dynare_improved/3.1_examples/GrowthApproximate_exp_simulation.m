%simulation code for all endogenous variables for the log-level variable
%mod file in GrowthApproximate_exp.mod.
%this file should be run after  dynare GrowthApproximate_exp;
%this file is adapted at April,3,2015@SuZhou

n=4; %number of controls, y,i,c,lab
m=2;%number of states,k,z

%list variables in DR-order
 M_.endo_names(oo_.dr.order_var,:)
% y  
% i  
% k  
% z  
% c  
% lab

 %retriving model solution representation from oo_ object:
psi=oo_.dr.ghx; %policy function
omega=oo_.dr.ghu; %coefficient of shocks
ss=oo_.dr.ys(oo_.dr.order_var); %declaration order to DR-order

T=200; %number of periods for simulation
es=oo_.exo_simul(:,1); %1000*1

Xsim=zeros(n+m,T);%the simulation results

%the first entry, Xsim is percentage deviation form
%endogenous vars starting from steady states, Xsim from zeros
Xsim(:,1)=omega*es(1,1);

%Xsim(3:,) = capital
%Xsim(4:,) = technology shock
for j=2:T
    Xsim(:,j)=psi*[Xsim(3:4,j-1)]+omega*es(j,1);
end

figure;
subplot(3,2,6)
plot(1:1:T,100*Xsim(1,:))
title('output')
ylabel('% deviation from ss')
axis tight

subplot(3,2,1)
plot(1:1:T,100*Xsim(2,:))
title('investment')
ylabel('% deviation from ss')
axis tight

subplot(3,2,2)
plot(1:1:T,100*Xsim(3,:))
title('capital stock')
ylabel('% deviation from ss')
axis tight

subplot(3,2,3)
plot(1:1:T,100*Xsim(4,:))
title('technology shock')
ylabel('% deviation from ss')
axis tight

subplot(3,2,4)
plot(1:1:T,100*Xsim(5,:))
title('consumption')
ylabel('% deviation from ss')
axis tight

subplot(3,2,5)
plot(1:1:T,100*Xsim(6,:))
title('labor')
ylabel('% deviation from ss')
axis tight

%converting to levels from log-level
%log-difference = log-level -   log-level-steady-states;
%level = exp(log-difference + log-level-steady-states)
for j=1:T
    Xsim(:,j)=exp(Xsim(:,j)+ss);
end
figure;
subplot(3,2,1)
plot(1:1:T,Xsim(1,:))
title('output')
ylabel('levels')
axis tight

subplot(3,2,2)
plot(1:1:T,Xsim(2,:))
title('investment')
ylabel('log-levels')
axis tight

subplot(3,2,3)
plot(1:1:T,Xsim(3,:))
title('capital stock')
ylabel('levels')
axis tight

subplot(3,2,4)
plot(1:1:T,Xsim(4,:))
title('technology')
ylabel('levels')
axis tight

subplot(3,2,5)
plot(1:1:T,Xsim(5,:))
title('consumption')
ylabel('levels')
axis tight

subplot(3,2,6)
plot(1:1:T,Xsim(6,:))
title('labor')
ylabel('levels')
axis tight

%comparing the our simulation results with Dynare results
y_dynare = oo_.endo_simul(1,1:20); %log-level
y_simul  = log(Xsim(1,1:20)); % log(level)
%the results are the same
[y_dynare',y_simul']
%% simulation for technology shock itself
ee=oo_.exo_simul(1:100,1); 
zz=zeros(100,1);
zz(1,1)=rho*0+ee(1,1); %the first value; 
for t=2:100
   zz(t,1)=rho*zz(t-1,1)+ee(t,1); 
end
%the 4th row is the technology variable z

%
[zz oo_.endo_simul(6,1:100)' Xsim(4,1:100)']


%% impulse response for level variable mod file
%please load the dynare estimation results after running level variable mod file:
%load GrowthApproximate_results.mat first before running this cell

%use save command to save desired result into mat file;
load GrowthApproximate_results.mat

%retrive the solution from oo_ object in the mat file;
matE=oo_.dr.ghx;
matF=oo_.dr.ghu;

%number of IRF to be calculated;
T=40;

%where to store the IRF
matIRFS=zeros(size(matE,1),T);

%s is standard deviation of technology shock; initial setup;
matIRFS(:,1)=matF*s;

%the DR order of the model
% k  
% z  
% c  
% lab

%calculating the IRF,we do not have matF here since shocks after period 2
%is zero.
for ii=2:T
    matIRFS(:,ii)=matE*matIRFS(1:2,ii-1);
end

%absolute deviation = level - level_steady_states
figure;
subplot(2,2,1)
plot(1:1:40,matIRFS(1,:),'b-o')
title('capital stock')
ylabel('absolute deviation')
axis tight

subplot(2,2,2)
plot(1:1:40,matIRFS(2,:),'b-o')
title('technology')
ylabel('absolute deviation')
axis tight

subplot(2,2,3)
plot(1:1:40,matIRFS(3,:),'b-o')
title('consumption')
ylabel('absolute deviation')
axis tight

subplot(2,2,4)
plot(1:1:40,matIRFS(4,:),'b-o')
title('labor')
ylabel('absolute deviation')
axis tight

%percentage deviations
%converting to DR order
ys=oo_.dr.ys(oo_.dr.order_var);

figure;
subplot(2,2,1)
plot(1:1:40,100*matIRFS(1,:)/ys(1,1),'b-o')
title('capital stock')
ylabel('%,percentage deviation')
axis tight

%technology has steady state value zero, it no longer need divided it
%steady state, it is all ready in percentage deviation
subplot(2,2,2)
plot(1:1:40,100*matIRFS(2,:),'b-o')
title('technology')
ylabel('%,percentage deviation')
axis tight

subplot(2,2,3)
plot(1:1:40,100*matIRFS(3,:)/ys(3,1),'b-o')
title('consumption')
ylabel('%,percentage deviation')
axis tight

subplot(2,2,4)
plot(1:1:40,100*matIRFS(4,:)/ys(4,1),'b-o')
title('labor')
ylabel('%,percentage deviation')
axis tight

%% comparison of steady state values and irf of level and log-level
%level results
load GrowthApproximate_results.mat
ys1=oo_.dr.ys % steady states in level
varname1=M_.endo_names 
%converting absolute deviation to percentage deviation by dividing steady
%states
irf1=100*oo_.irfs.c_e/ys1(1);

%log-level results
load GrowthApproximate_exp_results.mat
ys2=oo_.dr.ys
ys2exp=exp(oo_.dr.ys) %converting to level
varname2=M_.endo_names
irf2=100*oo_.irfs.c_e; % already in percentage deviation;

%they are the same
[irf1' irf2']