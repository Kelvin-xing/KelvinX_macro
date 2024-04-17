%This version includes the programs to find the stationary states,i.e., the
%steady state values of endogenous variables.
%The scripts to find the steady states of labor and capital are hbarfind.m
%and kbarfind.m respectively.
%This script was written by Xiangyang Li, 2015-7-20
%find the steady states
hbarfind; %the steady states of labor
kbarfind; %the steady states of capital
rbar=1/beta-(1-delta);
ybar=kbar^alpha*hbar^(1-alpha);
cbar=ybar-delta*kbar;

%find the solution of the system
A=[0 -kbar 0 0]';
B=[0 (1-delta)*kbar alpha -1]';
C=[1 -1 -1/(1-hbar) 0
    ybar -cbar 0 0
    -1 0 1-alpha 0
    1 0 0 -1];
D=[0 0 1 0]';
F=[0];
G=F;
H=F;
J=[0 -1 0 beta*rbar];
K=[0 1 0 0];
L=F;
M=F;
N=[.95];

a=F-J*(C\A);
b=-(J*(C\B)-G+K*(C\A));
c=-K*(C\B)+H;
P1=(-b+sqrt(b^2-4*a*c))/(2*a);
P2=(-b-sqrt(b^2-4*a*c))/(2*a);

%we need a stable P because it is the coefficient of AR(1)
if abs(P1)<1
    P=P1;
else
    P=P2;
end
R=-C\(A*P+B);
Q=(J*(C\D)-L)*N+K*(C\D)-M;
QD=kron(N',(F-J*(C\A)))+(J*R+F*P+G-K*(C\A));
Q=Q/QD;
S=-C\(A*Q+D);
