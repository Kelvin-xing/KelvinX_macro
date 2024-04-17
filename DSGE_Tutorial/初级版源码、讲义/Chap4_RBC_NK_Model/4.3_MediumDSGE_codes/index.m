clear all
close all
%addpath G:\dynare\4.4.1\matlab

beta = 0.99;
alpha =1/3;
delta = 0.025;
b = 0.65;
tau = 2;
eta = 1;
omegags = 0.2;


epsw = 10;
epsp = 10;
phiw =0.75;
phip=0.75;
phipi = 1.5;
phiy = 0.5;
chi1 = 1/beta - (1-delta);
chi2 = 0.01;

rhog = 0.9;
rhoa =0.9;
rhoi =0.9;
rhoz =0.9;
sdg =0.01;
sda =0.01;
sdi =0.01;
sdz =0.01;


vps = 1;
us  =1;
zs =1;
As =1;
pis = 1; %zeros inflation steady state
pistars =1;
is = pis/beta; %nominal interest rate

%the following algorithm is only valid under the assumption of zero
%inflation steady states. For a non-zeros inflation s.s, you need
%re-calculate, this is left for homework.

%index.m +mediumdsge.mod will lay out a framework within which to
%recursively run mod file for sensitivity test, or different parameter values.

%xs means the steady state value of x;
Rs= 1/beta - (1-delta); %capital rate
mcs = (epsw-1)/epsw;
ws = (1-alpha)*(mcs*alpha^alpha/Rs^alpha)^(1/(1-alpha));
kn= ws/Rs*alpha/(1-alpha);
yn = kn^alpha;
yk = yn/kn;
iy = delta/yk;
cy = 1 - iy -omegags;
Ns = 1/3;
c1 = (1-beta*b)/(1-b);
Ks = kn*Ns;
Ys = Ks*yk;
invs = delta*Ks;
Gs = Ys*omegags;
Cs = Ys - invs - Gs;
psi = ws*(epsw-1)/epsw/Cs/Ns^eta*c1;

Kbars = Ks;
lams = c1/Cs;
mus = lams;
qs = mus/lams;
wstars = ws;
x1s = lams*mcs*Ys/(1-phip*beta);
x2s = lams*Ys/(1-phip*beta);
h1s = psi*Ns^(1+eta)/(1-phiw*beta);
h2s = lams*Ns/(1-phiw*beta);

%save all the variables into mat file and then load in mod file
save mediumdsge;

%run the dynare mod file or the compiled mod m file
dynare mediumdsge