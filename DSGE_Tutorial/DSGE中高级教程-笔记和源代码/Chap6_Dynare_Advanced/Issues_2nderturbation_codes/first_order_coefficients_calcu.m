%2015-10-14@Beijing
%written by Xiangyang  Li

beta = 0.99;
gamma = 2; %20
alpha = 0.36;
delta =0.02;
rho =0.95;
Veps = 0.01^2;

Kbar = (alpha*beta/(1-(1-delta)*beta))^(1/(1-alpha));
kbar = log(Kbar);

%find g_k, reads as g sub k;
ca = 1;  % coefficient of quad. equation for gk;
as = 0 ; %steady state of technology
fbar = exp(alpha*kbar+as)+(1-delta)*exp(kbar);
cbar = fbar - exp(kbar);

%the 1st and 2nd derivatives of utility
uprime = cbar^(-gamma);
udprime = -gamma*cbar^(-1-gamma);

fK = alpha*Kbar^(alpha-1)*exp(as)+(1-delta);
fk = Kbar*fK;
fKk= alpha*(alpha-1)*exp((alpha-1)*kbar+as);
cb = -(1+1/beta + uprime/udprime*fKk/fk);
cc = 1/beta;

gk1 = (-cb + sqrt(cb^2 - 4*ca*cc))/2/ca; 
gk2 = (-cb - sqrt(cb^2 - 4*ca*cc))/2/ca;
if gk1 >1 
    gk = gk2;
else
    gk =gk1;
end

%find g_a conditional on g_k and E_a = 0;
fa = fbar - (1-delta)*exp(kbar);
fKa = fK - (1-delta);
ga1 = udprime*fa-beta*uprime*fKa*rho-beta*udprime*fa*rho*fK;
ga2 = udprime*exp(kbar)+beta*uprime*fKk+beta*udprime...
          *(fk - exp(kbar)*(rho+gk))*fK;
ga = ga1/ga2;

%display
gk
ga