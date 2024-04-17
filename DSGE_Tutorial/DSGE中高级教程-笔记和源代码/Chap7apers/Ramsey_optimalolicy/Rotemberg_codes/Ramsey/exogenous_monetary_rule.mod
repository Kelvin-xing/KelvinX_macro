/*
 * This file is prepared for exogenous monetary policy rule
 */

// Endogenous variables
var R h pie C Util Welf;

// Exogenous variables
var A;

// Innovations
varexo eps_A;

// PARAMETERS
parameters nbeta chi beta epsil phi rho alpha nu pietarget A_SS h_SS R_SS pie_SS C_SS;

beta = 0.99;
epsil = 5;
phi = 100;
rho = 0.9;
alpha = 1.5;
nbeta = beta;
nu=1/(epsil-1);
chi=1;
pietarget=1;

%[A_SS,h_SS,R_SS,pie_SS,C_SS] = steadys(pietarget,beta,phi,nu,epsil,chi);
A_SS=1;

%from the intertemporal Euler equation:
R_SS=pietarget/beta-1;

%from the monetary policy rule:
pie_SS=pietarget;    


h_SS=(((1+nu)*(epsil-1)+phi*(pietarget-1)*pietarget*(1-beta)/(1+(phi/2)*(pietarget-1)^2)))/(epsil*chi);


h_SS=sqrt(h_SS);

C_SS=h_SS/(1+(phi/2)*(pie_SS-1)^2);


model;

Util=log(C)-chi*h^2/2;

Welf=Util + nbeta*Welf(+1);

// Monetary Policy Rule
R = (pietarget/beta)-1+alpha*(pie(+1)-pietarget);

//Euler
1/(1+R) = beta*C/(C(+1)*pie(+1));

//new Keynesian Phillips Curve
((1+nu)*(1-epsil)+epsil*(chi*h*C/A))*(1+phi*(pie-1)^2/2)-phi*(pie-1)*pie
   + beta*phi*(pie(+1)-1)*pie(+1);

// resource constraint
C*(1+(phi/2)*(pie-1)^2)-A*h;

//Technology Shocks
log(A) = rho*log(A(-1))+eps_A;


end;

initval;
R=R_SS; 
h=h_SS; 
pie=pie_SS; 
C=C_SS;
A=1; 
Util=log(C)-chi*h^2/2;
Welf=Util/(1-nbeta);
end;

shocks;
var eps_A;
stderr .01;
end;

//steady;

stoch_simul(order=1,irf=20) pie h R C A Util Welf;
