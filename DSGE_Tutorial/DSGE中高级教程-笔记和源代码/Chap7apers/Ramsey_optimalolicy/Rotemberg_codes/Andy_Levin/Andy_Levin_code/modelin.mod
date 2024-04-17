/*
 * This file is prepared for Andy levin's code to produce Ramsey mod file
 * There are some comments lines to be added in order to produce correct mod file
 * // Endogenous variables
 * // Exogenous variables
 * // Monetary Policy Rule
 * period utility and welfare function must be declared as follows
 * More info., please refer to readme.m in Andy_levin directory;
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

[A_SS,h_SS,R_SS,pie_SS,C_SS] = steadys(pietarget,beta,phi,nu,epsil,chi);


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

steady;

stoch_simul(order=1,irf=15,nograph) pie h R C A;
