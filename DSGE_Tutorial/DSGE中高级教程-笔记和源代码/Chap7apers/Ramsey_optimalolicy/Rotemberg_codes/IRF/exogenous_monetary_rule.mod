/*
 * This file is prepared for exogenous monetary policy rule
 * All variables in levels
 */

// Endogenous variables
var R //nominal interest rate, not gross
h      // hours worked
pie   // gross inflation
C     // consumption
Util // period utility
Welf //welfare metric
;

// Exogenous variables
var A;

// Innovations
varexo eps_A;

// PARAMETERS
parameters nbeta chi beta epsilon phi rho alpha nu pietarget A_SS h_SS R_SS pie_SS C_SS;

beta = 0.99;
epsilon = 5;
phi = 100;
rho = 0.9;
alpha = 1.5;
nbeta = beta;

%recursively running mechansim
%nu=1/(epsilon-1);
%nu=0; %tax subsidy does not satisfiy: nu=1/(epsilon-1);
load parameterfile_irf;
set_param_value('nu',nu);

chi=1;
pietarget=1;

A_SS=1;

%from the intertemporal Euler equation:
R_SS=pietarget/beta-1;

%from the monetary policy rule:
pie_SS=pietarget;    

h_SS=(((1+nu)*(epsilon-1)+phi*(pietarget-1)*pietarget*(1-beta)/(1+(phi/2)*(pietarget-1)^2)))/(epsilon*chi);

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
((1+nu)*(1-epsilon)+epsilon*(chi*h*C/A))*(1+phi*(pie-1)^2/2)-phi*(pie-1)*pie
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
A=A_SS; 
Util=log(C)-chi*h^2/2;
Welf=Util/(1-nbeta);
end;

shocks;
var eps_A;
stderr .01;
end;

//steady;

stoch_simul(order=1,irf=40,nograph) pie h R C A Util Welf;
