
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
 
 
var lmult1,  lmult2,  lmult3; 
  
parameters lmult1_SS   lmult2_SS   lmult3_SS ; 
  lmult1_SS=0;
   lmult2_SS=0;
   lmult3_SS=-1;
model; 
 
Util=log(C)-chi*h^2/2; 
 
Welf=Util + nbeta*Welf(+1); 
 
// Policymaker's First-Order Conditions 
lmult3*((phi*(pie - 1)^2)/2 + 1) + 1/C - (beta*lmult1)/(C(+1)*pie(+1)) + (C(-1)*beta*lmult1(-1))/(C^2*nbeta*pie) + (chi*epsil*h*lmult2*((phi*(pie - 1)^2)/2 + 1))/A; 
(C*chi*epsil*lmult2*((phi*(pie - 1)^2)/2 + 1))/A - chi*h - A*lmult3; 
(lmult2(-1)*(beta*phi*(pie - 1) + beta*phi*pie))/nbeta - lmult2*(phi*pie + phi*(pie - 1) + (phi*(2*pie - 2)*((epsil - 1)*(nu + 1) - (C*chi*epsil*h)/A))/2) + (C*lmult3*phi*(2*pie - 2))/2 + (C(-1)*beta*lmult1(-1))/(C*nbeta*pie^2); 
-lmult1/(R + 1)^2; 
 
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
lmult1 = lmult1_SS; 
lmult2 = lmult2_SS; 
lmult3 = lmult3_SS; 
end; 
 
shocks; 
var eps_A; 
stderr .01; 
end; 
 
steady; 
 
stoch_simul(order=1,irf=15,nograph) pie h R C A; 
