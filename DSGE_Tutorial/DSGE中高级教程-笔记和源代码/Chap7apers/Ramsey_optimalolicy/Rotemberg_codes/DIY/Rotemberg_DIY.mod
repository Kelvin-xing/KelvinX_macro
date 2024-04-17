//Written By Xiangyang Li 2015-10-21@BJ
// Endogenous variables 
var R h pie C A Util Welf lmult1 lmult2 lmult3; 

// Innovations 
varexo eps_A; 
   
// PARAMETERS 
parameters nbeta chi beta epsil phi rho alpha nu pietarget; 
parameters as rs pis hs cs Utils;
parameters lmult1_SS  lmult2_SS  lmult3_SS; 
beta = 0.99; 
epsil = 5; 
phi = 100; 
rho = 0.9; 
alpha = 1.5; 
nbeta = beta; 
//nu=1/(epsil-1); //lmult2 will equal to zero
nu =0;    //lmult2 will not equal to zero
chi=1; 
pietarget=1.; 
 
%technology steady state
as=1;
%nominal rate ss from the intertemporal Euler equation
rs=pietarget/beta-1;

%from the monetary policy rule:
pis=pietarget;    

hs=(((1+nu)*(epsil-1)+phi*(pietarget-1)*pietarget*(1-beta)/
(1+(phi/2)*(pietarget-1)^2)))/(epsil*chi);

hs=sqrt(hs);
cs=hs/(1+(phi/2)*(pis-1)^2);
Utils = log(cs) -chi*hs^2/2;

% inital values, actually ss under the current parameterization;
lmult1_SS = 0;
lmult2_SS = 0;
lmult3_SS  = -1; 

model; 
//(1) period utility 
Util=log(C)-chi*h^2/2; 

//(2) policymakers' social welfare 
Welf=Util + nbeta*Welf(+1); 
 
// Ramsey optimal FOCs
// produced by Andy Levin's codes;
//(3) FOCs w.r.t nomial rate;
-lmult1/(R + 1)^2=0; 

//(4) FOCs w.r.t consumption;
lmult3*((phi*(pie - 1)^2)/2 + 1) + 1/C - 
(beta*lmult1)/(C(+1)*pie(+1)) + (C(-1)*beta*lmult1(-1))/(C^2*nbeta*pie) +
 (chi*epsil*h*lmult2*((phi*(pie - 1)^2)/2 + 1))/A; 

//(5) FOCs w.r.t. labor;
(C*chi*epsil*lmult2*((phi*(pie - 1)^2)/2 + 1))/A - chi*h - A*lmult3; 

//(6) FOCs w.r.t. inflation
(lmult2(-1)*(beta*phi*(pie - 1) + beta*phi*pie))/nbeta 
- lmult2*(phi*pie + phi*(pie - 1) + (phi*(2*pie - 2)*((epsil - 1)*(nu + 1) 
- (C*chi*epsil*h)/A))/2) + (C*lmult3*phi*(2*pie - 2))/2 
+ (C(-1)*beta*lmult1(-1))/(C*nbeta*pie^2); 
 
//(7)Original FOCs - Home Euler Equation
1/(1+R) = beta*C/(C(+1)*pie(+1)); 
 
//(8)Original FOCs -new Keynesian Phillips Curve 
((1+nu)*(1-epsil)+epsil*(chi*h*C/A))*(1+phi*(pie-1)^2/2)-phi*(pie-1)*pie
   + beta*phi*(pie(+1)-1)*pie(+1);
 
//(9)Original FOCs -resource constraint 
C*(1+(phi/2)*(pie-1)^2)-A*h; 
 
//(10)Technology Shocks 
log(A) = rho*log(A(-1))+eps_A; 
 
end; 
 
initval; 
R=rs;  
h=hs;  
pie=pis;  
C=cs; 
A=1;  
lmult1 = lmult1_SS; 
lmult2 = lmult2_SS; 
lmult3 = lmult3_SS; 
Util=log(cs)-chi*hs^2/2;  
Welf=log(cs)-chi*hs^2/2/(1-nbeta);
end; 
 
shocks; 
var eps_A; 
stderr .01; 
end; 
 
%steady; 
 
stoch_simul(order=1,irf=20,periods=200) pie h R C A lmult1 lmult2 lmult3 Welf Util; 
