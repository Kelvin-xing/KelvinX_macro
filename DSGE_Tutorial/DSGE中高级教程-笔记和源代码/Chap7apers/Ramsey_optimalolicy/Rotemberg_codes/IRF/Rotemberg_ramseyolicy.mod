//2015-10-21@Beijing
// 6 endogenous variables including period utility
// all variables in levels
var R //Nominal interest rate, not gross 
h       //labor
pie     //gross inflation
C     //consumption
A     //technology shock
Util   //period utility, planner's objective
;

// Innovations
varexo eps_A;

// PARAMETERS
parameters nbeta chi beta epsilon phi rho nu pietarget;
parameters as rs pis hs cs Utils;
beta = 0.99;
epsilon = 5;
phi = 100;
rho = 0.9;
nbeta = beta;

%recursively running mechansim
%nu=1/(epsilon-1);
%nu=0; %tax subsidy does not satisfiy: nu=1/(epsilon-1);
load parameterfile_irf;
set_param_value('nu',nu);

chi=1;
pietarget=1.;

%technology steady state
as=1;

%nominal rate ss from the intertemporal Euler equation
rs=pietarget/beta-1;

%the inflation target if we set it a monetary policy rule
pis=pietarget;    

%steady state labor
hs=((1+nu)*(epsilon-1)+phi*(pis-1)*pis*(1-beta)/
(1+(phi/2)*(pis-1)^2))/(epsilon*chi);
hs=sqrt(hs);

%consumption and utility s.s.
cs=hs/(1+(phi/2)*(pis-1)^2);
Utils = log(cs) -chi*hs^2/2;

model;
//(1) home Euler equation--intertemporal condition
1/(1+R) = beta*C/(C(+1)*pie(+1));

//(2) new Keynesian Phillips Curve
((1+nu)*(1-epsilon)+epsilon*(chi*h*C/A))*(1+phi*(pie-1)^2/2)-phi*(pie-1)*pie
   + beta*phi*(pie(+1)-1)*pie(+1);

//(3) resource constraint
C*(1+(phi/2)*(pie-1)^2)-A*h;

//(4) Technology Shocks
log(A) = rho*log(A(-1))+eps_A;

//(5) utlility function
Util=log(C)-chi*h^2/2;

end;

initval;
	Util=Utils;
	R=rs; 
	h=hs; 
	pie=pis; 
	C=cs;
	A=as; 
end;

shocks;
var eps_A;
stderr .01;
end;

//steady;

planner_objective Util;
ramsey_policy(order=1,irf=40, planner_discount=0.99,nograph);
