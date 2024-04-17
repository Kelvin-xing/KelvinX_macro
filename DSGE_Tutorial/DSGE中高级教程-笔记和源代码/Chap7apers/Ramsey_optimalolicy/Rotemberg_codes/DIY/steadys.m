function [A_SS,h_SS,R_SS,pie_SS,C_SS] = steadys(pietarget,beta,phi,nu,epsil,chi)


A_SS=1;

%from the intertemporal Euler equation:
R_SS=pietarget/beta-1;

%from the monetary policy rule:
pie_SS=pietarget;    


h_SS=(((1+nu)*(epsil-1)+phi*(pietarget-1)*pietarget*(1-beta)/(1+(phi/2)*(pietarget-1)^2)))/(epsil*chi);

if h_SS < 0
    error('fatal (steadys) complex steady state hours!')
end

h_SS=sqrt(h_SS);

C_SS=h_SS/(1+(phi/2)*(pie_SS-1)^2);

err(1)=(1/(1+R_SS))-beta/pie_SS;

err(2)=((1+nu) *(1-epsil)+epsil*chi*h_SS*C_SS/A_SS)*(1+phi*(pie_SS-1)^2/2) ...
   + beta*(pie_SS-1)*pie_SS*phi-phi*(pie_SS-1)*pie_SS;

err(3)=C_SS*(1+(phi/2)*(pie_SS-1)^2)-h_SS;
    
if max(abs(err)) > .1e-10
    error('fatal (steadys) state state not satisfied')
end