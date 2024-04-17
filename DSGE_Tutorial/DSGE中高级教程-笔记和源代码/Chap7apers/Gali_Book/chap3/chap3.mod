/*
 * This file implements the baseline New Keynesian model of Jordi Gali(2008):
 * Monetary Policy, Inflation,
 * and the Business Cycle, Princeton University Press, Chapter 3
 *
 * Note that all variables are expressed in log-deviations; 
 */

var pi          //inflation
    y_gap       //output gap
    y_nat       //natural output 
    y           //output
    r_nat       //natural interest rate
    r      //real interest rate
    R           //nominal interst rate
    n           //hours worked
    m_growth_ann//money growth
    nu          //AR(1) monetary policy shock process
    a           //AR(1) technology shock process
    r_ann  //annualized real interest rate
    R_ann       //annualized nominal interest rate
    r_nat_ann   //annualized natural interest rate
    pi_ann;     //annualized inflation rate

varexo eps_a    //technology shock
       eps_nu;  //monetary policy shock

parameters beta sigma psi_n_ya rho_nu rho_a phi; 
parameters phi_pi phi_y kappa alpha epsilon eta;

%----------------------------------------------------------------
% Parametrization, p. 52
%----------------------------------------------------------------
sigma = 1;      //log utility
phi=1;          //unitary Frisch elasticity
phi_pi = 1.5;   //inflation feedback Taylor Rule
phi_y  = .5/4;  //output feedback Taylor Rule
theta=2/3;      //Calvo parameter
//rho_nu = 0.9;   //high persistent monetary policy shock
rho_nu =0.5; // moderately persistent mon. pol. shock
rho_a  = 0.9;   //autocorrelation technology shock
beta = 0.99;    //discount factor
eta  =4;        // semi-elasticity of money demand
alpha=1/3;      //capital share
epsilon=6;      //demand elasticity

//Composite parameters
Omega=(1-alpha)/(1-alpha+alpha*epsilon);  //defined on page 47
psi_n_ya=(1+phi)/(sigma*(1-alpha)+phi+alpha); //defined on page 48
lambda=(1-theta)*(1-beta*theta)/theta*Omega; //defined on page 47
kappa=lambda*(sigma+(phi+alpha)/(1-alpha));  //defined on page 49

%----------------------------------------------------------------
% First Order Conditions
%----------------------------------------------------------------

model(linear); 
//1. New Keynesian Phillips Curve eq. (21)
pi=beta*pi(+1)+kappa*y_gap;

//2. Dynamic IS Curve eq. (22)
y_gap=-1/sigma*(R-pi(+1)-r_nat)+y_gap(+1);

//3. Interest Rate Rule eq. (25)
R=phi_pi*pi+phi_y*y_gap+nu;

//4. Definition natural rate of interest eq. (23)
r_nat=sigma*psi_n_ya*(a(+1)-a);

//5. Definition real interest rate
r=R-pi(+1);

//6. Definition of natural output
y_nat = psi_n_ya*a;

//7. Definition output gap
y_gap=y-y_nat;

//8. Monetary policy shock
nu=rho_nu*nu(-1)+eps_nu;

//9. TFP shock
a=rho_a*a(-1)+eps_a;

//10. Production function (eq. 13)
y=a+(1-alpha)*n;

//11. Nominal Money growth (derived from eq. (4))
m_growth_ann=4*(y-y(-1)-eta*(R-R(-1))+pi);

//12. Annualized nominal interest rate
R_ann=4*R;

//13. Annualized real interest rate
r_ann=4*r;

//14. Annualized natural interest rate
r_nat_ann=4*r_nat;

//15. Annualized inflation
pi_ann=4*pi;
end;

%----------------------------------------------------------------
%  define shock variances
%---------------------------------------------------------------

shocks;
var eps_nu = 0.25^2; //1 standard deviation shock of 25 basis points
end;

%----------------------------------------------------------------
%  steady states: all 0 due to linear model, no initval block anymore
%---------------------------------------------------------------
resid(1);
steady;
check;

%----------------------------------------------------------------
% generate IRFs, replicates Figures 3.1, p. 53
%----------------------------------------------------------------
stoch_simul(order = 1,irf=15) y_gap pi_ann R_ann r_ann m_growth_ann nu;


shocks;
var eps_nu = 0;   //shut off monetary policy shock
var eps_a  = 1^2; //unit shock to technology
end;

%----------------------------------------------------------------
% generate IRFs, replicates Figures 3.2, p. 55
%----------------------------------------------------------------
stoch_simul(order = 1,irf=15) y_gap pi_ann y n R_ann r_ann m_growth_ann a ;

