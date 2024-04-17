/*
 * This mod file implements the baseline Classical Monetary Economy model of Jordi Gal?(2008): Monetary Policy, Inflation,
 * and the Business Cycle, Princeton University Press, Chapter 2
 *
 * Note that this mod-file implements the log-level of non-linear first order conditions 
 * and that the IRFs show the percentage deviations (not linear deviation) from steady state.
 *
 * It demonstrate the neutrality of money by showing that real variables do not move after a monetary policy shock
 *
 * This implementation was written by Xiangyang Li. 2016-3@Shanghai
 *
 */

var C           //Consumption
    w      //Real Wage
    pi         //inflation
    A           //AR(1) technology process
    N           //Hours worked
    R           //Nominal Interest Rate
    r //Real Interest Rate
    Y           //Output
    m_growth_ann;//money growth
varexo eps_A    //technology shock
       eps_m;   //monetary policy shock

parameters alpha beta rho sigma phi phi_pi eta;
parameters Cs ws pis As Ns Rs rs Ys m_growth_anns;
%----------------------------------------------------------------
% Follows parametrization of Chapter 3, p. 52
%----------------------------------------------------------------

alpha=0.33; //capital share
beta=0.99;  //discount factor
rho=0.9;     //autocorrelation technology shock
sigma=1;    //log utility   
phi=1;       //unitary Frisch elasticity
phi_pi=1.5;  //inflation feedback Taylor Rule
eta  =4;     // semi-elasticity of money demand

As=1;
Rs=1/beta;
pis=1;
rs=Rs;
Ns=(1-alpha)^(1/((1-sigma)*alpha+phi+sigma));
Cs=As*(Ns^(1-alpha));
ws=((1-alpha)*As*Ns^(-alpha));
Ys=Cs;
m_growth_anns=0;

%----------------------------------------------------------------
% First Order Conditions
%----------------------------------------------------------------

model;
//1. labor demand, eq. (6) in Chap 2;
exp(w)=exp(sigma*C+phi*N);

//2. Euler equation eq. (7) in Chap 2;
1/exp(R)=beta*(exp(-sigma*(C(+1)- C)) - pi(+1));

//3. Production function eq. (8)
exp(A+(1-alpha)*N)= exp(C);

//4. Wage rate, eq. (13)
exp(w)=(1-alpha)*exp(A-alpha*N);

//5. Definition Real interest rate
exp(r)=exp(R-pi(+1));

//6. Monetary Policy Rule, eq. (22,24)
exp(R)=1/beta+phi_pi*(pi)+eps_m;

//7. Market Clearing, eq. (15)
exp(C)=exp(Y);

//8. Technology Shock
A=rho*A(-1)+eps_A;

//9. Nominal Money growth (derived from eq. (10))
// m - p = y -eta*R
m_growth_ann=4*(Y-Y(-1)-eta*(R-R(-1))+pi);

end;

%----------------------------------------------------------------
%  define shock variances
%---------------------------------------------------------------

shocks;
var eps_A; stderr 1;
var eps_m; stderr 1;
end;

%----------------------------------------------------------------
%  Initial Values for steady state
%---------------------------------------------------------------

initval;
A=log(As);
R=log(Rs);
pi=log(pis);
r=log(rs);
N=log(Ns);
C=log(Cs);
w=log(ws);
Y=log(Ys);
m_growth_ann=m_growth_anns;
end;


resid(1);
steady;
check;

%----------------------------------------------------------------
% generate IRFs to show neutrality of money
%----------------------------------------------------------------

stoch_simul(order=1) Y C pi R r m_growth_ann;