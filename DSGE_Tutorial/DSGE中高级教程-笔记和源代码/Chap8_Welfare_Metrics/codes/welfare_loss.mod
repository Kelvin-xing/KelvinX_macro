/*
 * This file tries to implements the Open Economy model of Jordi Gali(2008):
 * Monetary Policy, Inflation,
 * and the Business Cycle, Princeton University Press, Chapter 7
 * 
 * This file is written by Xiangyang Li@Shanghai, 2016-3-25
 * I have tried to use the model equilibrium conditions in textbook
 * but it seems that there are colinearity or indeterminacy problem.
 * Thus, I borrow some equlibrium conditions from 
 * Gali & Monacelli, NBER w.p. 8905, April 2002
 * Note that all model variables are expressed in log-deviations; 
 */
var y_star  //overseas output
  a   // domestic technology shock;
 ygap //domestic output gap
rnat  //natrual rate
R      //nominal interest rate;
y      //domestic output
ynat // natural level of output
pi    // CPI inflation
pi_h //domestic home infaltion
pi_star // foreign inflation
s	  //effective term of trade
q	//  effective real exchange rate
e	//  effective nominal exchange rate
p_h  //domestic price level
cpi_level // CPI price level
mc // marginal cost
nx // net export
r_star //foreign nominal interest rate
mc_star // foreign marginal cost
a_star // foreign technology shock;
n   //domestic labor 
w  //domestic real wage
;		
varexo  eps_a eps_y_star eps_a_star;

parameters beta sigma alpha eta epsilon phi theta rho_a rho_y_star 
      phi_pi phi_y phi_pi_star phi_a_star rho_a_star a_shock_correl
      kappa_a omega sigma_a lambda  BigGamma_a BigGamma_star BigTheta
     ;


// Calibrations as per p.174
beta = 0.99;									// Pure temporal discount factor
sigma = 1;										// Intertemporal consumption elastiticy
alpha = 0.4;									// Degree of 'openness' in the Home economy
eta = 1;										// Elast. of sub. between Home and Foreign goods
%epsilon = 6;									// Dixit-Stiglitz parameter for within-sector consumption 
%phi = 3;										// Labour disutility parameter 
theta = 0.75;									// Calvo probability
gamma =1; 

%mechansim on recursively running
load parameterfile_welfare;
set_param_value('epsilon',epsilon);
set_param_value('phi',phi);


// Coefficient on marginal cost in the Phillips Curve
lambda = (1-(beta*theta))*(1-theta)/theta;	
	
omega = sigma*gamma+(1-alpha)*(sigma*eta-1); 
sigma_a = sigma/(1+alpha*(omega-1));
kappa_a = lambda*(phi +sigma_a);		// Real rigidity; see eq(37) in textbook, henceforth
BigTheta = omega - 1;
BigGamma_a = (1+phi)/(sigma_a + phi);  //See eq(36);
BigGamma_star = -alpha*BigTheta*sigma_a/(sigma_a + phi); //See eq(36);

// Parameters of the productivity shocks (p.174 in textbook)
rho_a = 0.66;
rho_y_star = 0.86;
rho_a_star =0.9; //not present in textbook
phi_pi = 1.5;
phi_y = 0;

//// See (51) in the working paper and note 19;
phi_a_star = -(sigma*(1+phi)*(1-rho_a_star)) / (phi+sigma);
phi_pi_star = 1.01;
a_shock_correl = 0.3;



model(linear);
//(1) Home CPI inflation eq(15), P155 in textbook, the same afterward
pi = pi_h + alpha*(s - s(-1));

//(2) An identity to pin down the relative price of home goods,P155
p_h = p_h(-1) + pi_h;									

//(3)An identity to pin down the consumer price level
cpi_level = cpi_level(-1) + pi;		

//(4)Real exchange rate P156				
q = (1-alpha)*s;	

//(5) term of trade, eq(16), differenced version
s - s(-1) = e - e(-1) + pi_star - pi_h;	

//(6)Market clearing eq(29)			
y = y_star +s/sigma_a;		

//(7) Definition of Home output (p. 164)
y = ynat + ygap;										

//(8)Home's Phillips curve, eq(37)
pi_h = beta*pi_h(+1) + kappa_a*ygap;						

//(9)Home's IS curve, eq(38)
ygap = ygap(+1) - (1/sigma_a)*(R - pi_h(+1) - rnat);

//(10) Home's natural level output (39)
ynat = BigGamma_a*a + BigGamma_star*y_star;			

//(11)The definition of Home's Wicksellian interest rate, eq(39):										
rnat = -sigma_a*(1-rho_a)*BigGamma_a*a
		 - phi*BigGamma_star*(y_star(+1) - y_star);

//(12)The net export
nx = alpha*(omega/sigma-1)*s;

//(13)The home marginal cost
mc = (sigma_a + phi)*ygap;

//(14)The home monetary policy
//Home's monetary policy; 
//pi = 0;						// Strict inflation targeting (SCIT)
//R =  phi_pi*pi; //Domestic inflation Taylor rule (CIT)
e = 0;							// Exchange rate peg (PEG)
//pi_h = 0;						// Strict Domestic inflation targeting (SDIT)
//R = rnat +phi_pi*pi_h + phi_y*ygap;  //optimal policy, equi. to SDIT
// R = 0.5*R(-1) + phi_pi*pi + phi_y*ygap;		// Simple Taylor rule
//R = rnat; // indeterminacy problem arises, can not run;


//(15)Home technology shock
// you can turn off the correlation by setting a_shock_correl = 0;
a = rho_a*a(-1) + eps_a +a_shock_correl*eps_a_star;	

//(16)Foreign Euler condition,in working paper, eq(22):
y_star = y_star(+1) - (r_star - pi_star(+1))/sigma;	

//AR(1) for y_star in textbook(P174) like this also works
//y_star = rho_y_star*y_star(-1) + eps_y_star;

//(17) Foreign's marginal cost, ,in working paper, eq(32): 
mc_star = (sigma + phi)*y_star - (1+phi)*a_star;
					
// (18)Foreign's Phillips curve,in working paper, eq(31):
pi_star = beta*pi_star(+1) + lambda*mc_star;					

//(19)foreign interest rate rule, Taylor rule in working paper eq(51)
r_star = phi_pi_star*pi_star(+1) + phi_a_star*a_star;	

// (20)Foreign's technology process
a_star = rho_a_star*a_star(-1) + eps_a_star;		

//(21) domestic real wage,in P154, eq(8), log-linearized:
w = sigma*y +phi*n;

//(22) domestic production technology,in P162, eq(32):
y = a+n;
end;

//you can turn on all the shocks; I only turn on the technology shock here;
//using the standard deviation in P174, it is unrealistic to have 
// that sizes of IRFs in the textbook.
shocks;
var eps_a; stderr 0.01; 
var eps_y_star; stderr 0; 
var eps_a_star; stderr 0;
end;

stoch_simul(irf=16,nograph,noprint)	 ygap pi_h R	pi	e	q	p_h	cpi_level;
