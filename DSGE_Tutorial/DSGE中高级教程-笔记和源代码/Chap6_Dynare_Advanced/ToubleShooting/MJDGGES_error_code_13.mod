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
;		
varexo  eps_a eps_y_star;

parameters beta sigma alpha eta epsilon phi theta rho_a rho_y_star phi_pi phi_y
kappa_a	omega sigma_a lambda  BigGamma_a BigGamma_star BigTheta;


// Calibrations as per p.174
beta = 0.99;									// Pure temporal discount factor
sigma = 1;										// Intertemporal consumption elastiticy
alpha = 0.4;									// Degree of 'openness' in the Home economy
eta = 1;										// Elast. of sub. between Home and Foreign goods
epsilon = 6;									// Dixit-Stiglitz parameter for within-sector consumption 
phi = 3;										// Labour disutility parameter 
theta = 0.75;									// Calvo probability
gamma =1; 
lambda = (1-(beta*theta))*(1-theta)/theta;		// Coefficient on marginal cost in the Phillips Curve
omega = sigma*gamma+(1-alpha)*(sigma*eta-1); 
sigma_a = sigma/(1+alpha*(omega-1));
kappa_a = lambda*(phi +sigma_a);		// Real rigidity; see (37)
BigTheta = omega - 1;
BigGamma_a = (1+phi)/(sigma_a + phi);  //See(36);
BigGamma_star = -alpha*BigTheta*sigma_a/(sigma_a + phi); //See (36);

// Parameters of the productivity shocks (p.174)
rho_a = 0.66;
rho_y_star = 0.86;
phi_pi = 1.5;
phi_y = 0;

model(linear);
//(1) Home CPI inflation eq(15), P155
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
	//pi = 0;						// Strict inflation targeting (CIT)
//			e = 0;							// Exchange rate peg (PEG)
pi_h = 0;						// Domestic inflation targeting (DIT)
//				R = 0.5*R(-1) + phi_pi*pi + phi_y*ygap;		// Simple Taylor rule with coefficients plucked out of the air

//(15)Home technology shock
a = rho_a*a(-1) + eps_a;	

//(16) Foreign output 
y_star = rho_y_star*y_star(-1) + eps_y_star;

//(17)Euler condition
y = y(+1) -1/sigma*(R - pi(+1));

end;

shocks;
var eps_a; stderr 1; //0.0071;
var eps_y_star; stderr 1; //0.0078;
end;

stoch_simul(irf=16)	 ygap pi_h R	pi	e	q	p_h	cpi_level;
