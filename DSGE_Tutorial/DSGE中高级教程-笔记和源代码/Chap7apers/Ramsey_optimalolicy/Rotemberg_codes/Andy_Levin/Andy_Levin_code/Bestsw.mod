// LOWW Project
// March 2005

//DECLARATIONS.

// Variables of the Economy with frictions
// Endogenous variables
var lY, lN, lK, lC, linvest, lwr, lrk, Welf, Util, vbleq, lx0, lDisp, z1, z2,
lInt, lRint, linfl, lMC, lwx0, lDispw, zw1, zw2, lwinfl, lpindexa, lwindexa, 
au, logu, lCtil, lUc, lYgap, lCtilpm, lUcpm;

// Exogenous variables
var lA, lB, lG, lI, lL, letaP, letaQ, letaW, lTaxp, lTaxw, 
    zepsA, zepsB, zepsG, zepsI, zepsL, zepsTp, zepsTw, zepsP, zepsQ, zepsW, 
    zmpsA, zmpsB, zmpsG, zmpsI, zmpsL, zmpsTp, zmpsTw, zmpsP, zmpsQ, zmpsW,
    zlA, zlB, zlG, zlI, zlL, zlTp, zlTw, zlP, zlQ, zlW;

// Variables of the Frictionless Economy (Exogenous variables)
var lYf, lNf, lKf, lCf, linvestf, lwrf, lrkf, Welff, Utilf, vbleqf, 
auf, loguf, lCtilf, lUcf, lIntf, lRintf, lYgapf, linflf,lwinflf;

// Stochastic Exogenous Shocks

varexo epsA, epsB, epsG, epsI, epsL, etaP, etaQ, etaW, epsTp, epsTw;

//PARAMETERS

parameters alph, bet, tau, chi, chi0, sig, epp, epw, nkappa, xip, nbeta, 
xiw, gamp, gamw, theta, mps_flag, neta, 
nrhoa, nrhob, nrhog, nrhoi, nrhol, nrhoTp, nrhoTw, 
nirc, nuflag, na0, CapPhi, exthabit,
lTaxpBar, lTaxwBar, inflBar, winflBar, intBar, 
ABar, BBar, GBar, IBar, LBar, NBar, YBar,
mpr1, mpp0, mpp1, mpdp0, mpw0, mpygap0, mpygap1, mpdygap0, mpy0, mpdy0, mpdi0, mppf1,mpest,
sigpm, chipm, thetapm, netapm, xippm, xiwpm, gamppm, gamwpm, nkappapm, chipm0;

model;           

// Model Equations in the Economy with Frictions

// Welfare

Util = exp(lB)*((exp(lCtil)^(1-sig))/(1-sig) - exp(lDispw)*chi0*exp(lL)*(((exp(lN))^(1+chi))/(1+chi)));   
Welf = Util + bet*Welf(+1);

// Marginal Utility of Consumption 

exp(lCtil) = exp(lC) - theta*exp(lC(-1));
exp(lUc) = exp(lB)*((exp(lCtil))^(-sig)) - (1-exthabit)*exp(lB(+1))*bet*theta*((exp(lCtil(+1)))^(-sig));


// Marginal Utility of Consumption 

exp(lCtilpm) = exp(lC) - thetapm*exp(lC(-1));
exp(lUcpm) = exp(lB)*((exp(lCtilpm))^(-sigpm)) - (1-exthabit)*exp(lB(+1))*bet*thetapm*((exp(lCtilpm(+1)))^(-sigpm));


// Indexation Rules (checked)

lwindexa = gamw*linfl(-1);
lpindexa = gamp*linfl(-1);

// Variable Utilization Function (checked; na0=RK)

logu = log((1-nuflag) + nuflag*(exp(lrk)/na0)^(1/neta));
au = na0*((exp(logu)^(1+neta) - 1)/(1+neta));

//Consumers FOCs (checked)

exp(lUc) = bet*exp(lRint)*exp(lUc(+1));
exp(lUc) = bet*(exp(lInt)/exp(linfl(+1)))*exp(lUc(+1));
vbleq = exp(letaQ)*((bet*(exp(lUc(+1))/exp(lUc)))*(exp(lrk(+1))*(1 + nuflag*(exp(logu(+1))-1)) - nuflag*au(+1) + (1 - tau)*vbleq(+1)));
1 = (vbleq*(1 - (3/2)*nkappa*(((exp(lI+linvest)/exp(linvest(-1)))-1)^(2))- nkappa*((exp(lI+linvest)/exp(linvest(-1)))-1) ))
    + ((bet*(exp(lUc(+1))/exp(lUc)))*vbleq(+1)*(nkappa*((exp(lI(+1)+linvest(+1))/exp(linvest))-1)^(3) +
    2*nkappa*((exp(lI(+1)+linvest(+1))/exp(linvest))-1)^(2) + nkappa*((exp(lI(+1)+linvest(+1))/exp(linvest))-1) ) );
  
// Capital Accumulation with Adjustment Costs (check)

exp(lK) = (1-tau)*exp(lK(-1)) + exp(linvest)*(1 - 0.5*nkappa*(((exp(lI+linvest)/exp(linvest(-1)))-1)^2));

//Production (check)

lY = log(exp(-lDisp)*exp(lA)*((exp(lK(-1))*nuflag*exp(logu))^alph)*(exp(lN)^(1-alph))-CapPhi);               

// Price stickiness and Wage Stickiness

exp(lDisp) = (1-xip)*exp(lx0)^(-epp) + xip*exp(lDisp(-1))*((exp(linfl)/exp(lpindexa))^(epp));
lx0 = log(z2/z1) + log(epp/(epp-1));
z1 =  exp(lY)*(exp(lUc))/exp(lTaxp) + bet*xip*(exp(lpindexa(+1))*(exp(linfl(+1)-lpindexa(+1)))^(epp))*z1(+1);
z2 =  exp(lMC+letaP)*(exp(lY))*(exp(lUc)) +  bet*xip*(exp(lpindexa(+1))*(exp(linfl(+1)-lpindexa(+1)))^(1+epp))*z2(+1);
exp(linfl)  = exp(lpindexa)*(xip^(-1/(epp-1)))*(1-(1-xip)*(exp(lx0)^(1-epp)))^(1/(epp-1));   

exp(lDispw) = (1-xiw)*exp(lwx0)^(-epw*(1+chi)) + xiw*exp(lDispw(-1))*(exp(lwinfl)/exp(lwindexa))^(epw*(1+chi));
lwx0 = (log(zw2/zw1) + log(epw/(epw-1)))/(1+chi*epw);
zw1 =  exp(lN)*exp(lwr)*exp(lUc)/exp(lTaxw) + bet*xiw*((exp(lwinfl(+1)-lwindexa(+1)))^(epw-1))*zw1(+1) ;
zw2 = chi0*exp(lL+lB+letaW)*((exp(lN))^(1+chi)) +  bet*xiw*((exp(lwinfl(+1)-lwindexa(+1)))^(epw*(1+chi)))*zw2(+1);
exp(lwinfl)  = exp(lwindexa)*(xiw^(-1/(epw-1)))*(1-(1-xiw)*(exp(lwx0)^(1-epw)))^(1/(epw-1));   

// Equilibrium Capital-Labor Ratio (check)

lwr + nirc*lInt - lrk = log((1-alph)/alph) + nuflag*logu + lK(-1) - lN;

// Real Marginal Cost (check)

lMC = alph*(lrk - log(alph)) + (1-alph)*(lwr + nirc*lInt - log(1-alph)) - lA;

//GDP identity and government  (check)

linvest = log(exp(lY) - exp(lC) - nuflag*au*exp(lK(-1)) - exp(lG));

// Real Wages Growth (check)

lwinfl - linfl = lwr - lwr(-1);

// Equations in the Frictionless Economy

//Welfare

Utilf = exp(lB)*((exp(lCtilf)^(1-sig))/(1-sig) - chi0*exp(lL)*(((exp(lNf))^(1+chi))/(1+chi)));   
exp(lCtilf) = exp(lCf) - theta*exp(lCf(-1));
Welff = Utilf + bet*Welff(+1);

// Marginal Utility of Consumption

exp(lUcf) = exp(lB)*((exp(lCtilf))^(-sig)) - (1-exthabit)*exp(lB(+1))*bet*theta*((exp(lCtilf(+1)))^(-sig));

// Real rate
exp(lUcf) = bet*exp(lRintf)*exp(lUcf(+1));
exp(lUcf) = bet*(exp(lIntf)/inflBar)*exp(lUcf(+1));

// Variable Utilization Function

loguf = log((1-nuflag) + nuflag*(exp(lrkf)/na0)^(1/neta));
auf = na0*((exp(loguf)^(1+neta) - 1)/(1+neta));

//Consumers FOCs - Frictionless Economy

vbleqf = (bet*(exp(lUcf(+1))/exp(lUcf)))*(exp(lrkf(+1))*(1 + nuflag*(exp(loguf(+1))-1)) - nuflag*auf(+1) + (1 - tau)*vbleqf(+1));
1 = (vbleqf*(1 - (3/2)*nkappa*(((exp(lI+linvestf)/exp(linvestf(-1)))-1)^(2))- nkappa*((exp(lI+linvestf)/exp(linvestf(-1)))-1) ))
    + ((bet*(exp(lUcf(+1))/exp(lUcf)))*vbleqf(+1)*(nkappa*((exp(lI(+1)+linvestf(+1))/exp(linvestf))-1)^(3) +
    2*nkappa*((exp(lI(+1)+linvestf(+1))/exp(linvestf))-1)^(2) + nkappa*((exp(lI(+1)+linvestf(+1))/exp(linvestf))-1) ) );
  
// Capital Accumulation with Adjustment Costs

exp(lKf) = (1-tau)*exp(lKf(-1)) + exp(linvestf)*(1 - 0.5*nkappa*(((exp(lI+linvestf)/exp(linvestf(-1)))-1)^2));


//Production Frictionless Economy

lYf = log(exp(lA)*((exp(lKf(-1))*nuflag*exp(loguf))^alph)*(exp(lNf)^(1-alph))-CapPhi);               

// Input Prices Frictionless
 
lwrf + lUcf -lB = log(chi0*exp(lL)) + chi*lNf + log(epw/(epw-1)) + lTaxw;
 
lTaxp + log(epp/(epp-1)) + lrkf = log(alph) + log(exp(lYf)+1*CapPhi) - lKf(-1)- loguf;
lTaxp + log(epp/(epp-1)) + lwrf = log(1-alph) + log(exp(lYf)+1*CapPhi) - (lNf + nirc*log(intBar));

//GDP identity and government 
linvestf = log(exp(lYf) - exp(lCf) - nuflag*auf*exp(lKf(-1)) - exp(lG));

// Output Gap

lYgap = lY - lYf;

// Output Gap

lYgapf = lYf - lYf;

linflf = 0*lYgapf;

lwinflf = lwrf - lwrf(-1);

// Monetary Policy Rule

lInt = (1-mpest)*((1-mpr1)*(linfl -log(bet)) 
       +  mpr1*lInt(-1) 
       +  mpp0*(linfl - log(inflBar)) 
       +  mppf1*(linfl(+1) - log(inflBar)) 
       +  mpw0*(lwinfl - log(winflBar))
       +  mpygap0*lYgap 
       +  mpygap1*lYgap(-1) 
       +  mpy0*(lY - log(YBar))
       +  mpdy0*(lY-lY(-1))
       +  mpdi0*(linvest-linvest(-1)))
      + mpest*(mpr1*lInt(-1)
        + (1-mpr1)*(log(inflBar)-log(bet))
        + (1-mpr1)*(mpp1*(linfl(-1)-log(inflBar))+mpygap1*lYgap(-1))
        + mpdp0*(linfl-linfl(-1))
        + mpdy0*(lYgap-lYgap(-1))                 );


// Exogenous Disturbances
// Auxiliary variables used for determining mean-preserving spreads

zmpsA = (nrhoa^2)*zmpsA(-1) + 0.5*mps_flag*zepsA(+1)^2;
zmpsB = (nrhob^2)*zmpsB(-1) + 0.5*mps_flag*zepsB(+1)^2;
zmpsG = (nrhog^2)*zmpsG(-1) + 0.5*mps_flag*zepsG(+1)^2;
zmpsI = (nrhoi^2)*zmpsI(-1) + 0.5*mps_flag*zepsI(+1)^2;
zmpsL = (nrhol^2)*zmpsL(-1) + 0.5*mps_flag*zepsL(+1)^2;
zmpsTp = (nrhoTp^2)*zmpsTp(-1) + 0.5*mps_flag*zepsTp(+1)^2;
zmpsTw = (nrhoTw^2)*zmpsTw(-1) + 0.5*mps_flag*zepsTw(+1)^2;

zmpsP = (0)*zmpsP(-1) + 0.5*mps_flag*zepsP(+1)^2;
zmpsQ = (0)*zmpsQ(-1) + 0.5*mps_flag*zepsQ(+1)^2;
zmpsW = (0)*zmpsW(-1) + 0.5*mps_flag*zepsW(+1)^2;

zlA =  log(ABar) + nrhoa*(zlA(-1) - log(ABar)) + epsA ;
zlB =  log(BBar) + nrhob*(zlB(-1) - log(BBar)) + epsB ;
zlG =  log(GBar) + nrhog*(zlG(-1) - log(GBar)) + epsG ;
zlI =  log(IBar) + nrhoi*(zlI(-1) - log(IBar)) + epsI ;
zlL =  log(LBar) + nrhol*(zlL(-1) - log(LBar)) + epsL ;
zlTp = lTaxpBar + nrhoTp*(zlTp(-1)-lTaxpBar) + epsTp;
zlTw = lTaxwBar + nrhoTw*(zlTw(-1)-lTaxwBar) + epsTw;

zlP = etaP;
zlQ = etaQ;
zlW = etaW;

lA = zlA - zmpsA;
lB = zlB - zmpsB;
lG = zlG - zmpsG;
lI = zlI - zmpsI;
lL = zlL - zmpsL;
lTaxp = zlTp - zmpsTp;
lTaxw = zlTw - zmpsTw;

letaP = zlP - zmpsP;
letaQ = zlQ - zmpsQ;
letaW = zlW - zmpsW;

zepsA =  epsA; 
zepsB =  epsB;
zepsG =  epsG;
zepsI =  epsI;
zepsL =  epsL;
zepsTp = epsTp;
zepsTw = epsTw;

zepsP = etaP;
zepsQ = etaQ;
zepsW = etaW;

end;


//INITIAL VALUES.
initval;

lA = log(ABar);
lB = log(BBar);
lG = log(GBar);
lI = log(IBar);
lL = log(LBar);

lTaxp = lTaxpBar;
lTaxw = lTaxwBar;

lY = log(YBar);

linfl = log(inflBar);
lwinfl = log(inflBar);
vbleq = 1;
logu=0;
au = 0;

zepsA =  0; 
zepsB =  0;
zepsG =  0;
zepsI =  0;
zepsL =  0;
zepsTp = 0;
zepsTw = 0;
zepsP = 0;
zepsQ = 0;
zepsW = 0;

etaP=0;
letaP=0;
etaQ=0;
letaQ=0;
etaW=0;
letaW=0;

lwindexa = gamw*log(inflBar);
lpindexa = gamp*log(inflBar);

lx0 = log(((1 - xip*(inflBar^((gamp-1)*(1-epp))))/(1-xip))^(1/(1-epp)));
lwx0 = log(((1 - xiw*(inflBar^((gamp-1)*(1-epw))))/(1-xiw))^(1/(1-epw)));
lDisp = log(((1-xip)*(exp(lx0)^(-epp)))/(1-xip*inflBar^((1-gamp)*epp)));
lDispw = log(((1-xiw)*(exp(lwx0)^(-epw*(1+chi))))/(1-xiw*inflBar^((1-gamw)*epp*(1+chi))));

lMC = lx0 + log(1-bet*xip*inflBar^(1+epp)) + log((epp-1)/epp) - (lTaxp + log(1-bet*xip*inflBar^(epp+gamp)));
lrk = log((bet^(-1)) - (1.0 - tau));

lK = log(alph) + lY + lMC - lrk;
lN = (lY - lA - alph*(lK) + lDisp)/(1-alph);

lG = log(GBar);
lC = lY + log(1 - (tau*alph/exp(lrk))*exp(lMC) - exp(lG-lY) - au*exp(lK));
lCtil = log((1-theta)*exp(lC));
lUc = log((exp(lCtil)^(-sig))*(1 - (1-exthabit)*bet*theta));
lwr = log((((1-alph)*exp(lY-lN)*exp(lMC))/(intBar)^(nirc)));
z2 =  exp(lMC)*exp(lY)*(exp(lUc))/(1-bet*xip*(inflBar^(1+epp)));
z1 = z2 *(epp/(epp-1))*exp(lTaxp)/exp(lx0);

zw1 = exp(lN)*exp(lwr)*(exp(lUc))/(1-bet*xiw*exp(lwindexa)*(inflBar^(epw-1))) ;
zw2 = chi0*exp(lL)*((exp(lN))^(1+chi))/(1- bet*xiw*(inflBar^(epw*(1+chi))));
linvest = lK + log(tau);
lRint = -log(bet);
lInt = log(inflBar) - log(bet); 
Util = (exp(lCtil)^(1-sig))/(1-sig) - exp(lDispw)*chi0*(((exp(lN))^(1+chi))/(1+chi));   
Welf = Util/(1-bet);


lCtilpm = lCtil;
lUcpm = lUc;

//periodloss=0;
//adhocloss=0;

// Frictionless Economy

lYf = log(YBar);
vbleqf = 1;
loguf=0;
auf = 0;
lrkf = log((bet^(-1)) - (1.0 - tau));

lKf = log(alph*(exp(lYf)*exp(lMC))/exp(lrkf));
lNf = (lYf - lA - alph*(lKf))/(1-alph);

lCf = lYf + log(1 - (tau*alph/exp(lrkf))*exp(lMC) - exp(lG-lYf) - auf*exp(lKf));
lCtilf = log((1-theta)*exp(lCf));
lUcf = log((exp(lCtilf)^(-sig))*(1 - (1-exthabit)*bet*theta));
lwrf = log((((1-alph)*exp(lYf-lNf)*exp(lMC))/(intBar)^(nirc)));
linvestf = lKf + log(tau);
lRint = -log(bet);
Utilf = (exp(lCtilf)^(1-sig))/(1-sig) - chi0*(((exp(lNf))^(1+chi))/(1+chi));   
Welff = Utilf/(1-bet);


lIntf=lInt;
lRintf=lRint;
lYgap = lY - lYf;
lYgapf = 0;

linflf = log(inflBar);
lwinflf = log(inflBar);

epsA=0;
epsB=0;
epsG=0;
epsI=0;
epsL=0;
epsTp=0;
epsTw=0;

zlA=0; 
zlB=0; 
zlG=0; 
zlI=0; 
zlL=0; 
zlTp=0; 
zlTw=0; 

zlP=0; 
zlQ=0; 
zlW=0; 

zmpsA=0; 
zmpsB=0; 
zmpsG=0; 
zmpsI=0; 
zmpsL=0; 
zmpsTp=0; 
zmpsTw=0;
zmpsP=0;
zmpsQ=0;
zmpsW=0;

end;

