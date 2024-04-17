%this file is written by Xiangyang Li @2015.8.8
%
var lam C R z i u mu inv wstar w x1 x2 Y A N vp pi pistar h1 h2 mc K Kbar G q omegag;
varexo eg ei ea ez;
parameters beta alpha delta b tau eta omegags psi epsw epsp phiw phip chi1 chi2;
parameters rhog rhoi rhoa rhoz sdg sdi sda sdz phipi phiy;
parameters lams Cs Rs zs is us mus invs wstars ws x1s x2s Ys As Ns;
parameters vps pis pistars h1s h2s mcs Ks Kbars Gs qs;

load mediumdsge
set_param_value('beta',beta);
set_param_value('alpha',alpha);
set_param_value('delta',delta);
set_param_value('b',b);
set_param_value('tau',tau);
set_param_value('eta',eta);
set_param_value('omegags',omegags);
set_param_value('psi',psi);
set_param_value('epsw',epsw);
set_param_value('epsp',epsp);
set_param_value('phip',phip);
set_param_value('phiw',phiw);
set_param_value('chi1',chi1);
set_param_value('chi2',chi2);
set_param_value('rhog',rhog);
set_param_value('rhoa',rhoa);
set_param_value('rhoz',rhoz);
set_param_value('rhoi',rhoi);
set_param_value('sdg',sdg);
set_param_value('sda',sda);
set_param_value('sdz',sdz);
set_param_value('sdi',sdi);
set_param_value('phipi',phipi);
set_param_value('phiy',phiy);
set_param_value('lams',lams);
set_param_value('Cs',Cs);
set_param_value('Rs',Rs);
set_param_value('zs',zs);
set_param_value('is',is);
set_param_value('us',us);
set_param_value('mus',mus);
set_param_value('invs',invs);
set_param_value('wstars',wstars);
set_param_value('ws',ws);
set_param_value('x1s',x1s);
set_param_value('x2s',x2s);
set_param_value('Ys',Ys);
set_param_value('As',As);
set_param_value('Ns',Ns);
set_param_value('vps',vps);
set_param_value('pis',pis);
set_param_value('pistars',pistars);
set_param_value('h1s',h1s);
set_param_value('h2s',h2s);
set_param_value('mcs',mcs);
set_param_value('Ks',Ks);
set_param_value('Kbars',Kbars);
set_param_value('Gs',Gs);
set_param_value('qs',qs);


model;
# invt  = exp(inv)/exp(inv(-1)) -1;
# invtp= exp(inv(+1))/exp(inv) -1;
# invr  =exp(inv)/exp(inv(-1));
# invrp =exp(inv(+1))/exp(inv);
# rc     = chi1*(exp(u)-1) + chi2/2*(exp(u)-1)^2;
# rcp     = chi1*(exp(u(+1))-1) + chi2/2*(exp(u(+1))-1)^2;

%(1) home Euler equation 1
exp(lam) = 1/(exp(C) - b*exp(C(-1))) - beta*b/(exp(C(+1)) - b*exp(C));

%(2) utilization cost, already in ratio
R= (chi1 +chi2*(exp(u)-1))/exp(z);

%(3) home Euler equation 2
exp(lam) = beta*exp(lam(+1))*exp(i)/exp(pi(+1));

%(4)investment decision
exp(lam) = exp(mu)*exp(z)*(1-tau/2*invt^2 - tau*invt*invr)
                 +beta*exp(mu(+1))*exp(z(+1))*tau*invtp*invrp^2;

%(5)capital stock desicion
exp(mu) = beta*(exp(lam(+1)) *(R(+1)*exp(u(+1))- rcp/exp(z(+1))) 
                +exp(mu(+1))*(1-delta));

%(6) optimal wage decision
exp(wstar) = epsw/(epsw-1) * exp(h1 - h2);

%(7) auxiliary h1
exp(h1) = psi*exp(epsw*(1+eta)*(w-wstar))*exp((1+eta)*N)
                +phiw*beta*exp(epsw*(1+eta)*pi(+1))*exp(epsw*(1+eta)*(wstar(+1)-wstar))*exp(h1(+1));

%(8) auxiliary h2
exp(h2) = exp(lam)*exp(epsw*(w-wstar))*exp(N) 
                +phiw*beta*exp((epsw-1)*pi(+1))*exp(epsw*(wstar(+1)-wstar))*exp(h2(+1));

%(9) wage index
exp((1-epsw)*w) = (1-phiw)*exp((1-epsw)*wstar)+exp((epsw-1)*pi)*phiw*exp((1-epsw)*w(-1));

%(10) production technology
exp(Y) = exp(A)*exp(alpha*Kbar)*exp((1-alpha)*N)/exp(vp);

%(11) the price dispersion
exp(vp) = exp(epsw*pi)*((1-phip)*exp(-epsw*pistar)+phip*exp(vp(-1)));

%(12) the CPI 
exp((1-epsp)*pi) = (1-phip)*exp((1-epsp)*pistar) +phip;

%(13) optimal price decision
exp(pistar) = epsp/(epsp-1)*exp(pi)*exp(x1-x2);

%(14) auxiliary variable x1
exp(x1) = exp(lam+mc+Y)+phip*beta*exp(epsp*pi(+1))*exp(x1(+1));

%(15) auxiliary variable x2
exp(x2) = exp(lam+Y)+phip*beta*exp((epsp-1)*pi(+1))*exp(x2(+1));

%(16)capital-labor ratio
exp(w)/R=(1-alpha)/alpha*exp(Kbar-N);

%(17)wage decision
exp(w) = exp(mc+A)*(1-alpha)*exp(alpha*(Kbar-N));

%(18)resource constraint
exp(Y) = exp(C) + exp(inv) + exp(G) + rc*exp(K(-1))/exp(z);

%(19)capital evolution 
exp(K) = exp(z)*(1-tau/2*invt^2)*exp(inv)+(1-delta)*exp(K(-1));

%(20)capital service
exp(Kbar) = exp(u)*exp(K(-1));

%(21) government spending
exp(G) = omegag*exp(Y);

%(22) law of government spending share,omegag could be zero or negative in simulation
omegag = (1-rhog)*omegags + rhog*omegag(-1) + eg;

%(23)technology shock
A = rhoa*A(-1) + ea;

%(24)investment-specific shock
z = rhoz*z(-1) + ez; 

%(25)Hayashi Q
exp(q) = exp(mu - lam);

%(26) Taylor rule
i = (1-rhoi)*log(is) + rhoi*i(-1) + (1-rhoi)*(phipi*(pi - log(pis))+phiy*(Y-log(Ys)))+ei;

end;

initval;
lam = log(lams);
C = log(Cs);
R =Rs;
z  =log(zs);
i   = log(is);
u  = log(us);
mu =log(mus);
inv  =log(invs);
wstar =log(wstars);
w =log(ws);
x1 =log(x1s);
x2 =log(x2s);
Y  =log(Ys);
A =log(As);
N =log(Ns);
vp =log(vps);
pi =log(pis);
pistar =log(pistars);
h1 =log(h1s);
h2 =log(h2s);
mc =log(mcs);
K   =log(Ks);
Kbar =log(Kbars);
G =log(Gs);
q =log(qs);
omegag = omegags;
end;

shocks;
var eg = sdg^2;
var ea = sda^2;
var ei  = sdi^2;
var ez = sdz^2;
end;

resid(1);
steady;
check;

stoch_simul(order=1)  Y C inv N i pi w u z A G;