clear;

% Parameters

alpha = 1/3;
beta = 1.05^(-1/4);
delta = 0.025;
gamma = 2;
mu = 0.5;             
phi = 0;  
rho = 1-beta;
Gamma = Inf;

Na = 300;

a_max = 500;        
a_min = phi;
agrd = linspace(a_min,a_max,Na)';

daf = agrd(2:Na)-agrd(1:Na-1);
dab = -(agrd(1:Na-1)-agrd(2:Na));

% Transition matrix in discrete time: g_{t+1}=T*g_{t}

T = [0.95, 1-0.95; 0.8, 1-0.8];

% Convert to continuous time such that dg_t=T*g_t

T = [-T(1,2) T(2,1);T(1,2) -T(2,1)];

% The steady state values of g satisfy 0=T*g. That is g is an eigenvector
% associated with a zero eigenvalue normalized to sum to one. But since T
% has a zero eigenvalue it is also singular. Find the steady state value of
% employment and unemployment.

% Compute the B matrix used in HJB equation (see notes)

B = [-eye(Na)*T(2,1),eye(Na)*T(2,1);eye(Na)*T(1,2),-eye(Na)*T(1,2)];

% Tax rate to balance budget: (1-e)*mu*w=e*tau*w.

tau = 0.05;

% Utility and consumption functions

u = @(c) (c.^(1-gamma))./(1-gamma);
up = @(c) c.^(-gamma);
c = @(dv) dv.^(-1/gamma);

% Matrices for derivatives

Df = -diag([ones(Na-1,1)./daf;1]);
Db = diag([1;ones(Na-1,1)./dab]);
Df(1:end-1,2:end)=Df(1:end-1,2:end)+diag(ones(Na-1,1)./daf);
Db(2:end,1:end-1)=Db(2:end,1:end-1)-diag(ones(Na-1,1)./dab);

r = 0.0115;

w = 2;

% Ok, problem set up. Let's solve it

% Guess for value function

ve = u((r*agrd+w*(1-tau)))./rho;
ve(end) = u(r*a_max+w*(1-tau))./rho;
vu = u((r*agrd+w*mu))./rho;
vu(1) = u(r*a_min+w*mu)./rho;

metric = 1;

tic

while metric>1e-8

% Calculate derivatives

dvef = Df*ve;
dveb = Db*ve;

dvuf = Df*vu;
dvub = Db*vu;

dvef(end) = up(w*(1-tau)+r*a_max);
dvuf(end) = up(w*mu+r*a_max);

dveb(1) = up(w*(1-tau)+r*a_min);
dvub(1) = up(w*mu+r*a_min);

cef = c(dvef);
cuf = c(dvuf);
ceb = c(dveb);
cub = c(dvub);

sef = r*agrd+w*(1-tau)-cef;
seb = r*agrd+w*(1-tau)-ceb;
suf = r*agrd+w*mu-cuf;
sub = r*agrd+w*mu-cub;

Ief = sef>0;
Ieb = seb<0;
Iuf = suf>0;
Iub = sub<0;

I0e = (1-Ief-Ieb);
I0u = (1-Iuf-Iub);

ce = Ief.*cef+Ieb.*ceb+I0e.*(w*(1-tau)+r*agrd);
cu = Iuf.*cuf+Iub.*cub+I0u.*(w*mu+r*agrd);

se = Ief.*sef+Ieb.*seb;
su = Iuf.*suf+Iub.*sub;

Dse = diag(Ief.*sef)*Df+diag(Ieb.*seb)*Db;
Dsu = diag(Iuf.*suf)*Df+diag(Iub.*sub)*Db;

P = [Dse,zeros(Na,Na);zeros(Na,Na),Dsu]+B;

A = (1/Gamma+rho)*speye(Na*2)-P;
b = [u(ce);u(cu)]+[ve;vu]./Gamma;

v = A\b;

ven = v(1:Na);
vun = v(Na+1:end);

metric = max(max(abs([ven-ve vun-vu])));

ve = ven;
vu = vun;

end

toc






