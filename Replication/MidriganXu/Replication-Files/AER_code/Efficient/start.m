clear all;
clc;

mu     = 1.08;               % growth rate
beta   = 0.92*mu;            

rho    = 0.2457;             % persistence transitory productivity
se     = 0.4958;

kappau  = 1.1930;            % fixed cost of joining u
kappap  = kappau + 1;        % fixed cost of joining p 

alpha   = 2/3;
eta     = 0.85;              % decreasing returns
delta   = 0.06; 
Lbar    = 1;

phiu    = 1/(1-eta)*0.2;     % productivity gap unproductive
phip    = phiu;              % productivity gap productive

delta   = delta + mu/beta - 1; 

% Discretize using Rouwenhorst method (Adapted a code from Karen Kopecky)

k = 9;   % number of nodes for e

q = (rho+1)/2;
nu = sqrt((k-1)/(1-rho^2))*se;

P = [q 1-q; 1-q q];

for i = 2:k-1
    
   P =     q*[P zeros(i,1); zeros(1,i+1)] + (1-q)*[zeros(i,1) P; zeros(1,i+1)] + ...
       (1-q)*[zeros(1,i+1); P zeros(i,1)] +     q*[zeros(1,i+1); zeros(i,1) P];
   
   P(2:i,:) = P(2:i,:)/2;
   
end

Perg  = P^5000; Perg = Perg(1,:)';
egrid = linspace(-nu, +nu ,k)';

% Guess a cutoff rule ebar: send everyone above ebart to unproductive
% sector, above ebaru to productive sector

ebargrid = [egrid; egrid(end) + .1];
ebargrid = gridmake(ebargrid, ebargrid);  % cutoff for t is first, cutoff for u is second

Ne  = size(ebargrid,1);


Ce  = zeros(Ne,1);
Ke  = zeros(Ne,1);
Ye  = zeros(Ne,1);
Yte = zeros(Ne,1);
Yme = zeros(Ne,1);
Lte = zeros(Ne,1);
Lme = zeros(Ne,1);
Ate = zeros(Ne,1);
Ame = zeros(Ne,1);
nte = cell(Ne,1);
npe = cell(Ne,1);
nue = cell(Ne,1);

for jt = 1 : Ne

ebart = ebargrid(jt,1);
ebaru = ebargrid(jt,2);

% Solve ergodic distributions nm and nt

ergodic;    % gives egrid (nodes for e) and nt and nm -- measures of producers by type over egrid

At = nt'*exp(egrid);
Am = nu'*exp(egrid+phiu) + np'*exp(egrid+phip);

Lm = 0.95*Lbar;

for it = 1 : 1500

    Lmold = Lm;
    
    Lt  = alpha^(-1/(1-eta))*((1-alpha)*eta/delta)^(-(1-alpha)*eta/(1-eta)/(1-(1-alpha)*eta))*...
           At*Am^(-1/(1-(1-alpha)*eta))*Lm^(1/(1-(1-alpha)*eta));
       
    Lm  = 1/100*(Lbar - Lt)+99/100*Lmold;   

    %fprintf('%4i %6.2e \n',[it, norm(Lm - Lmold)]);

    if norm(Lm-Lmold) < 1e-10, break, end

end

K = ((1-alpha)*eta/delta)^(1/(1-(1-alpha)*eta))*Am^((1-eta)/(1-(1-alpha)*eta))*Lm.^(alpha*eta/(1-(1-alpha)*eta));

Y = At^(1-eta)*Lt.^eta + Am^(1-eta)*(Lm.^alpha.*K.^(1-alpha)).^eta;

finvest = sum(np)*(mu-1)*kappap + sum(nu)*(mu-1)*kappau;

Yt = At^(1-eta)*Lt.^eta;
Ym = Am^(1-eta)*(Lm.^alpha.*K.^(1-alpha)).^eta;

if sum(nu) + sum(np) <1e-6
    
    Lm = 0;
    Lt = Lbar;
    Am = 0;
    K  = 0;
    Y  = At^(1-eta)*Lt.^eta;
    Ym = 0;
    
end

C = Y - delta*K - finvest/beta;

Ye(jt)  = Y;
Yte(jt) = Yt;
Yme(jt) = Ym;
Lte(jt) = Lt;
Lme(jt) = Lm;
Ce(jt)  = C;
nte{jt} = nt;
nue{jt} = nu;
npe{jt} = np;
Ate(jt) = At;
Ame(jt) = Am;
Ke(jt)  = K;

end

[ii,jj] = max(Ce); 

At = Ate(jj);
Lt = Lte(jj);
Yt = Yte(jj);
Am = Ame(jj);
K  = Ke(jj);
Lm = Lme(jj);
Ym = Yme(jj);
Y  = Ye(jj);
C  = Ce(jj);
nt = nte{jj};
np = npe{jj};
nu = nue{jj};

ebart = ebargrid(jj,1);
ebaru = ebargrid(jj,2);

delta = delta + 1 - mu/beta; % back to original delta


finvest = sum(np)*(mu-1)*kappap + sum(nu)*(mu-1)*kappau;

C = Y + (1 - delta - mu)*K - finvest;


fprintf('Productivity cutoff t    =  %9.3f \n',       ebart);
fprintf('Productivity cutoff u    =  %9.3f \n',       ebaru);
fprintf('Measure productive       =  %9.3f \n',       sum(np));
fprintf('Measure unproductive     =  %9.3f \n',       sum(nu));
fprintf('Measure traditional      =  %9.3f \n',       sum(nt));
fprintf('Output                   =  %9.3f \n',       Y);
fprintf('C                        =  %9.3f \n',       C);
fprintf('Share empl. traditional  =  %9.3f \n',       Lt/Lbar);
fprintf('Share outp. traditional  =  %9.3f \n',       Yt/Y);
fprintf('Labor prod. traditional  =  %9.3f \n',       Yt/Lt);
fprintf('Labor prod. modern       =  %9.3f \n',       Ym/Lm);
fprintf('TFP modern               =  %9.3f \n',       Ym/(Lm.^alpha.*K.^(1-alpha)).^eta);
fprintf('TFP traditional          =  %9.3f \n',       Yt/(Lt).^eta);




