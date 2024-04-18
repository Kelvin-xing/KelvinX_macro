T = 500;  % periods we wait for weights to converge

% Guess initial measures

np = 0.25*1/k*ones(k,1);
nu = 0.50*1/k*ones(k,1);
nt = 0.25*1/k*ones(k,1);

for t = 1 : T
    
ntold = nt;
nuold = nu;
npold = np;

% traditional sector

nt = bsxfun(@times, ntold(egrid<ebart), P(egrid<ebart,:)); nt = sum(nt, 1)';
nt = 1/mu*nt + (1-1/mu)*Perg;

% unproductive sector

ntu = bsxfun(@times, ntold(egrid>=ebart), P(egrid>=ebart,:)); ntu = sum(ntu, 1)';
nuu = bsxfun(@times, nuold(egrid <ebaru), P(egrid <ebaru,:)); nuu = sum(nuu, 1)';

nu = 1/mu*ntu + 1/mu*nuu;

% productive sector

nup = bsxfun(@times, nuold(egrid>=ebaru), P(egrid>=ebaru,:)); nup = sum(nup, 1)';
npp = bsxfun(@times, npold              , P                ); npp = sum(npp, 1)';

np = 1/mu*nup + 1/mu*npp;

%fprintf('%4i %6.2e %6.2e %6.2e \n',[t, norm(nt-ntold), norm(nu-nuold), norm(np-npold)]);

if norm(nt-ntold) < 1e-10 && norm(nu-nuold) < 1e-10 && norm(np-npold) < 1e-10, break, end

end

