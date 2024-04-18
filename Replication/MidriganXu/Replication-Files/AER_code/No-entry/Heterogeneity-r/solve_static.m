% compute factor choices and output given (A)

A = s(:,1);
E = exp(egrid(s(:,2)));

% 1: evaluate f'(a)

fpa = (1-alpha)*eta*(alpha*eta/W)^(alpha*eta/(1-alpha*eta))*E.^((1-eta)/(1-alpha*eta)).*A.^((1-alpha)*eta/(1-alpha*eta)-1);

rr = r2.*(fpa >= r2 + delta) + r1.*(fpa < r1 + delta) + (fpa-delta).*(fpa>=r1+delta & fpa < r2 + delta);

      
cr = funfitxy(fspace, s, rr);

