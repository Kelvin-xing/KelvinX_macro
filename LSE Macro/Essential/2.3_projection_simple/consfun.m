function c = consfun(k,z,coef)

c_log = coef(1) + coef(2)*log(k) + coef(3)*log(z) + coef(4)*(log(k)).^2 + coef(5)*(log(z)).^2 + coef(6)*log(k).*log(z);

c = exp(c_log);

