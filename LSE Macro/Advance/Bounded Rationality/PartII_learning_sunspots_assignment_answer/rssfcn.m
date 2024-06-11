function rss = rssfcn(coef,Y,X,X_dzeta,eta_dzeta)

rss = (Y-exp(X*coef+X_dzeta*eta_dzeta))'*(Y-exp(X*coef+X_dzeta*eta_dzeta));
% disp(coef')
% disp(rss)
% pause
