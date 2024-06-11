function rss = rssfcn_sinus(coef,eta_dzeta,Y,X1,X2,N_s1,N_s2)

coef1 = coef(     1:N_s1     );
coef2 = coef(N_s1+1:N_s1+N_s2);

rss = (Y-exp(X1*coef1)-eta_dzeta*sin(X2*coef2))' ...
     *(Y-exp(X1*coef1)-eta_dzeta*sin(X2*coef2));
 
% disp(coef')
% disp([rss rss])
% pause
