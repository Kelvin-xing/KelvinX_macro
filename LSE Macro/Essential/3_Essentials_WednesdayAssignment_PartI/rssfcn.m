function rss = rssfcn(coef,Y,X)

rss = (Y-exp(X*coef))'*(Y-exp(X*coef));

