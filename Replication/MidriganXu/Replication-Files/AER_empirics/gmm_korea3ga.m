function mom = gmm_korea3ga(par) 

global lnv lnv_1 lnv_2 lnv_3 lnl lnl_1 lnl_2 lnl_3 lnk lnk_1 lnk_2 lnk_3
eta=exp(par(1,1))/(1+exp(par(1,1)));
alpha=exp(par(1,2))/(1+exp(par(1,2)));
rho=exp(par(1,3))/(1+exp(par(1,3)));

dify1=lnv-lnv_1;
dify2=lnv_1-lnv_2;
difk1=lnk - lnk_1;
difk2=lnk_1 - lnk_2;
difl1=lnl - lnl_1;
difl2=lnl_1 - lnl_2;

res = (dify1-(eta*(1-alpha))*difk1-eta*alpha*difl1-rho*(dify2-eta*(1-alpha)*difk2-eta*alpha*difl2));
nw = lnv-eta*(1-alpha)*lnk-rho*eta*(1-alpha)*lnk_1+eta*alpha*lnl-rho*eta*alpha*lnl_1 - rho*lnv_1-1+rho;

difl3=lnl_2-lnl_3;
difk3=lnk_2-lnk_3;
dify3=lnv_2-lnv_3;

mom1 = res.*lnl_3;
mom2 = res.*lnv_3;
mom3 = res.*lnk_3;
mom4 = nw.*difl3;
mom5 = nw.*difk3;
mom6 = nw.*dify3;

mom = norm([mom1,mom2,mom3, mom4, mom5, mom6]);
end
