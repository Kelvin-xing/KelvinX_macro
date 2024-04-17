aa=oo_.dr.ghx;
bb=oo_.dr.ghu;
cc=oo_.dr.ghxx;
dd=oo_.dr.ghuu;
ee=oo_.dr.ghxu;
ys=oo_.steady_state;
ys=ys(oo_.dr.order_var);
delta2=oo_.dr.ghs2;
cs=oo_.endo_simul(1,1:10)';
et=oo_.exo_simul(1:10,1);
ccs=zeros(10,1);
ccs(1,1)=ys(1,1)+bb(3,1)*et(1,1)+0.5*dd(3,1)*et(1,1)^2+0.5*delta2(3,1);

yh1=[oo_.endo_simul(2,1:10)-ys(1,1);oo_.endo_simul(4,1:10)];
for t=1:10
    yh2(:,t)=kron(yh1(:,t),yh1(:,t)); 
end
u=oo_.exo_simul(1:10,1);
uu=u.^2;
yhu=yh1;
yhu(1,:)=yhu(1,:).*oo_.exo_simul(1:10,1)';
yhu(2,:)=yhu(2,:).*oo_.exo_simul(1:10,1)';

yt=aa*yh1+bb*u'+0.5*cc*yh2+0.5*dd*uu'+ee*yhu;
yt(3,:)+ys(3,1)+0.5*delta2(1,1)