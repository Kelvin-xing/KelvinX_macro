function [residual, g1, g2] = mediumdsge_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 26, 1);

%
% Model equations
%

invt = 0;
invtp = 0;
invr = 1;
invrp = 1;
rc = params(13)*(exp(y(6))-1)+params(14)/2*(exp(y(6))-1)^2;
rcp = params(13)*(exp(y(6))-1)+params(14)/2*(exp(y(6))-1)^2;
T46 = exp(y(1))*params(1)*exp(y(5))/exp(y(17));
T56 = 1-params(5)/2*invt^2;
T70 = exp(y(4))*exp(y(7))*(T56-params(5)*invt*invr)+params(5)*exp(y(4))*params(1)*exp(y(7))*invtp*invrp^2;
T154 = exp(y(14))*exp(params(2)*y(23))*exp(y(15)*(1-params(2)))/exp(y(16));
T341 = params(13)*exp(y(6))+params(14)/2*exp(y(6))*2*(exp(y(6))-1);
lhs =exp(y(1));
rhs =1/(exp(y(2))-exp(y(2))*params(4))-params(4)*params(1)/(exp(y(2))-exp(y(2))*params(4));
residual(1)= lhs-rhs;
lhs =y(3);
rhs =(params(13)+(exp(y(6))-1)*params(14))/exp(y(4));
residual(2)= lhs-rhs;
lhs =exp(y(1));
rhs =T46;
residual(3)= lhs-rhs;
lhs =exp(y(1));
rhs =T70;
residual(4)= lhs-rhs;
lhs =exp(y(7));
rhs =params(1)*(exp(y(1))*(exp(y(6))*y(3)-rcp/exp(y(4)))+exp(y(7))*(1-params(3)));
residual(5)= lhs-rhs;
lhs =exp(y(9));
rhs =params(9)/(params(9)-1)*exp(y(19)-y(20));
residual(6)= lhs-rhs;
lhs =exp(y(19));
rhs =params(8)*exp(params(9)*(1+params(6))*(y(10)-y(9)))*exp((1+params(6))*y(15))+exp(y(19))*params(1)*params(11)*exp(y(17)*params(9)*(1+params(6)));
residual(7)= lhs-rhs;
lhs =exp(y(20));
rhs =exp(y(1))*exp(params(9)*(y(10)-y(9)))*exp(y(15))+exp(y(20))*params(1)*params(11)*exp(y(17)*(params(9)-1));
residual(8)= lhs-rhs;
lhs =exp(y(10)*(1-params(9)));
rhs =(1-params(11))*exp(y(9)*(1-params(9)))+exp(y(10)*(1-params(9)))*params(11)*exp(y(17)*(params(9)-1));
residual(9)= lhs-rhs;
lhs =exp(y(13));
rhs =T154;
residual(10)= lhs-rhs;
lhs =exp(y(16));
rhs =exp(y(17)*params(9))*((1-params(12))*exp((-params(9))*y(18))+exp(y(16))*params(12));
residual(11)= lhs-rhs;
lhs =exp(y(17)*(1-params(10)));
rhs =params(12)+(1-params(12))*exp(y(18)*(1-params(10)));
residual(12)= lhs-rhs;
lhs =exp(y(18));
rhs =exp(y(17))*params(10)/(params(10)-1)*exp(y(11)-y(12));
residual(13)= lhs-rhs;
lhs =exp(y(11));
rhs =exp(y(13)+y(1)+y(21))+exp(y(11))*params(1)*params(12)*exp(y(17)*params(10));
residual(14)= lhs-rhs;
lhs =exp(y(12));
rhs =exp(y(1)+y(13))+exp(y(12))*params(1)*params(12)*exp(y(17)*(params(10)-1));
residual(15)= lhs-rhs;
lhs =exp(y(10))/y(3);
rhs =(1-params(2))/params(2)*exp(y(23)-y(15));
residual(16)= lhs-rhs;
lhs =exp(y(10));
rhs =(1-params(2))*exp(y(14)+y(21))*exp(params(2)*(y(23)-y(15)));
residual(17)= lhs-rhs;
lhs =exp(y(13));
rhs =exp(y(8))+exp(y(2))+exp(y(24))+rc*exp(y(22))/exp(y(4));
residual(18)= lhs-rhs;
lhs =exp(y(22));
rhs =exp(y(8))*exp(y(4))*T56+(1-params(3))*exp(y(22));
residual(19)= lhs-rhs;
lhs =exp(y(23));
rhs =exp(y(6))*exp(y(22));
residual(20)= lhs-rhs;
lhs =exp(y(24));
rhs =exp(y(13))*y(26);
residual(21)= lhs-rhs;
lhs =y(26);
rhs =(1-params(15))*params(7)+y(26)*params(15)+x(1);
residual(22)= lhs-rhs;
lhs =y(14);
rhs =y(14)*params(17)+x(3);
residual(23)= lhs-rhs;
lhs =y(4);
rhs =y(4)*params(18)+x(4);
residual(24)= lhs-rhs;
lhs =exp(y(25));
rhs =exp(y(7)-y(1));
residual(25)= lhs-rhs;
lhs =y(5);
rhs =(1-params(16))*log(params(29))+y(5)*params(16)+(1-params(16))*(params(23)*(y(17)-log(params(41)))+params(24)*(y(13)-log(params(37))))+x(2);
residual(26)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(26, 26);

%
% Jacobian matrix
%

  g1(1,1)=exp(y(1));
  g1(1,2)=(-((-(exp(y(2))-exp(y(2))*params(4)))/((exp(y(2))-exp(y(2))*params(4))*(exp(y(2))-exp(y(2))*params(4)))-(-((exp(y(2))-exp(y(2))*params(4))*params(4)*params(1)))/((exp(y(2))-exp(y(2))*params(4))*(exp(y(2))-exp(y(2))*params(4)))));
  g1(2,3)=1;
  g1(2,4)=(-((-((params(13)+(exp(y(6))-1)*params(14))*exp(y(4))))/(exp(y(4))*exp(y(4)))));
  g1(2,6)=(-(exp(y(6))*params(14)/exp(y(4))));
  g1(3,1)=exp(y(1))-T46;
  g1(3,5)=(-T46);
  g1(3,17)=(-((-(exp(y(1))*params(1)*exp(y(5))*exp(y(17))))/(exp(y(17))*exp(y(17)))));
  g1(4,1)=exp(y(1));
  g1(4,4)=(-T70);
  g1(4,7)=(-T70);
  g1(5,1)=(-(params(1)*exp(y(1))*(exp(y(6))*y(3)-rcp/exp(y(4)))));
  g1(5,3)=(-(params(1)*exp(y(6))*exp(y(1))));
  g1(5,4)=(-(params(1)*exp(y(1))*(-((-(exp(y(4))*rcp))/(exp(y(4))*exp(y(4)))))));
  g1(5,6)=(-(params(1)*exp(y(1))*(exp(y(6))*y(3)-T341/exp(y(4)))));
  g1(5,7)=exp(y(7))-params(1)*exp(y(7))*(1-params(3));
  g1(6,9)=exp(y(9));
  g1(6,19)=(-(params(9)/(params(9)-1)*exp(y(19)-y(20))));
  g1(6,20)=(-(params(9)/(params(9)-1)*(-exp(y(19)-y(20)))));
  g1(7,9)=(-(exp((1+params(6))*y(15))*params(8)*exp(params(9)*(1+params(6))*(y(10)-y(9)))*(-(params(9)*(1+params(6))))));
  g1(7,10)=(-(exp((1+params(6))*y(15))*params(8)*params(9)*(1+params(6))*exp(params(9)*(1+params(6))*(y(10)-y(9)))));
  g1(7,15)=(-(params(8)*exp(params(9)*(1+params(6))*(y(10)-y(9)))*(1+params(6))*exp((1+params(6))*y(15))));
  g1(7,17)=(-(exp(y(19))*params(1)*params(11)*params(9)*(1+params(6))*exp(y(17)*params(9)*(1+params(6)))));
  g1(7,19)=exp(y(19))-exp(y(19))*params(1)*params(11)*exp(y(17)*params(9)*(1+params(6)));
  g1(8,1)=(-(exp(y(1))*exp(params(9)*(y(10)-y(9)))*exp(y(15))));
  g1(8,9)=(-(exp(y(15))*exp(y(1))*exp(params(9)*(y(10)-y(9)))*(-params(9))));
  g1(8,10)=(-(exp(y(15))*exp(y(1))*params(9)*exp(params(9)*(y(10)-y(9)))));
  g1(8,15)=(-(exp(y(1))*exp(params(9)*(y(10)-y(9)))*exp(y(15))));
  g1(8,17)=(-(exp(y(20))*params(1)*params(11)*(params(9)-1)*exp(y(17)*(params(9)-1))));
  g1(8,20)=exp(y(20))-exp(y(20))*params(1)*params(11)*exp(y(17)*(params(9)-1));
  g1(9,9)=(-((1-params(11))*(1-params(9))*exp(y(9)*(1-params(9)))));
  g1(9,10)=(1-params(9))*exp(y(10)*(1-params(9)))-params(11)*exp(y(17)*(params(9)-1))*(1-params(9))*exp(y(10)*(1-params(9)));
  g1(9,17)=(-(exp(y(10)*(1-params(9)))*params(11)*(params(9)-1)*exp(y(17)*(params(9)-1))));
  g1(10,13)=exp(y(13));
  g1(10,14)=(-T154);
  g1(10,15)=(-(exp(y(14))*exp(params(2)*y(23))*(1-params(2))*exp(y(15)*(1-params(2)))/exp(y(16))));
  g1(10,16)=(-((-(exp(y(14))*exp(params(2)*y(23))*exp(y(15)*(1-params(2)))*exp(y(16))))/(exp(y(16))*exp(y(16)))));
  g1(10,23)=(-(exp(y(15)*(1-params(2)))*exp(y(14))*params(2)*exp(params(2)*y(23))/exp(y(16))));
  g1(11,16)=exp(y(16))-exp(y(17)*params(9))*exp(y(16))*params(12);
  g1(11,17)=(-(((1-params(12))*exp((-params(9))*y(18))+exp(y(16))*params(12))*params(9)*exp(y(17)*params(9))));
  g1(11,18)=(-(exp(y(17)*params(9))*(1-params(12))*(-params(9))*exp((-params(9))*y(18))));
  g1(12,17)=(1-params(10))*exp(y(17)*(1-params(10)));
  g1(12,18)=(-((1-params(12))*(1-params(10))*exp(y(18)*(1-params(10)))));
  g1(13,11)=(-(exp(y(17))*params(10)/(params(10)-1)*exp(y(11)-y(12))));
  g1(13,12)=(-(exp(y(17))*params(10)/(params(10)-1)*(-exp(y(11)-y(12)))));
  g1(13,17)=(-(exp(y(17))*params(10)/(params(10)-1)*exp(y(11)-y(12))));
  g1(13,18)=exp(y(18));
  g1(14,1)=(-exp(y(13)+y(1)+y(21)));
  g1(14,11)=exp(y(11))-exp(y(11))*params(1)*params(12)*exp(y(17)*params(10));
  g1(14,13)=(-exp(y(13)+y(1)+y(21)));
  g1(14,17)=(-(exp(y(11))*params(1)*params(12)*params(10)*exp(y(17)*params(10))));
  g1(14,21)=(-exp(y(13)+y(1)+y(21)));
  g1(15,1)=(-exp(y(1)+y(13)));
  g1(15,12)=exp(y(12))-exp(y(12))*params(1)*params(12)*exp(y(17)*(params(10)-1));
  g1(15,13)=(-exp(y(1)+y(13)));
  g1(15,17)=(-(exp(y(12))*params(1)*params(12)*(params(10)-1)*exp(y(17)*(params(10)-1))));
  g1(16,3)=(-exp(y(10)))/(y(3)*y(3));
  g1(16,10)=exp(y(10))/y(3);
  g1(16,15)=(-((1-params(2))/params(2)*(-exp(y(23)-y(15)))));
  g1(16,23)=(-((1-params(2))/params(2)*exp(y(23)-y(15))));
  g1(17,10)=exp(y(10));
  g1(17,14)=(-((1-params(2))*exp(y(14)+y(21))*exp(params(2)*(y(23)-y(15)))));
  g1(17,15)=(-((1-params(2))*exp(y(14)+y(21))*exp(params(2)*(y(23)-y(15)))*(-params(2))));
  g1(17,21)=(-((1-params(2))*exp(y(14)+y(21))*exp(params(2)*(y(23)-y(15)))));
  g1(17,23)=(-((1-params(2))*exp(y(14)+y(21))*params(2)*exp(params(2)*(y(23)-y(15)))));
  g1(18,2)=(-exp(y(2)));
  g1(18,4)=(-((-(exp(y(4))*rc*exp(y(22))))/(exp(y(4))*exp(y(4)))));
  g1(18,6)=(-(exp(y(22))*T341/exp(y(4))));
  g1(18,8)=(-exp(y(8)));
  g1(18,13)=exp(y(13));
  g1(18,22)=(-(rc*exp(y(22))/exp(y(4))));
  g1(18,24)=(-exp(y(24)));
  g1(19,4)=(-(exp(y(8))*exp(y(4))*T56));
  g1(19,8)=(-(exp(y(8))*exp(y(4))*T56));
  g1(19,22)=exp(y(22))-(1-params(3))*exp(y(22));
  g1(20,6)=(-(exp(y(6))*exp(y(22))));
  g1(20,22)=(-(exp(y(6))*exp(y(22))));
  g1(20,23)=exp(y(23));
  g1(21,13)=(-(exp(y(13))*y(26)));
  g1(21,24)=exp(y(24));
  g1(21,26)=(-exp(y(13)));
  g1(22,26)=1-params(15);
  g1(23,14)=1-params(17);
  g1(24,4)=1-params(18);
  g1(25,1)=exp(y(7)-y(1));
  g1(25,7)=(-exp(y(7)-y(1)));
  g1(25,25)=exp(y(25));
  g1(26,5)=1-params(16);
  g1(26,13)=(-((1-params(16))*params(24)));
  g1(26,17)=(-((1-params(16))*params(23)));
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],26,676);
end
end
