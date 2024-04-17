function [residual, g1, g2] = ZLB_four_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 33, 1);

%
% Model equations
%

T20 = exp((-params(1))*y(1))*params(2)*exp(y(2))/exp(y(3));
T86 = exp(y(3))*params(6)/(params(6)-1)*exp(y(13))/exp(y(14));
T142 = ((params(6)-1)/params(6)/params(3)/(1-params(9)-params(10))^params(1))^(1/(params(1)+params(4)));
lhs =exp((-params(1))*y(1));
rhs =T20;
residual(1)= lhs-rhs;
lhs =params(3)*exp(params(4)*y(5));
rhs =exp((-params(1))*y(1))*exp(y(6));
residual(2)= lhs-rhs;
lhs =exp(y(7));
rhs =exp(y(6))/exp(y(8));
residual(3)= lhs-rhs;
lhs =exp(y(1))+exp(y(15))+exp(y(16));
rhs =exp(y(9));
residual(4)= lhs-rhs;
lhs =exp(y(9));
rhs =exp(y(8))*exp(y(5))/exp(y(11));
residual(5)= lhs-rhs;
lhs =exp(y(11));
rhs =(1-params(5))*exp((-params(6))*y(12))*exp(y(3)*params(6))+exp(y(11))*params(5)*exp(y(3)*params(6));
residual(6)= lhs-rhs;
lhs =exp(y(3)*(1-params(6)));
rhs =params(5)+(1-params(5))*exp(y(12)*(1-params(6)));
residual(7)= lhs-rhs;
lhs =exp(y(12));
rhs =T86;
residual(8)= lhs-rhs;
lhs =exp(y(13));
rhs =exp(y(7))*exp((-params(1))*y(1))*exp(y(9))+exp(y(13))*exp(y(3)*params(6))*params(2)*params(5);
residual(9)= lhs-rhs;
lhs =exp(y(14));
rhs =exp((-params(1))*y(1))*exp(y(9))+exp(y(14))*params(2)*params(5)*exp(y(3)*(params(6)-1));
residual(10)= lhs-rhs;
lhs =y(8);
rhs =y(8)*params(13)+x(1);
residual(11)= lhs-rhs;
lhs =y(15);
rhs =(1-params(11))*log(params(32))+y(15)*params(11)+x(3);
residual(12)= lhs-rhs;
lhs =y(16);
rhs =(1-params(12))*log(params(33))+y(16)*params(12)+x(4);
residual(13)= lhs-rhs;
lhs =exp(y(4));
rhs =exp(y(2))/exp(y(3));
residual(14)= lhs-rhs;
lhs =exp(y(10));
rhs =T142*exp(y(8)*(1+params(4))/(params(1)+params(4)));
residual(15)= lhs-rhs;
lhs =y(2);
rhs =y(17);
residual(16)= lhs-rhs;
lhs =y(17);
rhs =y(18);
residual(17)= lhs-rhs;
lhs =y(18);
rhs =y(19);
residual(18)= lhs-rhs;
lhs =y(19);
rhs =y(20);
residual(19)= lhs-rhs;
lhs =y(20);
rhs =(1-params(14))*log(params(20))+y(20)*params(14)+(1-params(14))*(params(7)*(y(23)-log(params(21)))+params(8)*(y(26)-y(29)))+y(33);
residual(20)= lhs-rhs;
lhs =y(21);
rhs =y(3);
residual(21)= lhs-rhs;
lhs =y(22);
rhs =y(21);
residual(22)= lhs-rhs;
lhs =y(23);
rhs =y(22);
residual(23)= lhs-rhs;
lhs =y(24);
rhs =y(9);
residual(24)= lhs-rhs;
lhs =y(25);
rhs =y(24);
residual(25)= lhs-rhs;
lhs =y(26);
rhs =y(25);
residual(26)= lhs-rhs;
lhs =y(27);
rhs =y(10);
residual(27)= lhs-rhs;
lhs =y(28);
rhs =y(27);
residual(28)= lhs-rhs;
lhs =y(29);
rhs =y(28);
residual(29)= lhs-rhs;
lhs =y(30);
rhs =x(2);
residual(30)= lhs-rhs;
lhs =y(31);
rhs =y(30);
residual(31)= lhs-rhs;
lhs =y(32);
rhs =y(31);
residual(32)= lhs-rhs;
lhs =y(33);
rhs =y(32);
residual(33)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(33, 33);

%
% Jacobian matrix
%

  g1(1,1)=(-params(1))*exp((-params(1))*y(1))-exp(y(2))*params(2)*(-params(1))*exp((-params(1))*y(1))/exp(y(3));
  g1(1,2)=(-T20);
  g1(1,3)=(-((-(exp((-params(1))*y(1))*params(2)*exp(y(2))*exp(y(3))))/(exp(y(3))*exp(y(3)))));
  g1(2,1)=(-(exp(y(6))*(-params(1))*exp((-params(1))*y(1))));
  g1(2,5)=params(3)*params(4)*exp(params(4)*y(5));
  g1(2,6)=(-(exp((-params(1))*y(1))*exp(y(6))));
  g1(3,6)=(-(exp(y(6))/exp(y(8))));
  g1(3,7)=exp(y(7));
  g1(3,8)=(-((-(exp(y(6))*exp(y(8))))/(exp(y(8))*exp(y(8)))));
  g1(4,1)=exp(y(1));
  g1(4,9)=(-exp(y(9)));
  g1(4,15)=exp(y(15));
  g1(4,16)=exp(y(16));
  g1(5,5)=(-(exp(y(8))*exp(y(5))/exp(y(11))));
  g1(5,8)=(-(exp(y(8))*exp(y(5))/exp(y(11))));
  g1(5,9)=exp(y(9));
  g1(5,11)=(-((-(exp(y(8))*exp(y(5))*exp(y(11))))/(exp(y(11))*exp(y(11)))));
  g1(6,3)=(-((1-params(5))*exp((-params(6))*y(12))*params(6)*exp(y(3)*params(6))+exp(y(11))*params(5)*params(6)*exp(y(3)*params(6))));
  g1(6,11)=exp(y(11))-exp(y(11))*params(5)*exp(y(3)*params(6));
  g1(6,12)=(-(exp(y(3)*params(6))*(1-params(5))*(-params(6))*exp((-params(6))*y(12))));
  g1(7,3)=(1-params(6))*exp(y(3)*(1-params(6)));
  g1(7,12)=(-((1-params(5))*(1-params(6))*exp(y(12)*(1-params(6)))));
  g1(8,3)=(-T86);
  g1(8,12)=exp(y(12));
  g1(8,13)=(-T86);
  g1(8,14)=(-((-(exp(y(3))*params(6)/(params(6)-1)*exp(y(13))*exp(y(14))))/(exp(y(14))*exp(y(14)))));
  g1(9,1)=(-(exp(y(7))*exp(y(9))*(-params(1))*exp((-params(1))*y(1))));
  g1(9,3)=(-(exp(y(13))*params(2)*params(5)*params(6)*exp(y(3)*params(6))));
  g1(9,7)=(-(exp(y(7))*exp((-params(1))*y(1))*exp(y(9))));
  g1(9,9)=(-(exp(y(7))*exp((-params(1))*y(1))*exp(y(9))));
  g1(9,13)=exp(y(13))-exp(y(13))*exp(y(3)*params(6))*params(2)*params(5);
  g1(10,1)=(-(exp(y(9))*(-params(1))*exp((-params(1))*y(1))));
  g1(10,3)=(-(exp(y(14))*params(2)*params(5)*(params(6)-1)*exp(y(3)*(params(6)-1))));
  g1(10,9)=(-(exp((-params(1))*y(1))*exp(y(9))));
  g1(10,14)=exp(y(14))-exp(y(14))*params(2)*params(5)*exp(y(3)*(params(6)-1));
  g1(11,8)=1-params(13);
  g1(12,15)=1-params(11);
  g1(13,16)=1-params(12);
  g1(14,2)=(-(exp(y(2))/exp(y(3))));
  g1(14,3)=(-((-(exp(y(2))*exp(y(3))))/(exp(y(3))*exp(y(3)))));
  g1(14,4)=exp(y(4));
  g1(15,8)=(-(T142*exp(y(8)*(1+params(4))/(params(1)+params(4)))*(1+params(4))/(params(1)+params(4))));
  g1(15,10)=exp(y(10));
  g1(16,2)=1;
  g1(16,17)=(-1);
  g1(17,17)=1;
  g1(17,18)=(-1);
  g1(18,18)=1;
  g1(18,19)=(-1);
  g1(19,19)=1;
  g1(19,20)=(-1);
  g1(20,20)=1-params(14);
  g1(20,23)=(-((1-params(14))*params(7)));
  g1(20,26)=(-((1-params(14))*params(8)));
  g1(20,29)=(-((1-params(14))*(-params(8))));
  g1(20,33)=(-1);
  g1(21,3)=(-1);
  g1(21,21)=1;
  g1(22,21)=(-1);
  g1(22,22)=1;
  g1(23,22)=(-1);
  g1(23,23)=1;
  g1(24,9)=(-1);
  g1(24,24)=1;
  g1(25,24)=(-1);
  g1(25,25)=1;
  g1(26,25)=(-1);
  g1(26,26)=1;
  g1(27,10)=(-1);
  g1(27,27)=1;
  g1(28,27)=(-1);
  g1(28,28)=1;
  g1(29,28)=(-1);
  g1(29,29)=1;
  g1(30,30)=1;
  g1(31,30)=(-1);
  g1(31,31)=1;
  g1(32,31)=(-1);
  g1(32,32)=1;
  g1(33,32)=(-1);
  g1(33,33)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],33,1089);
end
end
