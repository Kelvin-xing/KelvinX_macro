function [residual, g1, g2, g3] = welfare_loss_dynamic(y, x, params, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(22, 1);
lhs =y(14);
rhs =y(15)+params(3)*(y(17)-y(2));
residual(1)= lhs-rhs;
lhs =y(20);
rhs =y(15)+y(4);
residual(2)= lhs-rhs;
lhs =y(21);
rhs =y(14)+y(5);
residual(3)= lhs-rhs;
lhs =y(18);
rhs =y(17)*(1-params(3));
residual(4)= lhs-rhs;
lhs =y(17)-y(2);
rhs =y(19)-y(3)+y(16)-y(15);
residual(5)= lhs-rhs;
lhs =y(12);
rhs =y(7)+y(17)/params(18);
residual(6)= lhs-rhs;
lhs =y(12);
rhs =y(13)+y(9);
residual(7)= lhs-rhs;
lhs =y(15);
rhs =params(1)*y(31)+y(9)*params(16);
residual(8)= lhs-rhs;
lhs =y(9);
rhs =y(30)-1/params(18)*(y(11)-y(31)-y(10));
residual(9)= lhs-rhs;
lhs =y(13);
rhs =params(20)*y(8)+y(7)*params(21);
residual(10)= lhs-rhs;
lhs =y(10);
rhs =y(8)*params(20)*(-params(18))*(1-params(8))-params(21)*params(6)*(y(29)-y(7));
residual(11)= lhs-rhs;
lhs =y(23);
rhs =y(17)*params(3)*(params(17)/params(2)-1);
residual(12)= lhs-rhs;
lhs =y(22);
rhs =y(9)*(params(18)+params(6));
residual(13)= lhs-rhs;
residual(14) = y(19);
lhs =y(8);
rhs =params(8)*y(1)+x(it_, 1)+params(15)*x(it_, 3);
residual(15)= lhs-rhs;
lhs =y(7);
rhs =y(29)-(y(24)-y(32))/params(2);
residual(16)= lhs-rhs;
lhs =y(25);
rhs =y(7)*(params(6)+params(2))-(1+params(6))*y(26);
residual(17)= lhs-rhs;
lhs =y(16);
rhs =params(1)*y(32)+y(25)*params(19);
residual(18)= lhs-rhs;
lhs =y(24);
rhs =y(32)*params(12)+y(26)*params(13);
residual(19)= lhs-rhs;
lhs =y(26);
rhs =x(it_, 3)+params(14)*y(6);
residual(20)= lhs-rhs;
lhs =y(28);
rhs =y(12)*params(2)+params(6)*y(27);
residual(21)= lhs-rhs;
lhs =y(12);
rhs =y(8)+y(27);
residual(22)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(22, 35);

%
% Jacobian matrix
%

g1(1,14)=1;
g1(1,15)=(-1);
g1(1,2)=params(3);
g1(1,17)=(-params(3));
g1(2,15)=(-1);
g1(2,4)=(-1);
g1(2,20)=1;
g1(3,14)=(-1);
g1(3,5)=(-1);
g1(3,21)=1;
g1(4,17)=(-(1-params(3)));
g1(4,18)=1;
g1(5,15)=1;
g1(5,16)=(-1);
g1(5,2)=(-1);
g1(5,17)=1;
g1(5,3)=1;
g1(5,19)=(-1);
g1(6,7)=(-1);
g1(6,12)=1;
g1(6,17)=(-(1/params(18)));
g1(7,9)=(-1);
g1(7,12)=1;
g1(7,13)=(-1);
g1(8,9)=(-params(16));
g1(8,15)=1;
g1(8,31)=(-params(1));
g1(9,9)=1;
g1(9,30)=(-1);
g1(9,10)=(-(1/params(18)));
g1(9,11)=1/params(18);
g1(9,31)=(-(1/params(18)));
g1(10,7)=(-params(21));
g1(10,8)=(-params(20));
g1(10,13)=1;
g1(11,7)=(-(params(21)*params(6)));
g1(11,29)=params(21)*params(6);
g1(11,8)=(-(params(20)*(-params(18))*(1-params(8))));
g1(11,10)=1;
g1(12,17)=(-(params(3)*(params(17)/params(2)-1)));
g1(12,23)=1;
g1(13,9)=(-(params(18)+params(6)));
g1(13,22)=1;
g1(14,19)=1;
g1(15,1)=(-params(8));
g1(15,8)=1;
g1(15,33)=(-1);
g1(15,35)=(-params(15));
g1(16,7)=1;
g1(16,29)=(-1);
g1(16,32)=(-1)/params(2);
g1(16,24)=1/params(2);
g1(17,7)=(-(params(6)+params(2)));
g1(17,25)=1;
g1(17,26)=1+params(6);
g1(18,16)=1;
g1(18,32)=(-params(1));
g1(18,25)=(-params(19));
g1(19,32)=(-params(12));
g1(19,24)=1;
g1(19,26)=(-params(13));
g1(20,6)=(-params(14));
g1(20,26)=1;
g1(20,35)=(-1);
g1(21,12)=(-params(2));
g1(21,27)=(-params(6));
g1(21,28)=1;
g1(22,8)=(-1);
g1(22,12)=1;
g1(22,27)=(-1);
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],22,1225);
end
if nargout >= 4,
%
% Third order derivatives
%

  g3 = sparse([],[],[],22,42875);
end
end
