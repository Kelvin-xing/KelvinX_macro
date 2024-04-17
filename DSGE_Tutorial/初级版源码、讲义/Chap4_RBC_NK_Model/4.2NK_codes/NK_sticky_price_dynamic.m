function [residual, g1, g2, g3] = NK_sticky_price_dynamic(y, x, params, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(17, 1);
T23 = params(2)*exp((-params(1))*y(23))*exp(y(7))/exp(y(24));
T116 = params(6)/(params(6)-1)*exp(y(8))*exp(y(17))/exp(y(18));
T149 = ((params(6)-1)/params(6)/params(3))^(1/(params(1)+params(4)));
lhs =exp((-params(1))*y(6));
rhs =T23;
residual(1)= lhs-rhs;
lhs =params(3)*exp(params(4)*y(9));
rhs =exp((-params(1))*y(6))*exp(y(10));
residual(2)= lhs-rhs;
lhs =exp(y(12));
rhs =exp(y(10))/exp(y(13));
residual(3)= lhs-rhs;
lhs =exp(y(11));
rhs =exp(y(7))*params(7)*exp(params(1)*y(6))/(exp(y(7))-1);
residual(4)= lhs-rhs;
lhs =y(19);
rhs =(1-params(9))*(params(14)-1)-y(8)+params(9)*y(5)+params(9)*y(1)+x(it_, 2);
residual(5)= lhs-rhs;
lhs =y(19);
rhs =y(11)-y(2);
residual(6)= lhs-rhs;
lhs =exp(y(6));
rhs =exp(y(14));
residual(7)= lhs-rhs;
lhs =exp(y(14));
rhs =exp(y(13))*exp(y(9))/exp(y(15));
residual(8)= lhs-rhs;
lhs =exp(y(15));
rhs =(1-params(5))*exp((-params(6))*y(16))*exp(y(8)*params(6))+params(5)*exp(y(8)*params(6))*exp(y(4));
residual(9)= lhs-rhs;
lhs =exp(y(8)*(1-params(6)));
rhs =params(5)+(1-params(5))*exp(y(16)*(1-params(6)));
residual(10)= lhs-rhs;
lhs =exp(y(16));
rhs =T116;
residual(11)= lhs-rhs;
lhs =exp(y(17));
rhs =exp(y(12))*exp((-params(1))*y(6))*exp(y(14))+params(2)*params(5)*exp(y(24)*params(6))*exp(y(25));
residual(12)= lhs-rhs;
lhs =exp(y(18));
rhs =exp((-params(1))*y(6))*exp(y(14))+params(2)*params(5)*exp(y(24)*(params(6)-1))*exp(y(26));
residual(13)= lhs-rhs;
lhs =y(13);
rhs =params(8)*y(3)+x(it_, 1);
residual(14)= lhs-rhs;
lhs =exp(y(20));
rhs =T149*exp(y(13)*(1+params(4))/(params(1)+params(4)));
residual(15)= lhs-rhs;
lhs =y(21);
rhs =y(14)-y(20);
residual(16)= lhs-rhs;
lhs =y(22);
rhs =exp(y(7))-exp(y(24));
residual(17)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(17, 28);

%
% Jacobian matrix
%

g1(1,6)=(-params(1))*exp((-params(1))*y(6));
g1(1,23)=(-(exp(y(7))*params(2)*(-params(1))*exp((-params(1))*y(23))/exp(y(24))));
g1(1,7)=(-T23);
g1(1,24)=(-((-(params(2)*exp((-params(1))*y(23))*exp(y(7))*exp(y(24))))/(exp(y(24))*exp(y(24)))));
g1(2,6)=(-(exp(y(10))*(-params(1))*exp((-params(1))*y(6))));
g1(2,9)=params(3)*params(4)*exp(params(4)*y(9));
g1(2,10)=(-(exp((-params(1))*y(6))*exp(y(10))));
g1(3,10)=(-(exp(y(10))/exp(y(13))));
g1(3,12)=exp(y(12));
g1(3,13)=(-((-(exp(y(10))*exp(y(13))))/(exp(y(13))*exp(y(13)))));
g1(4,6)=(-(exp(y(7))*params(7)*params(1)*exp(params(1)*y(6))/(exp(y(7))-1)));
g1(4,7)=(-((exp(y(7))*params(7)*exp(params(1)*y(6))*(exp(y(7))-1)-exp(y(7))*exp(y(7))*params(7)*exp(params(1)*y(6)))/((exp(y(7))-1)*(exp(y(7))-1))));
g1(4,11)=exp(y(11));
g1(5,1)=(-params(9));
g1(5,8)=1;
g1(5,5)=(-params(9));
g1(5,19)=1;
g1(5,28)=(-1);
g1(6,2)=1;
g1(6,11)=(-1);
g1(6,19)=1;
g1(7,6)=exp(y(6));
g1(7,14)=(-exp(y(14)));
g1(8,9)=(-(exp(y(13))*exp(y(9))/exp(y(15))));
g1(8,13)=(-(exp(y(13))*exp(y(9))/exp(y(15))));
g1(8,14)=exp(y(14));
g1(8,15)=(-((-(exp(y(13))*exp(y(9))*exp(y(15))))/(exp(y(15))*exp(y(15)))));
g1(9,8)=(-((1-params(5))*exp((-params(6))*y(16))*params(6)*exp(y(8)*params(6))+exp(y(4))*params(5)*params(6)*exp(y(8)*params(6))));
g1(9,4)=(-(params(5)*exp(y(8)*params(6))*exp(y(4))));
g1(9,15)=exp(y(15));
g1(9,16)=(-(exp(y(8)*params(6))*(1-params(5))*(-params(6))*exp((-params(6))*y(16))));
g1(10,8)=(1-params(6))*exp(y(8)*(1-params(6)));
g1(10,16)=(-((1-params(5))*(1-params(6))*exp(y(16)*(1-params(6)))));
g1(11,8)=(-T116);
g1(11,16)=exp(y(16));
g1(11,17)=(-T116);
g1(11,18)=(-((-(params(6)/(params(6)-1)*exp(y(8))*exp(y(17))*exp(y(18))))/(exp(y(18))*exp(y(18)))));
g1(12,6)=(-(exp(y(12))*exp(y(14))*(-params(1))*exp((-params(1))*y(6))));
g1(12,24)=(-(exp(y(25))*params(2)*params(5)*params(6)*exp(y(24)*params(6))));
g1(12,12)=(-(exp(y(12))*exp((-params(1))*y(6))*exp(y(14))));
g1(12,14)=(-(exp(y(12))*exp((-params(1))*y(6))*exp(y(14))));
g1(12,17)=exp(y(17));
g1(12,25)=(-(params(2)*params(5)*exp(y(24)*params(6))*exp(y(25))));
g1(13,6)=(-(exp(y(14))*(-params(1))*exp((-params(1))*y(6))));
g1(13,24)=(-(exp(y(26))*params(2)*params(5)*(params(6)-1)*exp(y(24)*(params(6)-1))));
g1(13,14)=(-(exp((-params(1))*y(6))*exp(y(14))));
g1(13,18)=exp(y(18));
g1(13,26)=(-(params(2)*params(5)*exp(y(24)*(params(6)-1))*exp(y(26))));
g1(14,3)=(-params(8));
g1(14,13)=1;
g1(14,27)=(-1);
g1(15,13)=(-(T149*(1+params(4))/(params(1)+params(4))*exp(y(13)*(1+params(4))/(params(1)+params(4)))));
g1(15,20)=exp(y(20));
g1(16,14)=(-1);
g1(16,20)=1;
g1(16,21)=1;
g1(17,7)=(-exp(y(7)));
g1(17,24)=exp(y(24));
g1(17,22)=1;
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],17,784);
end
if nargout >= 4,
%
% Third order derivatives
%

  g3 = sparse([],[],[],17,21952);
end
end
