function [residual, g1, g2, g3] = NK_sticky_taylor_dynamic(y, x, params, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(16, 1);
T23 = params(2)*exp((-params(1))*y(20))*exp(y(5))/exp(y(21));
T109 = params(6)/(params(6)-1)*exp(y(6))*exp(y(15))/exp(y(16));
T142 = ((params(6)-1)/params(6)/params(3))^(1/(params(1)+params(4)));
lhs =exp((-params(1))*y(4));
rhs =T23;
residual(1)= lhs-rhs;
lhs =params(3)*exp(params(4)*y(8));
rhs =exp((-params(1))*y(4))*exp(y(9));
residual(2)= lhs-rhs;
lhs =exp(y(10));
rhs =exp(y(9))/exp(y(11));
residual(3)= lhs-rhs;
lhs =y(5);
rhs =(1-params(11))*log(params(15))+params(11)*y(1)+(1-params(11))*(params(8)*(y(6)-log(params(16)))+params(9)*y(18))+x(it_, 2);
residual(4)= lhs-rhs;
lhs =exp(y(4));
rhs =exp(y(12));
residual(5)= lhs-rhs;
lhs =exp(y(12));
rhs =exp(y(11))*exp(y(8))/exp(y(13));
residual(6)= lhs-rhs;
lhs =exp(y(13));
rhs =(1-params(5))*exp((-params(6))*y(14))*exp(y(6)*params(6))+params(5)*exp(y(6)*params(6))*exp(y(3));
residual(7)= lhs-rhs;
lhs =exp(y(6)*(1-params(6)));
rhs =params(5)+(1-params(5))*exp(y(14)*(1-params(6)));
residual(8)= lhs-rhs;
lhs =exp(y(14));
rhs =T109;
residual(9)= lhs-rhs;
lhs =exp(y(15));
rhs =exp(y(10))*exp((-params(1))*y(4))*exp(y(12))+params(2)*params(5)*exp(y(21)*params(6))*exp(y(22));
residual(10)= lhs-rhs;
lhs =exp(y(16));
rhs =exp((-params(1))*y(4))*exp(y(12))+params(2)*params(5)*exp(y(21)*(params(6)-1))*exp(y(23));
residual(11)= lhs-rhs;
lhs =y(11);
rhs =params(10)*y(2)+x(it_, 1);
residual(12)= lhs-rhs;
lhs =exp(y(17));
rhs =T142*exp(y(11)*(1+params(4))/(params(1)+params(4)));
residual(13)= lhs-rhs;
lhs =y(18);
rhs =y(12)-y(17);
residual(14)= lhs-rhs;
lhs =exp(y(7));
rhs =exp(y(5))/exp(y(21));
residual(15)= lhs-rhs;
lhs =exp(y(19));
rhs =exp(y(5))*params(7)/(exp(y(5))-1)*exp(params(1)*y(4));
residual(16)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(16, 25);

%
% Jacobian matrix
%

g1(1,4)=(-params(1))*exp((-params(1))*y(4));
g1(1,20)=(-(exp(y(5))*params(2)*(-params(1))*exp((-params(1))*y(20))/exp(y(21))));
g1(1,5)=(-T23);
g1(1,21)=(-((-(params(2)*exp((-params(1))*y(20))*exp(y(5))*exp(y(21))))/(exp(y(21))*exp(y(21)))));
g1(2,4)=(-(exp(y(9))*(-params(1))*exp((-params(1))*y(4))));
g1(2,8)=params(3)*params(4)*exp(params(4)*y(8));
g1(2,9)=(-(exp((-params(1))*y(4))*exp(y(9))));
g1(3,9)=(-(exp(y(9))/exp(y(11))));
g1(3,10)=exp(y(10));
g1(3,11)=(-((-(exp(y(9))*exp(y(11))))/(exp(y(11))*exp(y(11)))));
g1(4,1)=(-params(11));
g1(4,5)=1;
g1(4,6)=(-((1-params(11))*params(8)));
g1(4,18)=(-((1-params(11))*params(9)));
g1(4,25)=(-1);
g1(5,4)=exp(y(4));
g1(5,12)=(-exp(y(12)));
g1(6,8)=(-(exp(y(11))*exp(y(8))/exp(y(13))));
g1(6,11)=(-(exp(y(11))*exp(y(8))/exp(y(13))));
g1(6,12)=exp(y(12));
g1(6,13)=(-((-(exp(y(11))*exp(y(8))*exp(y(13))))/(exp(y(13))*exp(y(13)))));
g1(7,6)=(-((1-params(5))*exp((-params(6))*y(14))*params(6)*exp(y(6)*params(6))+exp(y(3))*params(5)*params(6)*exp(y(6)*params(6))));
g1(7,3)=(-(params(5)*exp(y(6)*params(6))*exp(y(3))));
g1(7,13)=exp(y(13));
g1(7,14)=(-(exp(y(6)*params(6))*(1-params(5))*(-params(6))*exp((-params(6))*y(14))));
g1(8,6)=(1-params(6))*exp(y(6)*(1-params(6)));
g1(8,14)=(-((1-params(5))*(1-params(6))*exp(y(14)*(1-params(6)))));
g1(9,6)=(-T109);
g1(9,14)=exp(y(14));
g1(9,15)=(-T109);
g1(9,16)=(-((-(params(6)/(params(6)-1)*exp(y(6))*exp(y(15))*exp(y(16))))/(exp(y(16))*exp(y(16)))));
g1(10,4)=(-(exp(y(10))*exp(y(12))*(-params(1))*exp((-params(1))*y(4))));
g1(10,21)=(-(exp(y(22))*params(2)*params(5)*params(6)*exp(y(21)*params(6))));
g1(10,10)=(-(exp(y(10))*exp((-params(1))*y(4))*exp(y(12))));
g1(10,12)=(-(exp(y(10))*exp((-params(1))*y(4))*exp(y(12))));
g1(10,15)=exp(y(15));
g1(10,22)=(-(params(2)*params(5)*exp(y(21)*params(6))*exp(y(22))));
g1(11,4)=(-(exp(y(12))*(-params(1))*exp((-params(1))*y(4))));
g1(11,21)=(-(exp(y(23))*params(2)*params(5)*(params(6)-1)*exp(y(21)*(params(6)-1))));
g1(11,12)=(-(exp((-params(1))*y(4))*exp(y(12))));
g1(11,16)=exp(y(16));
g1(11,23)=(-(params(2)*params(5)*exp(y(21)*(params(6)-1))*exp(y(23))));
g1(12,2)=(-params(10));
g1(12,11)=1;
g1(12,24)=(-1);
g1(13,11)=(-(T142*(1+params(4))/(params(1)+params(4))*exp(y(11)*(1+params(4))/(params(1)+params(4)))));
g1(13,17)=exp(y(17));
g1(14,12)=(-1);
g1(14,17)=1;
g1(14,18)=1;
g1(15,5)=(-(exp(y(5))/exp(y(21))));
g1(15,21)=(-((-(exp(y(5))*exp(y(21))))/(exp(y(21))*exp(y(21)))));
g1(15,7)=exp(y(7));
g1(16,4)=(-(exp(y(5))*params(7)/(exp(y(5))-1)*params(1)*exp(params(1)*y(4))));
g1(16,5)=(-(exp(params(1)*y(4))*(exp(y(5))*params(7)*(exp(y(5))-1)-exp(y(5))*exp(y(5))*params(7))/((exp(y(5))-1)*(exp(y(5))-1))));
g1(16,19)=exp(y(19));
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],16,625);
end
if nargout >= 4,
%
% Third order derivatives
%

  g3 = sparse([],[],[],16,15625);
end
end
