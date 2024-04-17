function [residual, g1, g2, g3] = modelout_dynamic(y, x, params, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(10, 1);
T31 = 1+params(5)*(y(7)-1)^2/2;
T49 = y(7)*params(1)*y(8)^2;
T94 = params(5)*y(7)+params(5)*(y(7)-1)+params(5)*(2*y(7)-2)*((params(4)-1)*(1+params(8))-y(6)*params(4)*y(8)*params(2)/y(11))/2;
T104 = y(8)*params(1)*y(7)^2;
T111 = (1+y(5))^2;
T134 = 1+(y(7)-1)^2*params(5)/2;
T174 = params(5)*2*(y(7)-1)/2;
lhs =y(9);
rhs =log(y(8))-params(2)*y(6)^2/2;
residual(1)= lhs-rhs;
lhs =y(10);
rhs =y(9)+params(1)*y(17);
residual(2)= lhs-rhs;
residual(3) = y(14)*T31+1/y(8)-params(3)*y(12)/(y(16)*y(15))+params(3)*y(1)*y(3)/T49+T31*y(6)*params(2)*params(4)*y(13)/y(11);
residual(4) = T31*y(13)*params(4)*y(8)*params(2)/y(11)-params(2)*y(6)-y(14)*y(11);
residual(5) = y(4)*((y(7)-1)*params(5)*params(3)+y(7)*params(5)*params(3))/params(1)-y(13)*T94+(2*y(7)-2)*params(5)*y(8)*y(14)/2+params(3)*y(1)*y(3)/T104;
residual(6) = (-y(12))/T111;
lhs =1/(1+y(5));
rhs =y(8)*params(3)/(y(16)*y(15));
residual(7)= lhs-rhs;
residual(8) = T31*((1+params(8))*(1-params(4))+params(4)*y(8)*params(2)*y(6)/y(11))-y(7)*params(5)*(y(7)-1)+y(15)*params(5)*params(3)*(y(15)-1);
residual(9) = y(8)*T134-y(6)*y(11);
lhs =log(y(11));
rhs =params(6)*log(y(2))+x(it_, 1);
residual(10)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(10, 18);

%
% Jacobian matrix
%

g1(1,6)=params(2)*2*y(6)/2;
g1(1,8)=(-(1/y(8)));
g1(1,9)=1;
g1(2,9)=(-1);
g1(2,10)=1;
g1(2,17)=(-params(1));
g1(3,6)=T31*params(2)*params(4)*y(13)/y(11);
g1(3,7)=y(14)*T174+(-(params(3)*y(1)*y(3)*params(1)*y(8)^2))/(T49*T49)+y(6)*params(2)*params(4)*y(13)*T174/y(11);
g1(3,15)=(-((-(params(3)*y(12)*y(16)))/(y(16)*y(15)*y(16)*y(15))));
g1(3,1)=params(3)*y(3)/T49;
g1(3,8)=(-1)/(y(8)*y(8))+(-(params(3)*y(1)*y(3)*y(7)*params(1)*2*y(8)))/(T49*T49);
g1(3,16)=(-((-(params(3)*y(12)*y(15)))/(y(16)*y(15)*y(16)*y(15))));
g1(3,11)=(-(T31*y(6)*params(2)*params(4)*y(13)))/(y(11)*y(11));
g1(3,3)=params(3)*y(1)/T49;
g1(3,12)=(-(params(3)/(y(16)*y(15))));
g1(3,13)=T31*y(6)*params(2)*params(4)/y(11);
g1(3,14)=T31;
g1(4,6)=(-params(2));
g1(4,7)=y(13)*params(4)*y(8)*params(2)*T174/y(11);
g1(4,8)=T31*params(2)*params(4)*y(13)/y(11);
g1(4,11)=(-(T31*y(13)*params(4)*y(8)*params(2)))/(y(11)*y(11))-y(14);
g1(4,13)=T31*params(4)*y(8)*params(2)/y(11);
g1(4,14)=(-y(11));
g1(5,6)=(-(y(13)*params(5)*(2*y(7)-2)*(-(params(4)*y(8)*params(2)/y(11)))/2));
g1(5,7)=y(4)*(params(5)*params(3)+params(5)*params(3))/params(1)-y(13)*(params(5)+params(5)+((params(4)-1)*(1+params(8))-y(6)*params(4)*y(8)*params(2)/y(11))*2*params(5)/2)+2*params(5)*y(8)*y(14)/2+(-(params(3)*y(1)*y(3)*2*y(7)*y(8)*params(1)))/(T104*T104);
g1(5,1)=params(3)*y(3)/T104;
g1(5,8)=(-(y(13)*params(5)*(2*y(7)-2)*(-(y(6)*params(2)*params(4)/y(11)))/2))+(2*y(7)-2)*y(14)*params(5)/2+(-(params(3)*y(1)*y(3)*params(1)*y(7)^2))/(T104*T104);
g1(5,11)=(-(y(13)*params(5)*(2*y(7)-2)*(-((-(y(6)*params(4)*y(8)*params(2)))/(y(11)*y(11))))/2));
g1(5,3)=params(3)*y(1)/T104;
g1(5,4)=((y(7)-1)*params(5)*params(3)+y(7)*params(5)*params(3))/params(1);
g1(5,13)=(-T94);
g1(5,14)=(2*y(7)-2)*y(8)*params(5)/2;
g1(6,5)=(-((-y(12))*2*(1+y(5))))/(T111*T111);
g1(6,12)=(-1)/T111;
g1(7,5)=(-1)/((1+y(5))*(1+y(5)));
g1(7,15)=(-((-(y(16)*y(8)*params(3)))/(y(16)*y(15)*y(16)*y(15))));
g1(7,8)=(-(params(3)/(y(16)*y(15))));
g1(7,16)=(-((-(y(15)*y(8)*params(3)))/(y(16)*y(15)*y(16)*y(15))));
g1(8,6)=T31*params(4)*y(8)*params(2)/y(11);
g1(8,7)=((1+params(8))*(1-params(4))+params(4)*y(8)*params(2)*y(6)/y(11))*T174-(params(5)*y(7)+params(5)*(y(7)-1));
g1(8,15)=params(5)*params(3)*(y(15)-1)+y(15)*params(5)*params(3);
g1(8,8)=T31*params(4)*params(2)*y(6)/y(11);
g1(8,11)=T31*params(4)*(-(y(8)*params(2)*y(6)))/(y(11)*y(11));
g1(9,6)=(-y(11));
g1(9,7)=y(8)*params(5)/2*2*(y(7)-1);
g1(9,8)=T134;
g1(9,11)=(-y(6));
g1(10,2)=(-(params(6)*1/y(2)));
g1(10,11)=1/y(11);
g1(10,18)=(-1);
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],10,324);
end
if nargout >= 4,
%
% Third order derivatives
%

  g3 = sparse([],[],[],10,5832);
end
end
