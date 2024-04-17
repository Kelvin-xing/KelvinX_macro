function [residual, g1, g2] = NK_sticky_price_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 17, 1);

%
% Model equations
%

T20 = exp((-params(1))*y(1))*params(2)*exp(y(2))/exp(y(3));
T105 = exp(y(3))*params(6)/(params(6)-1)*exp(y(12))/exp(y(13));
T131 = ((params(6)-1)/params(6)/params(3))^(1/(params(1)+params(4)));
lhs =exp((-params(1))*y(1));
rhs =T20;
residual(1)= lhs-rhs;
lhs =params(3)*exp(params(4)*y(4));
rhs =exp((-params(1))*y(1))*exp(y(5));
residual(2)= lhs-rhs;
lhs =exp(y(7));
rhs =exp(y(5))/exp(y(8));
residual(3)= lhs-rhs;
lhs =exp(y(6));
rhs =exp(y(2))*params(7)*exp(params(1)*y(1))/(exp(y(2))-1);
residual(4)= lhs-rhs;
lhs =y(14);
rhs =(1-params(9))*(params(14)-1)-y(3)+y(14)*params(9)+y(3)*params(9)+x(2);
residual(5)= lhs-rhs;
residual(6) = y(14);
lhs =exp(y(1));
rhs =exp(y(9));
residual(7)= lhs-rhs;
lhs =exp(y(9));
rhs =exp(y(8))*exp(y(4))/exp(y(10));
residual(8)= lhs-rhs;
lhs =exp(y(10));
rhs =(1-params(5))*exp((-params(6))*y(11))*exp(y(3)*params(6))+exp(y(10))*params(5)*exp(y(3)*params(6));
residual(9)= lhs-rhs;
lhs =exp(y(3)*(1-params(6)));
rhs =params(5)+(1-params(5))*exp(y(11)*(1-params(6)));
residual(10)= lhs-rhs;
lhs =exp(y(11));
rhs =T105;
residual(11)= lhs-rhs;
lhs =exp(y(12));
rhs =exp(y(7))*exp((-params(1))*y(1))*exp(y(9))+exp(y(12))*exp(y(3)*params(6))*params(2)*params(5);
residual(12)= lhs-rhs;
lhs =exp(y(13));
rhs =exp((-params(1))*y(1))*exp(y(9))+exp(y(13))*params(2)*params(5)*exp(y(3)*(params(6)-1));
residual(13)= lhs-rhs;
lhs =y(8);
rhs =y(8)*params(8)+x(1);
residual(14)= lhs-rhs;
lhs =exp(y(15));
rhs =T131*exp(y(8)*(1+params(4))/(params(1)+params(4)));
residual(15)= lhs-rhs;
lhs =y(16);
rhs =y(9)-y(15);
residual(16)= lhs-rhs;
lhs =y(17);
rhs =exp(y(2))-exp(y(3));
residual(17)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(17, 17);

%
% Jacobian matrix
%

  g1(1,1)=(-params(1))*exp((-params(1))*y(1))-exp(y(2))*params(2)*(-params(1))*exp((-params(1))*y(1))/exp(y(3));
  g1(1,2)=(-T20);
  g1(1,3)=(-((-(exp((-params(1))*y(1))*params(2)*exp(y(2))*exp(y(3))))/(exp(y(3))*exp(y(3)))));
  g1(2,1)=(-(exp(y(5))*(-params(1))*exp((-params(1))*y(1))));
  g1(2,4)=params(3)*params(4)*exp(params(4)*y(4));
  g1(2,5)=(-(exp((-params(1))*y(1))*exp(y(5))));
  g1(3,5)=(-(exp(y(5))/exp(y(8))));
  g1(3,7)=exp(y(7));
  g1(3,8)=(-((-(exp(y(5))*exp(y(8))))/(exp(y(8))*exp(y(8)))));
  g1(4,1)=(-(exp(y(2))*params(7)*params(1)*exp(params(1)*y(1))/(exp(y(2))-1)));
  g1(4,2)=(-((exp(y(2))*params(7)*exp(params(1)*y(1))*(exp(y(2))-1)-exp(y(2))*exp(y(2))*params(7)*exp(params(1)*y(1)))/((exp(y(2))-1)*(exp(y(2))-1))));
  g1(4,6)=exp(y(6));
  g1(5,3)=(-((-1)+params(9)));
  g1(5,14)=1-params(9);
  g1(6,14)=1;
  g1(7,1)=exp(y(1));
  g1(7,9)=(-exp(y(9)));
  g1(8,4)=(-(exp(y(8))*exp(y(4))/exp(y(10))));
  g1(8,8)=(-(exp(y(8))*exp(y(4))/exp(y(10))));
  g1(8,9)=exp(y(9));
  g1(8,10)=(-((-(exp(y(8))*exp(y(4))*exp(y(10))))/(exp(y(10))*exp(y(10)))));
  g1(9,3)=(-((1-params(5))*exp((-params(6))*y(11))*params(6)*exp(y(3)*params(6))+exp(y(10))*params(5)*params(6)*exp(y(3)*params(6))));
  g1(9,10)=exp(y(10))-exp(y(10))*params(5)*exp(y(3)*params(6));
  g1(9,11)=(-(exp(y(3)*params(6))*(1-params(5))*(-params(6))*exp((-params(6))*y(11))));
  g1(10,3)=(1-params(6))*exp(y(3)*(1-params(6)));
  g1(10,11)=(-((1-params(5))*(1-params(6))*exp(y(11)*(1-params(6)))));
  g1(11,3)=(-T105);
  g1(11,11)=exp(y(11));
  g1(11,12)=(-T105);
  g1(11,13)=(-((-(exp(y(3))*params(6)/(params(6)-1)*exp(y(12))*exp(y(13))))/(exp(y(13))*exp(y(13)))));
  g1(12,1)=(-(exp(y(7))*exp(y(9))*(-params(1))*exp((-params(1))*y(1))));
  g1(12,3)=(-(exp(y(12))*params(2)*params(5)*params(6)*exp(y(3)*params(6))));
  g1(12,7)=(-(exp(y(7))*exp((-params(1))*y(1))*exp(y(9))));
  g1(12,9)=(-(exp(y(7))*exp((-params(1))*y(1))*exp(y(9))));
  g1(12,12)=exp(y(12))-exp(y(12))*exp(y(3)*params(6))*params(2)*params(5);
  g1(13,1)=(-(exp(y(9))*(-params(1))*exp((-params(1))*y(1))));
  g1(13,3)=(-(exp(y(13))*params(2)*params(5)*(params(6)-1)*exp(y(3)*(params(6)-1))));
  g1(13,9)=(-(exp((-params(1))*y(1))*exp(y(9))));
  g1(13,13)=exp(y(13))-exp(y(13))*params(2)*params(5)*exp(y(3)*(params(6)-1));
  g1(14,8)=1-params(8);
  g1(15,8)=(-(T131*(1+params(4))/(params(1)+params(4))*exp(y(8)*(1+params(4))/(params(1)+params(4)))));
  g1(15,15)=exp(y(15));
  g1(16,9)=(-1);
  g1(16,15)=1;
  g1(16,16)=1;
  g1(17,2)=(-exp(y(2)));
  g1(17,3)=exp(y(3));
  g1(17,17)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],17,289);
end
end
