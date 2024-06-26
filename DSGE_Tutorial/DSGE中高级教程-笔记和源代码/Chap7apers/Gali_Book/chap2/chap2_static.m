function [residual, g1, g2] = chap2_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 9, 1);

%
% Model equations
%

lhs =exp(y(2));
rhs =exp(params(4)*y(1)+params(5)*y(5));
residual(1)= lhs-rhs;
lhs =1/exp(y(6));
rhs =params(2)*(1-y(3));
residual(2)= lhs-rhs;
lhs =exp(y(4)+y(5)*(1-params(1)));
rhs =exp(y(1));
residual(3)= lhs-rhs;
lhs =exp(y(2));
rhs =(1-params(1))*exp(y(4)-y(5)*params(1));
residual(4)= lhs-rhs;
lhs =exp(y(7));
rhs =exp(y(6)-y(3));
residual(5)= lhs-rhs;
lhs =exp(y(6));
rhs =1/params(2)+y(3)*params(6)+x(2);
residual(6)= lhs-rhs;
lhs =exp(y(1));
rhs =exp(y(8));
residual(7)= lhs-rhs;
lhs =y(4);
rhs =y(4)*params(3)+x(1);
residual(8)= lhs-rhs;
lhs =y(9);
rhs =y(3)*4;
residual(9)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(9, 9);

%
% Jacobian matrix
%

  g1(1,1)=(-(params(4)*exp(params(4)*y(1)+params(5)*y(5))));
  g1(1,2)=exp(y(2));
  g1(1,5)=(-(params(5)*exp(params(4)*y(1)+params(5)*y(5))));
  g1(2,3)=params(2);
  g1(2,6)=(-exp(y(6)))/(exp(y(6))*exp(y(6)));
  g1(3,1)=(-exp(y(1)));
  g1(3,4)=exp(y(4)+y(5)*(1-params(1)));
  g1(3,5)=(1-params(1))*exp(y(4)+y(5)*(1-params(1)));
  g1(4,2)=exp(y(2));
  g1(4,4)=(-((1-params(1))*exp(y(4)-y(5)*params(1))));
  g1(4,5)=(-((1-params(1))*exp(y(4)-y(5)*params(1))*(-params(1))));
  g1(5,3)=exp(y(6)-y(3));
  g1(5,6)=(-exp(y(6)-y(3)));
  g1(5,7)=exp(y(7));
  g1(6,3)=(-params(6));
  g1(6,6)=exp(y(6));
  g1(7,1)=exp(y(1));
  g1(7,8)=(-exp(y(8)));
  g1(8,4)=1-params(3);
  g1(9,3)=(-4);
  g1(9,9)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],9,81);
end
end
