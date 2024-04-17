function [residual, g1, g2] = chap3_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                     columns: variables in declaration order
%                                                     rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 15, 1);

%
% Model equations
%

T17 = (-1)/params(2);
lhs =y(1);
rhs =y(1)*params(1)+params(9)*y(2);
residual(1)= lhs-rhs;
lhs =y(2);
rhs =y(2)+T17*(y(7)-y(1)-y(5));
residual(2)= lhs-rhs;
lhs =y(7);
rhs =y(1)*params(7)+y(2)*params(8)+y(10);
residual(3)= lhs-rhs;
residual(4) = y(5);
lhs =y(6);
rhs =y(7)-y(1);
residual(5)= lhs-rhs;
lhs =y(3);
rhs =params(3)*y(11);
residual(6)= lhs-rhs;
lhs =y(2);
rhs =y(4)-y(3);
residual(7)= lhs-rhs;
lhs =y(10);
rhs =y(10)*params(4)+x(2);
residual(8)= lhs-rhs;
lhs =y(11);
rhs =y(11)*params(5)+x(1);
residual(9)= lhs-rhs;
lhs =y(4);
rhs =y(11)+(1-params(10))*y(8);
residual(10)= lhs-rhs;
lhs =y(9);
rhs =y(1)*4;
residual(11)= lhs-rhs;
lhs =y(13);
rhs =y(7)*4;
residual(12)= lhs-rhs;
lhs =y(12);
rhs =y(6)*4;
residual(13)= lhs-rhs;
lhs =y(14);
rhs =y(5)*4;
residual(14)= lhs-rhs;
lhs =y(15);
rhs =y(1)*4;
residual(15)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(15, 15);

  %
  % Jacobian matrix
  %

  g1(1,1)=1-params(1);
  g1(1,2)=(-params(9));
  g1(2,1)=T17;
  g1(2,5)=T17;
  g1(2,7)=(-T17);
  g1(3,1)=(-params(7));
  g1(3,2)=(-params(8));
  g1(3,7)=1;
  g1(3,10)=(-1);
  g1(4,5)=1;
  g1(5,1)=1;
  g1(5,6)=1;
  g1(5,7)=(-1);
  g1(6,3)=1;
  g1(6,11)=(-params(3));
  g1(7,2)=1;
  g1(7,3)=1;
  g1(7,4)=(-1);
  g1(8,10)=1-params(4);
  g1(9,11)=1-params(5);
  g1(10,4)=1;
  g1(10,8)=(-(1-params(10)));
  g1(10,11)=(-1);
  g1(11,1)=(-4);
  g1(11,9)=1;
  g1(12,7)=(-4);
  g1(12,13)=1;
  g1(13,6)=(-4);
  g1(13,12)=1;
  g1(14,5)=(-4);
  g1(14,14)=1;
  g1(15,1)=(-4);
  g1(15,15)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],15,225);
end
end
