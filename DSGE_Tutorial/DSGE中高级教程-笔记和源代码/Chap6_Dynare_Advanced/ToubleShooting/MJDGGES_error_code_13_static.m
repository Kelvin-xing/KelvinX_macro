function [residual, g1, g2] = MJDGGES_error_code_13_static(y, x, params)
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

residual = zeros( 17, 1);

%
% Model equations
%

lhs =y(8);
rhs =y(9);
residual(1)= lhs-rhs;
lhs =y(14);
rhs =y(9)+y(14);
residual(2)= lhs-rhs;
lhs =y(15);
rhs =y(8)+y(15);
residual(3)= lhs-rhs;
lhs =y(12);
rhs =y(11)*(1-params(3));
residual(4)= lhs-rhs;
lhs =0;
rhs =y(10)-y(9);
residual(5)= lhs-rhs;
lhs =y(6);
rhs =y(1)+y(11)/params(14);
residual(6)= lhs-rhs;
lhs =y(6);
rhs =y(7)+y(3);
residual(7)= lhs-rhs;
lhs =y(9);
rhs =y(9)*params(1)+y(3)*params(12);
residual(8)= lhs-rhs;
lhs =y(3);
rhs =y(3)-1/params(14)*(y(5)-y(9)-y(4));
residual(9)= lhs-rhs;
lhs =y(7);
rhs =params(16)*y(2)+y(1)*params(17);
residual(10)= lhs-rhs;
lhs =y(4);
rhs =y(2)*params(16)*(-params(14))*(1-params(8));
residual(11)= lhs-rhs;
lhs =y(17);
rhs =y(11)*params(3)*(params(13)/params(2)-1);
residual(12)= lhs-rhs;
lhs =y(16);
rhs =y(3)*(params(14)+params(6));
residual(13)= lhs-rhs;
residual(14) = y(9);
lhs =y(2);
rhs =y(2)*params(8)+x(1);
residual(15)= lhs-rhs;
lhs =y(1);
rhs =y(1)*params(9)+x(2);
residual(16)= lhs-rhs;
lhs =y(6);
rhs =y(6)-1/params(2)*(y(5)-y(8));
residual(17)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(17, 17);

  %
  % Jacobian matrix
  %

  g1(1,8)=1;
  g1(1,9)=(-1);
  g1(2,9)=(-1);
  g1(3,8)=(-1);
  g1(4,11)=(-(1-params(3)));
  g1(4,12)=1;
  g1(5,9)=1;
  g1(5,10)=(-1);
  g1(6,1)=(-1);
  g1(6,6)=1;
  g1(6,11)=(-(1/params(14)));
  g1(7,3)=(-1);
  g1(7,6)=1;
  g1(7,7)=(-1);
  g1(8,3)=(-params(12));
  g1(8,9)=1-params(1);
  g1(9,4)=(-(1/params(14)));
  g1(9,5)=1/params(14);
  g1(9,9)=(-(1/params(14)));
  g1(10,1)=(-params(17));
  g1(10,2)=(-params(16));
  g1(10,7)=1;
  g1(11,2)=(-(params(16)*(-params(14))*(1-params(8))));
  g1(11,4)=1;
  g1(12,11)=(-(params(3)*(params(13)/params(2)-1)));
  g1(12,17)=1;
  g1(13,3)=(-(params(14)+params(6)));
  g1(13,16)=1;
  g1(14,9)=1;
  g1(15,2)=1-params(8);
  g1(16,1)=1-params(9);
  g1(17,5)=1/params(2);
  g1(17,8)=(-(1/params(2)));
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
