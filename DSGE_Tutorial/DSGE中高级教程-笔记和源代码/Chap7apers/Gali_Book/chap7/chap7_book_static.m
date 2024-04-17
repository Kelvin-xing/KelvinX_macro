function [residual, g1, g2] = chap7_book_static(y, x, params)
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

residual = zeros( 23, 1);

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
rhs =y(1)+y(11)/params(18);
residual(6)= lhs-rhs;
lhs =y(6);
rhs =y(7)+y(3);
residual(7)= lhs-rhs;
lhs =y(9);
rhs =y(9)*params(1)+y(3)*params(16);
residual(8)= lhs-rhs;
lhs =y(3);
rhs =y(3)-1/params(18)*(y(5)-y(9)-y(4));
residual(9)= lhs-rhs;
lhs =y(7);
rhs =params(20)*y(2)+y(1)*params(21);
residual(10)= lhs-rhs;
lhs =y(4);
rhs =y(2)*params(20)*(-params(18))*(1-params(8));
residual(11)= lhs-rhs;
lhs =y(17);
rhs =y(11)*params(3)*(params(17)/params(2)-1);
residual(12)= lhs-rhs;
lhs =y(16);
rhs =y(3)*(params(18)+params(6));
residual(13)= lhs-rhs;
lhs =y(5);
rhs =y(4)+y(9)*params(10)+y(3)*params(11)+y(23);
residual(14)= lhs-rhs;
lhs =y(2);
rhs =y(2)*params(8)+x(1)+params(15)*x(3);
residual(15)= lhs-rhs;
lhs =y(1);
rhs =y(1)-(y(18)-y(10))/params(2);
residual(16)= lhs-rhs;
lhs =y(19);
rhs =y(1)*(params(6)+params(2))-(1+params(6))*y(20);
residual(17)= lhs-rhs;
lhs =y(10);
rhs =y(10)*params(1)+y(19)*params(19);
residual(18)= lhs-rhs;
lhs =y(18);
rhs =y(10)*params(12)+y(20)*params(13);
residual(19)= lhs-rhs;
lhs =y(20);
rhs =x(3)+y(20)*params(14);
residual(20)= lhs-rhs;
lhs =y(22);
rhs =y(6)*params(2)+params(6)*y(21);
residual(21)= lhs-rhs;
lhs =y(6);
rhs =y(2)+y(21);
residual(22)= lhs-rhs;
lhs =y(23);
rhs =y(23)*params(23)+x(4);
residual(23)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(23, 23);

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
  g1(6,11)=(-(1/params(18)));
  g1(7,3)=(-1);
  g1(7,6)=1;
  g1(7,7)=(-1);
  g1(8,3)=(-params(16));
  g1(8,9)=1-params(1);
  g1(9,4)=(-(1/params(18)));
  g1(9,5)=1/params(18);
  g1(9,9)=(-(1/params(18)));
  g1(10,1)=(-params(21));
  g1(10,2)=(-params(20));
  g1(10,7)=1;
  g1(11,2)=(-(params(20)*(-params(18))*(1-params(8))));
  g1(11,4)=1;
  g1(12,11)=(-(params(3)*(params(17)/params(2)-1)));
  g1(12,17)=1;
  g1(13,3)=(-(params(18)+params(6)));
  g1(13,16)=1;
  g1(14,3)=(-params(11));
  g1(14,4)=(-1);
  g1(14,5)=1;
  g1(14,9)=(-params(10));
  g1(14,23)=(-1);
  g1(15,2)=1-params(8);
  g1(16,10)=(-1)/params(2);
  g1(16,18)=1/params(2);
  g1(17,1)=(-(params(6)+params(2)));
  g1(17,19)=1;
  g1(17,20)=1+params(6);
  g1(18,10)=1-params(1);
  g1(18,19)=(-params(19));
  g1(19,10)=(-params(12));
  g1(19,18)=1;
  g1(19,20)=(-params(13));
  g1(20,20)=1-params(14);
  g1(21,6)=(-params(2));
  g1(21,21)=(-params(6));
  g1(21,22)=1;
  g1(22,2)=(-1);
  g1(22,6)=1;
  g1(22,21)=(-1);
  g1(23,23)=1-params(23);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],23,529);
end
end
