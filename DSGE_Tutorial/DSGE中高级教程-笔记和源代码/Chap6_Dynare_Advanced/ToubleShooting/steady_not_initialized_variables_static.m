function [residual, g1, g2] = steady_not_initialized_variables_static(y, x, params)
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

residual = zeros( 7, 1);

%
% Model equations
%

T56 = 1+params(5)*(y(3)-1)^2/2;
T68 = 1+(y(3)-1)^2*params(5)/2;
lhs =y(5);
rhs =log(y(4))-params(2)*y(2)^2/2;
residual(1)= lhs-rhs;
lhs =y(6);
rhs =y(5)+y(6)*params(1);
residual(2)= lhs-rhs;
lhs =y(1);
rhs =params(9)/params(3)-1+params(7)*(y(3)-params(9));
residual(3)= lhs-rhs;
lhs =1/(1+y(1));
rhs =y(4)*params(3)/(y(4)*y(3));
residual(4)= lhs-rhs;
residual(5) = ((1+params(8))*(1-params(4))+params(4)*y(4)*params(2)*y(2)/y(7))*T56-y(3)*params(5)*(y(3)-1)+y(3)*(y(3)-1)*params(3)*params(5);
residual(6) = y(4)*T68-y(2)*y(7);
lhs =log(y(7));
rhs =log(y(7))*params(6)+x(1);
residual(7)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(7, 7);

  %
  % Jacobian matrix
  %

  g1(1,2)=params(2)*2*y(2)/2;
  g1(1,4)=(-(1/y(4)));
  g1(1,5)=1;
  g1(2,5)=(-1);
  g1(2,6)=1-params(1);
  g1(3,1)=1;
  g1(3,3)=(-params(7));
  g1(4,1)=(-1)/((1+y(1))*(1+y(1)));
  g1(4,3)=(-((-(y(4)*y(4)*params(3)))/(y(4)*y(3)*y(4)*y(3))));
  g1(4,4)=(-((params(3)*y(4)*y(3)-y(3)*y(4)*params(3))/(y(4)*y(3)*y(4)*y(3))));
  g1(5,2)=T56*params(4)*y(4)*params(2)/y(7);
  g1(5,3)=((1+params(8))*(1-params(4))+params(4)*y(4)*params(2)*y(2)/y(7))*params(5)*2*(y(3)-1)/2-(params(5)*(y(3)-1)+y(3)*params(5))+(y(3)-1)*params(3)*params(5)+y(3)*params(3)*params(5);
  g1(5,4)=T56*params(4)*params(2)*y(2)/y(7);
  g1(5,7)=T56*params(4)*(-(y(4)*params(2)*y(2)))/(y(7)*y(7));
  g1(6,2)=(-y(7));
  g1(6,3)=y(4)*params(5)/2*2*(y(3)-1);
  g1(6,4)=T68;
  g1(6,7)=(-y(2));
  g1(7,7)=1/y(7)-params(6)*1/y(7);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],7,49);
end
end
