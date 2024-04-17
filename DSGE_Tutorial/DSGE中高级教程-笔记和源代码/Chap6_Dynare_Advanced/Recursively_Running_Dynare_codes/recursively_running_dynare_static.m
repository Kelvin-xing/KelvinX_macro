function [residual, g1, g2] = recursively_running_dynare_static(y, x, params)
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

residual = zeros( 9, 1);

%
% Model equations
%

T12 = exp(y(9))^(-params(5));
T47 = exp(y(4))^(1-params(1));
T54 = T47*exp(y(7))*params(1)*exp(y(6))^(params(1)-1);
T60 = exp(y(6))^params(1)*exp(y(7))*(1-params(1))*exp(y(4))^(-params(1));
T111 = exp(y(9))*getPowerDeriv(exp(y(9)),(-params(5)),1);
lhs =T12;
rhs =T12*params(2)*(exp(y(8))+1-params(3));
residual(1)= lhs-rhs;
lhs =T12;
rhs =T12*params(2)*(1+y(3));
residual(2)= lhs-rhs;
lhs =params(6)*exp(y(4))^params(7);
rhs =T12*exp(y(5));
residual(3)= lhs-rhs;
lhs =exp(y(1));
rhs =exp(y(7))*exp(y(6))^params(1)*T47;
residual(4)= lhs-rhs;
lhs =exp(y(8));
rhs =T54;
residual(5)= lhs-rhs;
lhs =exp(y(5));
rhs =T60;
residual(6)= lhs-rhs;
lhs =exp(y(1));
rhs =exp(y(9))+exp(y(2));
residual(7)= lhs-rhs;
lhs =exp(y(6));
rhs =exp(y(2))+(1-params(3))*exp(y(6));
residual(8)= lhs-rhs;
lhs =y(7);
rhs =y(7)*params(4)+x(1);
residual(9)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(9, 9);

  %
  % Jacobian matrix
  %

  g1(1,8)=(-(T12*params(2)*exp(y(8))));
  g1(1,9)=T111-(exp(y(8))+1-params(3))*params(2)*T111;
  g1(2,3)=(-(T12*params(2)));
  g1(2,9)=T111-params(2)*(1+y(3))*T111;
  g1(3,4)=params(6)*exp(y(4))*getPowerDeriv(exp(y(4)),params(7),1);
  g1(3,5)=(-(T12*exp(y(5))));
  g1(3,9)=(-(exp(y(5))*T111));
  g1(4,1)=exp(y(1));
  g1(4,4)=(-(exp(y(7))*exp(y(6))^params(1)*exp(y(4))*getPowerDeriv(exp(y(4)),1-params(1),1)));
  g1(4,6)=(-(T47*exp(y(7))*exp(y(6))*getPowerDeriv(exp(y(6)),params(1),1)));
  g1(4,7)=(-(exp(y(7))*exp(y(6))^params(1)*T47));
  g1(5,4)=(-(exp(y(7))*params(1)*exp(y(6))^(params(1)-1)*exp(y(4))*getPowerDeriv(exp(y(4)),1-params(1),1)));
  g1(5,6)=(-(T47*exp(y(7))*params(1)*exp(y(6))*getPowerDeriv(exp(y(6)),params(1)-1,1)));
  g1(5,7)=(-T54);
  g1(5,8)=exp(y(8));
  g1(6,4)=(-(exp(y(6))^params(1)*exp(y(7))*(1-params(1))*exp(y(4))*getPowerDeriv(exp(y(4)),(-params(1)),1)));
  g1(6,5)=exp(y(5));
  g1(6,6)=(-(exp(y(4))^(-params(1))*exp(y(7))*(1-params(1))*exp(y(6))*getPowerDeriv(exp(y(6)),params(1),1)));
  g1(6,7)=(-T60);
  g1(7,1)=exp(y(1));
  g1(7,2)=(-exp(y(2)));
  g1(7,9)=(-exp(y(9)));
  g1(8,2)=(-exp(y(2)));
  g1(8,6)=exp(y(6))-(1-params(3))*exp(y(6));
  g1(9,7)=1-params(4);
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
