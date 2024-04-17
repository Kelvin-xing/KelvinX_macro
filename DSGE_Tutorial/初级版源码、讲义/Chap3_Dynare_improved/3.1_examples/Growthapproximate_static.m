function [residual, g1, g2] = Growthapproximate_static(y, x, params)
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
%                                                    columns: equations in order of declaration
%                                                    rows: variables in declaration order
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: equations in order of declaration
%                                                       rows: variables in declaration order
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 4, 1);

%
% Model equations
%

mar_c__ = (y(1)^params(2)*(1-y(3))^(1-params(2)))^(1-params(5));
mar_c1__ = (y(1)^params(2)*(1-y(3))^(1-params(2)))^(1-params(5));
T24 = params(1)*mar_c1__/y(1);
T34 = y(3)^(1-params(4));
T38 = 1+params(4)*exp(y(4))*y(2)^(params(4)-1)*T34-params(3);
T45 = exp(y(4))*(1-params(4))*params(2)/(1-params(2))*y(2)^params(4);
T48 = T45*y(3)^(-params(4));
T65 = getPowerDeriv(y(1)^params(2)*(1-y(3))^(1-params(2)),1-params(5),1);
T67 = y(1)*(1-y(3))^(1-params(2))*getPowerDeriv(y(1),params(2),1)*T65;
T94 = T65*y(1)^params(2)*(-(getPowerDeriv(1-y(3),1-params(2),1)))/y(1);
lhs =mar_c__/y(1);
rhs =T24*T38;
residual(1)= lhs-rhs;
lhs =y(1);
rhs =(1-y(3))*T48;
residual(2)= lhs-rhs;
lhs =y(2);
rhs =T34*exp(y(4))*y(2)^params(4)-y(1)+y(2)*(1-params(3));
residual(3)= lhs-rhs;
lhs =y(4);
rhs =y(4)*params(6)+x(1);
residual(4)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(4, 4);

  %
  % Jacobian matrix
  %

  g1(1,1)=(T67-mar_c__)/(y(1)*y(1))-T38*params(1)*(T67-mar_c1__)/(y(1)*y(1));
  g1(1,2)=(-(T24*T34*params(4)*exp(y(4))*getPowerDeriv(y(2),params(4)-1,1)));
  g1(1,3)=T94-(T38*params(1)*T94+T24*params(4)*exp(y(4))*y(2)^(params(4)-1)*getPowerDeriv(y(3),1-params(4),1));
  g1(1,4)=(-(T24*params(4)*exp(y(4))*y(2)^(params(4)-1)*T34));
  g1(2,1)=1;
  g1(2,2)=(-((1-y(3))*y(3)^(-params(4))*exp(y(4))*(1-params(4))*params(2)/(1-params(2))*getPowerDeriv(y(2),params(4),1)));
  g1(2,3)=(-((-T48)+(1-y(3))*T45*getPowerDeriv(y(3),(-params(4)),1)));
  g1(2,4)=(-((1-y(3))*T48));
  g1(3,1)=1;
  g1(3,2)=1-(1-params(3)+T34*exp(y(4))*getPowerDeriv(y(2),params(4),1));
  g1(3,3)=(-(exp(y(4))*y(2)^params(4)*getPowerDeriv(y(3),1-params(4),1)));
  g1(3,4)=(-(T34*exp(y(4))*y(2)^params(4)));
  g1(4,4)=1-params(6);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],4,16);
end
end
