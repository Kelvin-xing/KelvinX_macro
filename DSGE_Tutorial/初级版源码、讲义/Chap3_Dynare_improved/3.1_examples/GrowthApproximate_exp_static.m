function [residual, g1, g2] = GrowthApproximate_exp_static(y, x, params)
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

residual = zeros( 6, 1);

%
% Model equations
%

mar_c__ = (exp(y(2))^params(2)*(1-exp(y(5)))^(1-params(2)))^(1-params(5));
mar_c1__ = (exp(y(2))^params(2)*(1-exp(y(5)))^(1-params(2)))^(1-params(5));
T26 = params(1)*mar_c1__/exp(y(2));
T37 = exp(y(5))^(1-params(4));
T38 = params(4)*exp(y(6))*exp(y(3))^(params(4)-1)*T37;
T41 = 1+T38-params(3);
T48 = exp(y(6))*(1-params(4))*params(2)/(1-params(2))*exp(y(3))^params(4);
T51 = T48*exp(y(5))^(-params(4));
T75 = getPowerDeriv(exp(y(2))^params(2)*(1-exp(y(5)))^(1-params(2)),1-params(5),1);
T77 = exp(y(2))*(1-exp(y(5)))^(1-params(2))*exp(y(2))*getPowerDeriv(exp(y(2)),params(2),1)*T75;
T111 = T75*exp(y(2))^params(2)*(-exp(y(5)))*getPowerDeriv(1-exp(y(5)),1-params(2),1)/exp(y(2));
lhs =mar_c__/exp(y(2));
rhs =T26*T41;
residual(1)= lhs-rhs;
lhs =exp(y(2));
rhs =(1-exp(y(5)))*T51;
residual(2)= lhs-rhs;
lhs =exp(y(3));
rhs =exp(y(4))+exp(y(3))*(1-params(3));
residual(3)= lhs-rhs;
lhs =exp(y(1));
rhs =T37*exp(y(6))*exp(y(3))^params(4);
residual(4)= lhs-rhs;
lhs =exp(y(1));
rhs =exp(y(2))+exp(y(4));
residual(5)= lhs-rhs;
lhs =y(6);
rhs =y(6)*params(6)+x(1);
residual(6)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(6, 6);

  %
  % Jacobian matrix
  %

  g1(1,2)=(T77-exp(y(2))*mar_c__)/(exp(y(2))*exp(y(2)))-T41*params(1)*(T77-exp(y(2))*mar_c1__)/(exp(y(2))*exp(y(2)));
  g1(1,3)=(-(T26*T37*params(4)*exp(y(6))*exp(y(3))*getPowerDeriv(exp(y(3)),params(4)-1,1)));
  g1(1,5)=T111-(T41*params(1)*T111+T26*params(4)*exp(y(6))*exp(y(3))^(params(4)-1)*exp(y(5))*getPowerDeriv(exp(y(5)),1-params(4),1));
  g1(1,6)=(-(T26*T38));
  g1(2,2)=exp(y(2));
  g1(2,3)=(-((1-exp(y(5)))*exp(y(5))^(-params(4))*exp(y(6))*(1-params(4))*params(2)/(1-params(2))*exp(y(3))*getPowerDeriv(exp(y(3)),params(4),1)));
  g1(2,5)=(-(T51*(-exp(y(5)))+(1-exp(y(5)))*T48*exp(y(5))*getPowerDeriv(exp(y(5)),(-params(4)),1)));
  g1(2,6)=(-((1-exp(y(5)))*T51));
  g1(3,3)=exp(y(3))-exp(y(3))*(1-params(3));
  g1(3,4)=(-exp(y(4)));
  g1(4,1)=exp(y(1));
  g1(4,3)=(-(T37*exp(y(6))*exp(y(3))*getPowerDeriv(exp(y(3)),params(4),1)));
  g1(4,5)=(-(exp(y(6))*exp(y(3))^params(4)*exp(y(5))*getPowerDeriv(exp(y(5)),1-params(4),1)));
  g1(4,6)=(-(T37*exp(y(6))*exp(y(3))^params(4)));
  g1(5,1)=exp(y(1));
  g1(5,2)=(-exp(y(2)));
  g1(5,4)=(-exp(y(4)));
  g1(6,6)=1-params(6);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],6,36);
end
end
