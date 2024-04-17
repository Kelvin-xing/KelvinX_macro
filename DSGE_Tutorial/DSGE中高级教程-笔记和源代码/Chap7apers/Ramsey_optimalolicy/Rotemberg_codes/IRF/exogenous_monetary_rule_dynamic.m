function [residual, g1, g2, g3] = exogenous_monetary_rule_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [M_.exo_nbr by nperiods] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(7, 1);
T59 = 1+params(5)*(y(4)-1)^2/2;
T72 = 1+(y(4)-1)^2*params(5)/2;
lhs =y(6);
rhs =log(y(5))-params(2)*y(3)^2/2;
residual(1)= lhs-rhs;
lhs =y(7);
rhs =y(6)+params(1)*y(11);
residual(2)= lhs-rhs;
lhs =y(2);
rhs =params(9)/params(3)-1+params(7)*(y(9)-params(9));
residual(3)= lhs-rhs;
lhs =1/(1+y(2));
rhs =y(5)*params(3)/(y(9)*y(10));
residual(4)= lhs-rhs;
residual(5) = ((1+params(8))*(1-params(4))+params(4)*y(5)*params(2)*y(3)/y(8))*T59-y(4)*params(5)*(y(4)-1)+y(9)*params(3)*params(5)*(y(9)-1);
residual(6) = y(5)*T72-y(3)*y(8);
lhs =log(y(8));
rhs =params(6)*log(y(1))+x(it_, 1);
residual(7)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(7, 12);

  %
  % Jacobian matrix
  %

  g1(1,3)=params(2)*2*y(3)/2;
  g1(1,5)=(-(1/y(5)));
  g1(1,6)=1;
  g1(2,6)=(-1);
  g1(2,7)=1;
  g1(2,11)=(-params(1));
  g1(3,2)=1;
  g1(3,9)=(-params(7));
  g1(4,2)=(-1)/((1+y(2))*(1+y(2)));
  g1(4,9)=(-((-(y(5)*params(3)*y(10)))/(y(9)*y(10)*y(9)*y(10))));
  g1(4,5)=(-(params(3)/(y(9)*y(10))));
  g1(4,10)=(-((-(y(9)*y(5)*params(3)))/(y(9)*y(10)*y(9)*y(10))));
  g1(5,3)=T59*params(4)*y(5)*params(2)/y(8);
  g1(5,4)=((1+params(8))*(1-params(4))+params(4)*y(5)*params(2)*y(3)/y(8))*params(5)*2*(y(4)-1)/2-(params(5)*(y(4)-1)+params(5)*y(4));
  g1(5,9)=params(3)*params(5)*(y(9)-1)+y(9)*params(3)*params(5);
  g1(5,5)=T59*params(4)*params(2)*y(3)/y(8);
  g1(5,8)=T59*params(4)*(-(y(5)*params(2)*y(3)))/(y(8)*y(8));
  g1(6,3)=(-y(8));
  g1(6,4)=y(5)*params(5)/2*2*(y(4)-1);
  g1(6,5)=T72;
  g1(6,8)=(-y(3));
  g1(7,1)=(-(params(6)*1/y(1)));
  g1(7,8)=1/y(8);
  g1(7,12)=(-1);
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],7,144);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],7,1728);
end
end
