function [residual, g1, g2, g3] = recursively_running_dynare_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(9, 1);
T12 = exp(y(11))^(-params(5));
T17 = params(2)*exp(y(13))^(-params(5));
T50 = exp(y(6))^(1-params(1));
T59 = T50*exp(y(9))*params(1)*exp(y(1))^(params(1)-1);
T65 = exp(y(1))^params(1)*exp(y(9))*(1-params(1))*exp(y(6))^(-params(1));
T119 = exp(y(11))*getPowerDeriv(exp(y(11)),(-params(5)),1);
lhs =T12;
rhs =T17*(exp(y(12))+1-params(3));
residual(1)= lhs-rhs;
lhs =T12;
rhs =exp(y(13))^(-params(5))*params(2)*(1+y(5));
residual(2)= lhs-rhs;
lhs =params(6)*exp(y(6))^params(7);
rhs =T12*exp(y(7));
residual(3)= lhs-rhs;
lhs =exp(y(3));
rhs =exp(y(9))*exp(y(1))^params(1)*T50;
residual(4)= lhs-rhs;
lhs =exp(y(10));
rhs =T59;
residual(5)= lhs-rhs;
lhs =exp(y(7));
rhs =T65;
residual(6)= lhs-rhs;
lhs =exp(y(3));
rhs =exp(y(11))+exp(y(4));
residual(7)= lhs-rhs;
lhs =exp(y(8));
rhs =exp(y(4))+(1-params(3))*exp(y(1));
residual(8)= lhs-rhs;
lhs =y(9);
rhs =params(4)*y(2)+x(it_, 1);
residual(9)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(9, 14);

  %
  % Jacobian matrix
  %

  g1(1,12)=(-(T17*exp(y(12))));
  g1(1,11)=T119;
  g1(1,13)=(-((exp(y(12))+1-params(3))*params(2)*exp(y(13))*getPowerDeriv(exp(y(13)),(-params(5)),1)));
  g1(2,5)=(-T17);
  g1(2,11)=T119;
  g1(2,13)=(-(params(2)*(1+y(5))*exp(y(13))*getPowerDeriv(exp(y(13)),(-params(5)),1)));
  g1(3,6)=params(6)*exp(y(6))*getPowerDeriv(exp(y(6)),params(7),1);
  g1(3,7)=(-(T12*exp(y(7))));
  g1(3,11)=(-(exp(y(7))*T119));
  g1(4,3)=exp(y(3));
  g1(4,6)=(-(exp(y(9))*exp(y(1))^params(1)*exp(y(6))*getPowerDeriv(exp(y(6)),1-params(1),1)));
  g1(4,1)=(-(T50*exp(y(9))*exp(y(1))*getPowerDeriv(exp(y(1)),params(1),1)));
  g1(4,9)=(-(exp(y(9))*exp(y(1))^params(1)*T50));
  g1(5,6)=(-(exp(y(9))*params(1)*exp(y(1))^(params(1)-1)*exp(y(6))*getPowerDeriv(exp(y(6)),1-params(1),1)));
  g1(5,1)=(-(T50*exp(y(9))*params(1)*exp(y(1))*getPowerDeriv(exp(y(1)),params(1)-1,1)));
  g1(5,9)=(-T59);
  g1(5,10)=exp(y(10));
  g1(6,6)=(-(exp(y(1))^params(1)*exp(y(9))*(1-params(1))*exp(y(6))*getPowerDeriv(exp(y(6)),(-params(1)),1)));
  g1(6,7)=exp(y(7));
  g1(6,1)=(-(exp(y(6))^(-params(1))*exp(y(9))*(1-params(1))*exp(y(1))*getPowerDeriv(exp(y(1)),params(1),1)));
  g1(6,9)=(-T65);
  g1(7,3)=exp(y(3));
  g1(7,4)=(-exp(y(4)));
  g1(7,11)=(-exp(y(11)));
  g1(8,4)=(-exp(y(4)));
  g1(8,1)=(-((1-params(3))*exp(y(1))));
  g1(8,8)=exp(y(8));
  g1(9,2)=(-params(4));
  g1(9,9)=1;
  g1(9,14)=(-1);
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],9,196);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],9,2744);
end
end
