function [residual, g1, g2, g3] = chap2_m_real_dynamic(y, x, params, steady_state, it_)
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
lhs =exp(y(3));
rhs =exp(params(4)*y(2)+params(5)*y(6));
residual(1)= lhs-rhs;
lhs =1/exp(y(7));
rhs =params(2)*(exp((-params(4))*(y(11)-y(2)))-y(12));
residual(2)= lhs-rhs;
lhs =exp(y(5)+y(6)*(1-params(1)));
rhs =exp(y(2));
residual(3)= lhs-rhs;
lhs =exp(y(3));
rhs =(1-params(1))*exp(y(5)-y(6)*params(1));
residual(4)= lhs-rhs;
lhs =exp(y(8));
rhs =exp(y(7)-y(12));
residual(5)= lhs-rhs;
lhs =y(7);
rhs =log(1/params(2))+params(6)*y(4);
residual(6)= lhs-rhs;
lhs =exp(y(2));
rhs =exp(y(9));
residual(7)= lhs-rhs;
lhs =y(5);
rhs =params(3)*y(1)+x(it_, 1);
residual(8)= lhs-rhs;
lhs =y(10);
rhs =y(9)-y(7)*params(7)+x(it_, 2);
residual(9)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(9, 14);

  %
  % Jacobian matrix
  %

  g1(1,2)=(-(params(4)*exp(params(4)*y(2)+params(5)*y(6))));
  g1(1,3)=exp(y(3));
  g1(1,6)=(-(params(5)*exp(params(4)*y(2)+params(5)*y(6))));
  g1(2,2)=(-(params(2)*params(4)*exp((-params(4))*(y(11)-y(2)))));
  g1(2,11)=(-(params(2)*(-params(4))*exp((-params(4))*(y(11)-y(2)))));
  g1(2,12)=params(2);
  g1(2,7)=(-exp(y(7)))/(exp(y(7))*exp(y(7)));
  g1(3,2)=(-exp(y(2)));
  g1(3,5)=exp(y(5)+y(6)*(1-params(1)));
  g1(3,6)=(1-params(1))*exp(y(5)+y(6)*(1-params(1)));
  g1(4,3)=exp(y(3));
  g1(4,5)=(-((1-params(1))*exp(y(5)-y(6)*params(1))));
  g1(4,6)=(-((1-params(1))*exp(y(5)-y(6)*params(1))*(-params(1))));
  g1(5,12)=exp(y(7)-y(12));
  g1(5,7)=(-exp(y(7)-y(12)));
  g1(5,8)=exp(y(8));
  g1(6,4)=(-params(6));
  g1(6,7)=1;
  g1(7,2)=exp(y(2));
  g1(7,9)=(-exp(y(9)));
  g1(8,1)=(-params(3));
  g1(8,5)=1;
  g1(8,13)=(-1);
  g1(9,7)=params(7);
  g1(9,9)=(-1);
  g1(9,10)=1;
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
