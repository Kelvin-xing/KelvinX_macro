function [residual, g1, g2, g3] = chap2_m_growth_dynamic(y, x, params, steady_state, it_)
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
lhs =exp(y(5));
rhs =exp(params(4)*y(4)+params(5)*y(8));
residual(1)= lhs-rhs;
lhs =1/exp(y(9));
rhs =params(2)*(exp((-params(4))*(y(13)-y(4)))-y(14));
residual(2)= lhs-rhs;
lhs =exp(y(7)+y(8)*(1-params(1)));
rhs =exp(y(4));
residual(3)= lhs-rhs;
lhs =exp(y(5));
rhs =(1-params(1))*exp(y(7)-y(8)*params(1));
residual(4)= lhs-rhs;
lhs =exp(y(10));
rhs =exp(y(9)-y(14));
residual(5)= lhs-rhs;
lhs =y(9);
rhs =log(1/params(2))+params(6)*y(6)+x(it_, 2);
residual(6)= lhs-rhs;
lhs =exp(y(4));
rhs =exp(y(11));
residual(7)= lhs-rhs;
lhs =y(7);
rhs =params(3)*y(1)+x(it_, 1);
residual(8)= lhs-rhs;
lhs =y(12);
rhs =params(16)+4*(y(6)+y(11)-y(3)-params(7)*(y(9)-y(2)));
residual(9)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(9, 16);

  %
  % Jacobian matrix
  %

  g1(1,4)=(-(params(4)*exp(params(4)*y(4)+params(5)*y(8))));
  g1(1,5)=exp(y(5));
  g1(1,8)=(-(params(5)*exp(params(4)*y(4)+params(5)*y(8))));
  g1(2,4)=(-(params(2)*params(4)*exp((-params(4))*(y(13)-y(4)))));
  g1(2,13)=(-(params(2)*(-params(4))*exp((-params(4))*(y(13)-y(4)))));
  g1(2,14)=params(2);
  g1(2,9)=(-exp(y(9)))/(exp(y(9))*exp(y(9)));
  g1(3,4)=(-exp(y(4)));
  g1(3,7)=exp(y(7)+y(8)*(1-params(1)));
  g1(3,8)=(1-params(1))*exp(y(7)+y(8)*(1-params(1)));
  g1(4,5)=exp(y(5));
  g1(4,7)=(-((1-params(1))*exp(y(7)-y(8)*params(1))));
  g1(4,8)=(-((1-params(1))*exp(y(7)-y(8)*params(1))*(-params(1))));
  g1(5,14)=exp(y(9)-y(14));
  g1(5,9)=(-exp(y(9)-y(14)));
  g1(5,10)=exp(y(10));
  g1(6,6)=(-params(6));
  g1(6,9)=1;
  g1(6,16)=(-1);
  g1(7,4)=exp(y(4));
  g1(7,11)=(-exp(y(11)));
  g1(8,1)=(-params(3));
  g1(8,7)=1;
  g1(8,15)=(-1);
  g1(9,6)=(-4);
  g1(9,2)=(-(4*params(7)));
  g1(9,9)=(-(4*(-params(7))));
  g1(9,3)=4;
  g1(9,11)=(-4);
  g1(9,12)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],9,256);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],9,4096);
end
end
