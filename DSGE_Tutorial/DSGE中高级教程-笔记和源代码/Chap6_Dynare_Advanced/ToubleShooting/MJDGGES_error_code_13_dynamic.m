function [residual, g1, g2, g3] = MJDGGES_error_code_13_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(17, 1);
lhs =y(14);
rhs =y(15)+params(3)*(y(17)-y(3));
residual(1)= lhs-rhs;
lhs =y(20);
rhs =y(15)+y(5);
residual(2)= lhs-rhs;
lhs =y(21);
rhs =y(14)+y(6);
residual(3)= lhs-rhs;
lhs =y(18);
rhs =y(17)*(1-params(3));
residual(4)= lhs-rhs;
lhs =y(17)-y(3);
rhs =y(19)-y(4)+y(16)-y(15);
residual(5)= lhs-rhs;
lhs =y(12);
rhs =y(7)+y(17)/params(14);
residual(6)= lhs-rhs;
lhs =y(12);
rhs =y(13)+y(9);
residual(7)= lhs-rhs;
lhs =y(15);
rhs =params(1)*y(28)+y(9)*params(12);
residual(8)= lhs-rhs;
lhs =y(9);
rhs =y(25)-1/params(14)*(y(11)-y(28)-y(10));
residual(9)= lhs-rhs;
lhs =y(13);
rhs =params(16)*y(8)+y(7)*params(17);
residual(10)= lhs-rhs;
lhs =y(10);
rhs =y(8)*params(16)*(-params(14))*(1-params(8))-params(17)*params(6)*(y(24)-y(7));
residual(11)= lhs-rhs;
lhs =y(23);
rhs =y(17)*params(3)*(params(13)/params(2)-1);
residual(12)= lhs-rhs;
lhs =y(22);
rhs =y(9)*(params(14)+params(6));
residual(13)= lhs-rhs;
residual(14) = y(15);
lhs =y(8);
rhs =params(8)*y(2)+x(it_, 1);
residual(15)= lhs-rhs;
lhs =y(7);
rhs =params(9)*y(1)+x(it_, 2);
residual(16)= lhs-rhs;
lhs =y(12);
rhs =y(26)-1/params(2)*(y(11)-y(27));
residual(17)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(17, 30);

  %
  % Jacobian matrix
  %

  g1(1,14)=1;
  g1(1,15)=(-1);
  g1(1,3)=params(3);
  g1(1,17)=(-params(3));
  g1(2,15)=(-1);
  g1(2,5)=(-1);
  g1(2,20)=1;
  g1(3,14)=(-1);
  g1(3,6)=(-1);
  g1(3,21)=1;
  g1(4,17)=(-(1-params(3)));
  g1(4,18)=1;
  g1(5,15)=1;
  g1(5,16)=(-1);
  g1(5,3)=(-1);
  g1(5,17)=1;
  g1(5,4)=1;
  g1(5,19)=(-1);
  g1(6,7)=(-1);
  g1(6,12)=1;
  g1(6,17)=(-(1/params(14)));
  g1(7,9)=(-1);
  g1(7,12)=1;
  g1(7,13)=(-1);
  g1(8,9)=(-params(12));
  g1(8,15)=1;
  g1(8,28)=(-params(1));
  g1(9,9)=1;
  g1(9,25)=(-1);
  g1(9,10)=(-(1/params(14)));
  g1(9,11)=1/params(14);
  g1(9,28)=(-(1/params(14)));
  g1(10,7)=(-params(17));
  g1(10,8)=(-params(16));
  g1(10,13)=1;
  g1(11,7)=(-(params(17)*params(6)));
  g1(11,24)=params(17)*params(6);
  g1(11,8)=(-(params(16)*(-params(14))*(1-params(8))));
  g1(11,10)=1;
  g1(12,17)=(-(params(3)*(params(13)/params(2)-1)));
  g1(12,23)=1;
  g1(13,9)=(-(params(14)+params(6)));
  g1(13,22)=1;
  g1(14,15)=1;
  g1(15,2)=(-params(8));
  g1(15,8)=1;
  g1(15,29)=(-1);
  g1(16,1)=(-params(9));
  g1(16,7)=1;
  g1(16,30)=(-1);
  g1(17,11)=1/params(2);
  g1(17,12)=1;
  g1(17,26)=(-1);
  g1(17,27)=(-(1/params(2)));
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],17,900);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],17,27000);
end
end
