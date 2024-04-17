function [residual, g1, g2, g3] = chap7_book_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(23, 1);
lhs =y(15);
rhs =y(16)+params(3)*(y(18)-y(2));
residual(1)= lhs-rhs;
lhs =y(21);
rhs =y(16)+y(4);
residual(2)= lhs-rhs;
lhs =y(22);
rhs =y(15)+y(5);
residual(3)= lhs-rhs;
lhs =y(19);
rhs =y(18)*(1-params(3));
residual(4)= lhs-rhs;
lhs =y(18)-y(2);
rhs =y(20)-y(3)+y(17)-y(16);
residual(5)= lhs-rhs;
lhs =y(13);
rhs =y(8)+y(18)/params(18);
residual(6)= lhs-rhs;
lhs =y(13);
rhs =y(14)+y(10);
residual(7)= lhs-rhs;
lhs =y(16);
rhs =params(1)*y(33)+y(10)*params(16);
residual(8)= lhs-rhs;
lhs =y(10);
rhs =y(32)-1/params(18)*(y(12)-y(33)-y(11));
residual(9)= lhs-rhs;
lhs =y(14);
rhs =params(20)*y(9)+y(8)*params(21);
residual(10)= lhs-rhs;
lhs =y(11);
rhs =y(9)*params(20)*(-params(18))*(1-params(8))-params(21)*params(6)*(y(31)-y(8));
residual(11)= lhs-rhs;
lhs =y(24);
rhs =y(18)*params(3)*(params(17)/params(2)-1);
residual(12)= lhs-rhs;
lhs =y(23);
rhs =y(10)*(params(18)+params(6));
residual(13)= lhs-rhs;
lhs =y(12);
rhs =y(11)+y(16)*params(10)+y(10)*params(11)+y(30);
residual(14)= lhs-rhs;
lhs =y(9);
rhs =params(8)*y(1)+x(it_, 1)+params(15)*x(it_, 3);
residual(15)= lhs-rhs;
lhs =y(8);
rhs =y(31)-(y(25)-y(34))/params(2);
residual(16)= lhs-rhs;
lhs =y(26);
rhs =y(8)*(params(6)+params(2))-(1+params(6))*y(27);
residual(17)= lhs-rhs;
lhs =y(17);
rhs =params(1)*y(34)+y(26)*params(19);
residual(18)= lhs-rhs;
lhs =y(25);
rhs =y(34)*params(12)+y(27)*params(13);
residual(19)= lhs-rhs;
lhs =y(27);
rhs =x(it_, 3)+params(14)*y(6);
residual(20)= lhs-rhs;
lhs =y(29);
rhs =y(13)*params(2)+params(6)*y(28);
residual(21)= lhs-rhs;
lhs =y(13);
rhs =y(9)+y(28);
residual(22)= lhs-rhs;
lhs =y(30);
rhs =params(23)*y(7)+x(it_, 4);
residual(23)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(23, 38);

  %
  % Jacobian matrix
  %

  g1(1,15)=1;
  g1(1,16)=(-1);
  g1(1,2)=params(3);
  g1(1,18)=(-params(3));
  g1(2,16)=(-1);
  g1(2,4)=(-1);
  g1(2,21)=1;
  g1(3,15)=(-1);
  g1(3,5)=(-1);
  g1(3,22)=1;
  g1(4,18)=(-(1-params(3)));
  g1(4,19)=1;
  g1(5,16)=1;
  g1(5,17)=(-1);
  g1(5,2)=(-1);
  g1(5,18)=1;
  g1(5,3)=1;
  g1(5,20)=(-1);
  g1(6,8)=(-1);
  g1(6,13)=1;
  g1(6,18)=(-(1/params(18)));
  g1(7,10)=(-1);
  g1(7,13)=1;
  g1(7,14)=(-1);
  g1(8,10)=(-params(16));
  g1(8,16)=1;
  g1(8,33)=(-params(1));
  g1(9,10)=1;
  g1(9,32)=(-1);
  g1(9,11)=(-(1/params(18)));
  g1(9,12)=1/params(18);
  g1(9,33)=(-(1/params(18)));
  g1(10,8)=(-params(21));
  g1(10,9)=(-params(20));
  g1(10,14)=1;
  g1(11,8)=(-(params(21)*params(6)));
  g1(11,31)=params(21)*params(6);
  g1(11,9)=(-(params(20)*(-params(18))*(1-params(8))));
  g1(11,11)=1;
  g1(12,18)=(-(params(3)*(params(17)/params(2)-1)));
  g1(12,24)=1;
  g1(13,10)=(-(params(18)+params(6)));
  g1(13,23)=1;
  g1(14,10)=(-params(11));
  g1(14,11)=(-1);
  g1(14,12)=1;
  g1(14,16)=(-params(10));
  g1(14,30)=(-1);
  g1(15,1)=(-params(8));
  g1(15,9)=1;
  g1(15,35)=(-1);
  g1(15,37)=(-params(15));
  g1(16,8)=1;
  g1(16,31)=(-1);
  g1(16,34)=(-1)/params(2);
  g1(16,25)=1/params(2);
  g1(17,8)=(-(params(6)+params(2)));
  g1(17,26)=1;
  g1(17,27)=1+params(6);
  g1(18,17)=1;
  g1(18,34)=(-params(1));
  g1(18,26)=(-params(19));
  g1(19,34)=(-params(12));
  g1(19,25)=1;
  g1(19,27)=(-params(13));
  g1(20,6)=(-params(14));
  g1(20,27)=1;
  g1(20,37)=(-1);
  g1(21,13)=(-params(2));
  g1(21,28)=(-params(6));
  g1(21,29)=1;
  g1(22,9)=(-1);
  g1(22,13)=1;
  g1(22,28)=(-1);
  g1(23,7)=(-params(23));
  g1(23,30)=1;
  g1(23,38)=(-1);
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],23,1444);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],23,54872);
end
end
