function [residual, g1, g2, g3] = chap3_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(15, 1);
T18 = (-1)/params(2);
lhs =y(5);
rhs =params(1)*y(20)+params(9)*y(6);
residual(1)= lhs-rhs;
lhs =y(6);
rhs =T18*(y(11)-y(20)-y(9))+y(21);
residual(2)= lhs-rhs;
lhs =y(11);
rhs =y(5)*params(7)+y(6)*params(8)+y(14);
residual(3)= lhs-rhs;
lhs =y(9);
rhs =params(2)*params(3)*(y(22)-y(15));
residual(4)= lhs-rhs;
lhs =y(10);
rhs =y(11)-y(20);
residual(5)= lhs-rhs;
lhs =y(7);
rhs =params(3)*y(15);
residual(6)= lhs-rhs;
lhs =y(6);
rhs =y(8)-y(7);
residual(7)= lhs-rhs;
lhs =y(14);
rhs =params(4)*y(3)+x(it_, 2);
residual(8)= lhs-rhs;
lhs =y(15);
rhs =params(5)*y(4)+x(it_, 1);
residual(9)= lhs-rhs;
lhs =y(8);
rhs =y(15)+(1-params(10))*y(12);
residual(10)= lhs-rhs;
lhs =y(13);
rhs =4*(y(5)+y(8)-y(1)-params(12)*(y(11)-y(2)));
residual(11)= lhs-rhs;
lhs =y(17);
rhs =y(11)*4;
residual(12)= lhs-rhs;
lhs =y(16);
rhs =y(10)*4;
residual(13)= lhs-rhs;
lhs =y(18);
rhs =y(9)*4;
residual(14)= lhs-rhs;
lhs =y(19);
rhs =y(5)*4;
residual(15)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(15, 24);

  %
  % Jacobian matrix
  %

  g1(1,5)=1;
  g1(1,20)=(-params(1));
  g1(1,6)=(-params(9));
  g1(2,20)=T18;
  g1(2,6)=1;
  g1(2,21)=(-1);
  g1(2,9)=T18;
  g1(2,11)=(-T18);
  g1(3,5)=(-params(7));
  g1(3,6)=(-params(8));
  g1(3,11)=1;
  g1(3,14)=(-1);
  g1(4,9)=1;
  g1(4,15)=params(2)*params(3);
  g1(4,22)=(-(params(2)*params(3)));
  g1(5,20)=1;
  g1(5,10)=1;
  g1(5,11)=(-1);
  g1(6,7)=1;
  g1(6,15)=(-params(3));
  g1(7,6)=1;
  g1(7,7)=1;
  g1(7,8)=(-1);
  g1(8,3)=(-params(4));
  g1(8,14)=1;
  g1(8,24)=(-1);
  g1(9,4)=(-params(5));
  g1(9,15)=1;
  g1(9,23)=(-1);
  g1(10,8)=1;
  g1(10,12)=(-(1-params(10)));
  g1(10,15)=(-1);
  g1(11,5)=(-4);
  g1(11,1)=4;
  g1(11,8)=(-4);
  g1(11,2)=(-(4*params(12)));
  g1(11,11)=(-(4*(-params(12))));
  g1(11,13)=1;
  g1(12,11)=(-4);
  g1(12,17)=1;
  g1(13,10)=(-4);
  g1(13,16)=1;
  g1(14,9)=(-4);
  g1(14,18)=1;
  g1(15,5)=(-4);
  g1(15,19)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],15,576);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],15,13824);
end
end
