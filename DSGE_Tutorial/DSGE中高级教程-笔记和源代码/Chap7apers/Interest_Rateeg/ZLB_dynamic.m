function [residual, g1, g2, g3] = ZLB_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(16, 1);
T23 = params(2)*exp((-params(1))*y(22))*exp(y(7))/exp(y(23));
T93 = params(6)/(params(6)-1)*exp(y(8))*exp(y(18))/exp(y(19));
T158 = ((params(6)-1)/params(6)/params(3)/(1-params(9)-params(10))^params(1))^(1/(params(1)+params(4)));
lhs =exp((-params(1))*y(6));
rhs =T23;
residual(1)= lhs-rhs;
lhs =params(3)*exp(params(4)*y(10));
rhs =exp((-params(1))*y(6))*exp(y(11));
residual(2)= lhs-rhs;
lhs =exp(y(12));
rhs =exp(y(11))/exp(y(13));
residual(3)= lhs-rhs;
lhs =exp(y(6))+exp(y(20))+exp(y(21));
rhs =exp(y(14));
residual(4)= lhs-rhs;
lhs =exp(y(14));
rhs =exp(y(13))*exp(y(10))/exp(y(16));
residual(5)= lhs-rhs;
lhs =exp(y(16));
rhs =(1-params(5))*exp((-params(6))*y(17))*exp(params(6)*y(8))+params(5)*exp(params(6)*y(8))*exp(y(3));
residual(6)= lhs-rhs;
lhs =exp(y(8)*(1-params(6)));
rhs =params(5)+(1-params(5))*exp(y(17)*(1-params(6)));
residual(7)= lhs-rhs;
lhs =exp(y(17));
rhs =T93;
residual(8)= lhs-rhs;
lhs =exp(y(18));
rhs =exp(y(12))*exp((-params(1))*y(6))*exp(y(14))+params(2)*params(5)*exp(y(23)*params(6))*exp(y(24));
residual(9)= lhs-rhs;
lhs =exp(y(19));
rhs =exp((-params(1))*y(6))*exp(y(14))+params(2)*params(5)*exp(y(23)*(params(6)-1))*exp(y(25));
residual(10)= lhs-rhs;
lhs =y(13);
rhs =params(13)*y(2)+x(it_, 1);
residual(11)= lhs-rhs;
lhs =y(20);
rhs =(1-params(11))*log(params(32))+params(11)*y(4)+x(it_, 3);
residual(12)= lhs-rhs;
lhs =y(21);
rhs =(1-params(12))*log(params(33))+params(12)*y(5)+x(it_, 4);
residual(13)= lhs-rhs;
lhs =exp(y(9));
rhs =exp(y(7))/exp(y(23));
residual(14)= lhs-rhs;
lhs =exp(y(15));
rhs =T158*exp(y(13)*(1+params(4))/(params(1)+params(4)));
residual(15)= lhs-rhs;
lhs =y(7);
rhs =(1-params(14))*log(params(20))+params(14)*y(1)+(1-params(14))*(params(7)*(y(8)-log(params(21)))+params(8)*(y(14)-y(15)))+x(it_, 2);
residual(16)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(16, 29);

  %
  % Jacobian matrix
  %

  g1(1,6)=(-params(1))*exp((-params(1))*y(6));
  g1(1,22)=(-(exp(y(7))*params(2)*(-params(1))*exp((-params(1))*y(22))/exp(y(23))));
  g1(1,7)=(-T23);
  g1(1,23)=(-((-(params(2)*exp((-params(1))*y(22))*exp(y(7))*exp(y(23))))/(exp(y(23))*exp(y(23)))));
  g1(2,6)=(-(exp(y(11))*(-params(1))*exp((-params(1))*y(6))));
  g1(2,10)=params(3)*params(4)*exp(params(4)*y(10));
  g1(2,11)=(-(exp((-params(1))*y(6))*exp(y(11))));
  g1(3,11)=(-(exp(y(11))/exp(y(13))));
  g1(3,12)=exp(y(12));
  g1(3,13)=(-((-(exp(y(11))*exp(y(13))))/(exp(y(13))*exp(y(13)))));
  g1(4,6)=exp(y(6));
  g1(4,14)=(-exp(y(14)));
  g1(4,20)=exp(y(20));
  g1(4,21)=exp(y(21));
  g1(5,10)=(-(exp(y(13))*exp(y(10))/exp(y(16))));
  g1(5,13)=(-(exp(y(13))*exp(y(10))/exp(y(16))));
  g1(5,14)=exp(y(14));
  g1(5,16)=(-((-(exp(y(13))*exp(y(10))*exp(y(16))))/(exp(y(16))*exp(y(16)))));
  g1(6,8)=(-((1-params(5))*exp((-params(6))*y(17))*params(6)*exp(params(6)*y(8))+exp(y(3))*params(5)*params(6)*exp(params(6)*y(8))));
  g1(6,3)=(-(params(5)*exp(params(6)*y(8))*exp(y(3))));
  g1(6,16)=exp(y(16));
  g1(6,17)=(-(exp(params(6)*y(8))*(1-params(5))*(-params(6))*exp((-params(6))*y(17))));
  g1(7,8)=(1-params(6))*exp(y(8)*(1-params(6)));
  g1(7,17)=(-((1-params(5))*(1-params(6))*exp(y(17)*(1-params(6)))));
  g1(8,8)=(-T93);
  g1(8,17)=exp(y(17));
  g1(8,18)=(-T93);
  g1(8,19)=(-((-(params(6)/(params(6)-1)*exp(y(8))*exp(y(18))*exp(y(19))))/(exp(y(19))*exp(y(19)))));
  g1(9,6)=(-(exp(y(12))*exp(y(14))*(-params(1))*exp((-params(1))*y(6))));
  g1(9,23)=(-(exp(y(24))*params(2)*params(5)*params(6)*exp(y(23)*params(6))));
  g1(9,12)=(-(exp(y(12))*exp((-params(1))*y(6))*exp(y(14))));
  g1(9,14)=(-(exp(y(12))*exp((-params(1))*y(6))*exp(y(14))));
  g1(9,18)=exp(y(18));
  g1(9,24)=(-(params(2)*params(5)*exp(y(23)*params(6))*exp(y(24))));
  g1(10,6)=(-(exp(y(14))*(-params(1))*exp((-params(1))*y(6))));
  g1(10,23)=(-(exp(y(25))*params(2)*params(5)*(params(6)-1)*exp(y(23)*(params(6)-1))));
  g1(10,14)=(-(exp((-params(1))*y(6))*exp(y(14))));
  g1(10,19)=exp(y(19));
  g1(10,25)=(-(params(2)*params(5)*exp(y(23)*(params(6)-1))*exp(y(25))));
  g1(11,2)=(-params(13));
  g1(11,13)=1;
  g1(11,26)=(-1);
  g1(12,4)=(-params(11));
  g1(12,20)=1;
  g1(12,28)=(-1);
  g1(13,5)=(-params(12));
  g1(13,21)=1;
  g1(13,29)=(-1);
  g1(14,7)=(-(exp(y(7))/exp(y(23))));
  g1(14,23)=(-((-(exp(y(7))*exp(y(23))))/(exp(y(23))*exp(y(23)))));
  g1(14,9)=exp(y(9));
  g1(15,13)=(-(T158*exp(y(13)*(1+params(4))/(params(1)+params(4)))*(1+params(4))/(params(1)+params(4))));
  g1(15,15)=exp(y(15));
  g1(16,1)=(-params(14));
  g1(16,7)=1;
  g1(16,8)=(-((1-params(14))*params(7)));
  g1(16,14)=(-((1-params(14))*params(8)));
  g1(16,15)=(-((1-params(14))*(-params(8))));
  g1(16,27)=(-1);
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],16,841);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],16,24389);
end
end
