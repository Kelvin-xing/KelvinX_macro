function [residual, g1, g2, g3] = Rotemberg_DIY_periods_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(10, 1);
T28 = (1+y(5))^2;
T38 = 1+params(5)*(y(7)-1)^2/2;
T55 = y(7)*params(1)*y(8)^2;
T100 = params(5)*y(7)+params(5)*(y(7)-1)+params(5)*(2*y(7)-2)*((params(4)-1)*(1+params(8))-y(6)*params(4)*y(8)*params(2)/y(9))/2;
T110 = y(8)*params(1)*y(7)^2;
T134 = 1+(y(7)-1)^2*params(5)/2;
T174 = params(5)*2*(y(7)-1)/2;
lhs =y(10);
rhs =log(y(8))-params(2)*y(6)^2/2;
residual(1)= lhs-rhs;
lhs =y(11);
rhs =y(10)+params(1)*y(17);
residual(2)= lhs-rhs;
residual(3) = (-y(12))/T28;
residual(4) = y(14)*T38+1/y(8)-y(12)*params(3)/(y(16)*y(15))+params(3)*y(1)*y(3)/T55+T38*y(6)*params(2)*params(4)*y(13)/y(9);
residual(5) = T38*y(13)*params(4)*y(8)*params(2)/y(9)-params(2)*y(6)-y(14)*y(9);
residual(6) = y(4)*((y(7)-1)*params(5)*params(3)+y(7)*params(5)*params(3))/params(1)-y(13)*T100+(2*y(7)-2)*params(5)*y(8)*y(14)/2+params(3)*y(1)*y(3)/T110;
lhs =1/(1+y(5));
rhs =y(8)*params(3)/(y(16)*y(15));
residual(7)= lhs-rhs;
residual(8) = T38*((1+params(8))*(1-params(4))+params(4)*y(8)*params(2)*y(6)/y(9))-y(7)*params(5)*(y(7)-1)+y(15)*params(5)*params(3)*(y(15)-1);
residual(9) = y(8)*T134-y(6)*y(9);
lhs =log(y(9));
rhs =params(6)*log(y(2))+x(it_, 1);
residual(10)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(10, 18);

  %
  % Jacobian matrix
  %

  g1(1,6)=params(2)*2*y(6)/2;
  g1(1,8)=(-(1/y(8)));
  g1(1,10)=1;
  g1(2,10)=(-1);
  g1(2,11)=1;
  g1(2,17)=(-params(1));
  g1(3,5)=(-((-y(12))*2*(1+y(5))))/(T28*T28);
  g1(3,12)=(-1)/T28;
  g1(4,6)=T38*params(2)*params(4)*y(13)/y(9);
  g1(4,7)=y(14)*T174+(-(params(3)*y(1)*y(3)*params(1)*y(8)^2))/(T55*T55)+y(6)*params(2)*params(4)*y(13)*T174/y(9);
  g1(4,15)=(-((-(y(12)*params(3)*y(16)))/(y(16)*y(15)*y(16)*y(15))));
  g1(4,1)=params(3)*y(3)/T55;
  g1(4,8)=(-1)/(y(8)*y(8))+(-(params(3)*y(1)*y(3)*y(7)*params(1)*2*y(8)))/(T55*T55);
  g1(4,16)=(-((-(y(12)*params(3)*y(15)))/(y(16)*y(15)*y(16)*y(15))));
  g1(4,9)=(-(T38*y(6)*params(2)*params(4)*y(13)))/(y(9)*y(9));
  g1(4,3)=params(3)*y(1)/T55;
  g1(4,12)=(-(params(3)/(y(16)*y(15))));
  g1(4,13)=T38*y(6)*params(2)*params(4)/y(9);
  g1(4,14)=T38;
  g1(5,6)=(-params(2));
  g1(5,7)=y(13)*params(4)*y(8)*params(2)*T174/y(9);
  g1(5,8)=T38*params(2)*params(4)*y(13)/y(9);
  g1(5,9)=(-(T38*y(13)*params(4)*y(8)*params(2)))/(y(9)*y(9))-y(14);
  g1(5,13)=T38*params(4)*y(8)*params(2)/y(9);
  g1(5,14)=(-y(9));
  g1(6,6)=(-(y(13)*params(5)*(2*y(7)-2)*(-(params(4)*y(8)*params(2)/y(9)))/2));
  g1(6,7)=y(4)*(params(5)*params(3)+params(5)*params(3))/params(1)-y(13)*(params(5)+params(5)+((params(4)-1)*(1+params(8))-y(6)*params(4)*y(8)*params(2)/y(9))*2*params(5)/2)+2*params(5)*y(8)*y(14)/2+(-(params(3)*y(1)*y(3)*2*y(7)*y(8)*params(1)))/(T110*T110);
  g1(6,1)=params(3)*y(3)/T110;
  g1(6,8)=(-(y(13)*params(5)*(2*y(7)-2)*(-(y(6)*params(2)*params(4)/y(9)))/2))+(2*y(7)-2)*y(14)*params(5)/2+(-(params(3)*y(1)*y(3)*params(1)*y(7)^2))/(T110*T110);
  g1(6,9)=(-(y(13)*params(5)*(2*y(7)-2)*(-((-(y(6)*params(4)*y(8)*params(2)))/(y(9)*y(9))))/2));
  g1(6,3)=params(3)*y(1)/T110;
  g1(6,4)=((y(7)-1)*params(5)*params(3)+y(7)*params(5)*params(3))/params(1);
  g1(6,13)=(-T100);
  g1(6,14)=(2*y(7)-2)*y(8)*params(5)/2;
  g1(7,5)=(-1)/((1+y(5))*(1+y(5)));
  g1(7,15)=(-((-(y(16)*y(8)*params(3)))/(y(16)*y(15)*y(16)*y(15))));
  g1(7,8)=(-(params(3)/(y(16)*y(15))));
  g1(7,16)=(-((-(y(15)*y(8)*params(3)))/(y(16)*y(15)*y(16)*y(15))));
  g1(8,6)=T38*params(4)*y(8)*params(2)/y(9);
  g1(8,7)=((1+params(8))*(1-params(4))+params(4)*y(8)*params(2)*y(6)/y(9))*T174-(params(5)*y(7)+params(5)*(y(7)-1));
  g1(8,15)=params(5)*params(3)*(y(15)-1)+y(15)*params(5)*params(3);
  g1(8,8)=T38*params(4)*params(2)*y(6)/y(9);
  g1(8,9)=T38*params(4)*(-(y(8)*params(2)*y(6)))/(y(9)*y(9));
  g1(9,6)=(-y(9));
  g1(9,7)=y(8)*params(5)/2*2*(y(7)-1);
  g1(9,8)=T134;
  g1(9,9)=(-y(6));
  g1(10,2)=(-(params(6)*1/y(2)));
  g1(10,9)=1/y(9);
  g1(10,18)=(-1);
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],10,324);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],10,5832);
end
end
