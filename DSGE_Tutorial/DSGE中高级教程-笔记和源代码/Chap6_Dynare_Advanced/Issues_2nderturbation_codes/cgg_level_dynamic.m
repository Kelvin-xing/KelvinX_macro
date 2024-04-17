function [residual, g1, g2, g3] = cgg_level_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(4, 1);
lhs =exp(y(3))+exp(y(4));
rhs =exp(y(6));
residual(1)= lhs-rhs;
lhs =exp(y(6));
rhs =exp(y(5))*exp(params(2)*y(1))+(1-params(3))*exp(y(1));
residual(2)= lhs-rhs;
lhs =params(4)*exp((-params(1))*y(7))*(1-params(3)+params(2)*exp(y(4)*(params(2)-1))*exp(y(8)));
rhs =exp(y(3)*(-params(1)));
residual(3)= lhs-rhs;
lhs =y(5);
rhs =params(5)*y(2)+x(it_, 1);
residual(4)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(4, 9);

  %
  % Jacobian matrix
  %

  g1(1,3)=exp(y(3));
  g1(1,4)=exp(y(4));
  g1(1,6)=(-exp(y(6)));
  g1(2,1)=(-((1-params(3))*exp(y(1))+exp(y(5))*params(2)*exp(params(2)*y(1))));
  g1(2,5)=(-(exp(y(5))*exp(params(2)*y(1))));
  g1(2,6)=exp(y(6));
  g1(3,3)=(-((-params(1))*exp(y(3)*(-params(1)))));
  g1(3,7)=(1-params(3)+params(2)*exp(y(4)*(params(2)-1))*exp(y(8)))*params(4)*(-params(1))*exp((-params(1))*y(7));
  g1(3,4)=params(4)*exp((-params(1))*y(7))*exp(y(8))*params(2)*(params(2)-1)*exp(y(4)*(params(2)-1));
  g1(3,8)=params(4)*exp((-params(1))*y(7))*params(2)*exp(y(4)*(params(2)-1))*exp(y(8));
  g1(4,2)=(-params(5));
  g1(4,5)=1;
  g1(4,9)=(-1);
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  v2 = zeros(18,3);
  v2(1,1)=1;
  v2(1,2)=21;
  v2(1,3)=exp(y(3));
  v2(2,1)=1;
  v2(2,2)=31;
  v2(2,3)=exp(y(4));
  v2(3,1)=1;
  v2(3,2)=51;
  v2(3,3)=(-exp(y(6)));
  v2(4,1)=2;
  v2(4,2)=1;
  v2(4,3)=(-((1-params(3))*exp(y(1))+exp(y(5))*params(2)*params(2)*exp(params(2)*y(1))));
  v2(5,1)=2;
  v2(5,2)=37;
  v2(5,3)=(-(exp(y(5))*params(2)*exp(params(2)*y(1))));
  v2(6,1)=2;
  v2(6,2)=5;
  v2(6,3)=  v2(5,3);
  v2(7,1)=2;
  v2(7,2)=41;
  v2(7,3)=(-(exp(y(5))*exp(params(2)*y(1))));
  v2(8,1)=2;
  v2(8,2)=51;
  v2(8,3)=exp(y(6));
  v2(9,1)=3;
  v2(9,2)=21;
  v2(9,3)=(-((-params(1))*(-params(1))*exp(y(3)*(-params(1)))));
  v2(10,1)=3;
  v2(10,2)=61;
  v2(10,3)=(1-params(3)+params(2)*exp(y(4)*(params(2)-1))*exp(y(8)))*params(4)*(-params(1))*(-params(1))*exp((-params(1))*y(7));
  v2(11,1)=3;
  v2(11,2)=34;
  v2(11,3)=params(4)*(-params(1))*exp((-params(1))*y(7))*exp(y(8))*params(2)*(params(2)-1)*exp(y(4)*(params(2)-1));
  v2(12,1)=3;
  v2(12,2)=58;
  v2(12,3)=  v2(11,3);
  v2(13,1)=3;
  v2(13,2)=31;
  v2(13,3)=params(4)*exp((-params(1))*y(7))*exp(y(8))*params(2)*(params(2)-1)*(params(2)-1)*exp(y(4)*(params(2)-1));
  v2(14,1)=3;
  v2(14,2)=70;
  v2(14,3)=params(2)*exp(y(4)*(params(2)-1))*exp(y(8))*params(4)*(-params(1))*exp((-params(1))*y(7));
  v2(15,1)=3;
  v2(15,2)=62;
  v2(15,3)=  v2(14,3);
  v2(16,1)=3;
  v2(16,2)=67;
  v2(16,3)=params(4)*exp((-params(1))*y(7))*exp(y(8))*params(2)*(params(2)-1)*exp(y(4)*(params(2)-1));
  v2(17,1)=3;
  v2(17,2)=35;
  v2(17,3)=  v2(16,3);
  v2(18,1)=3;
  v2(18,2)=71;
  v2(18,3)=params(4)*exp((-params(1))*y(7))*params(2)*exp(y(4)*(params(2)-1))*exp(y(8));
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),4,81);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],4,729);
end
end
