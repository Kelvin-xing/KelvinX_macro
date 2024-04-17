function [residual, g1, g2, g3] = Growthestimate_dynamic(y, x, params, steady_state, it_)
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
%                                                           columns: equations in order of declaration
%                                                           rows: variables in order stored in M_.lead_lag_incidence
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              columns: equations in order of declaration
%                                                              rows: variables in order stored in M_.lead_lag_incidence
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              columns: equations in order of declaration
%                                                              rows: variables in order stored in M_.lead_lag_incidence
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(4, 1);
T14 = (1-y(5))^(1-params(5));
T15 = y(3)^params(5)*T14;
T25 = (1-y(8))^(1-params(5));
T26 = y(7)^params(5)*T25;
T29 = params(1)*T26^(1-params(6))/y(7);
T37 = params(3)*exp(y(6))*y(1)^(params(3)-1);
T39 = y(5)^(1-params(3));
T43 = 1+T37*T39-params(2);
T48 = exp(y(6))*(1-params(3))*params(5)/(1-params(5));
T50 = T48*y(1)^params(3);
T52 = y(5)^(-params(3));
T53 = T50*T52;
T57 = exp(y(6))*y(1)^params(3);
T73 = T14*getPowerDeriv(y(3),params(5),1);
T74 = getPowerDeriv(T15,1-params(6),1);
T75 = T73*T74;
T81 = T25*getPowerDeriv(y(7),params(5),1);
T82 = getPowerDeriv(T26,1-params(6),1);
T83 = T81*T82;
T88 = params(1)*(y(7)*T83-T26^(1-params(6)))/(y(7)*y(7));
T93 = T39*params(3)*exp(y(6))*getPowerDeriv(y(1),params(3)-1,1);
T107 = y(3)^params(5)*(-(getPowerDeriv(1-y(5),1-params(5),1)));
T110 = getPowerDeriv(y(5),1-params(3),1);
T111 = T37*T110;
T124 = y(7)^params(5)*(-(getPowerDeriv(1-y(8),1-params(5),1)));
T127 = params(1)*T82*T124/y(7);
lhs =T15^(1-params(6))/y(3);
rhs =T29*T43;
residual(1)= lhs-rhs;
lhs =y(3);
rhs =(1-y(5))*T53;
residual(2)= lhs-rhs;
lhs =y(4);
rhs =T39*T57-y(3)+y(1)*(1-params(2));
residual(3)= lhs-rhs;
lhs =y(6);
rhs =params(4)*y(2)+params(7)*x(it_, 1);
residual(4)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(4, 9);

  %
  % Jacobian matrix
  %

  g1(1,3)=(y(3)*T75-T15^(1-params(6)))/(y(3)*y(3));
  g1(1,7)=(-(T43*T88));
  g1(1,1)=(-(T29*T93));
  g1(1,5)=T74*T107/y(3)-T29*T111;
  g1(1,8)=(-(T43*T127));
  g1(1,6)=(-(T29*T37*T39));
  g1(2,3)=1;
  g1(2,1)=(-((1-y(5))*T52*T48*getPowerDeriv(y(1),params(3),1)));
  g1(2,5)=(-((-T53)+(1-y(5))*T50*getPowerDeriv(y(5),(-params(3)),1)));
  g1(2,6)=(-((1-y(5))*T53));
  g1(3,3)=1;
  g1(3,1)=(-(1-params(2)+T39*exp(y(6))*getPowerDeriv(y(1),params(3),1)));
  g1(3,4)=1;
  g1(3,5)=(-(T57*T110));
  g1(3,6)=(-(T39*T57));
  g1(4,2)=(-params(4));
  g1(4,6)=1;
  g1(4,9)=(-params(7));
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  v2 = zeros(46,3);
  v2(1,1)=1;
  v2(1,2)=21;
  v2(1,3)=(y(3)*y(3)*(T75+y(3)*(T74*T14*getPowerDeriv(y(3),params(5),2)+T73*T73*getPowerDeriv(T15,1-params(6),2))-T75)-(y(3)*T75-T15^(1-params(6)))*(y(3)+y(3)))/(y(3)*y(3)*y(3)*y(3));
  v2(2,1)=1;
  v2(2,2)=61;
  v2(2,3)=(-(T43*params(1)*(y(7)*y(7)*(T83+y(7)*(T82*T25*getPowerDeriv(y(7),params(5),2)+T81*T81*getPowerDeriv(T26,1-params(6),2))-T83)-(y(7)*T83-T26^(1-params(6)))*(y(7)+y(7)))/(y(7)*y(7)*y(7)*y(7))));
  v2(3,1)=1;
  v2(3,2)=7;
  v2(3,3)=(-(T88*T93));
  v2(4,1)=1;
  v2(4,2)=55;
  v2(4,3)=  v2(3,3);
  v2(5,1)=1;
  v2(5,2)=1;
  v2(5,3)=(-(T29*T39*params(3)*exp(y(6))*getPowerDeriv(y(1),params(3)-1,2)));
  v2(6,1)=1;
  v2(6,2)=39;
  v2(6,3)=(y(3)*(T107*T73*getPowerDeriv(T15,1-params(6),2)+T74*getPowerDeriv(y(3),params(5),1)*(-(getPowerDeriv(1-y(5),1-params(5),1))))-T74*T107)/(y(3)*y(3));
  v2(7,1)=1;
  v2(7,2)=23;
  v2(7,3)=  v2(6,3);
  v2(8,1)=1;
  v2(8,2)=43;
  v2(8,3)=(-(T88*T111));
  v2(9,1)=1;
  v2(9,2)=59;
  v2(9,3)=  v2(8,3);
  v2(10,1)=1;
  v2(10,2)=37;
  v2(10,3)=(-(T29*params(3)*exp(y(6))*getPowerDeriv(y(1),params(3)-1,1)*T110));
  v2(11,1)=1;
  v2(11,2)=5;
  v2(11,3)=  v2(10,3);
  v2(12,1)=1;
  v2(12,2)=41;
  v2(12,3)=(T107*T107*getPowerDeriv(T15,1-params(6),2)+T74*y(3)^params(5)*getPowerDeriv(1-y(5),1-params(5),2))/y(3)-T29*T37*getPowerDeriv(y(5),1-params(3),2);
  v2(13,1)=1;
  v2(13,2)=70;
  v2(13,3)=(-(T43*params(1)*(y(7)*(T124*T81*getPowerDeriv(T26,1-params(6),2)+T82*getPowerDeriv(y(7),params(5),1)*(-(getPowerDeriv(1-y(8),1-params(5),1))))-T82*T124)/(y(7)*y(7))));
  v2(14,1)=1;
  v2(14,2)=62;
  v2(14,3)=  v2(13,3);
  v2(15,1)=1;
  v2(15,2)=64;
  v2(15,3)=(-(T93*T127));
  v2(16,1)=1;
  v2(16,2)=8;
  v2(16,3)=  v2(15,3);
  v2(17,1)=1;
  v2(17,2)=68;
  v2(17,3)=(-(T111*T127));
  v2(18,1)=1;
  v2(18,2)=44;
  v2(18,3)=  v2(17,3);
  v2(19,1)=1;
  v2(19,2)=71;
  v2(19,3)=(-(T43*params(1)*(T124*T124*getPowerDeriv(T26,1-params(6),2)+T82*y(7)^params(5)*getPowerDeriv(1-y(8),1-params(5),2))/y(7)));
  v2(20,1)=1;
  v2(20,2)=52;
  v2(20,3)=(-(T37*T39*T88));
  v2(21,1)=1;
  v2(21,2)=60;
  v2(21,3)=  v2(20,3);
  v2(22,1)=1;
  v2(22,2)=46;
  v2(22,3)=(-(T29*T93));
  v2(23,1)=1;
  v2(23,2)=6;
  v2(23,3)=  v2(22,3);
  v2(24,1)=1;
  v2(24,2)=50;
  v2(24,3)=(-(T29*T111));
  v2(25,1)=1;
  v2(25,2)=42;
  v2(25,3)=  v2(24,3);
  v2(26,1)=1;
  v2(26,2)=53;
  v2(26,3)=(-(T37*T39*T127));
  v2(27,1)=1;
  v2(27,2)=69;
  v2(27,3)=  v2(26,3);
  v2(28,1)=1;
  v2(28,2)=51;
  v2(28,3)=(-(T29*T37*T39));
  v2(29,1)=2;
  v2(29,2)=1;
  v2(29,3)=(-((1-y(5))*T52*T48*getPowerDeriv(y(1),params(3),2)));
  v2(30,1)=2;
  v2(30,2)=37;
  v2(30,3)=(-((-(T52*T48*getPowerDeriv(y(1),params(3),1)))+(1-y(5))*T48*getPowerDeriv(y(1),params(3),1)*getPowerDeriv(y(5),(-params(3)),1)));
  v2(31,1)=2;
  v2(31,2)=5;
  v2(31,3)=  v2(30,3);
  v2(32,1)=2;
  v2(32,2)=41;
  v2(32,3)=(-((-(T50*getPowerDeriv(y(5),(-params(3)),1)))+(-(T50*getPowerDeriv(y(5),(-params(3)),1)))+(1-y(5))*T50*getPowerDeriv(y(5),(-params(3)),2)));
  v2(33,1)=2;
  v2(33,2)=46;
  v2(33,3)=(-((1-y(5))*T52*T48*getPowerDeriv(y(1),params(3),1)));
  v2(34,1)=2;
  v2(34,2)=6;
  v2(34,3)=  v2(33,3);
  v2(35,1)=2;
  v2(35,2)=50;
  v2(35,3)=(-((-T53)+(1-y(5))*T50*getPowerDeriv(y(5),(-params(3)),1)));
  v2(36,1)=2;
  v2(36,2)=42;
  v2(36,3)=  v2(35,3);
  v2(37,1)=2;
  v2(37,2)=51;
  v2(37,3)=(-((1-y(5))*T53));
  v2(38,1)=3;
  v2(38,2)=1;
  v2(38,3)=(-(T39*exp(y(6))*getPowerDeriv(y(1),params(3),2)));
  v2(39,1)=3;
  v2(39,2)=37;
  v2(39,3)=(-(exp(y(6))*getPowerDeriv(y(1),params(3),1)*T110));
  v2(40,1)=3;
  v2(40,2)=5;
  v2(40,3)=  v2(39,3);
  v2(41,1)=3;
  v2(41,2)=41;
  v2(41,3)=(-(T57*getPowerDeriv(y(5),1-params(3),2)));
  v2(42,1)=3;
  v2(42,2)=46;
  v2(42,3)=(-(T39*exp(y(6))*getPowerDeriv(y(1),params(3),1)));
  v2(43,1)=3;
  v2(43,2)=6;
  v2(43,3)=  v2(42,3);
  v2(44,1)=3;
  v2(44,2)=50;
  v2(44,3)=(-(T57*T110));
  v2(45,1)=3;
  v2(45,2)=42;
  v2(45,3)=  v2(44,3);
  v2(46,1)=3;
  v2(46,2)=51;
  v2(46,3)=(-(T39*T57));
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),4,81);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],4,729);
end
end
