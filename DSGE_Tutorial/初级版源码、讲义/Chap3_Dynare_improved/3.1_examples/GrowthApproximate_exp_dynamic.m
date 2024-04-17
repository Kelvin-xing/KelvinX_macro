function [residual, g1, g2, g3] = GrowthApproximate_exp_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(6, 1);
mar_c__ = (exp(y(4))^params(2)*(1-exp(y(7)))^(1-params(2)))^(1-params(5));
mar_c1__ = (exp(y(9))^params(2)*(1-exp(y(10)))^(1-params(2)))^(1-params(5));
T11 = exp(y(4))^params(2);
T16 = (1-exp(y(7)))^(1-params(2));
T23 = exp(y(9))^params(2);
T27 = (1-exp(y(10)))^(1-params(2));
T35 = params(1)*mar_c1__/exp(y(9));
T44 = params(4)*exp(y(8))*exp(y(1))^(params(4)-1);
T46 = exp(y(7))^(1-params(4));
T47 = T44*T46;
T50 = 1+T47-params(3);
T55 = exp(y(8))*(1-params(4))*params(2)/(1-params(2));
T57 = T55*exp(y(1))^params(4);
T59 = exp(y(7))^(-params(4));
T60 = T57*T59;
T73 = exp(y(8))*exp(y(1))^params(4);
T85 = exp(y(4))*getPowerDeriv(exp(y(4)),params(2),1);
T86 = T16*T85;
T87 = getPowerDeriv(T11*T16,1-params(5),1);
T89 = exp(y(4))*T86*T87;
T96 = exp(y(9))*getPowerDeriv(exp(y(9)),params(2),1);
T97 = T27*T96;
T98 = getPowerDeriv(T23*T27,1-params(5),1);
T100 = exp(y(9))*T97*T98;
T105 = params(1)*(T100-exp(y(9))*mar_c1__)/(exp(y(9))*exp(y(9)));
T110 = params(4)*exp(y(8))*exp(y(1))*getPowerDeriv(exp(y(1)),params(4)-1,1);
T111 = T46*T110;
T115 = exp(y(1))*getPowerDeriv(exp(y(1)),params(4),1);
T127 = (-exp(y(7)))*getPowerDeriv(1-exp(y(7)),1-params(2),1);
T128 = T11*T127;
T132 = exp(y(7))*getPowerDeriv(exp(y(7)),1-params(4),1);
T133 = T44*T132;
T137 = exp(y(7))*getPowerDeriv(exp(y(7)),(-params(4)),1);
T147 = (-exp(y(10)))*getPowerDeriv(1-exp(y(10)),1-params(2),1);
T148 = T23*T147;
T151 = params(1)*T98*T148/exp(y(9));
T164 = getPowerDeriv(T11*T16,1-params(5),2);
T165 = T86*T164;
T184 = getPowerDeriv(T23*T27,1-params(5),2);
T185 = T97*T184;
T238 = T132+exp(y(7))*exp(y(7))*getPowerDeriv(exp(y(7)),1-params(4),2);
T278 = T115+exp(y(1))*exp(y(1))*getPowerDeriv(exp(y(1)),params(4),2);
T288 = (-exp(y(7)))*T57*T137;
lhs =mar_c__/exp(y(4));
rhs =T35*T50;
residual(1)= lhs-rhs;
lhs =exp(y(4));
rhs =(1-exp(y(7)))*T60;
residual(2)= lhs-rhs;
lhs =exp(y(5));
rhs =exp(y(6))+exp(y(1))*(1-params(3));
residual(3)= lhs-rhs;
lhs =exp(y(3));
rhs =T46*T73;
residual(4)= lhs-rhs;
lhs =exp(y(3));
rhs =exp(y(4))+exp(y(6));
residual(5)= lhs-rhs;
lhs =y(8);
rhs =params(6)*y(2)+x(it_, 1);
residual(6)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(6, 11);

  %
  % Jacobian matrix
  %

  g1(1,4)=(T89-exp(y(4))*mar_c__)/(exp(y(4))*exp(y(4)));
  g1(1,9)=(-(T50*T105));
  g1(1,1)=(-(T35*T111));
  g1(1,7)=T87*T128/exp(y(4))-T35*T133;
  g1(1,10)=(-(T50*T151));
  g1(1,8)=(-(T35*T47));
  g1(2,4)=exp(y(4));
  g1(2,1)=(-((1-exp(y(7)))*T59*T55*T115));
  g1(2,7)=(-(T60*(-exp(y(7)))+(1-exp(y(7)))*T57*T137));
  g1(2,8)=(-((1-exp(y(7)))*T60));
  g1(3,1)=(-(exp(y(1))*(1-params(3))));
  g1(3,5)=exp(y(5));
  g1(3,6)=(-exp(y(6)));
  g1(4,3)=exp(y(3));
  g1(4,1)=(-(T46*exp(y(8))*T115));
  g1(4,7)=(-(T73*T132));
  g1(4,8)=(-(T46*T73));
  g1(5,3)=exp(y(3));
  g1(5,4)=(-exp(y(4)));
  g1(5,6)=(-exp(y(6)));
  g1(6,2)=(-params(6));
  g1(6,8)=1;
  g1(6,11)=(-1);
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  v2 = zeros(54,3);
  v2(1,1)=1;
  v2(1,2)=37;
  v2(1,3)=(exp(y(4))*exp(y(4))*(T89+exp(y(4))*(T87*T16*(T85+exp(y(4))*exp(y(4))*getPowerDeriv(exp(y(4)),params(2),2))+T86*T165)-(T89+exp(y(4))*mar_c__))-(T89-exp(y(4))*mar_c__)*(exp(y(4))*exp(y(4))+exp(y(4))*exp(y(4))))/(exp(y(4))*exp(y(4))*exp(y(4))*exp(y(4)));
  v2(2,1)=1;
  v2(2,2)=97;
  v2(2,3)=(-(T50*params(1)*(exp(y(9))*exp(y(9))*(T100+exp(y(9))*(T98*T27*(T96+exp(y(9))*exp(y(9))*getPowerDeriv(exp(y(9)),params(2),2))+T97*T185)-(T100+exp(y(9))*mar_c1__))-(T100-exp(y(9))*mar_c1__)*(exp(y(9))*exp(y(9))+exp(y(9))*exp(y(9))))/(exp(y(9))*exp(y(9))*exp(y(9))*exp(y(9)))));
  v2(3,1)=1;
  v2(3,2)=9;
  v2(3,3)=(-(T105*T111));
  v2(4,1)=1;
  v2(4,2)=89;
  v2(4,3)=  v2(3,3);
  v2(5,1)=1;
  v2(5,2)=1;
  v2(5,3)=(-(T35*T46*params(4)*exp(y(8))*(exp(y(1))*getPowerDeriv(exp(y(1)),params(4)-1,1)+exp(y(1))*exp(y(1))*getPowerDeriv(exp(y(1)),params(4)-1,2))));
  v2(6,1)=1;
  v2(6,2)=70;
  v2(6,3)=(exp(y(4))*(T128*T165+T87*T85*T127)-exp(y(4))*T87*T128)/(exp(y(4))*exp(y(4)));
  v2(7,1)=1;
  v2(7,2)=40;
  v2(7,3)=  v2(6,3);
  v2(8,1)=1;
  v2(8,2)=75;
  v2(8,3)=(-(T105*T133));
  v2(9,1)=1;
  v2(9,2)=95;
  v2(9,3)=  v2(8,3);
  v2(10,1)=1;
  v2(10,2)=67;
  v2(10,3)=(-(T35*T110*T132));
  v2(11,1)=1;
  v2(11,2)=7;
  v2(11,3)=  v2(10,3);
  v2(12,1)=1;
  v2(12,2)=73;
  v2(12,3)=(T128*T128*T164+T87*T11*(T127+(-exp(y(7)))*(-exp(y(7)))*getPowerDeriv(1-exp(y(7)),1-params(2),2)))/exp(y(4))-T35*T44*T238;
  v2(13,1)=1;
  v2(13,2)=108;
  v2(13,3)=(-(T50*params(1)*(exp(y(9))*(T148*T185+T98*T96*T147)-exp(y(9))*T98*T148)/(exp(y(9))*exp(y(9)))));
  v2(14,1)=1;
  v2(14,2)=98;
  v2(14,3)=  v2(13,3);
  v2(15,1)=1;
  v2(15,2)=100;
  v2(15,3)=(-(T111*T151));
  v2(16,1)=1;
  v2(16,2)=10;
  v2(16,3)=  v2(15,3);
  v2(17,1)=1;
  v2(17,2)=106;
  v2(17,3)=(-(T133*T151));
  v2(18,1)=1;
  v2(18,2)=76;
  v2(18,3)=  v2(17,3);
  v2(19,1)=1;
  v2(19,2)=109;
  v2(19,3)=(-(T50*params(1)*(T148*T148*T184+T98*T23*(T147+(-exp(y(10)))*(-exp(y(10)))*getPowerDeriv(1-exp(y(10)),1-params(2),2)))/exp(y(9))));
  v2(20,1)=1;
  v2(20,2)=86;
  v2(20,3)=(-(T47*T105));
  v2(21,1)=1;
  v2(21,2)=96;
  v2(21,3)=  v2(20,3);
  v2(22,1)=1;
  v2(22,2)=78;
  v2(22,3)=(-(T35*T111));
  v2(23,1)=1;
  v2(23,2)=8;
  v2(23,3)=  v2(22,3);
  v2(24,1)=1;
  v2(24,2)=84;
  v2(24,3)=(-(T35*T133));
  v2(25,1)=1;
  v2(25,2)=74;
  v2(25,3)=  v2(24,3);
  v2(26,1)=1;
  v2(26,2)=87;
  v2(26,3)=(-(T47*T151));
  v2(27,1)=1;
  v2(27,2)=107;
  v2(27,3)=  v2(26,3);
  v2(28,1)=1;
  v2(28,2)=85;
  v2(28,3)=(-(T35*T47));
  v2(29,1)=2;
  v2(29,2)=37;
  v2(29,3)=exp(y(4));
  v2(30,1)=2;
  v2(30,2)=1;
  v2(30,3)=(-((1-exp(y(7)))*T59*T55*T278));
  v2(31,1)=2;
  v2(31,2)=67;
  v2(31,3)=(-(T59*T55*T115*(-exp(y(7)))+(1-exp(y(7)))*T55*T115*T137));
  v2(32,1)=2;
  v2(32,2)=7;
  v2(32,3)=  v2(31,3);
  v2(33,1)=2;
  v2(33,2)=73;
  v2(33,3)=(-(T60*(-exp(y(7)))+T288+T288+(1-exp(y(7)))*T57*(T137+exp(y(7))*exp(y(7))*getPowerDeriv(exp(y(7)),(-params(4)),2))));
  v2(34,1)=2;
  v2(34,2)=78;
  v2(34,3)=(-((1-exp(y(7)))*T59*T55*T115));
  v2(35,1)=2;
  v2(35,2)=8;
  v2(35,3)=  v2(34,3);
  v2(36,1)=2;
  v2(36,2)=84;
  v2(36,3)=(-(T60*(-exp(y(7)))+(1-exp(y(7)))*T57*T137));
  v2(37,1)=2;
  v2(37,2)=74;
  v2(37,3)=  v2(36,3);
  v2(38,1)=2;
  v2(38,2)=85;
  v2(38,3)=(-((1-exp(y(7)))*T60));
  v2(39,1)=3;
  v2(39,2)=1;
  v2(39,3)=(-(exp(y(1))*(1-params(3))));
  v2(40,1)=3;
  v2(40,2)=49;
  v2(40,3)=exp(y(5));
  v2(41,1)=3;
  v2(41,2)=61;
  v2(41,3)=(-exp(y(6)));
  v2(42,1)=4;
  v2(42,2)=25;
  v2(42,3)=exp(y(3));
  v2(43,1)=4;
  v2(43,2)=1;
  v2(43,3)=(-(T46*exp(y(8))*T278));
  v2(44,1)=4;
  v2(44,2)=67;
  v2(44,3)=(-(exp(y(8))*T115*T132));
  v2(45,1)=4;
  v2(45,2)=7;
  v2(45,3)=  v2(44,3);
  v2(46,1)=4;
  v2(46,2)=73;
  v2(46,3)=(-(T73*T238));
  v2(47,1)=4;
  v2(47,2)=78;
  v2(47,3)=(-(T46*exp(y(8))*T115));
  v2(48,1)=4;
  v2(48,2)=8;
  v2(48,3)=  v2(47,3);
  v2(49,1)=4;
  v2(49,2)=84;
  v2(49,3)=(-(T73*T132));
  v2(50,1)=4;
  v2(50,2)=74;
  v2(50,3)=  v2(49,3);
  v2(51,1)=4;
  v2(51,2)=85;
  v2(51,3)=(-(T46*T73));
  v2(52,1)=5;
  v2(52,2)=25;
  v2(52,3)=exp(y(3));
  v2(53,1)=5;
  v2(53,2)=37;
  v2(53,3)=(-exp(y(4)));
  v2(54,1)=5;
  v2(54,2)=61;
  v2(54,3)=(-exp(y(6)));
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),6,121);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],6,1331);
end
end
