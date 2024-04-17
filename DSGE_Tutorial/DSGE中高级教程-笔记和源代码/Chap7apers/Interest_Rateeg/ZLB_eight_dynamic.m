function [residual, g1, g2, g3] = ZLB_eight_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(53, 1);
T23 = params(2)*exp((-params(1))*y(66))*exp(y(14))/exp(y(67));
T93 = params(6)/(params(6)-1)*exp(y(15))*exp(y(25))/exp(y(26));
T158 = ((params(6)-1)/params(6)/params(3)/(1-params(9)-params(10))^params(1))^(1/(params(1)+params(4)));
lhs =exp((-params(1))*y(13));
rhs =T23;
residual(1)= lhs-rhs;
lhs =params(3)*exp(params(4)*y(17));
rhs =exp((-params(1))*y(13))*exp(y(18));
residual(2)= lhs-rhs;
lhs =exp(y(19));
rhs =exp(y(18))/exp(y(20));
residual(3)= lhs-rhs;
lhs =exp(y(13))+exp(y(27))+exp(y(28));
rhs =exp(y(21));
residual(4)= lhs-rhs;
lhs =exp(y(21));
rhs =exp(y(20))*exp(y(17))/exp(y(22));
residual(5)= lhs-rhs;
lhs =exp(y(22));
rhs =(1-params(5))*exp((-params(6))*y(24))*exp(params(6)*y(15))+params(5)*exp(params(6)*y(15))*exp(y(2));
residual(6)= lhs-rhs;
lhs =exp(y(15)*(1-params(6)));
rhs =params(5)+(1-params(5))*exp(y(24)*(1-params(6)));
residual(7)= lhs-rhs;
lhs =exp(y(24));
rhs =T93;
residual(8)= lhs-rhs;
lhs =exp(y(25));
rhs =exp(y(19))*exp((-params(1))*y(13))*exp(y(21))+params(2)*params(5)*exp(y(67)*params(6))*exp(y(70));
residual(9)= lhs-rhs;
lhs =exp(y(26));
rhs =exp((-params(1))*y(13))*exp(y(21))+params(2)*params(5)*exp(y(67)*(params(6)-1))*exp(y(71));
residual(10)= lhs-rhs;
lhs =y(20);
rhs =params(13)*y(1)+x(it_, 1);
residual(11)= lhs-rhs;
lhs =y(27);
rhs =(1-params(11))*log(params(32))+params(11)*y(3)+x(it_, 3);
residual(12)= lhs-rhs;
lhs =y(28);
rhs =(1-params(12))*log(params(33))+params(12)*y(4)+x(it_, 4);
residual(13)= lhs-rhs;
lhs =exp(y(16));
rhs =exp(y(14))/exp(y(67));
residual(14)= lhs-rhs;
lhs =exp(y(23));
rhs =T158*exp(y(20)*(1+params(4))/(params(1)+params(4)));
residual(15)= lhs-rhs;
lhs =y(14);
rhs =y(5);
residual(16)= lhs-rhs;
lhs =y(29);
rhs =y(6);
residual(17)= lhs-rhs;
lhs =y(30);
rhs =y(7);
residual(18)= lhs-rhs;
lhs =y(31);
rhs =y(8);
residual(19)= lhs-rhs;
lhs =y(32);
rhs =y(9);
residual(20)= lhs-rhs;
lhs =y(33);
rhs =y(10);
residual(21)= lhs-rhs;
lhs =y(34);
rhs =y(11);
residual(22)= lhs-rhs;
lhs =y(35);
rhs =y(12);
residual(23)= lhs-rhs;
lhs =y(36);
rhs =(1-params(14))*log(params(20))+y(12)*params(14)+(1-params(14))*(params(7)*(y(78)-log(params(21)))+params(8)*(y(85)-y(92)))+y(100);
residual(24)= lhs-rhs;
lhs =y(37);
rhs =y(67);
residual(25)= lhs-rhs;
lhs =y(38);
rhs =y(72);
residual(26)= lhs-rhs;
lhs =y(39);
rhs =y(73);
residual(27)= lhs-rhs;
lhs =y(40);
rhs =y(74);
residual(28)= lhs-rhs;
lhs =y(41);
rhs =y(75);
residual(29)= lhs-rhs;
lhs =y(42);
rhs =y(76);
residual(30)= lhs-rhs;
lhs =y(43);
rhs =y(77);
residual(31)= lhs-rhs;
lhs =y(44);
rhs =y(68);
residual(32)= lhs-rhs;
lhs =y(45);
rhs =y(79);
residual(33)= lhs-rhs;
lhs =y(46);
rhs =y(80);
residual(34)= lhs-rhs;
lhs =y(47);
rhs =y(81);
residual(35)= lhs-rhs;
lhs =y(48);
rhs =y(82);
residual(36)= lhs-rhs;
lhs =y(49);
rhs =y(83);
residual(37)= lhs-rhs;
lhs =y(50);
rhs =y(84);
residual(38)= lhs-rhs;
lhs =y(51);
rhs =y(69);
residual(39)= lhs-rhs;
lhs =y(52);
rhs =y(86);
residual(40)= lhs-rhs;
lhs =y(53);
rhs =y(87);
residual(41)= lhs-rhs;
lhs =y(54);
rhs =y(88);
residual(42)= lhs-rhs;
lhs =y(55);
rhs =y(89);
residual(43)= lhs-rhs;
lhs =y(56);
rhs =y(90);
residual(44)= lhs-rhs;
lhs =y(57);
rhs =y(91);
residual(45)= lhs-rhs;
lhs =y(58);
rhs =x(it_, 2);
residual(46)= lhs-rhs;
lhs =y(59);
rhs =y(93);
residual(47)= lhs-rhs;
lhs =y(60);
rhs =y(94);
residual(48)= lhs-rhs;
lhs =y(61);
rhs =y(95);
residual(49)= lhs-rhs;
lhs =y(62);
rhs =y(96);
residual(50)= lhs-rhs;
lhs =y(63);
rhs =y(97);
residual(51)= lhs-rhs;
lhs =y(64);
rhs =y(98);
residual(52)= lhs-rhs;
lhs =y(65);
rhs =y(99);
residual(53)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(53, 104);

  %
  % Jacobian matrix
  %

  g1(1,13)=(-params(1))*exp((-params(1))*y(13));
  g1(1,66)=(-(exp(y(14))*params(2)*(-params(1))*exp((-params(1))*y(66))/exp(y(67))));
  g1(1,14)=(-T23);
  g1(1,67)=(-((-(params(2)*exp((-params(1))*y(66))*exp(y(14))*exp(y(67))))/(exp(y(67))*exp(y(67)))));
  g1(2,13)=(-(exp(y(18))*(-params(1))*exp((-params(1))*y(13))));
  g1(2,17)=params(3)*params(4)*exp(params(4)*y(17));
  g1(2,18)=(-(exp((-params(1))*y(13))*exp(y(18))));
  g1(3,18)=(-(exp(y(18))/exp(y(20))));
  g1(3,19)=exp(y(19));
  g1(3,20)=(-((-(exp(y(18))*exp(y(20))))/(exp(y(20))*exp(y(20)))));
  g1(4,13)=exp(y(13));
  g1(4,21)=(-exp(y(21)));
  g1(4,27)=exp(y(27));
  g1(4,28)=exp(y(28));
  g1(5,17)=(-(exp(y(20))*exp(y(17))/exp(y(22))));
  g1(5,20)=(-(exp(y(20))*exp(y(17))/exp(y(22))));
  g1(5,21)=exp(y(21));
  g1(5,22)=(-((-(exp(y(20))*exp(y(17))*exp(y(22))))/(exp(y(22))*exp(y(22)))));
  g1(6,15)=(-((1-params(5))*exp((-params(6))*y(24))*params(6)*exp(params(6)*y(15))+exp(y(2))*params(5)*params(6)*exp(params(6)*y(15))));
  g1(6,2)=(-(params(5)*exp(params(6)*y(15))*exp(y(2))));
  g1(6,22)=exp(y(22));
  g1(6,24)=(-(exp(params(6)*y(15))*(1-params(5))*(-params(6))*exp((-params(6))*y(24))));
  g1(7,15)=(1-params(6))*exp(y(15)*(1-params(6)));
  g1(7,24)=(-((1-params(5))*(1-params(6))*exp(y(24)*(1-params(6)))));
  g1(8,15)=(-T93);
  g1(8,24)=exp(y(24));
  g1(8,25)=(-T93);
  g1(8,26)=(-((-(params(6)/(params(6)-1)*exp(y(15))*exp(y(25))*exp(y(26))))/(exp(y(26))*exp(y(26)))));
  g1(9,13)=(-(exp(y(19))*exp(y(21))*(-params(1))*exp((-params(1))*y(13))));
  g1(9,67)=(-(exp(y(70))*params(2)*params(5)*params(6)*exp(y(67)*params(6))));
  g1(9,19)=(-(exp(y(19))*exp((-params(1))*y(13))*exp(y(21))));
  g1(9,21)=(-(exp(y(19))*exp((-params(1))*y(13))*exp(y(21))));
  g1(9,25)=exp(y(25));
  g1(9,70)=(-(params(2)*params(5)*exp(y(67)*params(6))*exp(y(70))));
  g1(10,13)=(-(exp(y(21))*(-params(1))*exp((-params(1))*y(13))));
  g1(10,67)=(-(exp(y(71))*params(2)*params(5)*(params(6)-1)*exp(y(67)*(params(6)-1))));
  g1(10,21)=(-(exp((-params(1))*y(13))*exp(y(21))));
  g1(10,26)=exp(y(26));
  g1(10,71)=(-(params(2)*params(5)*exp(y(67)*(params(6)-1))*exp(y(71))));
  g1(11,1)=(-params(13));
  g1(11,20)=1;
  g1(11,101)=(-1);
  g1(12,3)=(-params(11));
  g1(12,27)=1;
  g1(12,103)=(-1);
  g1(13,4)=(-params(12));
  g1(13,28)=1;
  g1(13,104)=(-1);
  g1(14,14)=(-(exp(y(14))/exp(y(67))));
  g1(14,67)=(-((-(exp(y(14))*exp(y(67))))/(exp(y(67))*exp(y(67)))));
  g1(14,16)=exp(y(16));
  g1(15,20)=(-(T158*exp(y(20)*(1+params(4))/(params(1)+params(4)))*(1+params(4))/(params(1)+params(4))));
  g1(15,23)=exp(y(23));
  g1(16,14)=1;
  g1(16,5)=(-1);
  g1(17,29)=1;
  g1(17,6)=(-1);
  g1(18,30)=1;
  g1(18,7)=(-1);
  g1(19,31)=1;
  g1(19,8)=(-1);
  g1(20,32)=1;
  g1(20,9)=(-1);
  g1(21,33)=1;
  g1(21,10)=(-1);
  g1(22,34)=1;
  g1(22,11)=(-1);
  g1(23,35)=1;
  g1(23,12)=(-1);
  g1(24,12)=(-params(14));
  g1(24,36)=1;
  g1(24,78)=(-((1-params(14))*params(7)));
  g1(24,85)=(-((1-params(14))*params(8)));
  g1(24,92)=(-((1-params(14))*(-params(8))));
  g1(24,100)=(-1);
  g1(25,67)=(-1);
  g1(25,37)=1;
  g1(26,72)=(-1);
  g1(26,38)=1;
  g1(27,73)=(-1);
  g1(27,39)=1;
  g1(28,74)=(-1);
  g1(28,40)=1;
  g1(29,75)=(-1);
  g1(29,41)=1;
  g1(30,76)=(-1);
  g1(30,42)=1;
  g1(31,77)=(-1);
  g1(31,43)=1;
  g1(32,68)=(-1);
  g1(32,44)=1;
  g1(33,79)=(-1);
  g1(33,45)=1;
  g1(34,80)=(-1);
  g1(34,46)=1;
  g1(35,81)=(-1);
  g1(35,47)=1;
  g1(36,82)=(-1);
  g1(36,48)=1;
  g1(37,83)=(-1);
  g1(37,49)=1;
  g1(38,84)=(-1);
  g1(38,50)=1;
  g1(39,69)=(-1);
  g1(39,51)=1;
  g1(40,86)=(-1);
  g1(40,52)=1;
  g1(41,87)=(-1);
  g1(41,53)=1;
  g1(42,88)=(-1);
  g1(42,54)=1;
  g1(43,89)=(-1);
  g1(43,55)=1;
  g1(44,90)=(-1);
  g1(44,56)=1;
  g1(45,91)=(-1);
  g1(45,57)=1;
  g1(46,102)=(-1);
  g1(46,58)=1;
  g1(47,93)=(-1);
  g1(47,59)=1;
  g1(48,94)=(-1);
  g1(48,60)=1;
  g1(49,95)=(-1);
  g1(49,61)=1;
  g1(50,96)=(-1);
  g1(50,62)=1;
  g1(51,97)=(-1);
  g1(51,63)=1;
  g1(52,98)=(-1);
  g1(52,64)=1;
  g1(53,99)=(-1);
  g1(53,65)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],53,10816);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],53,1124864);
end
end
