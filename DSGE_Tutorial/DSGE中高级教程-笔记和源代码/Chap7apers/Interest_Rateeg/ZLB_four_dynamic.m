function [residual, g1, g2, g3] = ZLB_four_dynamic(y, x, params, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(33, 1);
T23 = params(2)*exp((-params(1))*y(42))*exp(y(10))/exp(y(43));
T93 = params(6)/(params(6)-1)*exp(y(11))*exp(y(21))/exp(y(22));
T158 = ((params(6)-1)/params(6)/params(3)/(1-params(9)-params(10))^params(1))^(1/(params(1)+params(4)));
lhs =exp((-params(1))*y(9));
rhs =T23;
residual(1)= lhs-rhs;
lhs =params(3)*exp(params(4)*y(13));
rhs =exp((-params(1))*y(9))*exp(y(14));
residual(2)= lhs-rhs;
lhs =exp(y(15));
rhs =exp(y(14))/exp(y(16));
residual(3)= lhs-rhs;
lhs =exp(y(9))+exp(y(23))+exp(y(24));
rhs =exp(y(17));
residual(4)= lhs-rhs;
lhs =exp(y(17));
rhs =exp(y(16))*exp(y(13))/exp(y(19));
residual(5)= lhs-rhs;
lhs =exp(y(19));
rhs =(1-params(5))*exp((-params(6))*y(20))*exp(params(6)*y(11))+params(5)*exp(params(6)*y(11))*exp(y(2));
residual(6)= lhs-rhs;
lhs =exp(y(11)*(1-params(6)));
rhs =params(5)+(1-params(5))*exp(y(20)*(1-params(6)));
residual(7)= lhs-rhs;
lhs =exp(y(20));
rhs =T93;
residual(8)= lhs-rhs;
lhs =exp(y(21));
rhs =exp(y(15))*exp((-params(1))*y(9))*exp(y(17))+params(2)*params(5)*exp(y(43)*params(6))*exp(y(46));
residual(9)= lhs-rhs;
lhs =exp(y(22));
rhs =exp((-params(1))*y(9))*exp(y(17))+params(2)*params(5)*exp(y(43)*(params(6)-1))*exp(y(47));
residual(10)= lhs-rhs;
lhs =y(16);
rhs =params(13)*y(1)+x(it_, 1);
residual(11)= lhs-rhs;
lhs =y(23);
rhs =(1-params(11))*log(params(32))+params(11)*y(3)+x(it_, 3);
residual(12)= lhs-rhs;
lhs =y(24);
rhs =(1-params(12))*log(params(33))+params(12)*y(4)+x(it_, 4);
residual(13)= lhs-rhs;
lhs =exp(y(12));
rhs =exp(y(10))/exp(y(43));
residual(14)= lhs-rhs;
lhs =exp(y(18));
rhs =T158*exp(y(16)*(1+params(4))/(params(1)+params(4)));
residual(15)= lhs-rhs;
lhs =y(10);
rhs =y(5);
residual(16)= lhs-rhs;
lhs =y(25);
rhs =y(6);
residual(17)= lhs-rhs;
lhs =y(26);
rhs =y(7);
residual(18)= lhs-rhs;
lhs =y(27);
rhs =y(8);
residual(19)= lhs-rhs;
lhs =y(28);
rhs =(1-params(14))*log(params(20))+y(8)*params(14)+(1-params(14))*(params(7)*(y(50)-log(params(21)))+params(8)*(y(53)-y(56)))+y(60);
residual(20)= lhs-rhs;
lhs =y(29);
rhs =y(43);
residual(21)= lhs-rhs;
lhs =y(30);
rhs =y(48);
residual(22)= lhs-rhs;
lhs =y(31);
rhs =y(49);
residual(23)= lhs-rhs;
lhs =y(32);
rhs =y(44);
residual(24)= lhs-rhs;
lhs =y(33);
rhs =y(51);
residual(25)= lhs-rhs;
lhs =y(34);
rhs =y(52);
residual(26)= lhs-rhs;
lhs =y(35);
rhs =y(45);
residual(27)= lhs-rhs;
lhs =y(36);
rhs =y(54);
residual(28)= lhs-rhs;
lhs =y(37);
rhs =y(55);
residual(29)= lhs-rhs;
lhs =y(38);
rhs =x(it_, 2);
residual(30)= lhs-rhs;
lhs =y(39);
rhs =y(57);
residual(31)= lhs-rhs;
lhs =y(40);
rhs =y(58);
residual(32)= lhs-rhs;
lhs =y(41);
rhs =y(59);
residual(33)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(33, 64);

%
% Jacobian matrix
%

g1(1,9)=(-params(1))*exp((-params(1))*y(9));
g1(1,42)=(-(exp(y(10))*params(2)*(-params(1))*exp((-params(1))*y(42))/exp(y(43))));
g1(1,10)=(-T23);
g1(1,43)=(-((-(params(2)*exp((-params(1))*y(42))*exp(y(10))*exp(y(43))))/(exp(y(43))*exp(y(43)))));
g1(2,9)=(-(exp(y(14))*(-params(1))*exp((-params(1))*y(9))));
g1(2,13)=params(3)*params(4)*exp(params(4)*y(13));
g1(2,14)=(-(exp((-params(1))*y(9))*exp(y(14))));
g1(3,14)=(-(exp(y(14))/exp(y(16))));
g1(3,15)=exp(y(15));
g1(3,16)=(-((-(exp(y(14))*exp(y(16))))/(exp(y(16))*exp(y(16)))));
g1(4,9)=exp(y(9));
g1(4,17)=(-exp(y(17)));
g1(4,23)=exp(y(23));
g1(4,24)=exp(y(24));
g1(5,13)=(-(exp(y(16))*exp(y(13))/exp(y(19))));
g1(5,16)=(-(exp(y(16))*exp(y(13))/exp(y(19))));
g1(5,17)=exp(y(17));
g1(5,19)=(-((-(exp(y(16))*exp(y(13))*exp(y(19))))/(exp(y(19))*exp(y(19)))));
g1(6,11)=(-((1-params(5))*exp((-params(6))*y(20))*params(6)*exp(params(6)*y(11))+exp(y(2))*params(5)*params(6)*exp(params(6)*y(11))));
g1(6,2)=(-(params(5)*exp(params(6)*y(11))*exp(y(2))));
g1(6,19)=exp(y(19));
g1(6,20)=(-(exp(params(6)*y(11))*(1-params(5))*(-params(6))*exp((-params(6))*y(20))));
g1(7,11)=(1-params(6))*exp(y(11)*(1-params(6)));
g1(7,20)=(-((1-params(5))*(1-params(6))*exp(y(20)*(1-params(6)))));
g1(8,11)=(-T93);
g1(8,20)=exp(y(20));
g1(8,21)=(-T93);
g1(8,22)=(-((-(params(6)/(params(6)-1)*exp(y(11))*exp(y(21))*exp(y(22))))/(exp(y(22))*exp(y(22)))));
g1(9,9)=(-(exp(y(15))*exp(y(17))*(-params(1))*exp((-params(1))*y(9))));
g1(9,43)=(-(exp(y(46))*params(2)*params(5)*params(6)*exp(y(43)*params(6))));
g1(9,15)=(-(exp(y(15))*exp((-params(1))*y(9))*exp(y(17))));
g1(9,17)=(-(exp(y(15))*exp((-params(1))*y(9))*exp(y(17))));
g1(9,21)=exp(y(21));
g1(9,46)=(-(params(2)*params(5)*exp(y(43)*params(6))*exp(y(46))));
g1(10,9)=(-(exp(y(17))*(-params(1))*exp((-params(1))*y(9))));
g1(10,43)=(-(exp(y(47))*params(2)*params(5)*(params(6)-1)*exp(y(43)*(params(6)-1))));
g1(10,17)=(-(exp((-params(1))*y(9))*exp(y(17))));
g1(10,22)=exp(y(22));
g1(10,47)=(-(params(2)*params(5)*exp(y(43)*(params(6)-1))*exp(y(47))));
g1(11,1)=(-params(13));
g1(11,16)=1;
g1(11,61)=(-1);
g1(12,3)=(-params(11));
g1(12,23)=1;
g1(12,63)=(-1);
g1(13,4)=(-params(12));
g1(13,24)=1;
g1(13,64)=(-1);
g1(14,10)=(-(exp(y(10))/exp(y(43))));
g1(14,43)=(-((-(exp(y(10))*exp(y(43))))/(exp(y(43))*exp(y(43)))));
g1(14,12)=exp(y(12));
g1(15,16)=(-(T158*exp(y(16)*(1+params(4))/(params(1)+params(4)))*(1+params(4))/(params(1)+params(4))));
g1(15,18)=exp(y(18));
g1(16,10)=1;
g1(16,5)=(-1);
g1(17,25)=1;
g1(17,6)=(-1);
g1(18,26)=1;
g1(18,7)=(-1);
g1(19,27)=1;
g1(19,8)=(-1);
g1(20,8)=(-params(14));
g1(20,28)=1;
g1(20,50)=(-((1-params(14))*params(7)));
g1(20,53)=(-((1-params(14))*params(8)));
g1(20,56)=(-((1-params(14))*(-params(8))));
g1(20,60)=(-1);
g1(21,43)=(-1);
g1(21,29)=1;
g1(22,48)=(-1);
g1(22,30)=1;
g1(23,49)=(-1);
g1(23,31)=1;
g1(24,44)=(-1);
g1(24,32)=1;
g1(25,51)=(-1);
g1(25,33)=1;
g1(26,52)=(-1);
g1(26,34)=1;
g1(27,45)=(-1);
g1(27,35)=1;
g1(28,54)=(-1);
g1(28,36)=1;
g1(29,55)=(-1);
g1(29,37)=1;
g1(30,62)=(-1);
g1(30,38)=1;
g1(31,57)=(-1);
g1(31,39)=1;
g1(32,58)=(-1);
g1(32,40)=1;
g1(33,59)=(-1);
g1(33,41)=1;
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],33,4096);
end
if nargout >= 4,
%
% Third order derivatives
%

  g3 = sparse([],[],[],33,262144);
end
end
