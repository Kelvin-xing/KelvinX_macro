function [residual, g1, g2, g3] = mediumdsge_dynamic(y, x, params, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(26, 1);
invt = exp(y(17))/exp(y(4))-1;
invtp = exp(y(42))/exp(y(17))-1;
invr = exp(y(17))/exp(y(4));
invrp = exp(y(42))/exp(y(17));
rc = params(13)*(exp(y(15))-1)+params(14)/2*(exp(y(15))-1)^2;
rcp = params(13)*(exp(y(40))-1)+params(14)/2*(exp(y(40))-1)^2;
T12 = exp(y(17))/exp(y(4));
T79 = 1-params(5)/2*invt^2;
T96 = params(5)*params(1)*exp(y(41))*exp(y(39))*invtp*invrp^2;
T200 = exp(y(23))*exp(params(2)*y(32))*exp(y(24)*(1-params(2)))/exp(y(25));
T430 = (-(exp(y(17))*exp(y(4))))/(exp(y(4))*exp(y(4)));
T457 = (-(exp(y(17))*exp(y(42))))/(exp(y(17))*exp(y(17)));
lhs =exp(y(10));
rhs =1/(exp(y(11))-params(4)*exp(y(1)))-params(4)*params(1)/(exp(y(37))-exp(y(11))*params(4));
residual(1)= lhs-rhs;
lhs =y(12);
rhs =(params(13)+(exp(y(15))-1)*params(14))/exp(y(13));
residual(2)= lhs-rhs;
lhs =exp(y(10));
rhs =params(1)*exp(y(36))*exp(y(14))/exp(y(46));
residual(3)= lhs-rhs;
lhs =exp(y(10));
rhs =exp(y(13))*exp(y(16))*(T79-params(5)*invt*invr)+T96;
residual(4)= lhs-rhs;
lhs =exp(y(16));
rhs =params(1)*(exp(y(36))*(exp(y(40))*y(38)-rcp/exp(y(39)))+exp(y(41))*(1-params(3)));
residual(5)= lhs-rhs;
lhs =exp(y(18));
rhs =params(9)/(params(9)-1)*exp(y(28)-y(29));
residual(6)= lhs-rhs;
lhs =exp(y(28));
rhs =params(8)*exp(params(9)*(1+params(6))*(y(19)-y(18)))*exp((1+params(6))*y(24))+params(1)*params(11)*exp(y(46)*params(9)*(1+params(6)))*exp(params(9)*(1+params(6))*(y(43)-y(18)))*exp(y(47));
residual(7)= lhs-rhs;
lhs =exp(y(29));
rhs =exp(y(10))*exp(params(9)*(y(19)-y(18)))*exp(y(24))+params(1)*params(11)*exp(y(46)*(params(9)-1))*exp(params(9)*(y(43)-y(18)))*exp(y(48));
residual(8)= lhs-rhs;
lhs =exp(y(19)*(1-params(9)));
rhs =(1-params(11))*exp(y(18)*(1-params(9)))+params(11)*exp((params(9)-1)*y(26))*exp((1-params(9))*y(5));
residual(9)= lhs-rhs;
lhs =exp(y(22));
rhs =T200;
residual(10)= lhs-rhs;
lhs =exp(y(25));
rhs =exp(params(9)*y(26))*((1-params(12))*exp((-params(9))*y(27))+params(12)*exp(y(7)));
residual(11)= lhs-rhs;
lhs =exp(y(26)*(1-params(10)));
rhs =params(12)+(1-params(12))*exp(y(27)*(1-params(10)));
residual(12)= lhs-rhs;
lhs =exp(y(27));
rhs =params(10)/(params(10)-1)*exp(y(26))*exp(y(20)-y(21));
residual(13)= lhs-rhs;
lhs =exp(y(20));
rhs =exp(y(22)+y(10)+y(30))+params(1)*params(12)*exp(y(46)*params(10))*exp(y(44));
residual(14)= lhs-rhs;
lhs =exp(y(21));
rhs =exp(y(10)+y(22))+params(1)*params(12)*exp(y(46)*(params(10)-1))*exp(y(45));
residual(15)= lhs-rhs;
lhs =exp(y(19))/y(12);
rhs =(1-params(2))/params(2)*exp(y(32)-y(24));
residual(16)= lhs-rhs;
lhs =exp(y(19));
rhs =(1-params(2))*exp(y(23)+y(30))*exp(params(2)*(y(32)-y(24)));
residual(17)= lhs-rhs;
lhs =exp(y(22));
rhs =exp(y(17))+exp(y(11))+exp(y(33))+rc*exp(y(8))/exp(y(13));
residual(18)= lhs-rhs;
lhs =exp(y(31));
rhs =exp(y(17))*exp(y(13))*T79+(1-params(3))*exp(y(8));
residual(19)= lhs-rhs;
lhs =exp(y(32));
rhs =exp(y(15))*exp(y(8));
residual(20)= lhs-rhs;
lhs =exp(y(33));
rhs =exp(y(22))*y(35);
residual(21)= lhs-rhs;
lhs =y(35);
rhs =(1-params(15))*params(7)+params(15)*y(9)+x(it_, 1);
residual(22)= lhs-rhs;
lhs =y(23);
rhs =params(17)*y(6)+x(it_, 3);
residual(23)= lhs-rhs;
lhs =y(13);
rhs =params(18)*y(2)+x(it_, 4);
residual(24)= lhs-rhs;
lhs =exp(y(34));
rhs =exp(y(16)-y(10));
residual(25)= lhs-rhs;
lhs =y(14);
rhs =(1-params(16))*log(params(29))+params(16)*y(3)+(1-params(16))*(params(23)*(y(26)-log(params(41)))+params(24)*(y(22)-log(params(37))))+x(it_, 2);
residual(26)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(26, 52);

%
% Jacobian matrix
%

g1(1,10)=exp(y(10));
g1(1,1)=(-(params(4)*exp(y(1))/((exp(y(11))-params(4)*exp(y(1)))*(exp(y(11))-params(4)*exp(y(1))))));
g1(1,11)=(-((-exp(y(11)))/((exp(y(11))-params(4)*exp(y(1)))*(exp(y(11))-params(4)*exp(y(1))))-(-(params(4)*params(1)*(-(exp(y(11))*params(4)))))/((exp(y(37))-exp(y(11))*params(4))*(exp(y(37))-exp(y(11))*params(4)))));
g1(1,37)=(-(params(4)*params(1)*exp(y(37))))/((exp(y(37))-exp(y(11))*params(4))*(exp(y(37))-exp(y(11))*params(4)));
g1(2,12)=1;
g1(2,13)=(-((-((params(13)+(exp(y(15))-1)*params(14))*exp(y(13))))/(exp(y(13))*exp(y(13)))));
g1(2,15)=(-(exp(y(15))*params(14)/exp(y(13))));
g1(3,10)=exp(y(10));
g1(3,36)=(-(params(1)*exp(y(36))*exp(y(14))/exp(y(46))));
g1(3,14)=(-(params(1)*exp(y(36))*exp(y(14))/exp(y(46))));
g1(3,46)=(-((-(params(1)*exp(y(36))*exp(y(14))*exp(y(46))))/(exp(y(46))*exp(y(46)))));
g1(4,10)=exp(y(10));
g1(4,13)=(-(exp(y(13))*exp(y(16))*(T79-params(5)*invt*invr)));
g1(4,39)=(-T96);
g1(4,16)=(-(exp(y(13))*exp(y(16))*(T79-params(5)*invt*invr)));
g1(4,41)=(-T96);
g1(4,4)=(-(exp(y(13))*exp(y(16))*((-(params(5)/2*T430*2*invt))-(invr*params(5)*T430+params(5)*invt*T430))));
g1(4,17)=(-(exp(y(13))*exp(y(16))*((-(params(5)/2*T12*2*invt))-(invr*T12*params(5)+T12*params(5)*invt))+invrp^2*params(5)*params(1)*exp(y(41))*exp(y(39))*T457+params(5)*params(1)*exp(y(41))*exp(y(39))*invtp*T457*2*invrp));
g1(4,42)=(-(invrp^2*exp(y(42))/exp(y(17))*params(5)*params(1)*exp(y(41))*exp(y(39))+params(5)*params(1)*exp(y(41))*exp(y(39))*invtp*exp(y(42))/exp(y(17))*2*invrp));
g1(5,36)=(-(params(1)*exp(y(36))*(exp(y(40))*y(38)-rcp/exp(y(39)))));
g1(5,38)=(-(params(1)*exp(y(40))*exp(y(36))));
g1(5,39)=(-(params(1)*exp(y(36))*(-((-(exp(y(39))*rcp))/(exp(y(39))*exp(y(39)))))));
g1(5,40)=(-(params(1)*exp(y(36))*(exp(y(40))*y(38)-(params(13)*exp(y(40))+params(14)/2*exp(y(40))*2*(exp(y(40))-1))/exp(y(39)))));
g1(5,16)=exp(y(16));
g1(5,41)=(-(params(1)*exp(y(41))*(1-params(3))));
g1(6,18)=exp(y(18));
g1(6,28)=(-(params(9)/(params(9)-1)*exp(y(28)-y(29))));
g1(6,29)=(-(params(9)/(params(9)-1)*(-exp(y(28)-y(29)))));
g1(7,18)=(-(exp((1+params(6))*y(24))*params(8)*exp(params(9)*(1+params(6))*(y(19)-y(18)))*(-(params(9)*(1+params(6))))+exp(y(47))*params(1)*params(11)*exp(y(46)*params(9)*(1+params(6)))*exp(params(9)*(1+params(6))*(y(43)-y(18)))*(-(params(9)*(1+params(6))))));
g1(7,43)=(-(exp(y(47))*params(1)*params(11)*exp(y(46)*params(9)*(1+params(6)))*params(9)*(1+params(6))*exp(params(9)*(1+params(6))*(y(43)-y(18)))));
g1(7,19)=(-(exp((1+params(6))*y(24))*params(8)*params(9)*(1+params(6))*exp(params(9)*(1+params(6))*(y(19)-y(18)))));
g1(7,24)=(-(params(8)*exp(params(9)*(1+params(6))*(y(19)-y(18)))*(1+params(6))*exp((1+params(6))*y(24))));
g1(7,46)=(-(exp(y(47))*exp(params(9)*(1+params(6))*(y(43)-y(18)))*params(1)*params(11)*params(9)*(1+params(6))*exp(y(46)*params(9)*(1+params(6)))));
g1(7,28)=exp(y(28));
g1(7,47)=(-(params(1)*params(11)*exp(y(46)*params(9)*(1+params(6)))*exp(params(9)*(1+params(6))*(y(43)-y(18)))*exp(y(47))));
g1(8,10)=(-(exp(y(10))*exp(params(9)*(y(19)-y(18)))*exp(y(24))));
g1(8,18)=(-(exp(y(24))*exp(y(10))*exp(params(9)*(y(19)-y(18)))*(-params(9))+exp(y(48))*params(1)*params(11)*exp(y(46)*(params(9)-1))*exp(params(9)*(y(43)-y(18)))*(-params(9))));
g1(8,43)=(-(exp(y(48))*params(1)*params(11)*exp(y(46)*(params(9)-1))*params(9)*exp(params(9)*(y(43)-y(18)))));
g1(8,19)=(-(exp(y(24))*exp(y(10))*params(9)*exp(params(9)*(y(19)-y(18)))));
g1(8,24)=(-(exp(y(10))*exp(params(9)*(y(19)-y(18)))*exp(y(24))));
g1(8,46)=(-(exp(y(48))*exp(params(9)*(y(43)-y(18)))*params(1)*params(11)*(params(9)-1)*exp(y(46)*(params(9)-1))));
g1(8,29)=exp(y(29));
g1(8,48)=(-(params(1)*params(11)*exp(y(46)*(params(9)-1))*exp(params(9)*(y(43)-y(18)))*exp(y(48))));
g1(9,18)=(-((1-params(11))*(1-params(9))*exp(y(18)*(1-params(9)))));
g1(9,5)=(-(params(11)*exp((params(9)-1)*y(26))*(1-params(9))*exp((1-params(9))*y(5))));
g1(9,19)=(1-params(9))*exp(y(19)*(1-params(9)));
g1(9,26)=(-(exp((1-params(9))*y(5))*params(11)*(params(9)-1)*exp((params(9)-1)*y(26))));
g1(10,22)=exp(y(22));
g1(10,23)=(-T200);
g1(10,24)=(-(exp(y(23))*exp(params(2)*y(32))*(1-params(2))*exp(y(24)*(1-params(2)))/exp(y(25))));
g1(10,25)=(-((-(exp(y(23))*exp(params(2)*y(32))*exp(y(24)*(1-params(2)))*exp(y(25))))/(exp(y(25))*exp(y(25)))));
g1(10,32)=(-(exp(y(24)*(1-params(2)))*exp(y(23))*params(2)*exp(params(2)*y(32))/exp(y(25))));
g1(11,7)=(-(exp(params(9)*y(26))*params(12)*exp(y(7))));
g1(11,25)=exp(y(25));
g1(11,26)=(-(((1-params(12))*exp((-params(9))*y(27))+params(12)*exp(y(7)))*params(9)*exp(params(9)*y(26))));
g1(11,27)=(-(exp(params(9)*y(26))*(1-params(12))*(-params(9))*exp((-params(9))*y(27))));
g1(12,26)=(1-params(10))*exp(y(26)*(1-params(10)));
g1(12,27)=(-((1-params(12))*(1-params(10))*exp(y(27)*(1-params(10)))));
g1(13,20)=(-(params(10)/(params(10)-1)*exp(y(26))*exp(y(20)-y(21))));
g1(13,21)=(-(params(10)/(params(10)-1)*exp(y(26))*(-exp(y(20)-y(21)))));
g1(13,26)=(-(params(10)/(params(10)-1)*exp(y(26))*exp(y(20)-y(21))));
g1(13,27)=exp(y(27));
g1(14,10)=(-exp(y(22)+y(10)+y(30)));
g1(14,20)=exp(y(20));
g1(14,44)=(-(params(1)*params(12)*exp(y(46)*params(10))*exp(y(44))));
g1(14,22)=(-exp(y(22)+y(10)+y(30)));
g1(14,46)=(-(exp(y(44))*params(1)*params(12)*params(10)*exp(y(46)*params(10))));
g1(14,30)=(-exp(y(22)+y(10)+y(30)));
g1(15,10)=(-exp(y(10)+y(22)));
g1(15,21)=exp(y(21));
g1(15,45)=(-(params(1)*params(12)*exp(y(46)*(params(10)-1))*exp(y(45))));
g1(15,22)=(-exp(y(10)+y(22)));
g1(15,46)=(-(exp(y(45))*params(1)*params(12)*(params(10)-1)*exp(y(46)*(params(10)-1))));
g1(16,12)=(-exp(y(19)))/(y(12)*y(12));
g1(16,19)=exp(y(19))/y(12);
g1(16,24)=(-((1-params(2))/params(2)*(-exp(y(32)-y(24)))));
g1(16,32)=(-((1-params(2))/params(2)*exp(y(32)-y(24))));
g1(17,19)=exp(y(19));
g1(17,23)=(-((1-params(2))*exp(y(23)+y(30))*exp(params(2)*(y(32)-y(24)))));
g1(17,24)=(-((1-params(2))*exp(y(23)+y(30))*exp(params(2)*(y(32)-y(24)))*(-params(2))));
g1(17,30)=(-((1-params(2))*exp(y(23)+y(30))*exp(params(2)*(y(32)-y(24)))));
g1(17,32)=(-((1-params(2))*exp(y(23)+y(30))*params(2)*exp(params(2)*(y(32)-y(24)))));
g1(18,11)=(-exp(y(11)));
g1(18,13)=(-((-(exp(y(13))*rc*exp(y(8))))/(exp(y(13))*exp(y(13)))));
g1(18,15)=(-(exp(y(8))*(params(13)*exp(y(15))+params(14)/2*exp(y(15))*2*(exp(y(15))-1))/exp(y(13))));
g1(18,17)=(-exp(y(17)));
g1(18,22)=exp(y(22));
g1(18,8)=(-(rc*exp(y(8))/exp(y(13))));
g1(18,33)=(-exp(y(33)));
g1(19,13)=(-(exp(y(17))*exp(y(13))*T79));
g1(19,4)=(-(exp(y(17))*exp(y(13))*(-(params(5)/2*T430*2*invt))));
g1(19,17)=(-(exp(y(17))*exp(y(13))*T79+exp(y(17))*exp(y(13))*(-(params(5)/2*T12*2*invt))));
g1(19,8)=(-((1-params(3))*exp(y(8))));
g1(19,31)=exp(y(31));
g1(20,15)=(-(exp(y(15))*exp(y(8))));
g1(20,8)=(-(exp(y(15))*exp(y(8))));
g1(20,32)=exp(y(32));
g1(21,22)=(-(exp(y(22))*y(35)));
g1(21,33)=exp(y(33));
g1(21,35)=(-exp(y(22)));
g1(22,9)=(-params(15));
g1(22,35)=1;
g1(22,49)=(-1);
g1(23,6)=(-params(17));
g1(23,23)=1;
g1(23,51)=(-1);
g1(24,2)=(-params(18));
g1(24,13)=1;
g1(24,52)=(-1);
g1(25,10)=exp(y(16)-y(10));
g1(25,16)=(-exp(y(16)-y(10)));
g1(25,34)=exp(y(34));
g1(26,3)=(-params(16));
g1(26,14)=1;
g1(26,22)=(-((1-params(16))*params(24)));
g1(26,26)=(-((1-params(16))*params(23)));
g1(26,50)=(-1);
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],26,2704);
end
if nargout >= 4,
%
% Third order derivatives
%

  g3 = sparse([],[],[],26,140608);
end
end
