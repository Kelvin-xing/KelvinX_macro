function [residual, g1, g2] = ZLB_eight_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                     columns: variables in declaration order
%                                                     rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 53, 1);

%
% Model equations
%

T20 = exp((-params(1))*y(1))*params(2)*exp(y(2))/exp(y(3));
T86 = exp(y(3))*params(6)/(params(6)-1)*exp(y(13))/exp(y(14));
T142 = ((params(6)-1)/params(6)/params(3)/(1-params(9)-params(10))^params(1))^(1/(params(1)+params(4)));
lhs =exp((-params(1))*y(1));
rhs =T20;
residual(1)= lhs-rhs;
lhs =params(3)*exp(params(4)*y(5));
rhs =exp((-params(1))*y(1))*exp(y(6));
residual(2)= lhs-rhs;
lhs =exp(y(7));
rhs =exp(y(6))/exp(y(8));
residual(3)= lhs-rhs;
lhs =exp(y(1))+exp(y(15))+exp(y(16));
rhs =exp(y(9));
residual(4)= lhs-rhs;
lhs =exp(y(9));
rhs =exp(y(8))*exp(y(5))/exp(y(10));
residual(5)= lhs-rhs;
lhs =exp(y(10));
rhs =(1-params(5))*exp((-params(6))*y(12))*exp(y(3)*params(6))+exp(y(10))*params(5)*exp(y(3)*params(6));
residual(6)= lhs-rhs;
lhs =exp(y(3)*(1-params(6)));
rhs =params(5)+(1-params(5))*exp(y(12)*(1-params(6)));
residual(7)= lhs-rhs;
lhs =exp(y(12));
rhs =T86;
residual(8)= lhs-rhs;
lhs =exp(y(13));
rhs =exp(y(7))*exp((-params(1))*y(1))*exp(y(9))+exp(y(13))*exp(y(3)*params(6))*params(2)*params(5);
residual(9)= lhs-rhs;
lhs =exp(y(14));
rhs =exp((-params(1))*y(1))*exp(y(9))+exp(y(14))*params(2)*params(5)*exp(y(3)*(params(6)-1));
residual(10)= lhs-rhs;
lhs =y(8);
rhs =y(8)*params(13)+x(1);
residual(11)= lhs-rhs;
lhs =y(15);
rhs =(1-params(11))*log(params(32))+y(15)*params(11)+x(3);
residual(12)= lhs-rhs;
lhs =y(16);
rhs =(1-params(12))*log(params(33))+y(16)*params(12)+x(4);
residual(13)= lhs-rhs;
lhs =exp(y(4));
rhs =exp(y(2))/exp(y(3));
residual(14)= lhs-rhs;
lhs =exp(y(11));
rhs =T142*exp(y(8)*(1+params(4))/(params(1)+params(4)));
residual(15)= lhs-rhs;
lhs =y(2);
rhs =y(17);
residual(16)= lhs-rhs;
lhs =y(17);
rhs =y(18);
residual(17)= lhs-rhs;
lhs =y(18);
rhs =y(19);
residual(18)= lhs-rhs;
lhs =y(19);
rhs =y(20);
residual(19)= lhs-rhs;
lhs =y(20);
rhs =y(21);
residual(20)= lhs-rhs;
lhs =y(21);
rhs =y(22);
residual(21)= lhs-rhs;
lhs =y(22);
rhs =y(23);
residual(22)= lhs-rhs;
lhs =y(23);
rhs =y(24);
residual(23)= lhs-rhs;
lhs =y(24);
rhs =(1-params(14))*log(params(20))+y(24)*params(14)+(1-params(14))*(params(7)*(y(31)-log(params(21)))+params(8)*(y(38)-y(45)))+y(53);
residual(24)= lhs-rhs;
lhs =y(25);
rhs =y(3);
residual(25)= lhs-rhs;
lhs =y(26);
rhs =y(25);
residual(26)= lhs-rhs;
lhs =y(27);
rhs =y(26);
residual(27)= lhs-rhs;
lhs =y(28);
rhs =y(27);
residual(28)= lhs-rhs;
lhs =y(29);
rhs =y(28);
residual(29)= lhs-rhs;
lhs =y(30);
rhs =y(29);
residual(30)= lhs-rhs;
lhs =y(31);
rhs =y(30);
residual(31)= lhs-rhs;
lhs =y(32);
rhs =y(9);
residual(32)= lhs-rhs;
lhs =y(33);
rhs =y(32);
residual(33)= lhs-rhs;
lhs =y(34);
rhs =y(33);
residual(34)= lhs-rhs;
lhs =y(35);
rhs =y(34);
residual(35)= lhs-rhs;
lhs =y(36);
rhs =y(35);
residual(36)= lhs-rhs;
lhs =y(37);
rhs =y(36);
residual(37)= lhs-rhs;
lhs =y(38);
rhs =y(37);
residual(38)= lhs-rhs;
lhs =y(39);
rhs =y(11);
residual(39)= lhs-rhs;
lhs =y(40);
rhs =y(39);
residual(40)= lhs-rhs;
lhs =y(41);
rhs =y(40);
residual(41)= lhs-rhs;
lhs =y(42);
rhs =y(41);
residual(42)= lhs-rhs;
lhs =y(43);
rhs =y(42);
residual(43)= lhs-rhs;
lhs =y(44);
rhs =y(43);
residual(44)= lhs-rhs;
lhs =y(45);
rhs =y(44);
residual(45)= lhs-rhs;
lhs =y(46);
rhs =x(2);
residual(46)= lhs-rhs;
lhs =y(47);
rhs =y(46);
residual(47)= lhs-rhs;
lhs =y(48);
rhs =y(47);
residual(48)= lhs-rhs;
lhs =y(49);
rhs =y(48);
residual(49)= lhs-rhs;
lhs =y(50);
rhs =y(49);
residual(50)= lhs-rhs;
lhs =y(51);
rhs =y(50);
residual(51)= lhs-rhs;
lhs =y(52);
rhs =y(51);
residual(52)= lhs-rhs;
lhs =y(53);
rhs =y(52);
residual(53)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(53, 53);

  %
  % Jacobian matrix
  %

  g1(1,1)=(-params(1))*exp((-params(1))*y(1))-exp(y(2))*params(2)*(-params(1))*exp((-params(1))*y(1))/exp(y(3));
  g1(1,2)=(-T20);
  g1(1,3)=(-((-(exp((-params(1))*y(1))*params(2)*exp(y(2))*exp(y(3))))/(exp(y(3))*exp(y(3)))));
  g1(2,1)=(-(exp(y(6))*(-params(1))*exp((-params(1))*y(1))));
  g1(2,5)=params(3)*params(4)*exp(params(4)*y(5));
  g1(2,6)=(-(exp((-params(1))*y(1))*exp(y(6))));
  g1(3,6)=(-(exp(y(6))/exp(y(8))));
  g1(3,7)=exp(y(7));
  g1(3,8)=(-((-(exp(y(6))*exp(y(8))))/(exp(y(8))*exp(y(8)))));
  g1(4,1)=exp(y(1));
  g1(4,9)=(-exp(y(9)));
  g1(4,15)=exp(y(15));
  g1(4,16)=exp(y(16));
  g1(5,5)=(-(exp(y(8))*exp(y(5))/exp(y(10))));
  g1(5,8)=(-(exp(y(8))*exp(y(5))/exp(y(10))));
  g1(5,9)=exp(y(9));
  g1(5,10)=(-((-(exp(y(8))*exp(y(5))*exp(y(10))))/(exp(y(10))*exp(y(10)))));
  g1(6,3)=(-((1-params(5))*exp((-params(6))*y(12))*params(6)*exp(y(3)*params(6))+exp(y(10))*params(5)*params(6)*exp(y(3)*params(6))));
  g1(6,10)=exp(y(10))-exp(y(10))*params(5)*exp(y(3)*params(6));
  g1(6,12)=(-(exp(y(3)*params(6))*(1-params(5))*(-params(6))*exp((-params(6))*y(12))));
  g1(7,3)=(1-params(6))*exp(y(3)*(1-params(6)));
  g1(7,12)=(-((1-params(5))*(1-params(6))*exp(y(12)*(1-params(6)))));
  g1(8,3)=(-T86);
  g1(8,12)=exp(y(12));
  g1(8,13)=(-T86);
  g1(8,14)=(-((-(exp(y(3))*params(6)/(params(6)-1)*exp(y(13))*exp(y(14))))/(exp(y(14))*exp(y(14)))));
  g1(9,1)=(-(exp(y(7))*exp(y(9))*(-params(1))*exp((-params(1))*y(1))));
  g1(9,3)=(-(exp(y(13))*params(2)*params(5)*params(6)*exp(y(3)*params(6))));
  g1(9,7)=(-(exp(y(7))*exp((-params(1))*y(1))*exp(y(9))));
  g1(9,9)=(-(exp(y(7))*exp((-params(1))*y(1))*exp(y(9))));
  g1(9,13)=exp(y(13))-exp(y(13))*exp(y(3)*params(6))*params(2)*params(5);
  g1(10,1)=(-(exp(y(9))*(-params(1))*exp((-params(1))*y(1))));
  g1(10,3)=(-(exp(y(14))*params(2)*params(5)*(params(6)-1)*exp(y(3)*(params(6)-1))));
  g1(10,9)=(-(exp((-params(1))*y(1))*exp(y(9))));
  g1(10,14)=exp(y(14))-exp(y(14))*params(2)*params(5)*exp(y(3)*(params(6)-1));
  g1(11,8)=1-params(13);
  g1(12,15)=1-params(11);
  g1(13,16)=1-params(12);
  g1(14,2)=(-(exp(y(2))/exp(y(3))));
  g1(14,3)=(-((-(exp(y(2))*exp(y(3))))/(exp(y(3))*exp(y(3)))));
  g1(14,4)=exp(y(4));
  g1(15,8)=(-(T142*exp(y(8)*(1+params(4))/(params(1)+params(4)))*(1+params(4))/(params(1)+params(4))));
  g1(15,11)=exp(y(11));
  g1(16,2)=1;
  g1(16,17)=(-1);
  g1(17,17)=1;
  g1(17,18)=(-1);
  g1(18,18)=1;
  g1(18,19)=(-1);
  g1(19,19)=1;
  g1(19,20)=(-1);
  g1(20,20)=1;
  g1(20,21)=(-1);
  g1(21,21)=1;
  g1(21,22)=(-1);
  g1(22,22)=1;
  g1(22,23)=(-1);
  g1(23,23)=1;
  g1(23,24)=(-1);
  g1(24,24)=1-params(14);
  g1(24,31)=(-((1-params(14))*params(7)));
  g1(24,38)=(-((1-params(14))*params(8)));
  g1(24,45)=(-((1-params(14))*(-params(8))));
  g1(24,53)=(-1);
  g1(25,3)=(-1);
  g1(25,25)=1;
  g1(26,25)=(-1);
  g1(26,26)=1;
  g1(27,26)=(-1);
  g1(27,27)=1;
  g1(28,27)=(-1);
  g1(28,28)=1;
  g1(29,28)=(-1);
  g1(29,29)=1;
  g1(30,29)=(-1);
  g1(30,30)=1;
  g1(31,30)=(-1);
  g1(31,31)=1;
  g1(32,9)=(-1);
  g1(32,32)=1;
  g1(33,32)=(-1);
  g1(33,33)=1;
  g1(34,33)=(-1);
  g1(34,34)=1;
  g1(35,34)=(-1);
  g1(35,35)=1;
  g1(36,35)=(-1);
  g1(36,36)=1;
  g1(37,36)=(-1);
  g1(37,37)=1;
  g1(38,37)=(-1);
  g1(38,38)=1;
  g1(39,11)=(-1);
  g1(39,39)=1;
  g1(40,39)=(-1);
  g1(40,40)=1;
  g1(41,40)=(-1);
  g1(41,41)=1;
  g1(42,41)=(-1);
  g1(42,42)=1;
  g1(43,42)=(-1);
  g1(43,43)=1;
  g1(44,43)=(-1);
  g1(44,44)=1;
  g1(45,44)=(-1);
  g1(45,45)=1;
  g1(46,46)=1;
  g1(47,46)=(-1);
  g1(47,47)=1;
  g1(48,47)=(-1);
  g1(48,48)=1;
  g1(49,48)=(-1);
  g1(49,49)=1;
  g1(50,49)=(-1);
  g1(50,50)=1;
  g1(51,50)=(-1);
  g1(51,51)=1;
  g1(52,51)=(-1);
  g1(52,52)=1;
  g1(53,52)=(-1);
  g1(53,53)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],53,2809);
end
end
