function [residual, g1, g2] = ZLB_static(y, x, params)
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

residual = zeros( 16, 1);

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
rhs =exp(y(8))*exp(y(5))/exp(y(11));
residual(5)= lhs-rhs;
lhs =exp(y(11));
rhs =(1-params(5))*exp((-params(6))*y(12))*exp(y(3)*params(6))+exp(y(11))*params(5)*exp(y(3)*params(6));
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
lhs =exp(y(10));
rhs =T142*exp(y(8)*(1+params(4))/(params(1)+params(4)));
residual(15)= lhs-rhs;
lhs =y(2);
rhs =(1-params(14))*log(params(20))+y(2)*params(14)+(1-params(14))*(params(7)*(y(3)-log(params(21)))+params(8)*(y(9)-y(10)))+x(2);
residual(16)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(16, 16);

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
  g1(5,5)=(-(exp(y(8))*exp(y(5))/exp(y(11))));
  g1(5,8)=(-(exp(y(8))*exp(y(5))/exp(y(11))));
  g1(5,9)=exp(y(9));
  g1(5,11)=(-((-(exp(y(8))*exp(y(5))*exp(y(11))))/(exp(y(11))*exp(y(11)))));
  g1(6,3)=(-((1-params(5))*exp((-params(6))*y(12))*params(6)*exp(y(3)*params(6))+exp(y(11))*params(5)*params(6)*exp(y(3)*params(6))));
  g1(6,11)=exp(y(11))-exp(y(11))*params(5)*exp(y(3)*params(6));
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
  g1(15,10)=exp(y(10));
  g1(16,2)=1-params(14);
  g1(16,3)=(-((1-params(14))*params(7)));
  g1(16,9)=(-((1-params(14))*params(8)));
  g1(16,10)=(-((1-params(14))*(-params(8))));
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],16,256);
end
end
