function [residual, g1, g2] = Rotemberg_DIY_periods_static(y, x, params)
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

residual = zeros( 10, 1);

%
% Model equations
%

T27 = (1+y(1))^2;
T37 = 1+params(5)*(y(3)-1)^2/2;
T50 = y(3)*params(1)*y(4)^2;
T94 = params(5)*y(3)+params(5)*(y(3)-1)+params(5)*(2*y(3)-2)*((params(4)-1)*(1+params(8))-y(2)*params(4)*y(4)*params(2)/y(5))/2;
T104 = y(4)*params(1)*y(3)^2;
T125 = 1+(y(3)-1)^2*params(5)/2;
T163 = params(5)*2*(y(3)-1)/2;
lhs =y(6);
rhs =log(y(4))-params(2)*y(2)^2/2;
residual(1)= lhs-rhs;
lhs =y(7);
rhs =y(6)+y(7)*params(1);
residual(2)= lhs-rhs;
residual(3) = (-y(8))/T27;
residual(4) = y(10)*T37+1/y(4)-y(8)*params(3)/(y(4)*y(3))+y(8)*y(4)*params(3)/T50+T37*y(2)*params(2)*params(4)*y(9)/y(5);
residual(5) = T37*y(9)*params(4)*y(4)*params(2)/y(5)-params(2)*y(2)-y(10)*y(5);
residual(6) = y(9)*((y(3)-1)*params(5)*params(3)+y(3)*params(5)*params(3))/params(1)-y(9)*T94+(2*y(3)-2)*params(5)*y(4)*y(10)/2+y(8)*y(4)*params(3)/T104;
lhs =1/(1+y(1));
rhs =y(4)*params(3)/(y(4)*y(3));
residual(7)= lhs-rhs;
residual(8) = T37*((1+params(8))*(1-params(4))+params(4)*y(4)*params(2)*y(2)/y(5))-y(3)*params(5)*(y(3)-1)+y(3)*(y(3)-1)*params(5)*params(3);
residual(9) = y(4)*T125-y(2)*y(5);
lhs =log(y(5));
rhs =log(y(5))*params(6)+x(1);
residual(10)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(10, 10);

  %
  % Jacobian matrix
  %

  g1(1,2)=params(2)*2*y(2)/2;
  g1(1,4)=(-(1/y(4)));
  g1(1,6)=1;
  g1(2,6)=(-1);
  g1(2,7)=1-params(1);
  g1(3,1)=(-((-y(8))*2*(1+y(1))))/(T27*T27);
  g1(3,8)=(-1)/T27;
  g1(4,2)=T37*params(2)*params(4)*y(9)/y(5);
  g1(4,3)=y(10)*T163-(-(y(4)*y(8)*params(3)))/(y(4)*y(3)*y(4)*y(3))+(-(y(8)*y(4)*params(3)*params(1)*y(4)^2))/(T50*T50)+y(2)*params(2)*params(4)*y(9)*T163/y(5);
  g1(4,4)=(-1)/(y(4)*y(4))-(-(y(3)*y(8)*params(3)))/(y(4)*y(3)*y(4)*y(3))+(y(8)*params(3)*T50-y(8)*y(4)*params(3)*y(3)*params(1)*2*y(4))/(T50*T50);
  g1(4,5)=(-(T37*y(2)*params(2)*params(4)*y(9)))/(y(5)*y(5));
  g1(4,8)=(-(params(3)/(y(4)*y(3))))+y(4)*params(3)/T50;
  g1(4,9)=T37*y(2)*params(2)*params(4)/y(5);
  g1(4,10)=T37;
  g1(5,2)=(-params(2));
  g1(5,3)=y(9)*params(4)*y(4)*params(2)*T163/y(5);
  g1(5,4)=T37*params(2)*params(4)*y(9)/y(5);
  g1(5,5)=(-(T37*y(9)*params(4)*y(4)*params(2)))/(y(5)*y(5))-y(10);
  g1(5,9)=T37*params(4)*y(4)*params(2)/y(5);
  g1(5,10)=(-y(5));
  g1(6,2)=(-(y(9)*params(5)*(2*y(3)-2)*(-(params(4)*y(4)*params(2)/y(5)))/2));
  g1(6,3)=y(9)*(params(5)*params(3)+params(5)*params(3))/params(1)-y(9)*(params(5)+params(5)+((params(4)-1)*(1+params(8))-y(2)*params(4)*y(4)*params(2)/y(5))*2*params(5)/2)+2*params(5)*y(4)*y(10)/2+(-(y(8)*y(4)*params(3)*2*y(3)*y(4)*params(1)))/(T104*T104);
  g1(6,4)=(-(y(9)*params(5)*(2*y(3)-2)*(-(y(2)*params(2)*params(4)/y(5)))/2))+(2*y(3)-2)*y(10)*params(5)/2+(y(8)*params(3)*T104-y(8)*y(4)*params(3)*params(1)*y(3)^2)/(T104*T104);
  g1(6,5)=(-(y(9)*params(5)*(2*y(3)-2)*(-((-(y(2)*params(4)*y(4)*params(2)))/(y(5)*y(5))))/2));
  g1(6,8)=y(4)*params(3)/T104;
  g1(6,9)=((y(3)-1)*params(5)*params(3)+y(3)*params(5)*params(3))/params(1)-T94;
  g1(6,10)=(2*y(3)-2)*y(4)*params(5)/2;
  g1(7,1)=(-1)/((1+y(1))*(1+y(1)));
  g1(7,3)=(-((-(y(4)*y(4)*params(3)))/(y(4)*y(3)*y(4)*y(3))));
  g1(7,4)=(-((params(3)*y(4)*y(3)-y(3)*y(4)*params(3))/(y(4)*y(3)*y(4)*y(3))));
  g1(8,2)=T37*params(4)*y(4)*params(2)/y(5);
  g1(8,3)=(y(3)-1)*params(5)*params(3)+y(3)*params(5)*params(3)+((1+params(8))*(1-params(4))+params(4)*y(4)*params(2)*y(2)/y(5))*T163-(params(5)*y(3)+params(5)*(y(3)-1));
  g1(8,4)=T37*params(4)*params(2)*y(2)/y(5);
  g1(8,5)=T37*params(4)*(-(y(4)*params(2)*y(2)))/(y(5)*y(5));
  g1(9,2)=(-y(5));
  g1(9,3)=y(4)*params(5)/2*2*(y(3)-1);
  g1(9,4)=T125;
  g1(9,5)=(-y(2));
  g1(10,5)=1/y(5)-params(6)*1/y(5);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],10,100);
end
end
