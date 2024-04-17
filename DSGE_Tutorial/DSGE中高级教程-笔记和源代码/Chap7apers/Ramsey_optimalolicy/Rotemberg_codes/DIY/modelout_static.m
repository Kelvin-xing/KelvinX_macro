function [residual, g1, g2] = modelout_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 10, 1);

%
% Model equations
%

T30 = 1+params(5)*(y(3)-1)^2/2;
T44 = y(3)*params(1)*y(4)^2;
T88 = params(5)*y(3)+params(5)*(y(3)-1)+params(5)*(2*y(3)-2)*((params(4)-1)*(1+params(8))-y(2)*params(4)*y(4)*params(2)/y(7))/2;
T98 = y(4)*params(1)*y(3)^2;
T105 = (1+y(1))^2;
T125 = 1+(y(3)-1)^2*params(5)/2;
T163 = params(5)*2*(y(3)-1)/2;
lhs =y(5);
rhs =log(y(4))-params(2)*y(2)^2/2;
residual(1)= lhs-rhs;
lhs =y(6);
rhs =y(5)+y(6)*params(1);
residual(2)= lhs-rhs;
residual(3) = y(10)*T30+1/y(4)-params(3)*y(8)/(y(4)*y(3))+y(8)*y(4)*params(3)/T44+T30*y(2)*params(2)*params(4)*y(9)/y(7);
residual(4) = T30*y(9)*params(4)*y(4)*params(2)/y(7)-params(2)*y(2)-y(10)*y(7);
residual(5) = y(9)*((y(3)-1)*params(5)*params(3)+y(3)*params(5)*params(3))/params(1)-y(9)*T88+(2*y(3)-2)*params(5)*y(4)*y(10)/2+y(8)*y(4)*params(3)/T98;
residual(6) = (-y(8))/T105;
lhs =1/(1+y(1));
rhs =y(4)*params(3)/(y(4)*y(3));
residual(7)= lhs-rhs;
residual(8) = T30*((1+params(8))*(1-params(4))+params(4)*y(4)*params(2)*y(2)/y(7))-y(3)*params(5)*(y(3)-1)+y(3)*(y(3)-1)*params(5)*params(3);
residual(9) = y(4)*T125-y(2)*y(7);
lhs =log(y(7));
rhs =log(y(7))*params(6)+x(1);
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
  g1(1,5)=1;
  g1(2,5)=(-1);
  g1(2,6)=1-params(1);
  g1(3,2)=T30*params(2)*params(4)*y(9)/y(7);
  g1(3,3)=y(10)*T163-(-(y(4)*params(3)*y(8)))/(y(4)*y(3)*y(4)*y(3))+(-(y(8)*y(4)*params(3)*params(1)*y(4)^2))/(T44*T44)+y(2)*params(2)*params(4)*y(9)*T163/y(7);
  g1(3,4)=(-1)/(y(4)*y(4))-(-(y(3)*params(3)*y(8)))/(y(4)*y(3)*y(4)*y(3))+(params(3)*y(8)*T44-y(8)*y(4)*params(3)*y(3)*params(1)*2*y(4))/(T44*T44);
  g1(3,7)=(-(T30*y(2)*params(2)*params(4)*y(9)))/(y(7)*y(7));
  g1(3,8)=(-(params(3)/(y(4)*y(3))))+y(4)*params(3)/T44;
  g1(3,9)=T30*y(2)*params(2)*params(4)/y(7);
  g1(3,10)=T30;
  g1(4,2)=(-params(2));
  g1(4,3)=y(9)*params(4)*y(4)*params(2)*T163/y(7);
  g1(4,4)=T30*params(2)*params(4)*y(9)/y(7);
  g1(4,7)=(-(T30*y(9)*params(4)*y(4)*params(2)))/(y(7)*y(7))-y(10);
  g1(4,9)=T30*params(4)*y(4)*params(2)/y(7);
  g1(4,10)=(-y(7));
  g1(5,2)=(-(y(9)*params(5)*(2*y(3)-2)*(-(params(4)*y(4)*params(2)/y(7)))/2));
  g1(5,3)=y(9)*(params(5)*params(3)+params(5)*params(3))/params(1)-y(9)*(params(5)+params(5)+((params(4)-1)*(1+params(8))-y(2)*params(4)*y(4)*params(2)/y(7))*2*params(5)/2)+2*params(5)*y(4)*y(10)/2+(-(y(8)*y(4)*params(3)*2*y(3)*y(4)*params(1)))/(T98*T98);
  g1(5,4)=(-(y(9)*params(5)*(2*y(3)-2)*(-(y(2)*params(2)*params(4)/y(7)))/2))+(2*y(3)-2)*y(10)*params(5)/2+(params(3)*y(8)*T98-y(8)*y(4)*params(3)*params(1)*y(3)^2)/(T98*T98);
  g1(5,7)=(-(y(9)*params(5)*(2*y(3)-2)*(-((-(y(2)*params(4)*y(4)*params(2)))/(y(7)*y(7))))/2));
  g1(5,8)=y(4)*params(3)/T98;
  g1(5,9)=((y(3)-1)*params(5)*params(3)+y(3)*params(5)*params(3))/params(1)-T88;
  g1(5,10)=(2*y(3)-2)*y(4)*params(5)/2;
  g1(6,1)=(-((-y(8))*2*(1+y(1))))/(T105*T105);
  g1(6,8)=(-1)/T105;
  g1(7,1)=(-1)/((1+y(1))*(1+y(1)));
  g1(7,3)=(-((-(y(4)*y(4)*params(3)))/(y(4)*y(3)*y(4)*y(3))));
  g1(7,4)=(-((params(3)*y(4)*y(3)-y(3)*y(4)*params(3))/(y(4)*y(3)*y(4)*y(3))));
  g1(8,2)=T30*params(4)*y(4)*params(2)/y(7);
  g1(8,3)=(y(3)-1)*params(5)*params(3)+y(3)*params(5)*params(3)+((1+params(8))*(1-params(4))+params(4)*y(4)*params(2)*y(2)/y(7))*T163-(params(5)*y(3)+params(5)*(y(3)-1));
  g1(8,4)=T30*params(4)*params(2)*y(2)/y(7);
  g1(8,7)=T30*params(4)*(-(y(4)*params(2)*y(2)))/(y(7)*y(7));
  g1(9,2)=(-y(7));
  g1(9,3)=y(4)*params(5)/2*2*(y(3)-1);
  g1(9,4)=T125;
  g1(9,7)=(-y(2));
  g1(10,7)=1/y(7)-params(6)*1/y(7);
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
