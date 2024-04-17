function [residual, g1, g2] = Rotemberg_ramsey_policy_static(y, x, params)
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

residual = zeros( 11, 1);

%
% Model equations
%

T12 = (-1)/((1+y(1))*(1+y(1)));
T33 = 1+params(5)*(y(3)-1)^2/2;
T39 = T33*params(4)*params(2)*y(4)/y(5);
T47 = y(4)*params(5)/2*2*(y(3)-1);
T57 = (1+params(7))*(1-params(4))+params(4)*y(4)*params(2)*y(2)/y(5);
T59 = params(5)*2*(y(3)-1)/2;
T67 = params(15)^(-1);
T75 = (-((-(y(4)*y(4)*params(3)))/(y(3)*y(4)*y(3)*y(4))));
T92 = 1+(y(3)-1)^2*params(5)/2;
T96 = T33*params(4)*params(2)*y(2)/y(5);
T113 = 1/y(5);
T121 = T33*params(4)*(-(y(4)*params(2)*y(2)))/(y(5)*y(5));
T171 = y(8)*params(4)*params(2)*y(4)/y(5)*T59;
residual(1) = y(7)*T12;
residual(2) = y(11)*params(2)*2*y(2)/2+y(9)*(-y(5))+y(8)*T39;
residual(3) = y(9)*T47+y(8)*(T57*T59-(params(5)*(y(3)-1)+params(5)*y(3)))+T67*y(7)*T75+T67*y(8)*((y(3)-1)*params(5)*params(3)+y(3)*params(5)*params(3));
residual(4) = y(11)*(-(1/y(4)))+y(9)*T92+y(8)*T96+y(7)*(-(params(3)/(y(3)*y(4))))+T67*y(7)*(-((-(y(3)*y(4)*params(3)))/(y(3)*y(4)*y(3)*y(4))));
residual(5) = y(10)*T113+y(9)*(-y(2))+y(8)*T121+params(15)*y(10)*(-(T113*params(6)));
residual(6) = 1+y(11);
residual(7) = 1/(1+y(1))-y(4)*params(3)/(y(3)*y(4));
residual(8) = T33*T57-y(3)*params(5)*(y(3)-1)+y(3)*(y(3)-1)*params(5)*params(3);
residual(9) = y(4)*T92-y(2)*y(5);
residual(10) = log(y(5))-(params(6)*log(y(5))+x(1));
residual(11) = y(6)-(log(y(4))-params(2)*y(2)^2/2);
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(11, 11);

  %
  % Jacobian matrix
  %

  g1(1,1)=y(7)*(1+y(1)+1+y(1))/((1+y(1))*(1+y(1))*(1+y(1))*(1+y(1)));
  g1(1,7)=T12;
  g1(2,2)=y(11)*2*params(2)/2;
  g1(2,3)=T171;
  g1(2,4)=y(8)*T33*params(4)*params(2)/y(5);
  g1(2,5)=(-y(9))+y(8)*T33*params(4)*(-(params(2)*y(4)))/(y(5)*y(5));
  g1(2,8)=T39;
  g1(2,9)=(-y(5));
  g1(2,11)=params(2)*2*y(2)/2;
  g1(3,2)=T171;
  g1(3,3)=y(9)*y(4)*2*params(5)/2+y(8)*(T57*2*params(5)/2-(params(5)+params(5)))+T67*y(7)*(-((-((-(y(4)*y(4)*params(3)))*(y(4)*y(3)*y(4)+y(4)*y(3)*y(4))))/(y(3)*y(4)*y(3)*y(4)*y(3)*y(4)*y(3)*y(4))))+T67*y(8)*(params(5)*params(3)+params(5)*params(3));
  g1(3,4)=y(9)*params(5)/2*2*(y(3)-1)+y(8)*T59*params(4)*params(2)*y(2)/y(5)+T67*y(7)*(-((y(3)*y(4)*y(3)*y(4)*(-(y(4)*params(3)+y(4)*params(3)))-(-(y(4)*y(4)*params(3)))*(y(3)*y(3)*y(4)+y(3)*y(3)*y(4)))/(y(3)*y(4)*y(3)*y(4)*y(3)*y(4)*y(3)*y(4))));
  g1(3,5)=y(8)*T59*params(4)*(-(y(4)*params(2)*y(2)))/(y(5)*y(5));
  g1(3,7)=T67*T75;
  g1(3,8)=T57*T59-(params(5)*(y(3)-1)+params(5)*y(3))+T67*((y(3)-1)*params(5)*params(3)+y(3)*params(5)*params(3));
  g1(3,9)=T47;
  g1(4,2)=y(8)*T33*params(4)*params(2)/y(5);
  g1(4,3)=y(9)*params(5)/2*2*(y(3)-1)+y(8)*T59*params(4)*params(2)*y(2)/y(5)+y(7)*(-((-(y(4)*params(3)))/(y(3)*y(4)*y(3)*y(4))))+T67*y(7)*(-((y(3)*y(4)*y(3)*y(4)*(-(y(4)*params(3)))-(-(y(3)*y(4)*params(3)))*(y(4)*y(3)*y(4)+y(4)*y(3)*y(4)))/(y(3)*y(4)*y(3)*y(4)*y(3)*y(4)*y(3)*y(4))));
  g1(4,4)=y(11)*(-((-1)/(y(4)*y(4))))+y(7)*(-((-(y(3)*params(3)))/(y(3)*y(4)*y(3)*y(4))))+T67*y(7)*(-((y(3)*y(4)*y(3)*y(4)*(-(y(3)*params(3)))-(-(y(3)*y(4)*params(3)))*(y(3)*y(3)*y(4)+y(3)*y(3)*y(4)))/(y(3)*y(4)*y(3)*y(4)*y(3)*y(4)*y(3)*y(4))));
  g1(4,5)=y(8)*T33*params(4)*(-(params(2)*y(2)))/(y(5)*y(5));
  g1(4,7)=(-(params(3)/(y(3)*y(4))))+T67*(-((-(y(3)*y(4)*params(3)))/(y(3)*y(4)*y(3)*y(4))));
  g1(4,8)=T96;
  g1(4,9)=T92;
  g1(4,11)=(-(1/y(4)));
  g1(5,2)=(-y(9))+y(8)*T33*params(4)*(-(params(2)*y(4)))/(y(5)*y(5));
  g1(5,3)=y(8)*T59*params(4)*(-(y(4)*params(2)*y(2)))/(y(5)*y(5));
  g1(5,4)=y(8)*T33*params(4)*(-(params(2)*y(2)))/(y(5)*y(5));
  g1(5,5)=y(10)*(-1)/(y(5)*y(5))+y(8)*T33*params(4)*(-((-(y(4)*params(2)*y(2)))*(y(5)+y(5))))/(y(5)*y(5)*y(5)*y(5))+params(15)*y(10)*(-(params(6)*(-1)/(y(5)*y(5))));
  g1(5,8)=T121;
  g1(5,9)=(-y(2));
  g1(5,10)=T113+params(15)*(-(T113*params(6)));
  g1(6,11)=1;
  g1(7,1)=T12;
  g1(7,3)=T75;
  g1(7,4)=(-((params(3)*y(3)*y(4)-y(3)*y(4)*params(3))/(y(3)*y(4)*y(3)*y(4))));
  g1(8,2)=T39;
  g1(8,3)=T57*T59-(params(5)*(y(3)-1)+params(5)*y(3))+(y(3)-1)*params(5)*params(3)+y(3)*params(5)*params(3);
  g1(8,4)=T96;
  g1(8,5)=T121;
  g1(9,2)=(-y(5));
  g1(9,3)=T47;
  g1(9,4)=T92;
  g1(9,5)=(-y(2));
  g1(10,5)=T113-T113*params(6);
  g1(11,2)=params(2)*2*y(2)/2;
  g1(11,4)=(-(1/y(4)));
  g1(11,6)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],11,121);
end
end
