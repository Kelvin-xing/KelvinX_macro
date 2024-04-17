function [residual, g1, g2] = bgg_rbc_static(y, x, params)
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

residual = zeros( 13, 1);

%
% Model equations
%

z__ = (log(y(5))+y(8)^2/2)/y(8);
zplus__ = (log(y(5))+y(8)^2/2)/y(8);
F__ = normcdf(z__,0,1);
G__ = normcdf(z__-y(8),0,1);
d__ = params(2)*G__*(1+y(6))*y(1);
Fp1__ = normcdf(zplus__,0,1);
Gp1__ = normcdf(zplus__-y(8),0,1);
GAMMA__ = G__+y(5)*(1-F__);
GAMMAp1__ = y(5)*(1-Fp1__)+Gp1__;
dFp1__ = normpdf(zplus__,0,1)/y(5)/y(8);
T40 = normpdf(zplus__,0,1);
T66 = (1+y(6))/(1+y(4));
T137 = (-(1+y(6)))/((1+y(4))*(1+y(4)));
T156 = 1/y(5)/y(8);
T163 = exp((-((z__-y(8))*(z__-y(8))))/2)/2.506628274631;
T173 = exp((-(zplus__*zplus__))/2)/2.506628274631;
T175 = (-(T156*T173));
T182 = exp((-((zplus__-y(8))*(zplus__-y(8))))/2)/2.506628274631;
T183 = T156*T182;
T184 = 1-Fp1__+y(5)*T175+T183;
T217 = exp((-(z__*z__))/2)/2.506628274631;
T235 = 1/(1+y(4));
T261 = (y(8)*2*y(8)/2-(log(y(5))+y(8)^2/2))/(y(8)*y(8));
T311 = 1/params(1)/(y(8)/params(1));
lhs =y(1)^params(4);
rhs =y(3)+y(2)+d__;
residual(1)= lhs-rhs;
lhs =y(1);
rhs =y(2)+y(1)*(1-params(5));
residual(2)= lhs-rhs;
lhs =1;
rhs =params(6)*(1+y(4));
residual(3)= lhs-rhs;
residual(4) = (1-Fp1__)/(1-GAMMAp1__)-T66*(1-Fp1__-y(5)*params(2)*dFp1__)/(1-T66*(GAMMAp1__-params(2)*Gp1__));
lhs =1+y(6);
rhs =1+params(4)*y(1)^(params(4)-1)-params(5);
residual(5)= lhs-rhs;
residual(6) = y(1)*(1-T66*(GAMMA__-params(2)*G__))-y(7);
lhs =y(7);
rhs =y(1)*(1+y(6))*params(3)*(1-GAMMA__);
residual(7)= lhs-rhs;
lhs =y(10);
rhs =y(1)-y(7);
residual(8)= lhs-rhs;
lhs =y(9);
rhs =y(5)*T66*y(1)/y(10);
residual(9)= lhs-rhs;
lhs =log(y(8)/params(1));
rhs =log(y(8)/params(1))*params(7)+x(1);
residual(10)= lhs-rhs;
lhs =y(12);
rhs =y(3)+y(2);
residual(11)= lhs-rhs;
lhs =y(13);
rhs =1-(1+y(4))/(1+y(6));
residual(12)= lhs-rhs;
lhs =y(11);
rhs =F__;
residual(13)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(13, 13);

  %
  % Jacobian matrix
  %

  g1(1,1)=getPowerDeriv(y(1),params(4),1)-params(2)*G__*(1+y(6));
  g1(1,2)=(-1);
  g1(1,3)=(-1);
  g1(1,5)=(-(y(1)*(1+y(6))*params(2)*T156*T163));
  g1(1,6)=(-(params(2)*G__*y(1)));
  g1(1,8)=(-(y(1)*(1+y(6))*params(2)*T163*(T261-1)));
  g1(2,1)=1-(1-params(5));
  g1(2,2)=(-1);
  g1(3,4)=(-params(6));
  g1(4,4)=(-(((1-T66*(GAMMAp1__-params(2)*Gp1__))*(1-Fp1__-y(5)*params(2)*dFp1__)*T137-T66*(1-Fp1__-y(5)*params(2)*dFp1__)*(-((GAMMAp1__-params(2)*Gp1__)*T137)))/((1-T66*(GAMMAp1__-params(2)*Gp1__))*(1-T66*(GAMMAp1__-params(2)*Gp1__)))));
  g1(4,5)=((1-GAMMAp1__)*T175-(1-Fp1__)*(-T184))/((1-GAMMAp1__)*(1-GAMMAp1__))-((1-T66*(GAMMAp1__-params(2)*Gp1__))*T66*(T175-(params(2)*dFp1__+y(5)*params(2)*(y(5)*T40*zplus__*(-T156)-T40)/(y(5)*y(5))/y(8)))-T66*(1-Fp1__-y(5)*params(2)*dFp1__)*(-(T66*(T184-params(2)*T183))))/((1-T66*(GAMMAp1__-params(2)*Gp1__))*(1-T66*(GAMMAp1__-params(2)*Gp1__)));
  g1(4,6)=(-(((1-T66*(GAMMAp1__-params(2)*Gp1__))*(1-Fp1__-y(5)*params(2)*dFp1__)*T235-T66*(1-Fp1__-y(5)*params(2)*dFp1__)*(-((GAMMAp1__-params(2)*Gp1__)*T235)))/((1-T66*(GAMMAp1__-params(2)*Gp1__))*(1-T66*(GAMMAp1__-params(2)*Gp1__)))));
  g1(4,8)=((1-GAMMAp1__)*(-(T173*T261))-(1-Fp1__)*(-(y(5)*(-(T173*T261))+T182*(T261-1))))/((1-GAMMAp1__)*(1-GAMMAp1__))-((1-T66*(GAMMAp1__-params(2)*Gp1__))*T66*((-(T173*T261))-y(5)*params(2)*(y(8)*T40*zplus__*(-T261)/y(5)-T40/y(5))/(y(8)*y(8)))-T66*(1-Fp1__-y(5)*params(2)*dFp1__)*(-(T66*(y(5)*(-(T173*T261))+T182*(T261-1)-params(2)*T182*(T261-1)))))/((1-T66*(GAMMAp1__-params(2)*Gp1__))*(1-T66*(GAMMAp1__-params(2)*Gp1__)));
  g1(5,1)=(-(params(4)*getPowerDeriv(y(1),params(4)-1,1)));
  g1(5,6)=1;
  g1(6,1)=1-T66*(GAMMA__-params(2)*G__);
  g1(6,4)=y(1)*(-((GAMMA__-params(2)*G__)*T137));
  g1(6,5)=y(1)*(-(T66*(T156*T163+1-F__+y(5)*(-(T156*T217))-params(2)*T156*T163)));
  g1(6,6)=y(1)*(-((GAMMA__-params(2)*G__)*T235));
  g1(6,7)=(-1);
  g1(6,8)=y(1)*(-(T66*(T163*(T261-1)+y(5)*(-(T217*T261))-params(2)*T163*(T261-1))));
  g1(7,1)=(-((1+y(6))*params(3)*(1-GAMMA__)));
  g1(7,5)=(-(y(1)*(1+y(6))*params(3)*(-(T156*T163+1-F__+y(5)*(-(T156*T217))))));
  g1(7,6)=(-(y(1)*params(3)*(1-GAMMA__)));
  g1(7,7)=1;
  g1(7,8)=(-(y(1)*(1+y(6))*params(3)*(-(T163*(T261-1)+y(5)*(-(T217*T261))))));
  g1(8,1)=(-1);
  g1(8,7)=1;
  g1(8,10)=1;
  g1(9,1)=(-(y(5)*T66*1/y(10)));
  g1(9,4)=(-(y(5)*y(1)/y(10)*T137));
  g1(9,5)=(-(T66*y(1)/y(10)));
  g1(9,6)=(-(y(5)*y(1)/y(10)*T235));
  g1(9,9)=1;
  g1(9,10)=(-(y(5)*T66*(-y(1))/(y(10)*y(10))));
  g1(10,8)=T311-params(7)*T311;
  g1(11,2)=(-1);
  g1(11,3)=(-1);
  g1(11,12)=1;
  g1(12,4)=1/(1+y(6));
  g1(12,6)=(-(1+y(4)))/((1+y(6))*(1+y(6)));
  g1(12,13)=1;
  g1(13,5)=(-(T156*T217));
  g1(13,8)=(-(T217*T261));
  g1(13,11)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],13,169);
end
end
