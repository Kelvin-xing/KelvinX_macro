function [residual, g1, g2, g3] = bgg_rbc_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(13, 1);
z__ = (log(y(10))+y(4)^2/2)/y(4);
zplus__ = (log(y(20))+y(13)^2/2)/y(13);
F__ = normcdf(z__,0,1);
G__ = normcdf(z__-y(4),0,1);
d__ = params(2)*G__*(1+y(11))*y(1);
Fp1__ = normcdf(zplus__,0,1);
Gp1__ = normcdf(zplus__-y(13),0,1);
GAMMA__ = G__+y(10)*(1-F__);
GAMMAp1__ = y(20)*(1-Fp1__)+Gp1__;
dFp1__ = normpdf(zplus__,0,1)/y(20)/y(13);
T47 = normpdf(zplus__,0,1);
T70 = params(6)*(y(19)/y(8))^(-params(8));
T80 = (1+y(21))/(1+y(9));
T85 = T80*(1-Fp1__-y(20)*params(2)*dFp1__);
T89 = 1-T80*(GAMMAp1__-params(2)*Gp1__);
T101 = (1+y(11))/(1+y(2));
T158 = getPowerDeriv(y(19)/y(8),(-params(8)),1);
T193 = 1/y(10)/y(4);
T200 = exp((-((z__-y(4))*(z__-y(4))))/2)/2.506628274631;
T201 = T193*T200;
T210 = exp((-(z__*z__))/2)/2.506628274631;
T227 = 1/y(20)/y(13);
T232 = exp((-(zplus__*zplus__))/2)/2.506628274631;
T234 = (-(T227*T232));
T241 = exp((-((zplus__-y(13))*(zplus__-y(13))))/2)/2.506628274631;
T242 = T227*T241;
T243 = 1-Fp1__+y(20)*T234+T242;
T301 = (y(4)*2*y(4)/2-(log(y(10))+y(4)^2/2))/(y(4)*y(4));
T303 = T200*(T301-1);
T330 = (y(13)*2*y(13)/2-(log(y(20))+y(13)^2/2))/(y(13)*y(13));
T332 = (-(T232*T330));
lhs =y(1)^params(4);
rhs =y(8)+y(7)+d__;
residual(1)= lhs-rhs;
lhs =y(6);
rhs =y(7)+y(1)*(1-params(5));
residual(2)= lhs-rhs;
lhs =1;
rhs =T70*(1+y(9));
residual(3)= lhs-rhs;
residual(4) = (1-Fp1__)/(1-GAMMAp1__)-T85/T89;
lhs =1+y(11);
rhs =1+params(4)*y(1)^(params(4)-1)-params(5);
residual(5)= lhs-rhs;
residual(6) = y(1)*(1-T101*(GAMMA__-params(2)*G__))-y(3);
lhs =y(12);
rhs =y(1)*(1+y(11))*params(3)*(1-GAMMA__);
residual(7)= lhs-rhs;
lhs =y(15);
rhs =y(6)-y(12);
residual(8)= lhs-rhs;
lhs =y(14);
rhs =y(10)*T101*y(1)/y(5);
residual(9)= lhs-rhs;
lhs =log(y(13)/params(1));
rhs =params(7)*log(y(4)/params(1))+x(it_, 1);
residual(10)= lhs-rhs;
lhs =y(17);
rhs =y(8)+y(7);
residual(11)= lhs-rhs;
lhs =y(18);
rhs =1-(1+y(9))/(1+y(21));
residual(12)= lhs-rhs;
lhs =y(16);
rhs =F__;
residual(13)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(13, 22);

  %
  % Jacobian matrix
  %

  g1(1,1)=getPowerDeriv(y(1),params(4),1)-params(2)*G__*(1+y(11));
  g1(1,7)=(-1);
  g1(1,8)=(-1);
  g1(1,10)=(-(y(1)*(1+y(11))*params(2)*T201));
  g1(1,11)=(-(params(2)*G__*y(1)));
  g1(1,4)=(-(y(1)*(1+y(11))*params(2)*T303));
  g1(2,1)=(-(1-params(5)));
  g1(2,6)=1;
  g1(2,7)=(-1);
  g1(3,8)=(-((1+y(9))*params(6)*(-y(19))/(y(8)*y(8))*T158));
  g1(3,19)=(-((1+y(9))*params(6)*T158*1/y(8)));
  g1(3,9)=(-T70);
  g1(4,9)=(-((T89*(1-Fp1__-y(20)*params(2)*dFp1__)*(-(1+y(21)))/((1+y(9))*(1+y(9)))-T85*(-((GAMMAp1__-params(2)*Gp1__)*(-(1+y(21)))/((1+y(9))*(1+y(9))))))/(T89*T89)));
  g1(4,20)=((1-GAMMAp1__)*T234-(1-Fp1__)*(-T243))/((1-GAMMAp1__)*(1-GAMMAp1__))-(T89*T80*(T234-(params(2)*dFp1__+y(20)*params(2)*(y(20)*T47*zplus__*(-T227)-T47)/(y(20)*y(20))/y(13)))-T85*(-(T80*(T243-params(2)*T242))))/(T89*T89);
  g1(4,21)=(-((T89*(1-Fp1__-y(20)*params(2)*dFp1__)*1/(1+y(9))-T85*(-((GAMMAp1__-params(2)*Gp1__)*1/(1+y(9)))))/(T89*T89)));
  g1(4,13)=((1-GAMMAp1__)*T332-(1-Fp1__)*(-(y(20)*T332+T241*(T330-1))))/((1-GAMMAp1__)*(1-GAMMAp1__))-(T89*T80*(T332-y(20)*params(2)*(y(13)*T47*zplus__*(-T330)/y(20)-T47/y(20))/(y(13)*y(13)))-T85*(-(T80*(y(20)*T332+T241*(T330-1)-params(2)*T241*(T330-1)))))/(T89*T89);
  g1(5,1)=(-(params(4)*getPowerDeriv(y(1),params(4)-1,1)));
  g1(5,11)=1;
  g1(6,1)=1-T101*(GAMMA__-params(2)*G__);
  g1(6,2)=y(1)*(-((GAMMA__-params(2)*G__)*(-(1+y(11)))/((1+y(2))*(1+y(2)))));
  g1(6,10)=y(1)*(-(T101*(T201+1-F__+y(10)*(-(T193*T210))-params(2)*T201)));
  g1(6,11)=y(1)*(-((GAMMA__-params(2)*G__)*1/(1+y(2))));
  g1(6,3)=(-1);
  g1(6,4)=y(1)*(-(T101*(T303+y(10)*(-(T210*T301))-params(2)*T303)));
  g1(7,1)=(-((1+y(11))*params(3)*(1-GAMMA__)));
  g1(7,10)=(-(y(1)*(1+y(11))*params(3)*(-(T201+1-F__+y(10)*(-(T193*T210))))));
  g1(7,11)=(-(y(1)*params(3)*(1-GAMMA__)));
  g1(7,12)=1;
  g1(7,4)=(-(y(1)*(1+y(11))*params(3)*(-(T303+y(10)*(-(T210*T301))))));
  g1(8,6)=(-1);
  g1(8,12)=1;
  g1(8,15)=1;
  g1(9,1)=(-(y(10)*T101*1/y(5)));
  g1(9,2)=(-(y(10)*y(1)/y(5)*(-(1+y(11)))/((1+y(2))*(1+y(2)))));
  g1(9,10)=(-(T101*y(1)/y(5)));
  g1(9,11)=(-(y(10)*y(1)/y(5)*1/(1+y(2))));
  g1(9,14)=1;
  g1(9,5)=(-(y(10)*T101*(-y(1))/(y(5)*y(5))));
  g1(10,4)=(-(params(7)*1/params(1)/(y(4)/params(1))));
  g1(10,13)=1/params(1)/(y(13)/params(1));
  g1(10,22)=(-1);
  g1(11,7)=(-1);
  g1(11,8)=(-1);
  g1(11,17)=1;
  g1(12,9)=1/(1+y(21));
  g1(12,21)=(-(1+y(9)))/((1+y(21))*(1+y(21)));
  g1(12,18)=1;
  g1(13,10)=(-(T193*T210));
  g1(13,4)=(-(T210*T301));
  g1(13,16)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],13,484);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],13,10648);
end
end
