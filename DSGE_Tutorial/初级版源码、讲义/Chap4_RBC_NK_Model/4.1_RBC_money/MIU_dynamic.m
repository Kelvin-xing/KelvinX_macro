function [residual, g1, g2, g3] = MIU_dynamic(y, x, params, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(13, 1);
T23 = 1/exp(y(19));
T36 = exp(y(10))/exp(y(21));
T52 = exp(y(7))^(1-params(3));
T53 = params(3)*exp(y(12))*exp(y(1))^(params(3)-1)*T52;
T60 = exp(y(12))*(1-params(3))*exp(y(1))^params(3)*exp(y(7))^(-params(3));
T86 = params(5)^params(10)*exp(y(6))^params(10);
T89 = (exp(y(10))/(exp(y(10))-1))^params(10);
T192 = (-(exp(y(10))*exp(y(21))))/(exp(y(21))*exp(y(21)));
lhs =params(1)/(1-exp(y(7)));
rhs =exp(y(8))/exp(y(6));
residual(1)= lhs-rhs;
lhs =1/exp(y(6));
rhs =params(2)*T23*(1+exp(y(20))-params(4));
residual(2)= lhs-rhs;
lhs =1/exp(y(6));
rhs =params(2)*T23*T36;
residual(3)= lhs-rhs;
lhs =exp(y(9));
rhs =T53;
residual(4)= lhs-rhs;
lhs =exp(y(8));
rhs =T60;
residual(5)= lhs-rhs;
lhs =exp(y(13));
rhs =T52*exp(y(12))*exp(y(1))^params(3);
residual(6)= lhs-rhs;
lhs =exp(y(11));
rhs =exp(y(14))+exp(y(1))*(1-params(4));
residual(7)= lhs-rhs;
lhs =exp(y(13));
rhs =exp(y(6))+exp(y(14));
residual(8)= lhs-rhs;
lhs =exp(y(18));
rhs =T36;
residual(9)= lhs-rhs;
lhs =exp(y(17));
rhs =T86*T89;
residual(10)= lhs-rhs;
lhs =y(15);
rhs =(1-params(6))*log(params(7))-y(16)+params(6)*y(4)+params(6)*y(3)+x(it_, 2);
residual(11)= lhs-rhs;
lhs =y(12);
rhs =params(11)*y(2)+x(it_, 1);
residual(12)= lhs-rhs;
lhs =y(15);
rhs =y(17)-y(5);
residual(13)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(13, 23);

%
% Jacobian matrix
%

g1(1,6)=(-((-(exp(y(8))*exp(y(6))))/(exp(y(6))*exp(y(6)))));
g1(1,7)=(-(params(1)*(-exp(y(7)))))/((1-exp(y(7)))*(1-exp(y(7))));
g1(1,8)=(-(exp(y(8))/exp(y(6))));
g1(2,6)=(-exp(y(6)))/(exp(y(6))*exp(y(6)));
g1(2,19)=(-(params(2)*(1+exp(y(20))-params(4))*(-exp(y(19)))/(exp(y(19))*exp(y(19)))));
g1(2,20)=(-(params(2)*T23*exp(y(20))));
g1(3,6)=(-exp(y(6)))/(exp(y(6))*exp(y(6)));
g1(3,19)=(-(params(2)*T36*(-exp(y(19)))/(exp(y(19))*exp(y(19)))));
g1(3,10)=(-(params(2)*T23*T36));
g1(3,21)=(-(params(2)*T23*T192));
g1(4,7)=(-(params(3)*exp(y(12))*exp(y(1))^(params(3)-1)*exp(y(7))*getPowerDeriv(exp(y(7)),1-params(3),1)));
g1(4,9)=exp(y(9));
g1(4,1)=(-(T52*params(3)*exp(y(12))*exp(y(1))*getPowerDeriv(exp(y(1)),params(3)-1,1)));
g1(4,12)=(-T53);
g1(5,7)=(-(exp(y(12))*(1-params(3))*exp(y(1))^params(3)*exp(y(7))*getPowerDeriv(exp(y(7)),(-params(3)),1)));
g1(5,8)=exp(y(8));
g1(5,1)=(-(exp(y(7))^(-params(3))*exp(y(12))*(1-params(3))*exp(y(1))*getPowerDeriv(exp(y(1)),params(3),1)));
g1(5,12)=(-T60);
g1(6,7)=(-(exp(y(12))*exp(y(1))^params(3)*exp(y(7))*getPowerDeriv(exp(y(7)),1-params(3),1)));
g1(6,1)=(-(T52*exp(y(12))*exp(y(1))*getPowerDeriv(exp(y(1)),params(3),1)));
g1(6,12)=(-(T52*exp(y(12))*exp(y(1))^params(3)));
g1(6,13)=exp(y(13));
g1(7,1)=(-(exp(y(1))*(1-params(4))));
g1(7,11)=exp(y(11));
g1(7,14)=(-exp(y(14)));
g1(8,6)=(-exp(y(6)));
g1(8,13)=exp(y(13));
g1(8,14)=(-exp(y(14)));
g1(9,10)=(-T36);
g1(9,21)=(-T192);
g1(9,18)=exp(y(18));
g1(10,6)=(-(T89*params(5)^params(10)*exp(y(6))*getPowerDeriv(exp(y(6)),params(10),1)));
g1(10,10)=(-(T86*(exp(y(10))*(exp(y(10))-1)-exp(y(10))*exp(y(10)))/((exp(y(10))-1)*(exp(y(10))-1))*getPowerDeriv(exp(y(10))/(exp(y(10))-1),params(10),1)));
g1(10,17)=exp(y(17));
g1(11,3)=(-params(6));
g1(11,15)=1;
g1(11,4)=(-params(6));
g1(11,16)=1;
g1(11,23)=(-1);
g1(12,2)=(-params(11));
g1(12,12)=1;
g1(12,22)=(-1);
g1(13,15)=1;
g1(13,5)=1;
g1(13,17)=(-1);
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],13,529);
end
if nargout >= 4,
%
% Third order derivatives
%

  g3 = sparse([],[],[],13,12167);
end
end

%
% The k-th derivative of x^p
%
function dxp=getPowerDeriv(x,p,k)
    if (abs(x) < 1e-012) && (p > 0) && (k >= p) && (abs(p - round(p)) < 1e-012)
        dxp = 0;
    else
        dxp = x^(p-k);
        for i=0:k-1
            dxp = dxp*p;
            p = p-1;
        end
    end
end
