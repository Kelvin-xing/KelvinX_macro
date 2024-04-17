function [residual, g1, g2, g3] = CIA_dynamic(y, x, params, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(15, 1);
T34 = exp(y(10))/exp(y(22));
T50 = exp(y(7))^(1-params(3));
T51 = params(3)*exp(y(12))*exp(y(1))^(params(3)-1)*T50;
T58 = exp(y(12))*(1-params(3))*exp(y(1))^params(3)*exp(y(7))^(-params(3));
T170 = (-(exp(y(10))*exp(y(22))))/(exp(y(22))*exp(y(22)));
lhs =params(1)/(1-exp(y(7)));
rhs =exp(y(20))*exp(y(8));
residual(1)= lhs-rhs;
lhs =exp(y(20));
rhs =params(2)*exp(y(24))*(1+exp(y(21))-params(4));
residual(2)= lhs-rhs;
lhs =exp(y(20));
rhs =params(2)*exp(y(24))*T34;
residual(3)= lhs-rhs;
lhs =exp(y(9));
rhs =T51;
residual(4)= lhs-rhs;
lhs =exp(y(8));
rhs =T58;
residual(5)= lhs-rhs;
lhs =exp(y(13));
rhs =T50*exp(y(12))*exp(y(1))^params(3);
residual(6)= lhs-rhs;
lhs =exp(y(11));
rhs =exp(y(14))+exp(y(1))*(1-params(4));
residual(7)= lhs-rhs;
lhs =exp(y(13));
rhs =exp(y(14))+exp(y(6));
residual(8)= lhs-rhs;
lhs =exp(y(18));
rhs =T34;
residual(9)= lhs-rhs;
lhs =exp(y(20));
rhs =params(2)*(exp(y(24))+exp(y(23)))/exp(y(22));
residual(10)= lhs-rhs;
lhs =y(15);
rhs =(1-params(5))*log(params(6))-y(16)+params(5)*y(4)+params(5)*y(3)+x(it_, 2);
residual(11)= lhs-rhs;
lhs =y(12);
rhs =params(9)*y(2)+x(it_, 1);
residual(12)= lhs-rhs;
lhs =y(15);
rhs =y(17)-y(5);
residual(13)= lhs-rhs;
lhs =1/exp(y(6));
rhs =exp(y(20))+exp(y(19));
residual(14)= lhs-rhs;
lhs =exp(y(5));
rhs =exp(y(6))*exp(y(16));
residual(15)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(15, 26);

%
% Jacobian matrix
%

g1(1,7)=(-(params(1)*(-exp(y(7)))))/((1-exp(y(7)))*(1-exp(y(7))));
g1(1,8)=(-(exp(y(20))*exp(y(8))));
g1(1,20)=(-(exp(y(20))*exp(y(8))));
g1(2,21)=(-(params(2)*exp(y(24))*exp(y(21))));
g1(2,20)=exp(y(20));
g1(2,24)=(-(params(2)*exp(y(24))*(1+exp(y(21))-params(4))));
g1(3,10)=(-(params(2)*exp(y(24))*T34));
g1(3,22)=(-(params(2)*exp(y(24))*T170));
g1(3,20)=exp(y(20));
g1(3,24)=(-(params(2)*exp(y(24))*T34));
g1(4,7)=(-(params(3)*exp(y(12))*exp(y(1))^(params(3)-1)*exp(y(7))*getPowerDeriv(exp(y(7)),1-params(3),1)));
g1(4,9)=exp(y(9));
g1(4,1)=(-(T50*params(3)*exp(y(12))*exp(y(1))*getPowerDeriv(exp(y(1)),params(3)-1,1)));
g1(4,12)=(-T51);
g1(5,7)=(-(exp(y(12))*(1-params(3))*exp(y(1))^params(3)*exp(y(7))*getPowerDeriv(exp(y(7)),(-params(3)),1)));
g1(5,8)=exp(y(8));
g1(5,1)=(-(exp(y(7))^(-params(3))*exp(y(12))*(1-params(3))*exp(y(1))*getPowerDeriv(exp(y(1)),params(3),1)));
g1(5,12)=(-T58);
g1(6,7)=(-(exp(y(12))*exp(y(1))^params(3)*exp(y(7))*getPowerDeriv(exp(y(7)),1-params(3),1)));
g1(6,1)=(-(T50*exp(y(12))*exp(y(1))*getPowerDeriv(exp(y(1)),params(3),1)));
g1(6,12)=(-(T50*exp(y(12))*exp(y(1))^params(3)));
g1(6,13)=exp(y(13));
g1(7,1)=(-(exp(y(1))*(1-params(4))));
g1(7,11)=exp(y(11));
g1(7,14)=(-exp(y(14)));
g1(8,6)=(-exp(y(6)));
g1(8,13)=exp(y(13));
g1(8,14)=(-exp(y(14)));
g1(9,10)=(-T34);
g1(9,22)=(-T170);
g1(9,18)=exp(y(18));
g1(10,22)=(-((-(exp(y(22))*params(2)*(exp(y(24))+exp(y(23)))))/(exp(y(22))*exp(y(22)))));
g1(10,23)=(-(params(2)*exp(y(23))/exp(y(22))));
g1(10,20)=exp(y(20));
g1(10,24)=(-(params(2)*exp(y(24))/exp(y(22))));
g1(11,3)=(-params(5));
g1(11,15)=1;
g1(11,4)=(-params(5));
g1(11,16)=1;
g1(11,26)=(-1);
g1(12,2)=(-params(9));
g1(12,12)=1;
g1(12,25)=(-1);
g1(13,15)=1;
g1(13,5)=1;
g1(13,17)=(-1);
g1(14,6)=(-exp(y(6)))/(exp(y(6))*exp(y(6)));
g1(14,19)=(-exp(y(19)));
g1(14,20)=(-exp(y(20)));
g1(15,6)=(-(exp(y(6))*exp(y(16))));
g1(15,16)=(-(exp(y(6))*exp(y(16))));
g1(15,5)=exp(y(5));
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],15,676);
end
if nargout >= 4,
%
% Third order derivatives
%

  g3 = sparse([],[],[],15,17576);
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
