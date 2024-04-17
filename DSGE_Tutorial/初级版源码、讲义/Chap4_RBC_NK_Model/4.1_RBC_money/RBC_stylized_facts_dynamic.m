function [residual, g1, g2, g3] = RBC_stylized_facts_dynamic(y, x, params, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(10, 1);
T12 = exp(y(5))^(-params(5));
T48 = exp(y(4))^(1-params(2));
T49 = exp(y(3))*exp(y(2))^params(2)*T48;
T56 = T48*exp(y(3))*params(2)*exp(y(2))^(params(2)-1);
T62 = exp(y(2))^params(2)*exp(y(3))*(1-params(2))*exp(y(4))^(-params(2));
T108 = exp(y(5))*getPowerDeriv(exp(y(5)),(-params(5)),1);
lhs =T12;
rhs =params(4)*exp(y(13))^(-params(5))*(y(14)+1-params(3));
residual(1)= lhs-rhs;
lhs =T12;
rhs =exp(y(13))^(-params(5))*params(4)*exp(y(8));
residual(2)= lhs-rhs;
lhs =params(1)/(1-exp(y(4)));
rhs =T12*exp(y(9));
residual(3)= lhs-rhs;
lhs =exp(y(11));
rhs =T49;
residual(4)= lhs-rhs;
lhs =y(10);
rhs =T56;
residual(5)= lhs-rhs;
lhs =exp(y(9));
rhs =T62;
residual(6)= lhs-rhs;
lhs =exp(y(11));
rhs =exp(y(5))+exp(y(7));
residual(7)= lhs-rhs;
lhs =exp(y(6));
rhs =exp(y(7))+(1-params(3))*exp(y(2));
residual(8)= lhs-rhs;
lhs =y(3);
rhs =params(6)*y(1)+x(it_, 1);
residual(9)= lhs-rhs;
lhs =exp(y(12));
rhs =exp(y(11))/exp(y(4));
residual(10)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(10, 15);

%
% Jacobian matrix
%

g1(1,5)=T108;
g1(1,13)=(-((y(14)+1-params(3))*params(4)*exp(y(13))*getPowerDeriv(exp(y(13)),(-params(5)),1)));
g1(1,14)=(-(params(4)*exp(y(13))^(-params(5))));
g1(2,5)=T108;
g1(2,13)=(-(params(4)*exp(y(8))*exp(y(13))*getPowerDeriv(exp(y(13)),(-params(5)),1)));
g1(2,8)=(-(exp(y(13))^(-params(5))*params(4)*exp(y(8))));
g1(3,4)=(-(params(1)*(-exp(y(4)))))/((1-exp(y(4)))*(1-exp(y(4))));
g1(3,5)=(-(exp(y(9))*T108));
g1(3,9)=(-(T12*exp(y(9))));
g1(4,3)=(-T49);
g1(4,4)=(-(exp(y(3))*exp(y(2))^params(2)*exp(y(4))*getPowerDeriv(exp(y(4)),1-params(2),1)));
g1(4,2)=(-(T48*exp(y(3))*exp(y(2))*getPowerDeriv(exp(y(2)),params(2),1)));
g1(4,11)=exp(y(11));
g1(5,3)=(-T56);
g1(5,4)=(-(exp(y(3))*params(2)*exp(y(2))^(params(2)-1)*exp(y(4))*getPowerDeriv(exp(y(4)),1-params(2),1)));
g1(5,2)=(-(T48*exp(y(3))*params(2)*exp(y(2))*getPowerDeriv(exp(y(2)),params(2)-1,1)));
g1(5,10)=1;
g1(6,3)=(-T62);
g1(6,4)=(-(exp(y(2))^params(2)*exp(y(3))*(1-params(2))*exp(y(4))*getPowerDeriv(exp(y(4)),(-params(2)),1)));
g1(6,2)=(-(exp(y(4))^(-params(2))*exp(y(3))*(1-params(2))*exp(y(2))*getPowerDeriv(exp(y(2)),params(2),1)));
g1(6,9)=exp(y(9));
g1(7,5)=(-exp(y(5)));
g1(7,7)=(-exp(y(7)));
g1(7,11)=exp(y(11));
g1(8,2)=(-((1-params(3))*exp(y(2))));
g1(8,6)=exp(y(6));
g1(8,7)=(-exp(y(7)));
g1(9,1)=(-params(6));
g1(9,3)=1;
g1(9,15)=(-1);
g1(10,4)=(-((-(exp(y(4))*exp(y(11))))/(exp(y(4))*exp(y(4)))));
g1(10,11)=(-(exp(y(11))/exp(y(4))));
g1(10,12)=exp(y(12));
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],10,225);
end
if nargout >= 4,
%
% Third order derivatives
%

  g3 = sparse([],[],[],10,3375);
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
