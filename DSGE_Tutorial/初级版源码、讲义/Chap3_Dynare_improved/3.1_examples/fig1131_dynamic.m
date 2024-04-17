function [residual, g1, g2, g3] = Fig1131_dynamic(y, x, params, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(2, 1);
T29 = params(1)*y(4)^(-params(2));
T35 = T29*(1+x(it_-1, 2))/(1+x(it_, 2));
T48 = y(1)^(params(4)-1);
T50 = (1-params(3))*(1-x(it_, 1))/(1-x(it_-1, 1))+params(5)*params(4)*(1-x(it_, 3))/(1-x(it_-1, 1))*T48;
lhs =y(3);
rhs =params(5)*y(1)^params(4)+y(1)*(1-params(3))-y(2)-x(it_, 4);
residual(1)= lhs-rhs;
lhs =y(2)^(-params(2));
rhs =T35*T50;
residual(2)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(2, 8);

%
% Jacobian matrix
%

g1(1,2)=1;
g1(1,1)=(-(1-params(3)+params(5)*getPowerDeriv(y(1),params(4),1)));
g1(1,3)=1;
g1(1,8)=1;
g1(2,2)=getPowerDeriv(y(2),(-params(2)),1);
g1(2,4)=(-(T50*(1+x(it_-1, 2))/(1+x(it_, 2))*params(1)*getPowerDeriv(y(4),(-params(2)),1)));
g1(2,1)=(-(T35*params(5)*params(4)*(1-x(it_, 3))/(1-x(it_-1, 1))*getPowerDeriv(y(1),params(4)-1,1)));
g1(2,5)=(-(T35*((1-params(3))*(1-x(it_, 1))/((1-x(it_-1, 1))*(1-x(it_-1, 1)))+T48*params(5)*params(4)*(1-x(it_, 3))/((1-x(it_-1, 1))*(1-x(it_-1, 1))))));
g1(2,5)=(-(T35*(-(1-params(3)))/(1-x(it_-1, 1))));
g1(2,6)=(-(T50*T29*1/(1+x(it_, 2))));
g1(2,6)=(-(T50*T29*(-(1+x(it_-1, 2)))/((1+x(it_, 2))*(1+x(it_, 2)))));
g1(2,7)=(-(T35*T48*params(5)*params(4)*(-1)/(1-x(it_-1, 1))));
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],2,64);
end
if nargout >= 4,
%
% Third order derivatives
%

  g3 = sparse([],[],[],2,512);
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
