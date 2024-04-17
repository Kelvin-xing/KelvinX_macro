function [residual, g1, g2] = Fig1131_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 2, 1);

%
% Model equations
%

T41 = (1-params(3))*(1-x(1))/(1-x(1))+params(5)*params(4)*(1-x(3))/(1-x(1))*y(2)^(params(4)-1);
lhs =y(2);
rhs =params(5)*y(2)^params(4)+y(2)*(1-params(3))-y(1)-x(4);
residual(1)= lhs-rhs;
lhs =y(1)^(-params(2));
rhs =y(1)^(-params(2))*params(1)*T41;
residual(2)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(2, 2);

%
% Jacobian matrix
%

  g1(1,1)=1;
  g1(1,2)=1-(1-params(3)+params(5)*getPowerDeriv(y(2),params(4),1));
  g1(2,1)=getPowerDeriv(y(1),(-params(2)),1)-T41*params(1)*getPowerDeriv(y(1),(-params(2)),1);
  g1(2,2)=(-(y(1)^(-params(2))*params(1)*params(5)*params(4)*(1-x(3))/(1-x(1))*getPowerDeriv(y(2),params(4)-1,1)));
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],2,4);
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
