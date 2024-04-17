function [residual, g1, g2] = cv_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 12, 1);

%
% Model equations
%

T14 = y(4)^params(1);
T18 = y(7)^(-params(1));
T26 = y(7)^(1-params(1));
T43 = y(6)^(-params(6));
T110 = getPowerDeriv(y(6),(-params(6)),1);
lhs =y(8);
rhs =y(5)*(1-params(1))*T14*T18;
residual(1)= lhs-rhs;
lhs =y(11);
rhs =y(5)*params(1)*y(4)^(params(1)-1)*T26;
residual(2)= lhs-rhs;
lhs =y(5);
rhs =1-params(4)+y(5)*params(4)+x(1);
residual(3)= lhs-rhs;
lhs =params(8)*y(7)^params(7);
rhs =y(8)*T43;
residual(4)= lhs-rhs;
lhs =y(10);
rhs =y(4)-y(4)*(1-params(3));
residual(5)= lhs-rhs;
lhs =y(9);
rhs =T26*y(5)*T14;
residual(6)= lhs-rhs;
lhs =y(9);
rhs =y(6)+y(10);
residual(7)= lhs-rhs;
lhs =T43;
rhs =T43*params(2)*(1+y(12));
residual(8)= lhs-rhs;
lhs =1+y(12);
rhs =1+y(11)-params(3);
residual(9)= lhs-rhs;
lhs =y(2);
rhs =1/(1-params(6))*(y(6)^(1-params(6))-1)+params(2)*y(2);
residual(10)= lhs-rhs;
lhs =y(3);
rhs =(-params(8))*y(7)^(1+params(7))/(1+params(7))+params(2)*y(3);
residual(11)= lhs-rhs;
lhs =y(1);
rhs =y(2)+y(3);
residual(12)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(12, 12);

%
% Jacobian matrix
%

  g1(1,4)=(-(T18*y(5)*(1-params(1))*getPowerDeriv(y(4),params(1),1)));
  g1(1,5)=(-(T18*(1-params(1))*T14));
  g1(1,7)=(-(y(5)*(1-params(1))*T14*getPowerDeriv(y(7),(-params(1)),1)));
  g1(1,8)=1;
  g1(2,4)=(-(T26*y(5)*params(1)*getPowerDeriv(y(4),params(1)-1,1)));
  g1(2,5)=(-(T26*params(1)*y(4)^(params(1)-1)));
  g1(2,7)=(-(y(5)*params(1)*y(4)^(params(1)-1)*getPowerDeriv(y(7),1-params(1),1)));
  g1(2,11)=1;
  g1(3,5)=1-params(4);
  g1(4,6)=(-(y(8)*T110));
  g1(4,7)=params(8)*getPowerDeriv(y(7),params(7),1);
  g1(4,8)=(-T43);
  g1(5,4)=(-(1-(1-params(3))));
  g1(5,10)=1;
  g1(6,4)=(-(T26*y(5)*getPowerDeriv(y(4),params(1),1)));
  g1(6,5)=(-(T14*T26));
  g1(6,7)=(-(y(5)*T14*getPowerDeriv(y(7),1-params(1),1)));
  g1(6,9)=1;
  g1(7,6)=(-1);
  g1(7,9)=1;
  g1(7,10)=(-1);
  g1(8,6)=T110-params(2)*(1+y(12))*T110;
  g1(8,12)=(-(T43*params(2)));
  g1(9,11)=(-1);
  g1(9,12)=1;
  g1(10,2)=1-params(2);
  g1(10,6)=(-(1/(1-params(6))*getPowerDeriv(y(6),1-params(6),1)));
  g1(11,3)=1-params(2);
  g1(11,7)=(-((-params(8))*getPowerDeriv(y(7),1+params(7),1)/(1+params(7))));
  g1(12,1)=1;
  g1(12,2)=(-1);
  g1(12,3)=(-1);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],12,144);
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
