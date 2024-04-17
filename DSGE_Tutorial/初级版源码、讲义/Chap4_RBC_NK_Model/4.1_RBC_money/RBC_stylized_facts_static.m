function [residual, g1, g2] = RBC_stylized_facts_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 10, 1);

%
% Model equations
%

T12 = exp(y(3))^(-params(5));
T45 = exp(y(2))^(1-params(2));
T46 = exp(y(1))*exp(y(4))^params(2)*T45;
T52 = T45*exp(y(1))*params(2)*exp(y(4))^(params(2)-1);
T58 = exp(y(4))^params(2)*exp(y(1))*(1-params(2))*exp(y(2))^(-params(2));
T101 = exp(y(3))*getPowerDeriv(exp(y(3)),(-params(5)),1);
lhs =T12;
rhs =T12*params(4)*(y(8)+1-params(3));
residual(1)= lhs-rhs;
lhs =T12;
rhs =T12*params(4)*exp(y(6));
residual(2)= lhs-rhs;
lhs =params(1)/(1-exp(y(2)));
rhs =T12*exp(y(7));
residual(3)= lhs-rhs;
lhs =exp(y(9));
rhs =T46;
residual(4)= lhs-rhs;
lhs =y(8);
rhs =T52;
residual(5)= lhs-rhs;
lhs =exp(y(7));
rhs =T58;
residual(6)= lhs-rhs;
lhs =exp(y(9));
rhs =exp(y(3))+exp(y(5));
residual(7)= lhs-rhs;
lhs =exp(y(4));
rhs =exp(y(5))+(1-params(3))*exp(y(4));
residual(8)= lhs-rhs;
lhs =y(1);
rhs =y(1)*params(6)+x(1);
residual(9)= lhs-rhs;
lhs =exp(y(10));
rhs =exp(y(9))/exp(y(2));
residual(10)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(10, 10);

%
% Jacobian matrix
%

  g1(1,3)=T101-(y(8)+1-params(3))*params(4)*T101;
  g1(1,8)=(-(T12*params(4)));
  g1(2,3)=T101-params(4)*exp(y(6))*T101;
  g1(2,6)=(-(T12*params(4)*exp(y(6))));
  g1(3,2)=(-(params(1)*(-exp(y(2)))))/((1-exp(y(2)))*(1-exp(y(2))));
  g1(3,3)=(-(exp(y(7))*T101));
  g1(3,7)=(-(T12*exp(y(7))));
  g1(4,1)=(-T46);
  g1(4,2)=(-(exp(y(1))*exp(y(4))^params(2)*exp(y(2))*getPowerDeriv(exp(y(2)),1-params(2),1)));
  g1(4,4)=(-(T45*exp(y(1))*exp(y(4))*getPowerDeriv(exp(y(4)),params(2),1)));
  g1(4,9)=exp(y(9));
  g1(5,1)=(-T52);
  g1(5,2)=(-(exp(y(1))*params(2)*exp(y(4))^(params(2)-1)*exp(y(2))*getPowerDeriv(exp(y(2)),1-params(2),1)));
  g1(5,4)=(-(T45*exp(y(1))*params(2)*exp(y(4))*getPowerDeriv(exp(y(4)),params(2)-1,1)));
  g1(5,8)=1;
  g1(6,1)=(-T58);
  g1(6,2)=(-(exp(y(4))^params(2)*exp(y(1))*(1-params(2))*exp(y(2))*getPowerDeriv(exp(y(2)),(-params(2)),1)));
  g1(6,4)=(-(exp(y(2))^(-params(2))*exp(y(1))*(1-params(2))*exp(y(4))*getPowerDeriv(exp(y(4)),params(2),1)));
  g1(6,7)=exp(y(7));
  g1(7,3)=(-exp(y(3)));
  g1(7,5)=(-exp(y(5)));
  g1(7,9)=exp(y(9));
  g1(8,4)=exp(y(4))-(1-params(3))*exp(y(4));
  g1(8,5)=(-exp(y(5)));
  g1(9,1)=1-params(6);
  g1(10,2)=(-((-(exp(y(2))*exp(y(9))))/(exp(y(2))*exp(y(2)))));
  g1(10,9)=(-(exp(y(9))/exp(y(2))));
  g1(10,10)=exp(y(10));
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
