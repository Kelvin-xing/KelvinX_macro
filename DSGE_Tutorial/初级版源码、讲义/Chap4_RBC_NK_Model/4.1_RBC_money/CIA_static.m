function [residual, g1, g2] = CIA_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 15, 1);

%
% Model equations
%

T32 = exp(y(5))/exp(y(11));
T34 = params(2)*exp(y(15))*T32;
T46 = exp(y(2))^(1-params(3));
T47 = params(3)*exp(y(7))*exp(y(6))^(params(3)-1)*T46;
T54 = exp(y(7))*(1-params(3))*exp(y(6))^params(3)*exp(y(2))^(-params(3));
T153 = (-(exp(y(5))*exp(y(11))))/(exp(y(11))*exp(y(11)));
lhs =params(1)/(1-exp(y(2)));
rhs =exp(y(15))*exp(y(3));
residual(1)= lhs-rhs;
lhs =exp(y(15));
rhs =params(2)*exp(y(15))*(1+exp(y(4))-params(4));
residual(2)= lhs-rhs;
lhs =exp(y(15));
rhs =T34;
residual(3)= lhs-rhs;
lhs =exp(y(4));
rhs =T47;
residual(4)= lhs-rhs;
lhs =exp(y(3));
rhs =T54;
residual(5)= lhs-rhs;
lhs =exp(y(8));
rhs =T46*exp(y(7))*exp(y(6))^params(3);
residual(6)= lhs-rhs;
lhs =exp(y(6));
rhs =exp(y(9))+exp(y(6))*(1-params(4));
residual(7)= lhs-rhs;
lhs =exp(y(8));
rhs =exp(y(9))+exp(y(1));
residual(8)= lhs-rhs;
lhs =exp(y(13));
rhs =T32;
residual(9)= lhs-rhs;
lhs =exp(y(15));
rhs =params(2)*(exp(y(15))+exp(y(14)))/exp(y(11));
residual(10)= lhs-rhs;
lhs =y(10);
rhs =(1-params(5))*log(params(6))-y(11)+y(11)*params(5)+y(10)*params(5)+x(2);
residual(11)= lhs-rhs;
lhs =y(7);
rhs =y(7)*params(9)+x(1);
residual(12)= lhs-rhs;
residual(13) = y(10);
lhs =1/exp(y(1));
rhs =exp(y(15))+exp(y(14));
residual(14)= lhs-rhs;
lhs =exp(y(12));
rhs =exp(y(11))*exp(y(1));
residual(15)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(15, 15);

%
% Jacobian matrix
%

  g1(1,2)=(-(params(1)*(-exp(y(2)))))/((1-exp(y(2)))*(1-exp(y(2))));
  g1(1,3)=(-(exp(y(15))*exp(y(3))));
  g1(1,15)=(-(exp(y(15))*exp(y(3))));
  g1(2,4)=(-(params(2)*exp(y(15))*exp(y(4))));
  g1(2,15)=exp(y(15))-params(2)*exp(y(15))*(1+exp(y(4))-params(4));
  g1(3,5)=(-T34);
  g1(3,11)=(-(params(2)*exp(y(15))*T153));
  g1(3,15)=exp(y(15))-T34;
  g1(4,2)=(-(params(3)*exp(y(7))*exp(y(6))^(params(3)-1)*exp(y(2))*getPowerDeriv(exp(y(2)),1-params(3),1)));
  g1(4,4)=exp(y(4));
  g1(4,6)=(-(T46*params(3)*exp(y(7))*exp(y(6))*getPowerDeriv(exp(y(6)),params(3)-1,1)));
  g1(4,7)=(-T47);
  g1(5,2)=(-(exp(y(7))*(1-params(3))*exp(y(6))^params(3)*exp(y(2))*getPowerDeriv(exp(y(2)),(-params(3)),1)));
  g1(5,3)=exp(y(3));
  g1(5,6)=(-(exp(y(2))^(-params(3))*exp(y(7))*(1-params(3))*exp(y(6))*getPowerDeriv(exp(y(6)),params(3),1)));
  g1(5,7)=(-T54);
  g1(6,2)=(-(exp(y(7))*exp(y(6))^params(3)*exp(y(2))*getPowerDeriv(exp(y(2)),1-params(3),1)));
  g1(6,6)=(-(T46*exp(y(7))*exp(y(6))*getPowerDeriv(exp(y(6)),params(3),1)));
  g1(6,7)=(-(T46*exp(y(7))*exp(y(6))^params(3)));
  g1(6,8)=exp(y(8));
  g1(7,6)=exp(y(6))-exp(y(6))*(1-params(4));
  g1(7,9)=(-exp(y(9)));
  g1(8,1)=(-exp(y(1)));
  g1(8,8)=exp(y(8));
  g1(8,9)=(-exp(y(9)));
  g1(9,5)=(-T32);
  g1(9,11)=(-T153);
  g1(9,13)=exp(y(13));
  g1(10,11)=(-((-(exp(y(11))*params(2)*(exp(y(15))+exp(y(14)))))/(exp(y(11))*exp(y(11)))));
  g1(10,14)=(-(params(2)*exp(y(14))/exp(y(11))));
  g1(10,15)=exp(y(15))-exp(y(15))*params(2)/exp(y(11));
  g1(11,10)=1-params(5);
  g1(11,11)=(-((-1)+params(5)));
  g1(12,7)=1-params(9);
  g1(13,10)=1;
  g1(14,1)=(-exp(y(1)))/(exp(y(1))*exp(y(1)));
  g1(14,14)=(-exp(y(14)));
  g1(14,15)=(-exp(y(15)));
  g1(15,1)=(-(exp(y(11))*exp(y(1))));
  g1(15,11)=(-(exp(y(11))*exp(y(1))));
  g1(15,12)=exp(y(12));
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
%
% Hessian matrix
%

  g2 = sparse([],[],[],15,225);
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
