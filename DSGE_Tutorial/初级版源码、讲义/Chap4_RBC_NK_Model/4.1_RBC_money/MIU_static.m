function [residual, g1, g2] = MIU_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 13, 1);

%
% Model equations
%

T19 = 1/exp(y(1));
T33 = exp(y(5))/exp(y(11));
T47 = exp(y(2))^(1-params(3));
T48 = params(3)*exp(y(7))*exp(y(6))^(params(3)-1)*T47;
T55 = exp(y(7))*(1-params(3))*exp(y(6))^params(3)*exp(y(2))^(-params(3));
T79 = params(5)^params(10)*exp(y(1))^params(10);
T82 = (exp(y(5))/(exp(y(5))-1))^params(10);
T111 = (-exp(y(1)))/(exp(y(1))*exp(y(1)));
T175 = (-(exp(y(5))*exp(y(11))))/(exp(y(11))*exp(y(11)));
lhs =params(1)/(1-exp(y(2)));
rhs =exp(y(3))/exp(y(1));
residual(1)= lhs-rhs;
lhs =T19;
rhs =params(2)*T19*(1+exp(y(4))-params(4));
residual(2)= lhs-rhs;
lhs =T19;
rhs =params(2)*T19*T33;
residual(3)= lhs-rhs;
lhs =exp(y(4));
rhs =T48;
residual(4)= lhs-rhs;
lhs =exp(y(3));
rhs =T55;
residual(5)= lhs-rhs;
lhs =exp(y(8));
rhs =T47*exp(y(7))*exp(y(6))^params(3);
residual(6)= lhs-rhs;
lhs =exp(y(6));
rhs =exp(y(9))+exp(y(6))*(1-params(4));
residual(7)= lhs-rhs;
lhs =exp(y(8));
rhs =exp(y(1))+exp(y(9));
residual(8)= lhs-rhs;
lhs =exp(y(13));
rhs =T33;
residual(9)= lhs-rhs;
lhs =exp(y(12));
rhs =T79*T82;
residual(10)= lhs-rhs;
lhs =y(10);
rhs =(1-params(6))*log(params(7))-y(11)+y(11)*params(6)+y(10)*params(6)+x(2);
residual(11)= lhs-rhs;
lhs =y(7);
rhs =y(7)*params(11)+x(1);
residual(12)= lhs-rhs;
residual(13) = y(10);
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(13, 13);

%
% Jacobian matrix
%

  g1(1,1)=(-((-(exp(y(3))*exp(y(1))))/(exp(y(1))*exp(y(1)))));
  g1(1,2)=(-(params(1)*(-exp(y(2)))))/((1-exp(y(2)))*(1-exp(y(2))));
  g1(1,3)=(-(exp(y(3))/exp(y(1))));
  g1(2,1)=T111-params(2)*(1+exp(y(4))-params(4))*T111;
  g1(2,4)=(-(params(2)*T19*exp(y(4))));
  g1(3,1)=T111-params(2)*T33*T111;
  g1(3,5)=(-(params(2)*T19*T33));
  g1(3,11)=(-(params(2)*T19*T175));
  g1(4,2)=(-(params(3)*exp(y(7))*exp(y(6))^(params(3)-1)*exp(y(2))*getPowerDeriv(exp(y(2)),1-params(3),1)));
  g1(4,4)=exp(y(4));
  g1(4,6)=(-(T47*params(3)*exp(y(7))*exp(y(6))*getPowerDeriv(exp(y(6)),params(3)-1,1)));
  g1(4,7)=(-T48);
  g1(5,2)=(-(exp(y(7))*(1-params(3))*exp(y(6))^params(3)*exp(y(2))*getPowerDeriv(exp(y(2)),(-params(3)),1)));
  g1(5,3)=exp(y(3));
  g1(5,6)=(-(exp(y(2))^(-params(3))*exp(y(7))*(1-params(3))*exp(y(6))*getPowerDeriv(exp(y(6)),params(3),1)));
  g1(5,7)=(-T55);
  g1(6,2)=(-(exp(y(7))*exp(y(6))^params(3)*exp(y(2))*getPowerDeriv(exp(y(2)),1-params(3),1)));
  g1(6,6)=(-(T47*exp(y(7))*exp(y(6))*getPowerDeriv(exp(y(6)),params(3),1)));
  g1(6,7)=(-(T47*exp(y(7))*exp(y(6))^params(3)));
  g1(6,8)=exp(y(8));
  g1(7,6)=exp(y(6))-exp(y(6))*(1-params(4));
  g1(7,9)=(-exp(y(9)));
  g1(8,1)=(-exp(y(1)));
  g1(8,8)=exp(y(8));
  g1(8,9)=(-exp(y(9)));
  g1(9,5)=(-T33);
  g1(9,11)=(-T175);
  g1(9,13)=exp(y(13));
  g1(10,1)=(-(T82*params(5)^params(10)*exp(y(1))*getPowerDeriv(exp(y(1)),params(10),1)));
  g1(10,5)=(-(T79*(exp(y(5))*(exp(y(5))-1)-exp(y(5))*exp(y(5)))/((exp(y(5))-1)*(exp(y(5))-1))*getPowerDeriv(exp(y(5))/(exp(y(5))-1),params(10),1)));
  g1(10,12)=exp(y(12));
  g1(11,10)=1-params(6);
  g1(11,11)=(-((-1)+params(6)));
  g1(12,7)=1-params(11);
  g1(13,10)=1;
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
