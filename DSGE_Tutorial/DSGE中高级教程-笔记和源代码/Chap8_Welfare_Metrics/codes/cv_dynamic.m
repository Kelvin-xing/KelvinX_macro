function [residual, g1, g2, g3] = cv_dynamic(y, x, params, it_)
%
% Status : Computes dynamic model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(12, 1);
T14 = y(1)^params(1);
T18 = y(9)^(-params(1));
T25 = y(7)*params(1)*y(1)^(params(1)-1);
T26 = y(9)^(1-params(1));
T44 = y(8)^(-params(6));
T96 = getPowerDeriv(y(1),params(1),1);
T117 = getPowerDeriv(y(8),(-params(6)),1);
T126 = getPowerDeriv(y(9),(-params(1)),1);
T129 = getPowerDeriv(y(9),1-params(1),1);
lhs =y(10);
rhs =y(7)*(1-params(1))*T14*T18;
residual(1)= lhs-rhs;
lhs =y(13);
rhs =T25*T26;
residual(2)= lhs-rhs;
lhs =y(7);
rhs =1-params(4)+params(4)*y(2)+x(it_, 1);
residual(3)= lhs-rhs;
lhs =params(8)*y(9)^params(7);
rhs =y(10)*T44;
residual(4)= lhs-rhs;
lhs =y(12);
rhs =y(6)-y(1)*(1-params(3));
residual(5)= lhs-rhs;
lhs =y(11);
rhs =T26*y(7)*T14;
residual(6)= lhs-rhs;
lhs =y(11);
rhs =y(8)+y(12);
residual(7)= lhs-rhs;
lhs =T44;
rhs =params(2)*(1+y(14))*y(17)^(-params(6));
residual(8)= lhs-rhs;
lhs =1+y(14);
rhs =1+y(18)-params(3);
residual(9)= lhs-rhs;
lhs =y(4);
rhs =1/(1-params(6))*(y(8)^(1-params(6))-1)+params(2)*y(15);
residual(10)= lhs-rhs;
lhs =y(5);
rhs =(-params(8))*y(9)^(1+params(7))/(1+params(7))+params(2)*y(16);
residual(11)= lhs-rhs;
lhs =y(3);
rhs =y(4)+y(5);
residual(12)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(12, 19);

%
% Jacobian matrix
%

g1(1,1)=(-(T18*y(7)*(1-params(1))*T96));
g1(1,7)=(-(T18*(1-params(1))*T14));
g1(1,9)=(-(y(7)*(1-params(1))*T14*T126));
g1(1,10)=1;
g1(2,1)=(-(T26*y(7)*params(1)*getPowerDeriv(y(1),params(1)-1,1)));
g1(2,7)=(-(T26*params(1)*y(1)^(params(1)-1)));
g1(2,9)=(-(T25*T129));
g1(2,13)=1;
g1(3,2)=(-params(4));
g1(3,7)=1;
g1(3,19)=(-1);
g1(4,8)=(-(y(10)*T117));
g1(4,9)=params(8)*getPowerDeriv(y(9),params(7),1);
g1(4,10)=(-T44);
g1(5,1)=1-params(3);
g1(5,6)=(-1);
g1(5,12)=1;
g1(6,1)=(-(T26*y(7)*T96));
g1(6,7)=(-(T14*T26));
g1(6,9)=(-(y(7)*T14*T129));
g1(6,11)=1;
g1(7,8)=(-1);
g1(7,11)=1;
g1(7,12)=(-1);
g1(8,8)=T117;
g1(8,17)=(-(params(2)*(1+y(14))*getPowerDeriv(y(17),(-params(6)),1)));
g1(8,14)=(-(params(2)*y(17)^(-params(6))));
g1(9,18)=(-1);
g1(9,14)=1;
g1(10,4)=1;
g1(10,15)=(-params(2));
g1(10,8)=(-(1/(1-params(6))*getPowerDeriv(y(8),1-params(6),1)));
g1(11,5)=1;
g1(11,16)=(-params(2));
g1(11,9)=(-((-params(8))*getPowerDeriv(y(9),1+params(7),1)/(1+params(7))));
g1(12,3)=1;
g1(12,4)=(-1);
g1(12,5)=(-1);
end
if nargout >= 3,
%
% Hessian matrix
%

  v2 = zeros(34,3);
v2(1,1)=1;
v2(1,2)=1;
v2(1,3)=(-(T18*y(7)*(1-params(1))*getPowerDeriv(y(1),params(1),2)));
v2(2,1)=1;
v2(2,2)=115;
v2(2,3)=(-(T18*(1-params(1))*T96));
v2(3,1)=1;
v2(3,2)=7;
v2(3,3)=v2(2,3);
v2(4,1)=1;
v2(4,2)=153;
v2(4,3)=(-(y(7)*(1-params(1))*T96*T126));
v2(5,1)=1;
v2(5,2)=9;
v2(5,3)=v2(4,3);
v2(6,1)=1;
v2(6,2)=159;
v2(6,3)=(-((1-params(1))*T14*T126));
v2(7,1)=1;
v2(7,2)=123;
v2(7,3)=v2(6,3);
v2(8,1)=1;
v2(8,2)=161;
v2(8,3)=(-(y(7)*(1-params(1))*T14*getPowerDeriv(y(9),(-params(1)),2)));
v2(9,1)=2;
v2(9,2)=1;
v2(9,3)=(-(T26*y(7)*params(1)*getPowerDeriv(y(1),params(1)-1,2)));
v2(10,1)=2;
v2(10,2)=115;
v2(10,3)=(-(T26*params(1)*getPowerDeriv(y(1),params(1)-1,1)));
v2(11,1)=2;
v2(11,2)=7;
v2(11,3)=v2(10,3);
v2(12,1)=2;
v2(12,2)=153;
v2(12,3)=(-(y(7)*params(1)*getPowerDeriv(y(1),params(1)-1,1)*T129));
v2(13,1)=2;
v2(13,2)=9;
v2(13,3)=v2(12,3);
v2(14,1)=2;
v2(14,2)=159;
v2(14,3)=(-(params(1)*y(1)^(params(1)-1)*T129));
v2(15,1)=2;
v2(15,2)=123;
v2(15,3)=v2(14,3);
v2(16,1)=2;
v2(16,2)=161;
v2(16,3)=(-(T25*getPowerDeriv(y(9),1-params(1),2)));
v2(17,1)=4;
v2(17,2)=141;
v2(17,3)=(-(y(10)*getPowerDeriv(y(8),(-params(6)),2)));
v2(18,1)=4;
v2(18,2)=161;
v2(18,3)=params(8)*getPowerDeriv(y(9),params(7),2);
v2(19,1)=4;
v2(19,2)=179;
v2(19,3)=(-T117);
v2(20,1)=4;
v2(20,2)=143;
v2(20,3)=v2(19,3);
v2(21,1)=6;
v2(21,2)=1;
v2(21,3)=(-(T26*y(7)*getPowerDeriv(y(1),params(1),2)));
v2(22,1)=6;
v2(22,2)=115;
v2(22,3)=(-(T26*T96));
v2(23,1)=6;
v2(23,2)=7;
v2(23,3)=v2(22,3);
v2(24,1)=6;
v2(24,2)=153;
v2(24,3)=(-(y(7)*T96*T129));
v2(25,1)=6;
v2(25,2)=9;
v2(25,3)=v2(24,3);
v2(26,1)=6;
v2(26,2)=159;
v2(26,3)=(-(T14*T129));
v2(27,1)=6;
v2(27,2)=123;
v2(27,3)=v2(26,3);
v2(28,1)=6;
v2(28,2)=161;
v2(28,3)=(-(y(7)*T14*getPowerDeriv(y(9),1-params(1),2)));
v2(29,1)=8;
v2(29,2)=141;
v2(29,3)=getPowerDeriv(y(8),(-params(6)),2);
v2(30,1)=8;
v2(30,2)=321;
v2(30,3)=(-(params(2)*(1+y(14))*getPowerDeriv(y(17),(-params(6)),2)));
v2(31,1)=8;
v2(31,2)=264;
v2(31,3)=(-(params(2)*getPowerDeriv(y(17),(-params(6)),1)));
v2(32,1)=8;
v2(32,2)=318;
v2(32,3)=v2(31,3);
v2(33,1)=10;
v2(33,2)=141;
v2(33,3)=(-(1/(1-params(6))*getPowerDeriv(y(8),1-params(6),2)));
v2(34,1)=11;
v2(34,2)=161;
v2(34,3)=(-((-params(8))*getPowerDeriv(y(9),1+params(7),2)/(1+params(7))));
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),12,361);
end
if nargout >= 4,
%
% Third order derivatives
%

  g3 = sparse([],[],[],12,6859);
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
