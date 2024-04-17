
function z=BRA_SAMBA08_fff(y)
z=zeros(48,1);
global ex_ ex_det_ it_ recur_

global alpha beta byC bystarC coffispol cofintinf0 cofintinfb1 cofintinfb2  ...
cofintinfb3 cofintinfb4 
global cofintinff1 cofintinff2 cofintinff3 cofintinff4 cofintintb1  ...
cofintintb2 cofintintb3 cofintintb4 cofintout cofintoutb1 
global cofintoutb2 cofintoutb3 cofintoutb4 cofintoutf1 cofintoutf2  ...
cofintoutf3 cofintoutf4 del dela dels 
global fiistar gamab gamag gamapi gamar gamas gamay h iok kappa 
global kuu nuu omegabarb omegabarc omegabarn pessi rC rhoa rhoc rhofii 
global rhofiistar rhog rhoi rhomstar rhon rhopi rhopistar rhoq rhor  ...
rhorstar 
global rhosbar rstarC sc sd sg si sigma sm std_r_ std_r_quart 
global sva sx tet vi 
z(1) = y(15) -(y(32)*400);
z(2) = y(13) -((1/4)*(4*y(26)+4*y(26)+4*y(26)+4*y(26))*100);
z(3) = y(14) -(y(26)*400);
z(4) = y(25) -(y(41)*100);
z(5) = y(24) -(y(40)*100);
z(6) = y(8) -(y(10));
z(7) = y(15) -(cofintintb1*y(15)+cofintintb2*y(15)+cofintintb3*y(15)+ ...
cofintintb4*y(15)+cofintinf0*y(14)+cofintinfb1*y(14)+cofintinfb2*y(14)+ ...
cofintinfb3*y(14)+cofintinfb4*y(14)+cofintinff1*y(14)+cofintinff2*y(14)+ ...
cofintinff3*y(14)+cofintinff4*y(14)+cofintout*y(25)+cofintoutb1*y(25)+ ...
cofintoutb2*y(25)+cofintoutb3*y(25)+cofintoutb4*y(25)+cofintoutf1*y(25)+ ...
cofintoutf2*y(25)+cofintoutf3*y(25)+cofintoutf4*y(25)+std_r_*ex_(it_-4,8));
z(8) = y(8) -(coffispol*ex_(it_-4,5));
z(9) = y(11) -(gamag*y(11)+(1-gamag)*(gamas*y(36)-gamab*y(2))+y(45));
z(10) = y(5) -((1/(1+h))*y(5)+(h/(1+h))*y(5)-sigma^(-1)*((1-h)/(1+h))*( ...
y(32)-y(26))+sigma^(-1)*((1-h)/(1+h))*(1-rhoc)*y(42));
z(11) = y(6) -(y(38)+y(22));
z(12) = y(4) -((1-omegabarc)*y(5)+omegabarc*y(6));
z(13) = y(21) -(nuu^(-1)*(y(38)-(sigma/(1-h))*(y(5)-h*y(5))-y(47)));
z(14) = y(22) -(nuu^(-1)*(y(38)-(sigma/(1-h))*(y(6)-h*y(6))-y(47)));
z(15) = y(20) -((1-omegabarn)*y(21)+omegabarn*y(22));
z(16) = y(30) -(y(30)-((y(32)-y(26))-(y(34)+y(7)-y(28))));
z(17) = y(20) -(y(40)-(1-kuu)*y(1)-(kuu+alpha*(1-kuu))*y(38)+alpha*(1-kuu)* ...
y(33)+kuu*y(18));
z(18) = y(16)+y(37) -(y(40)-(1-kuu)*y(1)-(1-alpha*(1-kuu))*y(33)+(1-kuu)*(1 ...
-alpha)*y(38)+kuu*y(18));
z(19) = y(7) -(-pessi*y(3)+vi*y(44)+y(43));
z(20) = y(31) -(beta*(1-del)*y(31)+(1-beta*(1-del))*y(33)-(y(32)-y(26)));
z(21) = y(12) -((1/(dels*(1+beta)))*y(31)+(beta/(1+beta))*y(12)+(1/(1+beta) ...
)*y(12)+((1-rhoi*beta)/(1+beta))*y(46));
z(22) = y(16) -((1-del)*y(16)+iok*y(12));
z(23) = y(39) -(y(19)+kappa*y(30));
z(24) = y(17) -(y(40)-kuu*(y(30)-y(18)));
z(25) = y(33) -(dela*y(37));
z(26) = y(18) -(sd*(alpha*y(33)+(1-alpha)*y(38)-y(1))+(1-sd)*y(30));
z(27) = y(26) -(((1-tet*beta)*(1-omegabarb)*(1-tet)/(tet+omegabarb*(1-tet*( ...
1-beta))))*y(18)+(omegabarb/(tet+omegabarb*(1-tet*(1-beta))))*y(26)+(tet* ...
beta/(tet+omegabarb*(1-tet*(1-beta))))*y(26));
z(28) = y(3) -(fiistar*rstarC*(y(3)+y(23)+bystarC*(y(41)-y(41)+(1/sva)*( ...
y(30)-y(30))-y(28)))+bystarC*(y(7)+y(34)));
z(29) = y(23) -((sx/sva)*y(39)-(sm/sva)*y(17)-((sx-sm)/sva)*y(41)-(sm/sva)* ...
((1-sx)/sva)*y(30));
z(30) = y(36)+y(35) -(-y(11));
z(31) = y(2) -(rC*(y(2)+y(11)-byC*(y(41)-y(41)+y(29)))+byC*y(32));
z(32) = y(9) -(y(41)+(sva/sg)*y(11)-(sm/sva)*y(30));
z(33) = y(40) -(sc*y(4)+si*y(12)+sg*y(9)+sx*y(39));
z(34) = y(41) -((1/sva)*y(40)-(sm/sva)*y(17));
z(35) = y(29) -(y(26)-(sm/sva)*(y(30)-y(30)));
z(36) = y(27) -(rhopi*y(27)+ex_(it_-4,11));
z(37) = y(35) -(rhosbar*y(35)+ex_(it_-4,6));
z(38) = y(42) -(rhoc*y(42)+ex_(it_-4,2));
z(39) = y(47) -(rhon*y(47)+ex_(it_-4,10));
z(40) = y(46) -(rhoi*y(46)+ex_(it_-4,7));
z(41) = y(44) -(rhofiistar*y(44)+ex_(it_-4,4));
z(42) = y(43) -(rhofii*y(43)+ex_(it_-4,3));
z(43) = y(1) -(rhoa*y(1)+ex_(it_-4,1));
z(44) = y(48) -(rhor*y(48)+ex_(it_-4,13));
z(45) = y(45) -(rhog*y(45)+y(10));
z(46) = y(19) -(rhomstar*y(19)+ex_(it_-4,9));
z(47) = y(28) -(rhopistar*y(28)+ex_(it_-4,12));
z(48) = y(34) -(rhorstar*y(34)+ex_(it_-4,14));
if ~isreal(z)
  z = real(z)+imag(z).^2;
end
