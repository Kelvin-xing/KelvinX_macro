
function z=BRA_SAMBA08_ff(y)
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
z(1) = y(53) -(y(70)*400);
z(2) = y(51) -((1/4)*(4*y(64)+4*y(24)+4*y(11)+4*y(7))*100);
z(3) = y(52) -(y(64)*400);
z(4) = y(63) -(y(79)*100);
z(5) = y(62) -(y(78)*100);
z(6) = y(46) -(y(48));
z(7) = y(53) -(cofintintb1*y(20)+cofintintb2*y(9)+cofintintb3*y(5)+ ...
cofintintb4*y(2)+cofintinf0*y(52)+cofintinfb1*y(19)+cofintinfb2*y(8)+ ...
cofintinfb3*y(4)+cofintinfb4*y(1)+cofintinff1*y(89)+cofintinff2*y(96)+ ...
cofintinff3*y(98)+cofintinff4*y(100)+cofintout*y(63)+cofintoutb1*y(23)+ ...
cofintoutb2*y(10)+cofintoutb3*y(6)+cofintoutb4*y(3)+cofintoutf1*y(90)+ ...
cofintoutf2*y(97)+cofintoutf3*y(99)+cofintoutf4*y(101)+std_r_*ex_(it_-4,8));
z(8) = y(46) -(coffispol*ex_(it_-4,5));
z(9) = y(49) -(gamag*y(17)+(1-gamag)*(gamas*y(30)-gamab*y(13))+y(83));
z(10) = y(43) -((1/(1+h))*y(87)+(h/(1+h))*y(15)-sigma^(-1)*((1-h)/(1+h))*( ...
y(70)-y(91))+sigma^(-1)*((1-h)/(1+h))*(1-rhoc)*y(80));
z(11) = y(44) -(y(76)+y(60));
z(12) = y(42) -((1-omegabarc)*y(43)+omegabarc*y(44));
z(13) = y(59) -(nuu^(-1)*(y(76)-(sigma/(1-h))*(y(43)-h*y(15))-y(85)));
z(14) = y(60) -(nuu^(-1)*(y(76)-(sigma/(1-h))*(y(44)-h*y(16))-y(85)));
z(15) = y(58) -((1-omegabarn)*y(59)+omegabarn*y(60));
z(16) = y(68) -(y(93)-((y(70)-y(91))-(y(72)+y(45)-y(92))));
z(17) = y(58) -(y(78)-(1-kuu)*y(39)-(kuu+alpha*(1-kuu))*y(76)+alpha*(1-kuu) ...
*y(71)+kuu*y(56));
z(18) = y(21)+y(75) -(y(78)-(1-kuu)*y(39)-(1-alpha*(1-kuu))*y(71)+(1-kuu)*( ...
1-alpha)*y(76)+kuu*y(56));
z(19) = y(45) -(-pessi*y(41)+vi*y(82)+y(81));
z(20) = y(69) -(beta*(1-del)*y(94)+(1-beta*(1-del))*y(95)-(y(70)-y(91)));
z(21) = y(50) -((1/(dels*(1+beta)))*y(69)+(beta/(1+beta))*y(88)+(1/(1+beta) ...
)*y(18)+((1-rhoi*beta)/(1+beta))*y(84));
z(22) = y(54) -((1-del)*y(21)+iok*y(50));
z(23) = y(77) -(y(57)+kappa*y(68));
z(24) = y(55) -(y(78)-kuu*(y(68)-y(56)));
z(25) = y(71) -(dela*y(75));
z(26) = y(56) -(sd*(alpha*y(71)+(1-alpha)*y(76)-y(39))+(1-sd)*y(68));
z(27) = y(64) -(((1-tet*beta)*(1-omegabarb)*(1-tet)/(tet+omegabarb*(1-tet*( ...
1-beta))))*y(56)+(omegabarb/(tet+omegabarb*(1-tet*(1-beta))))*y(24)+(tet* ...
beta/(tet+omegabarb*(1-tet*(1-beta))))*y(91));
z(28) = y(41) -(fiistar*rstarC*(y(14)+y(61)+bystarC*(y(31)-y(79)+(1/sva)*( ...
y(68)-y(27))-y(66)))+bystarC*(y(45)+y(72)));
z(29) = y(61) -((sx/sva)*y(77)-(sm/sva)*y(55)-((sx-sm)/sva)*y(79)-(sm/sva)* ...
((1-sx)/sva)*y(68));
z(30) = y(74)+y(73) -(-y(49));
z(31) = y(40) -(rC*(y(13)+y(49)-byC*(y(79)-y(31)+y(67)))+byC*y(70));
z(32) = y(47) -(y(79)+(sva/sg)*y(49)-(sm/sva)*y(68));
z(33) = y(78) -(sc*y(42)+si*y(50)+sg*y(47)+sx*y(77));
z(34) = y(79) -((1/sva)*y(78)-(sm/sva)*y(55));
z(35) = y(67) -(y(64)-(sm/sva)*(y(68)-y(27)));
z(36) = y(65) -(rhopi*y(25)+ex_(it_-4,11));
z(37) = y(73) -(rhosbar*y(29)+ex_(it_-4,6));
z(38) = y(80) -(rhoc*y(32)+ex_(it_-4,2));
z(39) = y(85) -(rhon*y(37)+ex_(it_-4,10));
z(40) = y(84) -(rhoi*y(36)+ex_(it_-4,7));
z(41) = y(82) -(rhofiistar*y(34)+ex_(it_-4,4));
z(42) = y(81) -(rhofii*y(33)+ex_(it_-4,3));
z(43) = y(39) -(rhoa*y(12)+ex_(it_-4,1));
z(44) = y(86) -(rhor*y(38)+ex_(it_-4,13));
z(45) = y(83) -(rhog*y(35)+y(48));
z(46) = y(57) -(rhomstar*y(22)+ex_(it_-4,9));
z(47) = y(66) -(rhopistar*y(26)+ex_(it_-4,12));
z(48) = y(72) -(rhorstar*y(28)+ex_(it_-4,14));
