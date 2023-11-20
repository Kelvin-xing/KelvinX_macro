function [ys_, params, info] = US_PV15rep_steadystate2(ys_, exo_, params)
% Steady state generated by Dynare preprocessor
    info = 0;
    ys_(26)=1;
    ys_(4)=params(22)-1;
    ys_(5)=ys_(4);
    ys_(20)=params(24);
    ys_(31)=params(25)-1;
    params(30)=params(21);
    params(31)=1/params(21);
    ys_(35)=params(30)*(1-params(23))^((-1)/params(31));
    ys_(30)=(1+ys_(31))/ys_(35)*(1-params(21))-1;
    ys_(23)=1+ys_(30)-(1-params(3));
    ys_(37)=(params(30)/ys_(35))^params(31);
    ys_(38)=ys_(35)*params(31)/(params(31)-1);
    ys_(39)=(1-ys_(37)*ys_(38))/(1-ys_(37));
    ys_(36)=(1+ys_(4))/((params(30)/ys_(35))^params(31)+params(31)/(params(31)-1)*(1-params(37))*(1-(params(30)/ys_(35))^params(31))*(ys_(35)^(1-params(31))-params(30)^(1-params(31)))/(ys_(35)^(-params(31))-params(30)^(-params(31)))/ys_(35))-1;
    params(38)=ys_(31)/ys_(36);
    params(39)=params(10)/(params(10)-1);
    ys_(10)=ys_(20)/((1-(1-params(3))*params(1))/((1-params(3))*params(1)*(params(39)-1+ys_(26)*(1+ys_(30))*ys_(35)*ys_(37)*params(28)/100*params(2)/(ys_(23)*(params(31)-1))))*params(7)*(1+ys_(31)*params(29))*(1-params(2))+params(7)*params(3)/(1-params(3)));
    ys_(16)=ys_(10)^(1/(params(10)-1));
    ys_(13)=ys_(16)/params(39);
    ys_(11)=params(3)/(1-params(3))*ys_(10);
    ys_(21)=params(7)*ys_(11);
    ys_(22)=ys_(20)-ys_(21);
    ys_(3)=(1-params(2))*(ys_(13)*(params(2)/ys_(23))^params(2))^(1/(1-params(2)));
    ys_(24)=ys_(22)*ys_(3)*params(2)/(1-params(2))/ys_(23);
    ys_(18)=ys_(24)^params(2)*ys_(22)^(1-params(2))/ys_(10);
    ys_(25)=params(3)*ys_(24);
    ys_(19)=ys_(10)*ys_(18)*ys_(16)^params(10);
    ys_(17)=ys_(19)*(1-params(11))-ys_(25);
    ys_(8)=ys_(16)*ys_(18)*(1-1/params(39));
    ys_(2)=(ys_(17)-ys_(17)*params(12))^(-params(26));
    ys_(41)=ys_(3)/params(40);
    ys_(1)=ys_(2)*ys_(41);
    params(6)=ys_(1)*ys_(20)^(-params(5));
    ys_(32)=ys_(24)*ys_(26)*params(21);
    ys_(33)=ys_(26)*ys_(24)+ys_(21)*params(29)*ys_(3)-ys_(32);
    ys_(9)=params(30)^params(31)/(params(31)-1)*ys_(24)*ys_(26)*(1+ys_(30))*ys_(35)^(1-params(31))/ys_(10);
    params(33)=ys_(38)^(1-(1-params(34))/params(34)*params(35)/((1-params(34))/params(34)*params(35)-1));
    ys_(34)=(1+ys_(30))/(1+ys_(31));
    params(32)=ys_(33)-ys_(10)*ys_(9)*(1-params(3))*(1-params(28)/100);
    ys_(7)=(1-params(3))*params(1)/(1-(1-params(3))*params(1))*(ys_(8)+params(28)/100*ys_(9));
    ys_(40)=params(40);
    ys_(43)=params(39);
    ys_(44)=params(38);
    ys_(14)=0;
    ys_(12)=params(39);
    ys_(6)=1;
    ys_(15)=1;
    ys_(27)=1;
    ys_(28)=0;
    ys_(29)=ys_(24);
    ys_(42)=1;
    ys_(45)=0;
    ys_(46)=0;
    ys_(47)=0;
    ys_(48)=0;
    ys_(49)=0;
    ys_(50)=0;
    ys_(51)=0;
    ys_(52)=0;
    ys_(53)=0;
    ys_(54)=0;
    ys_(55)=0;
    ys_(56)=0;
    ys_(60)=0;
    ys_(61)=0;
    ys_(62)=0;
    ys_(63)=0;
    ys_(64)=0;
    ys_(65)=0;
    ys_(66)=0;
    ys_(57)=0;
    ys_(58)=0;
    ys_(59)=0;
    ys_(67)=0;
    ys_(68)=0;
    ys_(69)=0;
    ys_(70)=0;
    ys_(71)=0;
    ys_(72)=0;
    ys_(73)=0;
    ys_(74)=0;
    ys_(75)=0;
    ys_(76)=0;
    ys_(97)=1;
    params(58)=params(22)-1;
    ys_(80)=params(58);
    ys_(91)=params(24);
    ys_(102)=params(25)-1;
    ys_(106)=params(30)*(1-params(23))^((-1)/params(31));
    ys_(101)=(1-params(21))*(1+ys_(102))/ys_(106)-1;
    ys_(94)=1+ys_(101)-(1-params(3));
    ys_(108)=(params(30)/ys_(106))^params(31);
    ys_(109)=params(31)/(params(31)-1)*ys_(106);
    ys_(110)=(1-ys_(108)*ys_(109))/(1-ys_(108));
    ys_(107)=(1+params(58))/((params(30)/ys_(106))^params(31)+params(31)/(params(31)-1)*(1-params(37))*(1-(params(30)/ys_(106))^params(31))*(ys_(106)^(1-params(31))-params(30)^(1-params(31)))/(ys_(106)^(-params(31))-params(30)^(-params(31)))/ys_(106))-1;
    ys_(84)=ys_(91)/(params(7)*params(3)/(1-params(3))+(1-params(2))*params(7)*(1-(1-params(3))*params(1))/((1-params(3))*params(1)*(params(39)-1+params(2)*ys_(97)*(1+ys_(101))*ys_(106)*params(28)/100*ys_(108)/((params(31)-1)*ys_(94))))*(1+params(29)*ys_(102)));
    ys_(87)=ys_(84)^(1/(params(10)-1));
    ys_(86)=ys_(87)/params(39);
    ys_(85)=params(3)/(1-params(3))*ys_(84);
    ys_(92)=params(7)*ys_(85);
    ys_(93)=ys_(91)-ys_(92);
    ys_(79)=(1-params(2))*(ys_(86)*(params(2)/ys_(94))^params(2))^(1/(1-params(2)));
    ys_(95)=ys_(93)*params(2)/(1-params(2))*ys_(79)/ys_(94);
    ys_(89)=ys_(95)^params(2)*ys_(93)^(1-params(2))/ys_(84);
    ys_(96)=params(3)*ys_(95);
    ys_(90)=ys_(84)*ys_(89)*ys_(87)^params(10);
    ys_(88)=(1-params(11))*ys_(90)-ys_(96);
    ys_(82)=(1-1/params(39))*ys_(87)*ys_(89);
    ys_(78)=(ys_(88)-params(12)*ys_(88))^(-params(26));
    ys_(112)=ys_(79)/params(40);
    ys_(77)=ys_(78)*ys_(112);
    ys_(103)=ys_(95)*params(21)*ys_(97);
    ys_(104)=ys_(97)*ys_(95)+ys_(92)*params(29)*ys_(79)-ys_(103);
    ys_(83)=params(30)^params(31)/(params(31)-1)*ys_(95)*ys_(97)*(1+ys_(101))*ys_(106)^(1-params(31))/ys_(84);
    ys_(105)=(1+ys_(101))/(1+ys_(102));
    ys_(81)=(1-params(3))*params(1)/(1-(1-params(3))*params(1))*(ys_(82)+params(28)/100*ys_(83));
    ys_(111)=params(40);
    params(61)=params(39);
    ys_(114)=params(38);
    params(57)=0;
    params(56)=params(39);
    params(59)=1;
    params(60)=1;
    ys_(98)=1;
    ys_(99)=0;
    ys_(100)=ys_(95);
    ys_(113)=1;
    % Auxiliary equations
ys_(116)=exo_(7);
ys_(117)=exo_(8);
    check_=0;
end
