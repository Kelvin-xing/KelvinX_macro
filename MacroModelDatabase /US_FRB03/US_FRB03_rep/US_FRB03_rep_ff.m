
function z=US_FRB03_rep_ff(y)
z=zeros(279,1);
global ex_ ex_det_ it_ recur_

global anton bfi1 bfi2 datet dlpcr1 dlpxp1 dlpxp2 dmptay dmpvar dptr 
global ecd1 ecd2 ecnd1 ecnd2 ecnia1 eg1 eg2 egf1 egf2 egf3 
global egfo1 egs1 egs2 egs3 em1 em2 emn1 emn2 emp1 emp2 
global fcbn1 fcbn2 fcbn3 fcbn4 fnin1 fnin2 fpxr1 fynin1 fynin2 gfdbtn1 
global gfdbtn2 gfdbtn3 gfdefn1 gfdefn10 gfdefn11 gfdefn2 gfdefn3 gfdefn4  ...
gfdefn5 gfdefn6 
global gfdefn7 gfdefn8 gfdefn9 gfintn1 gfptn1 gsdbtn1 gsdbtn2 gsdbtn3  ...
gsector gshock4 
global gsintn1 gsintn2 gsptn1 gssrpn1 gssrpn10 gssrpn11 gssrpn2 gssrpn3  ...
gssrpn4 gssrpn5 
global gssrpn6 gssrpn7 gssrpn8 gssrpn9 iscurve jccan1 jccan2 jccan3 kcd1  ...
kcd2 
global kh1 kh2 ki1 ki2 kpd1 kpd2 kps1 kps2 laginfo lagpi1 
global lagpi2 lagpi3 lagpi4 leadinfo leadpi0 leadpi1 leadpi10 leadpi11  ...
leadpi12 leadpi13 
global leadpi14 leadpi15 leadpi16 leadpi17 leadpi18 leadpi19 leadpi1l  ...
leadpi2 leadpi20 leadpi3 
global leadpi4 leadpi5 leadpi6 leadpi7 leadpi8 leadpi9 leh1 leh2 leh3 leh4 
global ltfpt1 lur1 lurnat1 lwps1 pcrec pgasrec pitarg1 plminr1 poilr1  ...
poilr4 
global profit1 profit2 profit3 profit4 ptr1 pwfix qpxp1 qpxp10 qpxp2 qpxp3 
global qpxp4 qpxp5 qpxp6 qpxp7 qpxp8 qpxp9 qsector qyhibnr1 qyhibnr2  ...
qyhibnr3 
global qynid1 qynid2 qynid3 rcbe1 rg10e1 rg51 rg5e1 rstar1 rstarsh rtb1 
global rtr1 shortlag sprule stockoff taxon tayc0 taye0 tayfi0 tayg0 tayii0 
global taymu0 taymu1 tayp0 tayp1 tayp12 tayp16 tayp2 tayp3 tayp4 tayp8 
global taypl0 taypl1 taypl2 taypl3 taypl4 tayplm1 tayplm2 tayplm3 tayplm4  ...
taypm1 
global taypm10 taypm11 taypm12 taypm13 taypm14 taypm15 taypm16 taypm17  ...
taypm18 taypm19 
global taypm2 taypm20 taypm3 taypm4 taypm5 taypm6 taypm7 taypm8 taypm9  ...
tayr1 
global tayr2 tayr3 tayr4 tayre0 tayrl0 taysp0 tayu0 tayu1 tayul1 tayul2 
global tayul3 tayul4 tayx0 tayx1 tayx2 tayx3 tayx4 tayxm1 tayxm2 tayxm3 
global tayxm4 tayxm5 tayxm6 tayxm7 tayxm8 tfpn1 tfpn2 tfpn3 tfpn4 trfci1 
global trfp1 trfpt1 trfpt2 trsci1 trsib1 trsp1 trspt1 trspt2 trssi1 tryh1 
global tryh2 tryh3 tryh4 tspn1 tspn2 tspn3 tspn4 vpda1 vpda2 wpo1 
global wpo2 wpo3 wpo4 wpo5 wpo6 xb1 xb2 xb3 xbn1 xbn2 
global xbn3 xg1 xg2 xg3 xgaprho xgdp1 xgdp10 xgdp2 xgdp3 xgdp4 
global xgdp5 xgdp6 xgdp7 xgdp8 xgdp9 xgdpn1 xgdpn2 xgdpn3 xgdpn4 xgdpn5 
global xgdpn6 xgdpn7 xgdpn8 xgdpt1 xgdpt2 xgn1 xgn2 xgn3 xgn4 xgn5 
global xgpot1 xgpot2 xgpot3 xgpot4 xgv1 xgv2 xotht1 xotht2 xotht3 xp1 
global xp2 xp3 xp4 xp5 xp6 xp7 xp8 xp9 yh1 yh2 
global yh3 yhinr1 yhinr2 yhinr3 yhinr4 yhl1 yhln1 yhln2 yhln3 yhp1 
global yhp2 yhp3 yhpntfix yhpntn1 yhpntn2 yhpntn3 yhpntn4 yhpntn5 yhpntn6  ...
yhpntn7 
global yhpntn8 yhpntn9 yhptn1 yhptn2 yhptn3 yhtn1 yhtn2 ynicpn1 ynicpn2  ...
ynicpn3 
global ynicpn4 yniln1 yniln2 yniln3 ynin1 ynin2 ynin3 ynin4 ynin5 ynin6 
global ynin7 ypn1 ypn10 ypn2 ypn3 ypn4 ypn5 ypn6 ypn7 ypn8 
global ypn9 
z(1) = y(336) -(y(431));
z(2) = y(334) -((1/4)*(y(394)+y(195)+y(54)+y(2)));
z(3) = y(335) -(y(394));
z(4) = y(380) -(y(470));
z(5) = y(379) -(y(476)*100);
z(6) = y(336) -((1-tayr1)*y(442)+0.755226*y(144)+0.602691/4*(y(335)+y(143)+ ...
y(38)+y(1))+1.17616*y(380)-0.972390*y(185)+ex_(it_-3,22));
z(7) = y(279) -(y(431)-y(222));
z(8) = y(415) -(y(540)+0.0584068790054*(y(542)-y(540))-0.0655950765465*( ...
y(541)-y(540))+0.0325302560675*(1+iscurve*.5)*(y(354)+y(413)-y(389)-y(540))+ ...
0.1440652972670*(y(464)-y(540)));
z(9) = y(283) -(y(110)+0.154173479486*(y(209)-y(110))+0.207837076245*( ...
y(110)-y(15))+datet*y(502)+(1-datet)*y(163)+0.0080074377742*(datet*y(511)+(1 ...
-datet)*y(167))+0.0995119731721*(y(487)-y(261)-(datet*y(502)+(1-datet)* ...
y(163)))-0.0229678540295*(y(261)-y(91))+ex_(it_-3,6));
z(10) = y(424) -(-0.263228678243*0.25*(y(391)+y(193)+y(53)+y(194)));
z(11) = y(286) -(y(112)+0.339759596749*(y(216)+y(209)-y(112))-.257457530254 ...
*(y(112)-y(17))+3.03323716171*(datet*y(503)+(1-datet)*y(164))+3.03323716171* ...
(datet*y(535)+(1-datet)*y(180))+3.03323716171*0.008*(datet*y(512)+(1-datet)* ...
y(168))-0.029831192544*(1+iscurve/2.67)*(datet*y(530)+(1-datet)*y(178))+ ...
1.97908433908*(y(209)-y(63))+4.72167026462*(y(63)-y(210))+1.73447532682*( ...
y(210)-y(64))+2.41236583958*(y(64)-y(211))+ex_(it_-3,8));
z(12) = y(285) -(y(111)+0.0978812349793*(y(209)-y(111))+0.0762268913328*( ...
y(111)-y(16))+3.17577973362*(datet*y(504)+(1-datet)*y(165))+3.17577973362* ...
0.008*(datet*y(513)+(1-datet)*y(169))-0.018364159583*(1+iscurve/3.33)*(datet ...
*y(529)+(1-datet)*y(177))+0.762178056319*(y(209)-y(63))+0.483710408941*( ...
y(63)-y(210))-0.257987551449*(y(210)-y(64))+0.637934053208*(y(64)-y(211))+ ...
ex_(it_-3,7));
z(13) = y(305) -(y(126)+0.112177343577*(y(209)-y(126))+0.30940588957*( ...
y(126)-y(24))+14.1861305056*(datet*y(505)+(1-datet)*y(166))-0.116556799251*( ...
1+iscurve/2.58)*(datet*y(531)+(1-datet)*y(179))+5.33425802988*(y(209)-y(63)) ...
+2.61720944263*(y(63)-y(210))+1.26676734391*(y(210)-y(64))+2.37472643623*( ...
y(64)-y(211))+ex_(it_-3,14));
z(14) = y(338) -(kcd1*y(146)+kcd2*(ecd1*y(286)+ecd2*y(285)));
z(15) = y(339) -(kh1*y(147)+kh2*y(305));
z(16) = y(309) -(y(128)+0.0654821916769*(y(81)+y(76)-y(26))+0.0077433484769 ...
*(y(128)-y(26))+0.254415299022*(y(26)-y(129))+0.782307027592*(datet*(y(538)+ ...
(1+iscurve/13.635)*y(536))+(1-datet)*(y(182)+(1+iscurve/13.635)*y(181)))+ ...
0.151555911273*y(265)+0.0661370611345*y(95)+19.4767730841*(datet*y(500)+(1- ...
datet)*y(162))+ex_(it_-3,16));
z(17) = y(343) -(kpd1*y(151)+kpd2*y(309));
z(18) = y(457) -(-0.05055*(0.8*((1-.35*taxon)*y(437)-(datet*y(524)+(1-datet ...
)*y(174)))+0.2*2*y(430))-0.00897*(0.8*((1-.35*taxon)*y(437)-(datet*y(524)+(1 ...
-datet)*y(174)))+0.2*2*y(430)+taxon*(datet*y(524)+(1-datet)*y(174))));
z(19) = y(311) -(y(130)-0.0355252237239*(y(130)-y(249))+0.145094529937*( ...
y(130)-y(27))+0.0971793341412*(y(27)-y(131))+0.380504419177*y(280)+ ...
0.24901951649*y(107)+0.284404961332*y(13)+0.312474129243*y(108)+ ...
0.159040395765*y(14)-0.350082863562*y(109)+ex_(it_-3,17));
z(20) = y(344) -(kps1*y(152)+kps2*y(311));
z(21) = y(340) -(y(148)+.15032750476*(y(249)-y(148))+.23061958751*(y(148)- ...
y(41))+.118300982384*(y(41)-y(149))+.121631744703*(y(149)-y(42))+datet* ...
y(514)+(1-datet)*y(170)+ex_(it_-3,23));
z(22) = y(342) -(y(150)+.0387389192555*(y(249)-y(150))+.507651578953*( ...
y(150)-y(43))+.122173935168*(y(249)-y(81))+ex_(it_-3,24));
z(23) = y(313) -(y(132)+0.143161888685*((y(133)-y(135))-y(132))- ...
0.367783875925*(y(132)-y(28))+1.39918761877*(y(315)-y(133))-0.1*(y(320)- ...
y(135))-0.05*(y(135)-y(33))+ex_(it_-3,18));
z(24) = y(333) -(.775*y(142)+.225*(y(476)-y(253)));
z(25) = y(308) -(y(127)-0.426416797406*(y(127)-(y(253)+y(205)-y(200)))- ...
0.156189507787*(y(127)-y(25))+4.25725755003*y(333)-0.219943167475*(y(406)- ...
y(411)-(y(200)-y(205)))+ex_(it_-3,15));
z(26) = y(319) -(fnin1*y(134)+fnin2*y(314));
z(27) = y(320) -(0.04*(y(435)-y(522))+0.2*fpxr1*(y(319)-(y(476)+y(393)))+ ...
ex_(it_-3,3));
z(28) = y(271) -(y(320)-y(135)+y(269)-y(276)-(y(411)-y(205)));
z(29) = y(273) -(y(320)-y(135));
z(30) = y(476) -(xgdp1*y(288)+xgdp2*y(305)+xgdp3*y(309)+xgdp4*y(311)+xgdp5* ...
(ki1*y(340)+ki2*y(342))+xgdp6*(ki1*y(148)+ki2*y(150))+xgdp7*y(313)+xgdp8*( ...
em1*y(308)+em2*(emp1*y(468)+emp2*y(480)))+xgdp9*(egf1*y(295)+egf2*y(291)+ ...
egf3*y(293))+xgdp10*(egs1*y(303)+egs2*y(299)+egs3*y(301)));
z(31) = y(486) -(xp1*y(288)+xp2*y(305)+xp3*y(309)+xp4*y(311)+xp5*y(313)+xp6 ...
*y(291)+xp7*y(295)+xp8*y(299)+xp9*y(303));
z(32) = y(465) -(xb1*y(476)+xb2*(xgv1*y(293)+xgv2*y(301))+xb3*y(480));
z(33) = y(469) -(100*(y(468)-y(480)));
z(34) = y(470) -(100*(y(476)-y(478)));
z(35) = y(458) -(vpda1*y(243)+vpda2*(.1*(y(76)+y(245)+y(77)+y(246)+y(78)+ ...
y(247)+y(79)+y(248)+y(80)+y(244))));
z(36) = y(456) -(vpda1*y(241)+vpda2*(-.1*(y(46)+y(188)+y(47)+y(189)+y(48)+ ...
y(190)+y(49)+y(191)+y(50)+y(187))));
z(37) = y(351) -(ltfpt1*y(157));
z(38) = y(349) -(y(156)+0.151503615812*(y(212)-y(156))+0.37950418465*( ...
y(156)-y(45))+0.689067964742*(datet*y(516)+(1-datet)*y(171))+0.310932035258* ...
(y(418)-y(212))-0.118000008522*(y(212)-y(65))+ex_(it_-3,28));
z(39) = y(347) -(y(154)+0.11414*(y(156)-y(154))+0.17067*(y(154)-y(44))+ ...
0.63607*(y(349)-y(156))+ex_(it_-3,26));
z(40) = y(348) -(y(155)+0.250725221905*(-y(155))+0.0366539399533*(y(347)- ...
y(154))+0.0330173044951*y(266)+0.030274916672*y(96)+0.0284267764838*y(3)+ ...
0.0274728839308*y(97)+0.0274132390127*y(4)+0.0282478417297*y(98)+ ...
0.0299766920818*y(5)+ex_(it_-3,27));
z(41) = y(346) -(0.750172970581*y(153)-0.00035090025639*y(250)+ ...
ex_(it_-3,25));
z(42) = y(353) -(lurnat1*y(159));
z(43) = y(480) -(y(350)+xgpot1*y(345)+xgpot2*y(293)+xgpot3*y(301)+xgpot4* ...
y(353));
z(44) = y(478) -(xgdpt1*y(480)+(xgdpt2/12)*(xotht1*(xgv1*y(293)+xgv2*y(301) ...
+xgv1*y(115)+xgv2*y(122))+xotht2*(y(480)+y(255))+xotht3*(y(339)+y(147))+ ...
y(482)+y(257)+y(87)+y(258)+y(88)+y(259)+y(89)+y(260)+y(90)+y(256)));
z(45) = y(477) -(xgdpn1*y(486)+xgdpn2*(y(411)+ki1*y(340)+ki2*y(342))+xgdpn3 ...
*(y(411)+ki1*y(148)+ki2*y(150))+xgdpn4*(em1*(y(406)+y(308))+em2*(y(408)+ ...
y(411)+emp1*y(468)+emp2*y(480)))+xgdpn5*(y(405)+y(293))+xgdpn6*(y(405)+ ...
y(301))+xgdpn7*(y(288)+y(389))+xgdpn8*(trsib1*y(450)+y(288)+y(389)));
z(46) = y(337) -(y(411)+jccan1*y(147)+jccan2*y(151)+jccan3*y(152));
z(47) = y(488) -(y(262)-0.0569590277246*(y(262)-y(217))+0.671295777924*( ...
y(262)-y(92))+0.0829844820244*(y(389)+y(276)+ecd1*y(286)+ecd2*y(285)-(y(192) ...
+ecd1*y(112)+ecd2*y(111)))+ex_(it_-3,52));
z(48) = y(494) -(ynicpn1*y(497)+ynicpn2*(y(405)+yniln1*y(349)+yniln2*y(293) ...
+yniln3*y(301))+ynicpn3*y(496)+ynicpn4*(y(465)+y(411)));
z(49) = y(495) -(y(264)+0.0915371212941*(y(218)-y(264))+0.481011899739*( ...
y(264)-y(94)-(y(205)-y(61)))+datet*y(543)+(1-datet)*y(183)+y(411)-y(205)+ ...
ex_(it_-3,53));
z(50) = y(498) -(ypn1*y(497)+ypn2*y(494)+ypn3*y(495)+ypn4*(y(405)+yniln1* ...
y(349)+yniln2*y(293)+yniln3*y(301))+ypn5*(y(405)+trssi1*y(453)+yniln1*y(349) ...
+yniln2*y(293)+yniln3*y(301))+ypn6*(y(478)+gfptn1*y(326)+y(389))+ypn7*( ...
y(478)+gsptn1*y(332)+y(389))+ypn8*y(325)+ypn9*y(330)+ypn10*y(488));
z(51) = y(492) -(yhpntn1*(y(389)+y(146))+yhpntn2*y(488)+yhpntn3*y(494)+ ...
yhpntn4*(trfci1*y(446)+y(494))+yhpntn5*(trsci1*y(449)+y(494))+yhpntn6*y(495) ...
+yhpntn7*y(321)+yhpntn8*y(327)+(1-yhpntfix)*yhpntn9*(y(389)-y(52)+0.25*( ...
y(276)+y(105)+y(12)+y(106))));
z(52) = y(423) -(qpxp1*(y(422)+y(468))+qpxp2*(y(411)+y(339))+qpxp3*(y(405)+ ...
xgv1*y(293)+xgv2*y(301))+qpxp4*(y(411)+y(480))+qpxp5*(emp1*y(468)+emp2* ...
y(480)+y(408)+y(411))+qpxp6*(y(411)+ki1*y(340)+ki2*y(342))+qpxp7*(y(411)+ki1 ...
*y(148)+ki2*y(150))+qpxp8*(em1*(y(406)+y(308))+em2*(y(408)+y(411)+emp1* ...
y(468)+emp2*y(480)))+qpxp9*(y(405)+y(293))+qpxp10*(y(405)+y(301))-y(486));
z(53) = y(276) -(0.101323198072*y(215)+0.375610047295*y(105)+ ...
0.1906564348288402*y(12)+0.271258122587*dlpxp1*0.0424*(y(381)-y(186)+dlpxp2* ...
y(186))-0.0470874139451*dlpxp1*0.0424*(y(186)-y(46)+dlpxp2*y(46))+datet* ...
y(527)+(1-datet)*y(176)-(1+pwfix/1.58)*0.00307479678461*(datet*y(517)+(1- ...
datet)*y(172))+ex_(it_-3,5));
z(54) = y(419) -(y(350)+1.02040816327*(y(413))-0.02040816327*(y(381)+y(411) ...
));
z(55) = y(405) -(y(199)+y(274)-y(276));
z(56) = y(274) -(0.0298564306096*(y(213)-y(199))+0.231414207389*y(103)+ ...
0.210488499183*y(10)+0.2673366641778989*y(104)+datet*y(525)+(1-datet)*y(175) ...
-(1+pwfix/1.41)*0.0110869464785*(datet*y(519)+(1-datet)*y(173))+ ...
0.0226916552639*(y(404)-y(198)+y(411)-y(205)+y(276)-.25*(y(103)+y(10)+y(104) ...
+y(11)))+ex_(it_-3,4));
z(57) = y(406) -(y(200)+y(411)-y(205)+0.122770433843*(y(214)-y(200))+ ...
0.364224968853*(y(200)-y(205)-(y(59)-y(61)))+0.152572766338*(y(59)-y(61)-( ...
y(201)-y(206)))+0.820102224398*(y(269)-.25*(y(99)+y(6)+y(100)+y(7)))- ...
0.21894291106*(y(271)-.25*(y(101)+y(8)+y(102)+y(9)))+ex_(it_-3,33));
z(58) = y(389) -(y(192)-pcrec*y(192)+0.328778501444*(y(192)-y(51))- ...
0.671221498556*(dlpcr1*(y(450)-y(237)))+ex_(it_-3,31));
z(59) = y(394) -(400*(y(389)-y(192)+y(276)));
z(60) = y(381) -(y(186)-0.132995794168*(y(186)-y(46))-0.0791398065385* ...
y(186)+0.0709795052728*y(202)+0.73643556277*(y(408)-y(202))+ex_(it_-3,30));
z(61) = y(391) -(y(193)-pgasrec*y(193)+0.104506774375*(y(193)-y(53))+ ...
0.284297987947*(y(408)-y(202))+ex_(it_-3,32));
z(62) = y(411) -(y(467)-y(465));
z(63) = y(404) -(plminr1*y(198));
z(64) = y(414) -(y(208)+y(276));
z(65) = y(326) -(0.590828651494*y(138)-0.000220510973853*y(470)+ ...
ex_(it_-3,20));
z(66) = y(439) -(.95*(y(228)+0.33/100*(rtb1*y(443)-rtb1*y(231)))+0.05/100*( ...
0.67*rg51*y(437)+0.33*rtb1*y(443))+ex_(it_-3,42));
z(67) = y(324) -(gfdefn1*y(445)+gfdefn2*(trfci1*y(446)+y(494))+gfdefn3*( ...
y(288)+y(389))+gfdefn4*(y(405)+yniln1*y(349)+yniln2*y(293)+yniln3*y(301))+ ...
gfdefn5*(y(405)+y(293))+gfdefn6*(y(411)+y(295))+gfdefn7*(y(478)+gfptn1* ...
y(326)+y(389))+gfdefn8*(y(393)+y(480))+gfdefn9*y(325)+gfdefn10*(y(393)+ ...
y(480))+gfdefn11*(y(393)+y(480)));
z(68) = y(321) -(gfdbtn1*y(136)+gfdbtn2*y(324)+gfdbtn3*(y(411)+y(291)));
z(69) = y(332) -(0.868665070992*y(141)-0.0000208722728777*y(470)+ ...
ex_(it_-3,21));
z(70) = y(327) -(gsdbtn1*y(139)+gsdbtn2*y(331)+gsdbtn3*(y(411)+y(299)));
z(71) = y(449) -(0.689972923252*y(236)-0.000994932583632*y(470)+ ...
0.000775585371179*y(251)+0.151968645333*(y(446)-y(233))+ex_(it_-3,48));
z(72) = y(453) -(1.01974921383*y(240)-0.258873658191*y(75)- ...
0.0000175095693154*y(470)+ex_(it_-3,51));
z(73) = y(450) -(0.770695168907*y(237)-0.0000569254831607*y(470)+ ...
ex_(it_-3,49));
z(74) = y(446) -(0.567384382892*y(233)+0.00202396979478*y(470)+ ...
0.00521346564087*y(394)+ex_(it_-3,46));
z(75) = y(443) -(y(431)-0.221154289471*(y(431)-y(222))-0.139857827355*( ...
y(231)-y(70))+0.86006333688*(y(231)-y(222))+ex_(it_-3,44));
z(76) = y(441) -((1-0.23624262847)*y(229)+0.23624262847*((1-0.47896)*y(224) ...
+0.47896*y(231)+0.3262*(y(158)-y(159)))+(0.737093726763-0.0870651579706)*( ...
y(435)-y(224))+0.142601547171*(y(229)-y(69))+0.0870651579706*(y(443)-y(231)) ...
+ex_(it_-3,43));
z(77) = y(427) -((1-0.354780575489)*y(219)+0.354780575489*rg51*y(226)+ ...
0.208176943453*(rg51*y(437)-rg51*y(226))+0.0918891467196*(y(219)-y(66))+ ...
ex_(it_-3,36));
z(78) = y(438) -(rg5e1*y(227)+ex_(it_-3,41));
z(79) = y(436) -(rg10e1*y(225)+ex_(it_-3,40));
z(80) = y(429) -(rcbe1*y(220)+ex_(it_-3,37));
z(81) = y(430) -(0.39008463852*4.01*(y(495)-(y(354)+y(413)))-0.259329399006 ...
*4.01*(y(264)-(y(160)+y(207)))+0.869244760486*y(221)+ex_(it_-3,38));
z(82) = y(355) -(lwps1*y(161)+ex_(it_-3,29));
z(83) = y(315) -(qsector*(y(133)-0.0809494989281*y(133)+0.000970693267156* ...
y(470)+0.154305911514*(y(133)-y(29)))+ex_(it_-3,1));
z(84) = y(269) -(qsector*(0.25*(y(99)+y(6)+y(100)+y(7))+0.333258755084* ...
0.125*(y(315)-y(32))+0.0191493302666*(y(135)-y(33))+0.0117184776036*(y(408)- ...
y(202)))+ex_(it_-3,2));
z(85) = y(408) -(qsector*(y(202)-0.255883836726*y(202)+0.394001699237*( ...
y(202)-y(60)))+(1-qsector)*poilr1*y(202)+y(409)-poilr4*y(203));
z(86) = y(409) -(0*y(378)+ex_(it_-3,34));
z(87) = y(448) -(y(235)+trfpt1*(y(136)-y(254))+trfpt2*(y(136)-y(254)-y(137) ...
));
z(88) = y(447) -(trfp1*y(448)+0.531149977846*(y(234)-trfp1*y(235))+ ...
0.278195990818*(y(72)-trfp1*y(73))+0.000461604246182*y(251)+ex_(it_-3,47));
z(89) = y(452) -(y(239)+trspt1*(y(139)-y(254))+trspt2*(y(139)-y(254)-y(140) ...
));
z(90) = y(451) -(0.730539188266*y(238)+1.01604935879*trsp1*y(452)- ...
0.746588547058*trsp1*y(239)-0.0000421560050313*y(251)+0.0190800121832*( ...
y(447)-y(234))+ex_(it_-3,50));
z(91) = y(291) -(gsector*(y(113)-0.377061398783*y(113)-0.132195649068*( ...
y(113)-y(18))+0.0456246902781*(y(18)-y(114))+0.0103283511414*y(470)- ...
0.00861206503463*y(251))+ex_(it_-3,9));
z(92) = y(295) -(gsector*(y(117)-0.260684258075*y(117)-0.144159643547*( ...
y(117)-y(20))-0.220362386808*(y(20)-y(118))-0.0019609669422*y(470)+ ...
0.00326472002188*y(251))+(1-gsector)*egfo1*y(117)+y(297)-gshock4*y(119));
z(93) = y(297) -(0*y(378)+ex_(it_-3,19));
z(94) = y(293) -(gsector*(y(115)-0.111217107631*y(115)+0.329389469448*( ...
y(115)-y(19))-0.0631534124906*(y(19)-y(116))-0.0003987074422*y(470)+ ...
0.000639016169256*y(251))+ex_(it_-3,10));
z(95) = y(299) -(gsector*(y(120)-0.308315775954*y(120)+0.0480696288261*( ...
y(120)-y(21))-0.0952812353267*(y(21)-y(121))+0.00982878490768*y(470)- ...
0.00651765223808*y(251))+ex_(it_-3,11));
z(96) = y(303) -(gsector*(y(124)-0.155542033952*y(124)+0.347413651457*( ...
y(124)-y(23))+0.356797187808*(y(23)-y(125))-0.00063221659304*y(470)+ ...
0.000649024461471*y(251))+ex_(it_-3,13));
z(97) = y(301) -(gsector*(y(122)-0.225318312325*y(122)+0.258184643069*( ...
y(122)-y(22))-0.00328948799931*(y(22)-y(123))-0.000441199608342*y(470)+ ...
0.000764796265076*y(251))+ex_(it_-3,12));
z(98) = y(434) -(0.141521309909*(y(195)-y(197))-0.027351055434*(y(54)- ...
y(197))+0.053734862375*(y(196)-y(197))-0.137787267452*(y(55)-y(197))+ ...
0.992110641812*(y(222)-y(230)-y(204))-0.446416858945*(y(67)-y(230)-y(204))+ ...
0.488158237787*(y(223)-y(230)-y(204))-0.0877896495076*(y(68)-y(230)-y(204))+ ...
y(230)+y(204)+0.295954271439*y(250)-0.0837857095338*y(82)-0.136291571338* ...
y(252)+0.043304922611*y(84)+ex_(it_-3,39));
z(99) = y(431) -(dmptay*y(433)+dmpvar*y(434));
z(100) = y(442) -(rstar1*y(230)+(1-rstar1)*(y(431)-y(394))+rstarsh*(y(444)- ...
(y(410)+rtr1*(y(232)-y(204))+(1-rtr1)*(y(431)-y(394)))));
z(101) = y(444) -(y(410)+rtr1*(y(232)-y(204))+(1-rtr1)*(y(431)-y(394))+ ...
ex_(it_-3,45));
z(102) = y(403) -(pitarg1*y(197)+(1-pitarg1)*y(410));
z(103) = y(410) -(y(204)+dptr*ptr1*(y(195)-y(204))+ex_(it_-3,35));
z(104) = y(506) -(0*y(431));
z(105) = y(507) -(y(431)-(0.141521309909*(y(195)-y(204))-0.027351055434*( ...
y(54)-y(204))+0.053734862375*(y(196)-y(204))-0.137787267452*(y(55)-y(204))+ ...
0.992110641812*(y(222)-y(230)-y(204))-0.446416858945*(y(67)-y(230)-y(204))+ ...
0.488158237787*(y(223)-y(230)-y(204))-0.0877896495076*(y(68)-y(230)-y(204))+ ...
y(230)+y(204)+0.295954271439*y(250)-0.0837857095338*y(82)-0.136291571338* ...
y(252)+0.043304922611*y(84)));
z(106) = y(502) -(0.154173479486*(y(415)-y(209))-0.0307740637878*(y(547)- ...
y(415))+1.03259032482*y(558)-0.199606728022*y(614));
z(107) = y(357) -(y(558));
z(108) = y(511) -(0.0257495000128*y(470)+1.03259032482*y(565)- ...
0.199606728026*y(618));
z(109) = y(361) -(y(565));
z(110) = y(503) -(0.339759596749*(y(415)-y(209))+0.0840097094593*(y(547)- ...
y(415))+0.394727215537*y(559)+0.247262212056*y(615));
z(111) = y(358) -(y(559));
z(112) = y(512) -(0.121637527713*y(470)+0.394727215538*y(566)+ ...
0.247262212056*y(619));
z(113) = y(362) -(y(566));
z(114) = y(535) -(0.339759596749*(y(424)-y(216))+0.0840097094593*(y(551)- ...
y(424))+0.394727215537*y(589)+0.247262212056*y(633));
z(115) = y(374) -(y(589));
z(116) = y(530) -(0.121637527713*(y(427)-y(521))+0.394727215538*y(584)- ...
0.247262212056*y(631));
z(117) = y(372) -(y(584));
z(118) = y(504) -(0.178489847984*(y(415)-y(209))-0.0167789628047*(y(547)- ...
y(415))+0.901003559256*y(560)-0.0940051380753*y(616));
z(119) = y(359) -(y(560));
z(120) = y(513) -(0.0344488224641*y(470)+0.901003559255*y(567)- ...
0.094005138074*y(620));
z(121) = y(363) -(y(567));
z(122) = y(529) -(0.0344488224641*(y(427)-y(521))+0.901003559255*y(583)- ...
0.0940051380741*y(630));
z(123) = y(371) -(y(583));
z(124) = y(505) -(0.112177343577*(y(415)-y(209))-0.0333338808792*(y(547)- ...
y(415))+1.17328397508*y(561)-0.297153416346*y(617));
z(125) = y(360) -(y(561));
z(126) = y(531) -(0.0138953448721*((1-0.3*taxon)*y(441)-y(522))+ ...
1.1732839751*y(585)-0.297153416367*y(632));
z(127) = y(373) -(y(585));
z(128) = y(536) -(0.0654822044827*(y(457)-y(242))-0.000166173491039*(y(553) ...
-y(457))-0.0156805900825*(y(607)-y(553))+0.923403855288*y(590)+ ...
0.236924386007*y(634)-0.239462516617*y(591));
z(129) = y(375) -(y(590));
z(130) = y(538) -(0.0654822044827*(y(465)-y(249))-0.000166173491039*(y(554) ...
-y(465))-0.0156805900825*(y(608)-y(554))+0.923403855288*y(592)+ ...
0.236924386007*y(592)-0.239462516617*y(593));
z(131) = y(376) -(y(592));
z(132) = y(500) -(0.00518186709502*(y(465)-y(249))+0.923416081303*y(557)+ ...
0.236903312103*y(612)-0.239453408793*y(613));
z(133) = y(356) -(y(557));
z(134) = y(524) -(0.0588235294118*400*(y(411)-y(205)+y(276))+(1 ...
-0.0588235294118)*y(578));
z(135) = y(368) -(y(578));
z(136) = y(514) -(0.15032750476*(y(465)-y(249))-0.0326098165361*(y(554)- ...
y(465))-0.0163938700686*(y(608)-y(554))-0.016865129801*(y(555)-y(608))+ ...
1.05868624113*y(568)-0.107870788459*y(621)+0.00313488692336*y(569)- ...
0.112189248584*y(622));
z(137) = y(364) -(y(568));
z(138) = y(516) -(0.151503615812*(y(418)-y(212))-0.0552194044451*(y(548)- ...
y(418))+1.20344055746*y(570)-0.364475818937*y(623));
z(139) = y(365) -(y(570));
z(140) = y(540) -((1-0.930604859102)*y(487)+0.930604859102*y(594));
z(141) = y(541) -((1-0.930604859102)*y(491)+0.930604859102*y(595));
z(142) = y(542) -((1-0.930604859102)*y(493)+0.930604859102*y(596));
z(143) = y(543) -(0.0915371212941*(y(426)-y(218)-(y(411)-y(205)))- ...
0.0422868390037*(y(552)-y(426)-(y(546)-y(411)))+1.36168528288*y(597)- ...
0.461963828509*y(636));
z(144) = y(377) -(y(597));
z(145) = y(525) -(0.0298563916098*(y(419)-y(213)+y(276))-0.00635206734457*( ...
y(549)-y(419)+y(544))-0.00577131624186*(y(604)-y(549)+y(598))- ...
0.00736492314809*(y(605)-y(604)+y(545))+1.17710578781*y(579)-0.0194638751273 ...
*y(627)+0.0533774719806*y(580)-0.246677874267*y(628));
z(146) = y(369) -(y(579));
z(147) = y(527) -(0.101323198072*(y(423)-y(215)+y(276))-0.0361798554854*( ...
y(550)-y(423)+y(544))-0.018181865365*(y(606)-y(550)+y(598))+1.2488011121* ...
y(581)-0.177629510974*y(629)-0.179444250774*y(582));
z(148) = y(370) -(y(581));
z(149) = y(517) -(0.0109705311231*(y(352)-y(353))+1.24880111224*y(571)- ...
0.177629511165*y(624)-0.179444250696*y(572));
z(150) = y(366) -(y(571));
z(151) = y(519) -(0.00106432642675*(y(352)-y(353))+1.17752584206*y(573)- ...
0.0200926352835*y(625)+0.0535007491989*y(574)-0.246582176612*y(626));
z(152) = y(367) -(y(573));
z(153) = y(508) -(0.0588235294118*y(469)+(1-0.0588235294118)*y(562));
z(154) = y(509) -(0.0344827586207*y(469)+(1-0.0344827586207)*y(563));
z(155) = y(510) -(0.0243902439024*y(469)+(1-0.0243902439024)*y(564));
z(156) = y(521) -(0.0588235294118*y(394)+(1-0.0588235294118)*y(575));
z(157) = y(522) -(0.0344827586207*y(394)+(1-0.0344827586207)*y(576));
z(158) = y(523) -(0.0243902439024*y(394)+(1-0.0243902439024)*y(577));
z(159) = y(532) -(0.0588235294118*y(431)+(1-0.0588235294118)*y(586));
z(160) = y(533) -(0.0344827586207*y(431)+(1-0.0344827586207)*y(587));
z(161) = y(534) -(0.0243902439024*y(431)+(1-0.0243902439024)*y(588));
z(162) = y(499) -(y(495)-(y(413))-(y(264)-y(207))+0.98*y(556));
z(163) = y(378) -(y(184));
z(164) = y(416) -(y(63));
z(165) = y(417) -(y(64));
z(166) = y(396) -(y(600));
z(167) = y(397) -(y(601));
z(168) = y(398) -(y(602));
z(169) = y(399) -(y(54));
z(170) = y(400) -(y(55));
z(171) = y(401) -(y(56));
z(172) = y(402) -(y(57));
z(173) = y(395) -(y(58));
z(174) = y(432) -(y(67));
z(175) = y(388) -(y(51)+y(62));
z(176) = y(474) -(y(83));
z(177) = y(310) -(y(26));
z(178) = y(312) -(y(27));
z(179) = y(265) -(profit1*(y(145)-y(39))-(y(205)-y(61))+profit2*(y(263)- ...
y(93))+profit3*(trfci1*y(233)+y(263)-(trfci1*y(71)+y(93)))+profit4*(trsci1* ...
y(236)+y(263)-(trsci1*y(74)+y(93))));
z(180) = y(280) -(y(249)-y(81));
z(181) = y(281) -(y(13));
z(182) = y(282) -(y(14));
z(183) = y(341) -(y(41));
z(184) = y(460) -(y(76));
z(185) = y(461) -(y(77));
z(186) = y(462) -(y(78));
z(187) = y(463) -(y(79));
z(188) = y(459) -(y(80));
z(189) = y(383) -(y(46));
z(190) = y(384) -(y(47));
z(191) = y(385) -(y(48));
z(192) = y(386) -(y(49));
z(193) = y(382) -(y(50));
z(194) = y(266) -(y(154)-y(44));
z(195) = y(267) -(y(3));
z(196) = y(268) -(y(4));
z(197) = y(275) -(y(10));
z(198) = y(278) -(y(12));
z(199) = y(407) -(y(59));
z(200) = y(412) -(y(61));
z(201) = y(270) -(y(6));
z(202) = y(272) -(y(8));
z(203) = y(292) -(y(18));
z(204) = y(294) -(y(19));
z(205) = y(296) -(y(20));
z(206) = y(300) -(y(21));
z(207) = y(302) -(y(22));
z(208) = y(304) -(y(23));
z(209) = y(475) -(y(82));
z(210) = y(316) -(y(29));
z(211) = y(317) -(y(30));
z(212) = y(318) -(y(31));
z(213) = y(392) -(y(53));
z(214) = y(390) -(y(51));
z(215) = y(322) -(y(34)-y(85));
z(216) = y(323) -(y(35));
z(217) = y(328) -(y(36)-y(85));
z(218) = y(329) -(y(37));
z(219) = y(482) -(xotht1*(xgv1*y(19)+xgv2*y(22))+xotht2*y(86)+xotht3*y(40)) ...
;
z(220) = y(483) -(y(87));
z(221) = y(484) -(y(88));
z(222) = y(485) -(y(89));
z(223) = y(481) -(y(90));
z(224) = y(277) -(y(598));
z(225) = y(387) -(y(599)+y(603));
z(226) = y(471) -(y(609));
z(227) = y(472) -(y(610));
z(228) = y(473) -(y(611));
z(229) = y(466) -(y(608));
z(230) = y(420) -(y(604));
z(231) = y(515) -(y(621));
z(232) = y(526) -(y(627));
z(233) = y(520) -(y(625));
z(234) = y(537) -(y(634));
z(235) = y(539) -(y(635));
z(236) = y(501) -(y(612));
z(237) = y(528) -(y(629));
z(238) = y(518) -(y(624));
z(239) = y(413) -(y(479)-y(468));
z(240) = y(491) -(yhp1*y(454)+yhp2*(yhptn1*y(467)+yhptn2*y(489)+yhptn3* ...
y(495))+yhp3*y(492)-y(389));
z(241) = y(330) -(gsintn1*(y(139)+gfintn1*y(440))+gsintn2*y(467));
z(242) = y(445) -(tfpn4*y(447)+tfpn1*y(498)+tfpn2*(y(478)+gfptn1*y(326)+ ...
y(389))+tfpn3*(y(478)+gsptn1*y(332)+y(389)));
z(243) = y(287) -(ecnd1*y(283)+ecnd2*y(146));
z(244) = y(306) -(xgdp5*(ki1*y(340)+ki2*y(342))+xgdp6*(ki1*y(148)+ki2* ...
y(150)));
z(245) = y(314) -(fcbn1*(y(411)+y(313))+fcbn2*(emn1*(y(406)+y(308))+emn2*( ...
y(408)+y(411)+emp1*y(468)+emp2*y(480)))+fcbn3*(fynin1*y(134)+fynin2*(y(413)+ ...
y(480)))+fcbn4*(y(413)+y(480)));
z(246) = y(345) -((leh1*y(347)+y(346)+leh3*y(301)+leh4*y(293))/(1-leh2));
z(247) = y(352) -(lur1*(y(348)-y(345)));
z(248) = y(493) -(y(478)+yhtn1*gfptn1*y(326)+yhtn2*gsptn1*y(332));
z(249) = y(437) -(y(532)-anton*0.62346*y(508)+y(438));
z(250) = y(435) -(y(533)-anton*0.78598*y(509)+y(436));
z(251) = y(428) -(y(534)-anton*1.2078*y(510)+y(429));
z(252) = y(354) -(y(264)-y(207)+y(499)-(1/(1-0.98))*(y(428)-y(523))/400+ ...
y(355));
z(253) = y(290) -(egf1*y(295)+egf2*y(291)+egf3*y(293));
z(254) = y(298) -(egs1*y(303)+egs2*y(299)+egs3*y(301));
z(255) = y(284) -(ecd1*y(286)+ecd2*y(285));
z(256) = y(307) -(em1*y(308)+em2*(emp1*y(468)+emp2*y(480)));
z(257) = y(350) -(.257*y(458)+0.075*y(456));
z(258) = y(418) -(y(468)-y(350));
z(259) = y(425) -(qyhibnr1*y(427)+qyhibnr2*(y(389)+(ecd1*y(286)+ecd2*y(285) ...
))+qyhibnr3*(y(288)+y(389)));
z(260) = y(496) -(y(467));
z(261) = y(426) -(qynid1*y(494)+qynid2*(trfci1*y(446)+y(494))+qynid3*( ...
trsci1*y(449)+y(494)));
z(262) = y(422) -(0.98*(y(405)-y(350))+0.02*(y(381)+y(411)));
z(263) = y(421) -(y(411)-.786550*y(320));
z(264) = y(440) -(y(439));
z(265) = y(289) -(eg1*y(290)+eg2*y(298));
z(266) = y(468) -(xg1*y(465)+xg2*y(339)+xg3*(emp1*y(468)+emp2*y(480)));
z(267) = y(288) -(1.00*y(283)+ecnia1*(ecd1*y(286)+ecd2*y(285)-y(146)));
z(268) = y(467) -(xbn1*y(477)+xbn2*(y(405)+xgv1*y(293)+xgv2*y(301))+xbn3*( ...
y(411)+y(480)));
z(269) = y(479) -(xgn1*y(467)+xgn2*(emp1*y(468)+emp2*y(480)+y(408)+y(411))+ ...
xgn3*(y(411)+y(339))+xgn4*(y(288)+y(389))+xgn5*(trsib1*y(450)+y(288)+y(389)) ...
);
z(270) = y(497) -(ynin1*y(477)+ynin2*(fynin1*y(134)+fynin2*(y(413)+y(480))) ...
+ynin3*y(337)+ynin4*(y(288)+y(389))+ynin5*(trsib1*y(450)+y(288)+y(389))+ ...
ynin6*(y(393)+y(480))+ynin7*(y(393)+y(480)));
z(271) = y(489) -(yhinr1*y(496)+yhinr2*y(325)+yhinr3*y(330)+yhinr4*y(488));
z(272) = y(454) -(tryh1*y(445)+tryh2*y(455)-(tryh3*((y(405)+yniln1*y(349)+ ...
yniln2*y(293)+yniln3*y(301))+yhln3*trssi1*y(453))+tryh4*(yhptn1*y(467)+ ...
yhptn2*y(489)+yhptn3*y(495))));
z(273) = y(490) -(y(405)-y(389)+yniln1*y(349)+yniln2*y(293)+yniln3*y(301)+ ...
yhln3*trssi1*y(453)+yhl1*y(454));
z(274) = y(487) -(yh1*y(490)+yh2*y(493)+yh3*y(491));
z(275) = y(393) -(y(477)-y(476));
z(276) = y(325) -(y(136)+gfintn1*y(439));
z(277) = y(455) -(tspn4*y(451)+tspn1*y(498)+tspn2*(y(478)+gfptn1*y(326)+ ...
y(389))+tspn3*(y(478)+gsptn1*y(332)+y(389)));
z(278) = y(331) -(gssrpn1*y(455)+gssrpn2*(trsci1*y(449)+y(494))+gssrpn3*( ...
trsib1*y(450)+y(288)+y(389))+gssrpn4*(trssi1*y(453)+y(405)+yniln1*y(349)+ ...
yniln2*y(293)+yniln3*y(301))+gssrpn5*(y(393)+y(480))+gssrpn6*(y(393)+y(480)) ...
+gssrpn7*(y(405)+y(301))+gssrpn8*(y(411)+y(303))+gssrpn9*(y(478)+gsptn1* ...
y(332)+y(389))+gssrpn10*y(330)+gssrpn11*(y(393)+y(480)));
z(279) = y(464) -(wpo1*(y(321)-y(389))+wpo2*(y(327)-y(389))+wpo3*(y(319)- ...
y(389))+wpo4*(y(339)+y(411)-y(389))+wpo5*(y(338)+y(411)-y(389))+wpo6*y(480)) ...
;
