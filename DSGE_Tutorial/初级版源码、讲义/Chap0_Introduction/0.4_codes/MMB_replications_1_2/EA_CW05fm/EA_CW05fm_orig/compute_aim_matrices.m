% compute_aim _matrices()
%     This script will compute the G and H matrices.

  g = zeros(35, 70);
  h = zeros(35, 105);

  g(1226) = g(1226) + 1;
  g(1) = g(1) - (delta1*1);
  g(36) = g(36) - (delta2*1);
  g(806) = g(806) - (delta3*1);
  g(2346) = g(2346) - (delta3*(-1.00000000*((rlbar*0.00000000)*1)));
  g(2381) = g(2381) - (sigma_e_d*1);
  g(1262) = g(1262) + 1;
  g(2) = g(2) - 1;
  g(1298) = g(1298) + 1;
  g(1473) = g(1473) - (f0*1);
  g(1333) = g(1333) - (f1*1);
  g(1368) = g(1368) - (f2*1);
  g(1228) = g(1228) - (gamma*(f0*1));
  g(1403) = g(1403) - (gamma*(f1*1));
  g(1438) = g(1438) - (gamma*(f2*1));
  g(2418) = g(2418) - (sigma_e_cw*1);
  h(2593) = h(2593) - (f3*1);
  h(2663) = h(2663) - (gamma*(f3*1));
  g(1334) = g(1334) + 1;
  h(2699) = h(2699) - 1;
  g(1370) = g(1370) + 1;
  h(2560) = h(2560) - 1;
  g(1406) = g(1406) + 1;
  h(2456) = h(2456) - 1;
  g(1442) = g(1442) + 1;
  h(2632) = h(2632) - 1;
  g(1478) = g(1478) + 1;
  g(1303) = g(1303) - (f0*1);
  g(1513) = g(1513) - (f1*1);
  g(1548) = g(1548) - (f2*1);
  g(323) = g(323) - (f3*1);
  g(1514) = g(1514) + 1;
  g(79) = g(79) - 1;
  g(1550) = g(1550) + 1;
  g(290) = g(290) - 1;
  g(1586) = g(1586) + 1;
  g(1306) = g(1306) - ((f0*1)*(((f1+f2)+f3)^-1.00000000));
  g(1516) = g(1516) - ((f1*1)*(((f1+f2)+f3)^-1.00000000));
  g(1551) = g(1551) - ((f2*1)*(((f1+f2)+f3)^-1.00000000));
  g(326) = g(326) - ((f3*1)*(((f1+f2)+f3)^-1.00000000));
  g(1621) = g(1621) - ((-1.00000000*(f2*1))*(((f1+f2)+f3)^-1.00000000));
  g(1621) = g(1621) - ((-1.00000000*(f3*1))*(((f1+f2)+f3)^-1.00000000));
  g(1656) = g(1656) - ((-1.00000000*(f3*1))*(((f1+f2)+f3)^-1.00000000));
  g(1622) = g(1622) + 1;
  g(362) = g(362) - 1;
  g(1658) = g(1658) + 1;
  g(398) = g(398) - 1;
  g(1694) = g(1694) + 1;
  g(1589) = g(1589) - 1;
  g(1624) = g(1624) - 1;
  g(1659) = g(1659) - 1;
  g(434) = g(434) - 1;
  g(1730) = g(1730) + 1;
  g(2325) = g(2325) - 1;
  g(1766) = g(1766) + 1;
  g(1731) = g(1731) - 1;
  g(506) = g(506) - (-1.00000000*1);
  g(1802) = g(1802) + 1;
  g(1732) = g(1732) - (1*(8.00000000^-1.00000000));
  g(1837) = g(1837) - (1*(8.00000000^-1.00000000));
  g(1872) = g(1872) - (1*(8.00000000^-1.00000000));
  g(1907) = g(1907) - (1*(8.00000000^-1.00000000));
  g(1942) = g(1942) - (1*(8.00000000^-1.00000000));
  g(1977) = g(1977) - (1*(8.00000000^-1.00000000));
  g(2012) = g(2012) - (1*(8.00000000^-1.00000000));
  h(3237) = h(3237) - (1*(8.00000000^-1.00000000));
  g(1838) = g(1838) + 1;
  h(2958) = h(2958) - 1;
  g(1874) = g(1874) + 1;
  h(3064) = h(3064) - 1;
  g(1910) = g(1910) + 1;
  h(3100) = h(3100) - 1;
  g(1946) = g(1946) + 1;
  h(3136) = h(3136) - 1;
  g(1982) = g(1982) + 1;
  h(3172) = h(3172) - 1;
  g(2018) = g(2018) + 1;
  h(3208) = h(3208) - 1;
  g(2054) = g(2054) + 1;
  g(1809) = g(1809) - 1;
  g(2089) = g(2089) - (-1.00000000*(1*(2.00000000^-1.00000000)));
  g(2124) = g(2124) - (-1.00000000*(1*(2.00000000^-1.00000000)));
  g(2159) = g(2159) - (-1.00000000*(1*(2.00000000^-1.00000000)));
  g(2194) = g(2194) - (-1.00000000*(1*(2.00000000^-1.00000000)));
  g(2229) = g(2229) - (-1.00000000*(1*(2.00000000^-1.00000000)));
  g(2264) = g(2264) - (-1.00000000*(1*(2.00000000^-1.00000000)));
  g(2299) = g(2299) - (-1.00000000*(1*(2.00000000^-1.00000000)));
  h(3524) = h(3524) - (-1.00000000*(1*(2.00000000^-1.00000000)));
  g(2090) = g(2090) + 1;
  h(2825) = h(2825) - 1;
  g(2126) = g(2126) + 1;
  h(3316) = h(3316) - 1;
  g(2162) = g(2162) + 1;
  h(3352) = h(3352) - 1;
  g(2198) = g(2198) + 1;
  h(3388) = h(3388) - 1;
  g(2234) = g(2234) + 1;
  h(3424) = h(3424) - 1;
  g(2270) = g(2270) + 1;
  h(3460) = h(3460) - 1;
  g(2306) = g(2306) + 1;
  h(3496) = h(3496) - 1;
  g(2342) = g(2342) + 1;
  g(522) = g(522) - (rho*1);
  g(2377) = g(2377) - (((1.00000000*rlbar)*1)*0.00000000);
  g(2377) = g(2377) - ((((-1.00000000*rho)*rlbar)*1)*0.00000000);
  g(1712) = g(1712) - (alpha*1);
  g(2377) = g(2377) - (alpha*(-1.00000000*((pitarget*1)*0.00000000)));
  g(1257) = g(1257) - (beta*1);
  g(2378) = g(2378) + 1;
  g(1153) = g(1153) - (0.00000000*1);
  g(2414) = g(2414) + 1;
  g(2379) = g(2379) - (0.00000000*1);
  g(2450) = g(2450) + 1;
  g(2380) = g(2380) - (0.00000000*1);

  cofg = g;
  cofh = h;
