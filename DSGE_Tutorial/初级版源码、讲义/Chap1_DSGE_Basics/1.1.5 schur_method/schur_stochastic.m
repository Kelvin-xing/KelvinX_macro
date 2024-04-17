B=[12.6695 0 -1.2353 0 0;
0 1 0 0 0;
0 -1 .36 0 0;
0 0 1 0 0;
0 0 0 1 -.03475];

A=[12.353 0 0 -0.9186 0;
0 .95 0 0 0;
.36 0 0 -0.64 0;
1 0 0 0 1;
0 0 0 1 0];
G=[0 1 0 0 0 ]';

%the number forward-looking variables
%in McCandless(2008), the author put nx =2 who believe output is a
%predetermined variable which I think need further discussions.
nx=2;
% A,B is defined above
[N,L,C,D,alphabeta]=modelschurv2(A,B,G,nx);