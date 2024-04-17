function [N,L,C,D,alphabeta]=modelschurv2(A,B,G,nx)
%This program is modified by AHNULXY@2012-8-13,SUFE
%This program solve a model of the form
%     [ xt+1 ]  [xt]
%    B[      ]=A[  ]+G[et]
%     [Etyt+1]  [yt]
%using a schur decomposition of the matrices B and A
%nx is the number of expecational variables in Etyt+1
% if plotcode=1, impulse response are plotted
%solution is yt = -N xt - L et
%and xt+1 = C xt + D et
%Q*A*Z==AA,Q*B*Z==B
%QT*A*ZT=AAT,QT*B*ZT=BBT (transformed ordered version of QZ decomposition)
% qz is the built-in command; qzdiv is written by Sim C. 

% Q*A*Z==AA, Q*B*Z==BB
% V and W  whose columns are generalized eigenvectors.
% AA, BB are upper triangular matrix;
[AA,BB,Q,Z,V,W]=qz(A,B);

%Takes U.T. matrices A, B, orthonormal matrices Q,Z, rearranges them
% so that all cases of abs(B(i,i)/A(i,i))>stake are in lower right 
% corner, while preserving U.T. and orthonormal properties and Q'AZ' and
% Q'BZ'.  The columns of v are sorted correspondingly.
[BBT,AAT,QT,ZT,V]=qzdiv(1,BB,AA,Q,Z,V);

alpha=diag(AAT);
beta =diag(BBT);
alphabeta=[alpha beta];
Zp=ZT';
Qp=QT;
[a,b]=size(Z);
N=inv(Zp(a-nx+1:a,a-nx+1:a))*Zp(a-nx+1:a,1:a-nx);
L=inv(Zp(a-nx+1:a,a-nx+1:a))*inv(AAT(a-nx+1:a,a-nx+1:a));
L=L*(Qp(a-nx+1:a,1:a-nx)*G(1:a-nx,1)+Qp(a-nx+1:a,a-nx+1:a)*G(a-nx+1:a,1));
invBBN=inv(B(1:a-nx,1:a-nx)-B(1:a-nx,a-nx+1:a)*N);
C=invBBN*(A(1:a-nx,1:a-nx)-A(1:a-nx,a-nx+1:a)*N);
D=invBBN*(G(1:a-nx,1)-A(1:a-nx,a-nx+1:a)*L);