%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chebyshev polynomials
clear all; close all
N=120;
[D0,D1,D2,D4] = ChebMat(N);
vec = [0:N]';
y   = cos(pi*vec/N);
plot(y,'*');

close all
figure(1)
plot(y,D0(:,1),'k*-')
hold on
plot(y,D0(:,2),'r<-')
plot(y,D0(:,3),'g>-')
plot(y,D0(:,4),'ms-')
figure(2)
plot(y,D1(:,1),'k*-')
hold on
plot(y,D1(:,2),'r<-')
plot(y,D1(:,3),'g>-')
plot(y,D1(:,4),'ms-')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example of interpolation error for a sinusoidal function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 f=sin(pi*y);
 y_ex=-1:.01:1;
 f_ex=sin(pi*y_ex)
% 
% %D0*C=f
 C=inv(D0)*f;
 figure(4)
 plot(y,D0*C,'*-k')
 hold on
 plot(y_ex,f_ex,'r-')
 
 plot(y,D1*C,'*-k')
 plot(y_ex,cos(pi*y_ex)*pi,'r-')

 plot(y,D2*C,'*-b');plot(y_ex,-sin(pi*y_ex)*pi^2,'g-')




 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary value problem u''=1 with u(-1)=u(1)=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
er  = 1;
A = [er*D0(1,:); D2(2:N,:); er*D0(N+1,:)];

f=ones(N+1,1);
f=sin(pi*y);
f(1)=0;f(end)=0;%or 2


uh=inv(A)*f;
u=D0*uh;
figure(3)
plot(y,u,'k')
grid on
hold on
plot(y,D2*uh,'r')
plot(y,D1*uh,'g')
%aa=input('continue?')
%%
%close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eigenvalue problem u''= omega u with u(-1)=u(1)=0
% (D2- omega D0) uh=0
% solution eigenvalues= omega=(N*pi/2)^2
% eiegnsolution cos(pi/2*N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
er  = -2000*i;

A = [er*D0(1,:); D2(2:N,:); er*D0(N+1,:)];
B = -D0;
AB=inv(B)*A;
[S,L]=eig(AB);

Le=abs(diag(L));
[Lee,Iee]=sort(Le);

%Plot the eigenvalues against the analytical results (N*pi/2)^2
figure(4)
semilogy(Lee,'*')
hold on
semilogy(([1:50]*(pi/2)).^2,'ro')
axis([1 30 .1 Lee(50)])

%Plot the eigenfunctions
figure(5)
u=D0*S(:,Iee(1));
plot(y,u)
hold on
u=D0*S(:,Iee(2));
plot(y,u,'r')
u=D0*S(:,Iee(3));
plot(y,u,'g')
u=D0*S(:,Iee(4));
plot(y,u,'c')
u=D0*S(:,Iee(5));
plot(y,u,'m')
u=D0*S(:,Iee(6));
plot(y,u,'k')