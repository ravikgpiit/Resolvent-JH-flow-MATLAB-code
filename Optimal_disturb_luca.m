% Optimal Disturbance
%
% compute the Orr-Sommerfeld matrix for three-
% dimensional Poiseuille or Couette flows and 
% compute the energy weight matrix 
%
% INPUT 
%
% Re        = Reynolds number
% alpha     = alpha (streamwise wave number)
% beta      = beta  (spanwise wave number)
% iflow     = type of flow  (Poiseuille=1, Couette=2)  
% N         = total number of modes for normal velocity
% T         = time of optimal growth
%
    clear all; close all; clc;
    
    global D0 D1 D2 D4 
    global qb
    
    zi = sqrt(-1);
    %...input data
%     iflow  = input('Poiseuille (1) or Couette flow (2) ');
%     N      = input('Enter the number of Chebyshev polynomials: ');
%     Re     = input('Enter the Reynolds number: ');
%     alpha  = input('Enter alpha: ');
%     beta   = input('Enter beta: ');
%     T      = input('Enter T: '); 

    iflow=1;
    N=100;
    Re=250;
    alpha=0;
    beta=1;
                T=100;

    %...generate Chebyshev differentiation matrices
    [D0,D1,D2,D4] = ChebMat(N);

    %...set up Orr-Sommerfeld matrices A and B 
    if (iflow == 1)
      [A,B] = JefferyMatrix(N,alpha,beta,Re);
%       [A,B] = PoiseuilleMatrix(N,alpha,beta,Re);
    else
      [A,B] = CouetteMatrix(N,alpha,beta,Re);
    end

    %...generate energy weight matrix
    k2 = alpha^2 + beta^2;
    M  = EnergyMatrix(N+1,N+1,k2);

    %...compute the Orr-Sommerfeld matrix (by inverting B)
    OS = inv(B)*A;

    %...determine optimal initial condition and optimal output
    [flowin,flowot,gg] = Optimal(OS,T,M,k2,2);

    %...visualize the optimal perturbation 
    vin    = D0*flowin(1:N+1);
    etain  = D0*flowin(N+2:2*(N+1));
    vout   = D0*flowot(1:N+1); 
    etaout = D0*flowot(N+2:2*(N+1)); 
    ycoord = D0(:,2); 
    
    
    %u=i/k2*(alpha*Dv/Dy-beta*eta)
    uin=i/k2*(alpha*D1*flowin(1:N+1)-beta*etain);
    Duin=i/k2*(alpha*D2*flowin(1:N+1)-beta*D1*flowin(N+2:2*(N+1)));
    %w=i/k2*(beta*Dv/Dy+alpha*eta)
    win=i/k2*(beta*D1*flowin(1:N+1)+alpha*etain);
    %u=i/k2*(alpha*Dv/Dy-beta*eta)
    uout=i/k2*(alpha*D1*flowot(1:N+1)-beta*etaout);
    %w=i/k2*(beta*Dv/Dy+alpha*eta)
    wout=i/k2*(beta*D1*flowot(1:N+1)+alpha*etaout);
    
    
    nzz=80;ph=i*(pi/4+1.1);
    nxx = 80;
    z=linspace(0,pi,nzz);
    x = (linspace(0,pi,nxx))';
    for j=1:nzz
        fvv(:,j)=real(vin.*exp(i*2*z(j)+ph));
        fww(:,j)=real(win.*exp(i*2*z(j)+ph));
        fuu(:,j)=real(uout.*exp(i*2*z(j)+ph));
    end
    figure(10);subplot(1,1,1,'FontSize',16);
 quiver(z(1:2:end-1),ycoord(1:2:end),fww(1:2:end,1:2:end-1),fvv(1:2:end,1:2:end-1),'k');
    axis equal; axis([0 pi -1 1])
    
    
    [XX,YY] = meshgrid(x,ycoord);
    figure(101);
    contourf(XX,YY,fuu,50)
    shading flat
    
    figure(102);
    contourf(XX,YY,fvv,50)
    shading flat
   
    figure(103);
    contourf(XX,YY,fww,50)
    shading flat

   
    figure(11);subplot(1,1,1,'FontSize',16);
    contour3(z(1:2:end-1),ycoord(1:2:end),fuu(1:2:end,1:2:end-1),50)
%    surf(ycoord(1:2:end),z(1:2:end-1),fuu(1:2:end,1:2:end-1)')
[yy,zz]=meshgrid(ycoord(1:2:end),z(1:2:end-1));
    surf(zz,fuu(1:2:end,1:2:end-1)',yy);shading interp

    colormap('jet')
axis([0 pi -2 2])
    view(-70,40)
%    view(0,0)
   % break

    
    
    figure(1) 
    plot(real(vin),ycoord,'b',imag(vin),ycoord,'r',abs(vin),ycoord,'k');
    title('optimal initial condition (v)')
    figure(2)
    plot(real(win),ycoord,'b',imag(win),ycoord,'r',abs(win),ycoord,'k');
    title('optimal initial condition (w)')
    figure(3)
    plot(real(uin),ycoord,'b',imag(uin),ycoord,'r',abs(uin),ycoord,'k');
    title('optimal initial condition (u)')
    figure(4) 
    plot(real(vout),ycoord,'b',imag(vout),ycoord,'r',abs(vout),ycoord,'k');
    title('optimal output (v)')
    figure(5)
    plot(real(wout),ycoord,'b',imag(wout),ycoord,'r',abs(wout),ycoord,'k');
    title('optimal output (w)')
    figure(6)
    plot(real(uout),ycoord,'b',imag(uout),ycoord,'r',abs(uout),ycoord,'k');
    title('optimal output (u)')
    
    
