%=======================================
%
%  compute the resolvent norm for real and complex 
%  frequency omega 
%
%=======================================
% INPUT 
%
% Re        = Reynolds number
% alp       = alpha (streamwise wave number)
% beta      = beta  (spanwise wave number)
% iflow     = type of flow  (Poiseuille=1, Couette=2)  
% N         = total number of modes for normal velocity

    clear all; close all; clc;
    
    global D0 D1 D2 D4 

    zi=sqrt(-1);
    % input data

%     iflow  = input('Poiseuille (1) or Couette flow (2) ');
%     N      = input('Enter the number of Chebyshev polynomials: ');
%     Re     = input('Enter the Reynolds number: ');
%     alpha  = input('Enter alpha: ');
%     beta   = input('Enter beta: ');

    iflow  = 1;
    N      = 100;
    Re     = 250;
    alpha  = 0;
    beta   = 1;
    % generate Chebyshev differentiation matrices
    [D0,D1,D2,D4] = ChebMat(N);

    % set up Orr-Sommerfeld matrices A and B 

    if (iflow == 1)
      [A,B] = JefferyMatrix(N,alpha,beta,Re);
%       [A,B] = PoiseuilleMatrix(N,alpha,beta,Re);
      
    else
      [A,B] = CouetteMatrix(N,alpha,beta,Re);
    end

    % generate energy weight matrix
    k2 = alpha^2 + beta^2;
    M  = EnergyMatrix(N+1,N+1,k2);

    % compute the Orr-Sommerfeld matrix (by inverting B)
    OS = inv(B)*A;
    
    
    % compute the optimal
    eee = eig(OS); 
    %figure(1); 
    %plot(real(eee),imag(eee),'o')
    %axis([-2 2 -2 2]); 
    %axis image; 
  
    [F,e,invF] = GetMatrixParts(OS,M,k2);
    nreso = 50;
    for i=1:nreso
      for j=1:nreso
        zr = -0.5 + 2*(i-1)/(nreso-1);
        zi = -1 + 1.5*(j-1)/(nreso-1);
        zz = zr + sqrt(-1)*zi;
        dd = diag(1./(e-zz));
        Reso(i,j) = log(norm(F*dd*invF));
      end
    end
    for i=1:nreso
        zr = -0.5 + 2*(i-1)/(nreso-1);
        zz = zr;
        dd = diag(1./(e-zz));
        Reso_r(i) = (norm(F*dd*invF));
    end
    
%     figure(1);subplot(1,1,1,'Fontsize',14)
%     semilogy(linspace(-0.5,1.5,nreso),Reso_r)
    plot(linspace(-0.5,1.5,nreso),Reso_r,'b','LineWidth',2)
    subplot(1,1,1,'Fontsize',12)
    title('Resolvent norm')
%     legend('For converging angle \alpha = -0.005^\circ, Re = 11000, k_{x} = 0, k_{z} = 2');
    legend('For diverging angle \alpha = 1^\circ, Re = 250, k_{x} = 0, k_{z} = 1');
    ylabel('Resolvent Norm');xlabel('\omega Frequency')
%     grid on
    figure(2);
    subplot(1,1,1,'Fontsize',12)
    contour(linspace(-0.5,1.5,nreso),linspace(-1,0.5,nreso),Reso')
    hold on; 
    plot([-0.5 1.5],[0 0],'r')
    plot(real(e),imag(e),'*k');
    title('Pseudospectrum')
    hold off
    
%     figure(1)
%     hold on
%     num1 = xlsread('validation33'); 
%     om = num1(:,1); RR = num1(:,2);
%     plot(om, RR,'*r');
%     legend('Resolvant norm for JH flow with divering angle \alpha = 0^\circ for Re = 1000, k_{x} = 1, k_{z}', 'Resolvant norm for plane poiseuille flow for Re = 1000, k_{x} = 1, k_{z}')

    