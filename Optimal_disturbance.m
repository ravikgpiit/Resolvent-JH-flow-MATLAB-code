% Optimal Disturbance
%
% compute the Orr-Sommerfeld matrix for three-
% dimensional Poiseuille or Couette flows and 
% compute the energy weight matrix 
%
% compute and displays the optimal initial condition and
% the corresponding flow response: inpout/output
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
    clear 
    
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

iflow  = 1;
    N      = 100;
    Re     = 10000;
    alpha  = 1;
    beta   = 0;
    T   = 100; 
    
    %...generate Chebyshev differentiation matrices
    [D0,D1,D2,D4] = ChebMat(N);

    %...set up Orr-Sommerfeld matrices A and B 
    if (iflow == 1)
      [A,B] = PoiseuilleMatrix(N,alpha,beta,Re);
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
    
    figure(1) 
    plot(real(vin),ycoord,'b',imag(vin),ycoord,'r',abs(vin),ycoord,'k');
    title('optimal initial condition (v)')
    figure(2)
    plot(real(etain),ycoord,'b',imag(etain),ycoord,'r',abs(etain),ycoord,'k');
    title('optimal initial condition (eta)')
    figure(3) 
    plot(real(vout),ycoord,'b',imag(vout),ycoord,'r',abs(vout),ycoord,'k');
    title('optimal output (v)')
    figure(4)
    plot(real(etaout),ycoord,'b',imag(etaout),ycoord,'r',abs(etaout),ycoord,'k');
    title('optimal output (eta)')
    

    
    
