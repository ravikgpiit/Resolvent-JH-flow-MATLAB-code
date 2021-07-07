% Transient Growth
%
% compute the Orr-Sommerfeld matrix for three-
% dimensional Poiseuille,  
% compute the energy weight matrix and determines 
% maximum transient growth in the 
% alpha-Re-plane
% Plot the maximum optimal growth (fig 1) 
% and the least stable eigenvalue (fig 2)
% in the Reynolds-alpha plane
%
%
% INPUT 
%
% beta      = beta  (spanwise wave number)
% N         = total number of modes for normal velocity
% T         = compute maximum growth in time interval [0 T]
%

    clear 
    
    global D0 D1 D2 D4 
    global qb
    
    zi = sqrt(-1);

    %...input data
    N      = input('Enter the number of Chebyshev polynomials: ');
    beta   = input('Enter beta: ');
    Tmax   = input('Enter Tmax: ');
    T      = [0 Tmax];

    %...generate Chebyshev differentiation matrices
    [D0,D1,D2,D4] = ChebMat(N);

    nreso       = 8;
    Re_min      = 6000; Re_max    = 8000;
    alpha_min   = 0.8;  alpha_max = 1.2;
    Re_range    = linspace(Re_min,Re_max,nreso);
    alpha_range = linspace(alpha_min,alpha_max,nreso);
    
    for i=1:nreso
      for j=1:nreso
        alpha = alpha_range(i);
        Re    = Re_range(j);
        [A,B] = JefferyMatrix(N,alpha,beta,Re);

        %...generate energy weight matrix
        k2 = alpha^2 + beta^2;
        M  = EnergyMatrix(N+1,N+1,k2);
        
        %...compute the Orr-Sommerfeld matrix (by inverting B)
        OS = inv(B)*A;
        eOS = eig(OS);
        
        %...compute the optimal
        [flowin,flowot,gg] = Optimal(OS,T,M,k2,1);
        Gmax(i,j) = max(gg(:,2));
        emax(i,j) = max(imag(eOS));
      end
    end
    
    %...graphics
    figure(1) 
    contour(Re_range,alpha_range,Gmax,30);colorbar;
    contour(Re_range,alpha_range,log10(Gmax),30);colorbar;
    figure(2) 
    contour(Re_range,alpha_range,emax,30);colorbar
    hold on
    contour(Re_range,alpha_range,emax,[0 0],'k','LineWidth',3)
