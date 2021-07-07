function [flowin,flowot,gg] = Optimal(OS,T,M,k2,iflag);
%
% This function computes the initial flow structure which 
% achieves the maximum transient energy growth
%
% INPUT
% OS      = 3D Orr-Sommerfeld operator 
% T       = time 
% M       = energy weight matrix
% k2     = alpha^2+beta^2
% iflag   = flag
%           iflag = 1:  compute the maximum growth and 
%                       initial condition in time 
%                       interval [0,T]
%           iflag = 2:  compute the initial disturbance 
%                       yielding maximum growth at time T
% OUTPUT 
% flowin  = coefficients of optimal disturbance
%           flowin(1:Nos)         = normal velocity 
%                                   coefficients
%           flowin(Nos+1:Nos+Nsq) = normal vorticity 
%                                   coefficients
% flowot  = coefficients of field at optimal time
%           flowot(1:Nos)         = normal velocity 
%                                   coefficients
%           flowot(Nos+1:Nos+Nsq) = normal vorticity 
%                                   coefficients
% gg      = energy growth curve (t,G(t)) 

    global qb

    % Phase 1: Compute eigenvalues and eigenfunctions of 
    % Orr-Sommerfeld matrix and sort in order of descending 
    % imaginary part. The function nlize normalizes the 
    % eigenfunctions with respect to the weight matrix M.

    [xs,es] = OrderedEig(OS);
    xs      = ENormalize(xs,M);

    % Phase 2: Choose which modes are to be used in optimal 
    % calculation. Modes with imaginary part > 1 are neglected. 
    % Modes with imaginary part < imin are neglected as well.
    
    ishift = 1; 
    imin   = -1.5;

    while imag(es(ishift))>1, 
      ishift = ishift+1; 
    end

    [n1,n2] = ExtractEig(es,imin);

    cols    = (ishift:n2);
    xu      = xs(:,cols);
    eu      = es(cols); 
    ncols   = length(cols);
    fprintf('Number of modes used: %1.0f \n',ncols); 

    % Phase 3: Compute the reduced Orr-Sommerfeld operator

    [qb,invF] = TMatrix(M,xu,eu);
    
    % Phase 4: Compute the time for the maximum growth using 
    % the built-in Matlab routine 'fminbnd'

    if (iflag == 1)
      gcheck = NormMExp(1/100); 
      gcheck = gcheck^2;
      if (gcheck < 1)
        tformax = 0;
        mgrowth = 1;
      else
        ts = T(1); 
        tf = T(2); 
        options = [0 1e-3 1e-3];
        tformax = OptFunc(ts,tf,options);
        mgrowth = NormMExp(tformax);
        mgrowth = mgrowth^2;
      end
      fprintf('Time for maximum growth:  %e \n',tformax);
    else
      tformax = T;
    end 

    % Phase 5: Compute the initial condition that yields the 
    % maximum growth. This is obtained by
    % (1) computing the matrix exponential evaluated at the 
    %     optimal time;
    % (2) computing the SVD of the matrix exponential 
    %     exp(-i*A*t)=USV. 
    % The initial condition that yields the maximum growth is 
    % the first column of V. To convert the initial condition 
    % to a vector of coefficients in the eigenfunction basis 
    % multiply by the matrix of eigenfunctions and inv(F) 
    
    evol    = expm(tformax*qb);
    [U,S,V] = svd(evol);
    mgrowth = S(1,1)^2;
    fprintf('Maximum growth in energy:  %e \n',mgrowth);
    
    flowin = sqrt(2*k2)*xu*invF*V(:,1);
    flowot = sqrt(2*k2)*xu*invF*U(:,1);

    if (iflag==1)
        ts = 0; 
        tf = T(2);
    else
        ts=0;
        tf=T(1);
    end
    
    
     nit=100;
     for i = 1:nit
       tid     = ts + (tf-ts)/(nit-1)*(i-1);
       gg(i,1) = tid;
       gg(i,2) = norm(expm(tid*qb))^2;
     end 

%        nit=100;
%      gg(:,1) = logspace(-3,(log10(tf)),nit);
%     for i = 1:nit
%       gg(i,2) = norm(expm(gg(i,1)*qb))^2;
%     end