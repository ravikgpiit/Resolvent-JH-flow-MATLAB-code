function [A,B] = CouetteMatrix(N,alpha,beta,Re);

%
% Function to create Orr-Sommerfeld matrices using Chebyshev 
% pseudospectral discretization for plane Couette flow 
% profile
%
% N   = number of even or odd modes
% alpha     = alpha
% Re       = Reynolds number

    global D0 D1 D2 D4


    zi=sqrt(-1);

    % mean velocity

    k2  = alpha^2 + beta^2;
    Nos = N+1;
    Nsq = N+1;
    vec = [0:N]';
    u   = cos(pi*vec/N);
    du  = ones(length(u),1);

    % set up Orr-Sommerfeld matrix
    B11 = D2 - k2*D0;
    A11 = -(D4 - 2*k2*D2 + (k2*k2)*D0)/(zi*Re);
    A11 = A11 + alpha*(u*ones(1,length(u))).*B11;
    er  = -200*zi;
    A11 = [er*[D0(1,:); D1(1,:)]; A11(3:Nos-2,:); ... 
           er*[D1(Nos,:); D0(Nos,:)]];
    B11 = [D0(1,:); D1(1,:); B11(3:Nos-2,:); ... 
           D1(Nos,:); D0(Nos,:)];

    % set up Squire matrix and (cross-term) coupling matrix
    A21 = beta*(du*ones(1,length(u))).*D0;
    A22 = alpha*(u*ones(1,length(u))).*D0-(D2-k2*D0)/(zi*Re);
    B22 = D0;
    A22 = [er*D0(1,:); A22(2:Nsq-1,:); er*D0(Nsq,:)];
    
    % combine all the blocks 
    A = [A11 zeros(Nos,Nsq); A21 A22];
    B = [B11 zeros(Nos,Nsq); zeros(Nsq,Nos) B22];
