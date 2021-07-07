 function [A,B] = JefferyMatrix(N,alpha,beta,Re)
% clc; clear all; close all;
% create Orr-Sommerfeld matrices using Chebyshev 
% pseudospectral discretization for plane Poiseuille flow 
% profile 
%
% N  = 100; %number of modes
% alpha     = 1; %alpha
% beta    = 0; %beta 
% Re       = 1000;% Reynolds number

    global D0 D1 D2 D4
[D0,D1,D2,D4] = ChebMat(N);
    zi=sqrt(-1);

    % mean velocity
  
    k2 = alpha^2 + beta^2;
    Nos = N+1;
    Nsq = N+1;
    vec = [0:N]';
    
    [ U_final1,U_final2,U_final3 ] = Jeffrey_Hammel_Soln( );
%     u   = (ones(length(vec),1)-cos(pi*vec/N).^2);
%     du  = -2*cos(pi*vec/N);
    
    % set up Orr-Sommerfeld matrix
    B11 = D2 - k2*D0;
    A11 = -(D4 -2*k2*D2 + (k2*k2)*D0)/(zi*Re);
    A11 = A11 + alpha*(diag(U_final1))*B11 - alpha*diag(U_final3)*D0;
%       A11 = A11 + alpha*(u*ones(1,length(u))).*B11 + alpha*2*D0;
    er  = -200*zi;
    A11 = [er*D0(1,:); er*D1(1,:); A11(3:Nos-2,:); ... 
           er*D1(Nos,:); er*D0(Nos,:) ];
    B11 = [D0(1,:); D1(1,:); B11(3:Nos-2,:); ... 
           D1(Nos,:); D0(Nos,:)];
    
    % set up Squire matrix and (cross-term) coupling matrix
    A21 = beta*diag(U_final2)*D0(1:Nos,:); %JH flow
    A22 = alpha*diag(U_final1)*D0-(D2-k2*D0)/(zi*Re);

%     A21 = beta*(du*ones(1,length(u))).*D0(1:Nos,:);
%     A22 = alpha*(u*ones(1,length(u))).*D0-(D2-k2*D0)/(zi*Re);
    B22 = D0;
    A22 = [er*D0(1,:); A22(2:Nsq-1,:); er*D0(Nsq,:)];
    A21 = [zeros(1,Nos); A21(2:Nsq-1,:); zeros(1,Nos)];

    % combine all the blocks 
    A = [A11 zeros(Nos,Nsq); A21 A22];
    B = [B11 zeros(Nos,Nsq); zeros(Nsq,Nos) B22];

