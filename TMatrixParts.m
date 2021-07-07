function [F,e,invF] = TMatrixParts(M,xu,e)

%
% This function computes the matrix  Q=-i*F*diag(e)*inv(F) 
% which is used to compute the maximum transient growth 
% (see Reddy and Henningson, "Energy Growth in Viscous 
% Channel Flows", JFM 252, page 209, 1993). 
% 
% INPUT
% M     =  energy weight matrix
% xu    =  matrix of eigenfunctions (expansion coefficients)
% e     =  vector of eigenvalues of the stability matrix
%
% OUTPUT
% e     =  matrix of eigenvalues
% F     =  matrix defining the energy norm in the reduced space
%          based on eigenvectors E_norm(u)=2_Norm(Fu)
% invF  =  inverse of F  
%

    % Phase 1: compute inner product of the eigenfunctions 
    %          in energy norm
    work = M*xu;
    A    = work'*work;
    
    % Phase 2: compute decomposition A=F^*F
    [U,S,V] = svd(A);
    s       = sqrt(diag(S));
    F       = diag(s)*U';
    invF    = U*diag(ones(size(s))./s);
