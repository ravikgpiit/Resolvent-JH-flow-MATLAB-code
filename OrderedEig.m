function [xs,es] = OrderedEig(A)
%
% This function computes the eigenvalues of a matrix A and 
% orders the eigenvalues so that the imaginary parts are 
% decreasing.
%
% INPUT 
% A    = input matrix
%
% OUTPUT
% es   = ordered eigenvalues
% xs   = eigenvectors

  [v,e]      = eig(A);
  e          = diag(e);
  [eimag,is] = sort(-imag(e));
  xs         = v(:,is); 
  es         = e(is);
  