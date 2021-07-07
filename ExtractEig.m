function [n1,n2] = ExtractEig(e,a)

%
% This function computes the number of eigenvalues 
% satisfying
%
% a <= Imag(lambda) <= .5
%
% INPUT
% e   = eigenvalues ordered with decreasing imaginary 
%       part
%
% OUTPUT
% n1  = position of first eigenvalue in the interval
% n2  = position of last eigenvalue in the interval

   n1 = 1;
   while (imag(e(n1)) > 0.5)
     n1 = n1+1; 
   end
   
   n2 = n1;
   while (imag(e(n2)) > a)
     n2 = n2+1;
   end
   n2 = n2-1;
