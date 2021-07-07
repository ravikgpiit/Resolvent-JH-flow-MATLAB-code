function dcoef = DerivCoeff(N)

%
% Compute matrix which converts Chebyshev coefficients of a
% polynomial to coefficients of derivative
%
% Reference:
% Gottlieb and Orszag, Numerical Analysis of Spectral 
% Methods: Theory and Applications, SIAM, Philadelphia, 
% 1977.
%
% d1  = derivative matrix (N,N) 
% N   = number of coefficients 
%
   dcoef  = zeros(N,N);

   for i = 0:N-1
     for j = (i+1):2:(N-1)
       dcoef(i+1,j+1) = 2*j;
     end
   end
   dcoef(1,:) = dcoef(1,:)/2;
