function [D0,D1,D2,D4] = ChebMat(N)


%=========================================
%  Chebyshev differentiation matrices
%  (based on a hybrid pseudospectral method: 
%   these matrices act on the coefficients 
%   of a Chebyshev expansion, rather than 
%   on the function at the collocation 
%   points) 
%
%  input:    N  = resolution (# of coeff)
%  output:   D0 = 0th deriv. matrix
%            D1 = 1st deriv. matrix
%            D2 = 2nd deriv. matrix
%            D4 = 4th deriv. matrix
%=========================================

   %...create D0 (zeroth derivative)
   D0  = [];
   vec = (0:1:N)';
   for j = 0:N
     D0 = [D0 cos(j*pi*vec/N)];
   end

   %...create higher derivative matrices (using the 
   %   Chebyshev recursion relation)
   lv = length(vec);
   D1 = [zeros(lv,1) D0(:,1) 4*D0(:,2)];
   D2 = [zeros(lv,1) zeros(lv,1) 4*D0(:,1)];
   D3 = [zeros(lv,1) zeros(lv,1) zeros(lv,1)];
   D4 = [zeros(lv,1) zeros(lv,1) zeros(lv,1)];
   for j = 3:N
     D1 = [D1 2*j*D0(:,j)+j*D1(:,j-1)/(j-2)];
     D2 = [D2 2*j*D1(:,j)+j*D2(:,j-1)/(j-2)];
     D3 = [D3 2*j*D2(:,j)+j*D3(:,j-1)/(j-2)];
     D4 = [D4 2*j*D3(:,j)+j*D4(:,j-1)/(j-2)];
   end
   