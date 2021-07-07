function x = ENormalize(x,M)

%  
% This function normalizes the columns of x such that
%
%  || M x_i ||_2 = 1
%
    nc = size(x,2); 
    for i = 1:nc
      x(:,i) = x(:,i)/norm(M*x(:,i));
    end
