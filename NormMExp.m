
function a = NormMExp(t)
%
% This function computes the norm of the matrix exponential 
% of qb*t 
     
     global qb
     
     a = -norm(expm(t*qb));

