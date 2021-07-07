
function t=OptFunc(t1,t2,options);
%
% This function uses the built-in function 'fminbnd' to find 
% the maximum value of a function on the interval [t1,t2]. 
% The function is in the file NormMExp.m
% 
% INPUT: 
% t1,t2   = lower and upper bounds of interval
% options = input parameters for minimization routine 
%
% OUTPUT
% t       = value at which function NormMExp(t) is minimized

    f1     = NormMExp(t1);
    f2     = NormMExp(t2);
    tt     = fminbnd('NormMExp',t1,t2);
%    tt     = fminbnd('NormMExp',t1,t2,options);

    f3     = NormMExp(tt);
    f      = [f1 f2 f3]; 
    tm     = [t1 t2 tt];
    [y,is] = sort(f);
    t      = tm(is(1));





