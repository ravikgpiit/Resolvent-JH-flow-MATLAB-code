function dUdt = myJeffreyHammel(y,U,a,Re)
dUdt            = zeros(3,1);
dUdt(1)         = U(2);
dUdt(2)         = U(3);
dUdt(3)         = -2*a*Re*U(1).*U(2);
end