function M = EnergyMatrix(Nos,Nsq,k2);

%
% Program to compute the energy weight matrix for three- 
% dimensional Poiseuille and Couette flows
%
% INPUT 
% Nos   = Number of normal velocity modes
% Nsq   = Number of normal vorticity modes 
% k2    = alpha^2 + beta^2
%
% OUTPUT
% M     = weight matrix

   M   = eye(Nos+Nsq,Nos+Nsq); 
   Cos = TwoNormWeight(Nos); 
   Dos = DerivCoeff(Nos);
   Wos = Dos'*Cos*Dos + k2*Cos;
   Wsq = TwoNormWeight(Nsq);
   
   %...Orr-Sommerfeld part 
   [u,s,v] = svd(Wos); 
   s       = sqrt(diag(s));
   Mos     = diag(s)*u';
   
   %...Squire part 
   [u,s,v] = svd(Wsq); 
   s       = sqrt(diag(s));
   Msq     = diag(s)*u';
   
   %...patch together energy matrix 
   M       =[Mos zeros(Nos,Nsq); ...
             zeros(Nsq,Nos) Msq];
   