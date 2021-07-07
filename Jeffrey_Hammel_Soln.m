
function [ U_final1,U_final2,U_final3 ] = Jeffrey_Hammel_Soln( )
% clear all;clc;close all;
yspan       = [0 1];
lambda0     = 0.4;
dlambda     = 0.001;
%  S       = 0;
% a = 0;
a = 1*(pi/180);
% a = 0.5*(pi/180);
% a = -0.005*(pi/180);
Re = 250;
% S           = 1*(pi/180)*250;% diverging angle % S = alpha*Re; alpha = 1degree = pi/180 radian, Re = 250
% S           = 0.5*(pi/180)*250;
% S             = 1*(pi/180)*500;
% S             = 1*(pi/180)*100;
% S           = -0.01*(pi/180)*9000; % converging angle
error       = 10;
F           = zeros(2,1);
while error > 1e-4
    lambda      = [lambda0 lambda0+dlambda];
    for k = 1:length(lambda)
        init        = [1 0 lambda(k)];
%         [y,U]       = ode45(@(y,U)myJeffreyHammel(y,U,S),yspan,init);
         [y,U]       = ode45(@(y,U)myJeffreyHammel(y,U,a,Re),yspan,init);
        F(k)        = U(end,1);
    end
    dFdlamb     = (F(2)-F(1))/dlambda;
    lambda1     = lambda0 - F(1)/dFdlamb;
    error       = abs((lambda1 - lambda0)/(lambda0+eps));
    lambda0     = lambda1;
end
init        = [1 0 lambda0];
% [y,U]       = ode45(@(y,U)myJeffreyHammel(y,U,S),yspan,init);
[y,U]       = ode45(@(y,U)myJeffreyHammel(y,U,a,Re),yspan,init);
% F(k)        = U(end,1);
y_temp      = [-flipud(y(2:end))' y'];
U_temp      = zeros(length(y_temp),4);
% num_point   = 151;
% y1          = zeros(num_point,1);
% for k = 1:num_point
%     y1(k)           = cos(pi*(k-1)/(num_point-1));
% end
% F_int       = zeros(length(y1),4);
for k = 1:3
    if mod(k,2) == 1
        U_temp(1:length(y)-1,k)         = flipud(U(2:length(y),k));
        U_temp(length(y):end,k)         = U(:,k);
%         F_int(:,k)                      = interp1(y_temp,F_temp(:,k),y1);
    else
        U_temp(1:length(y)-1,k)         = -flipud(U(2:length(y),k));
        U_temp(length(y):end,k)         = U(:,k);
%         F_int(:,k)                      = interp1(y_temp,F_temp(:,k),y1);
    end
end
% % plot(y_temp,F_temp(:,4),y1,F_int(:,4));
% % plot(y1,F_int(:,2));

% figure(1);
% plot(U_temp(:,1),y_temp,'*b');
% figure(2)
% plot(U_temp(:,2),y_temp,'*r');
% figure(3)
% plot(U_temp(:,3),y_temp,'*g');
% figure(1);
% hold on;
% plot(U(:,1),y);
% xlabel('U');ylabel('y');
% hold on;
% plot(U(:,1), -y);
U1 = flipud(U(:,1));
U2 = -flipud(U(:,2));
U3 = flipud(U(:,3));
Uf = vertcat(U1,U(2:end,1));
Uf2 = vertcat(U2,U(2:end,2));
Uf3 = vertcat(U3,U(2:end,3));
y1 = flipud(y);
yf = vertcat(-y1,y(2:end));

%  figure(111)
hold on;
N = 100;
vec = [0:N]';
xq = cos(pi*vec/N);
U_final1 = interp1(yf,Uf,xq,'pchip');
U_final2 = interp1(yf,Uf2,xq,'pchip');
U_final3 = interp1(yf,Uf3,xq,'pchip');
% plot(U_final1,xq,'ob');
% figure(222)
% hold on;
% plot(U_final2,xq,'or');
% 
% figure(333)
% hold on;
% plot(U_final3,xq,'og');
end


