clc, clear all, close all;
yspan       = [0 1];
lambda0     = 0.4;
dlambda     = 0.001;
% a=-0.005*(pi/180);
a=1*(pi/180);
% Re=9000;
Re = 250;

% S           = a*Re;
error       = 10;
F           = zeros(2,1);
while error > 1e-4
    lambda      = [lambda0 lambda0+dlambda];
    for k = 1:length(lambda)
        init        = [1 0 lambda(k)];
        [y,U]       = ode45(@(y,U)myJeffreyHammel(y,U,a,Re),yspan,init);
        F(k)        = U(end,1);
    end
    dFdlamb     = (F(2)-F(1))/dlambda;
    lambda1     = lambda0 - F(1)/dFdlamb;
    error       = abs((lambda1 - lambda0)/(lambda0+eps));
    lambda0     = lambda1;
end
init        = [1 0 lambda0];
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
% plot(y_temp,F_temp(:,4),y1,F_int(:,4));
% plot(y1,F_int(:,2));
figure(99)
plot(U_temp(:,1),y_temp,'b','LineWidth',2);
hold on
plot(U_temp(:,2),y_temp,'r','LineWidth',2);
plot(U_temp(:,3),y_temp,'k','LineWidth',2);





