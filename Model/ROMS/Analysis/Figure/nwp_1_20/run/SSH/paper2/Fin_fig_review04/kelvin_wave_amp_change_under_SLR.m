clc; close all; clear all;

%% historical scenario

% x = 1000:1000:500000; % xgrid from the wall
x = 1000:1000:200000; % xgrid from the wall

% x = 100000:1000:500000; % xgrid from the wall
% x = 500000:1000:1000000; % xgrid from the wall
eta0=3.5 %  amplitude (m)
f= 2 * pi / (23*60*60 + 56*60 + 4.1)  %% coriolis frequency (rad/s)
g = 9.81 % gravitational acceleration (m/s)
% h= (x/15000+5) % depth (historical, m) 5 ~ 55
% h(1:100)=(x(1:100)/20000+5);
% 
% h(1:100)=5;
% h(101:200)= (x(1:100)/7000+5) % depth (historical, m) 5 ~ 20

h(1:50)=5;
h(51:200)= (x(1:150)/2800+5) % depth (historical, m) 5 ~ 20

% h= 5 % depth (historical, m)

Rd= sqrt(g*h)/f  % Rossby radius of deformation ~= 235km


eta = eta0*exp(-x ./ Rd)
plot(x/1000,eta*100)
xlabel('distance from the wall (km)')
ylabel('kelvin wave amplitude (cm)')

%% future scenario

eta0=3.5 % amplitude (m)
f= 2 * pi / (23*60*60 + 56*60 + 4.1)  %% coriolis frequency (rad/s)
g = 9.81 % gravitational acceleration (m/s)
h= h+0.7 % depth (historical, m)
Rd= sqrt(g*h)/f  % Rossby radius of deformation ~= 235km

% x = 0:1000:500000; % xgrid from the wall
eta2 = eta0*exp(-x ./ Rd)
hold on
plot(x/1000,eta2*100)
hold off

plot(x/1000, (eta2-eta)*100)
xlabel('distance from the wall (km)')
ylabel('difference of kelvin wave amplitude (cm)')



% a=3
% flag(1) = a<3;
% flag(2) = a==3;f
% flag(3) = a>3;
% flag(4) = a<4;
% for i=1:length(flag)
%     if flag(i)
%         b(i)=a
%     end
% end
% x=1:10
% datastats(x')