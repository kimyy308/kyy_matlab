
close all; clc; clear all;


%% white noise
% Y = rand(1,50);
Z=rand(1,5000);
Y=fft(Z);
% plot(real(Y(2:end)))
Y2=Y;
Y2(1)=[];
pxx=abs(Y2(1:floor(n/2))).^2 ./n .*2;
plot(pxx)
n = 5000;
X = ifft(Y,n,2, 'symmetric');
size(X)

plot(X)
% plot(X-Z)


%% red noise example
Z_low=lowpass(Z,100/5000,1);
plot(Z_low);
Y_low=fft(Z_low);
plot(real(Y_low(2:end/2)))
n=length(Y_low);
Y_low2=Y_low;
Y_low2(1)=[];
pxx_low=abs(Y_low2(1:floor(n/2))).^2 ./n .*2;
plot(pxx_low)

X_low= ifft(Y_low,n,2, 'symmetric');
plot(X_low);


%% red noise2 example
Z_low2=lowpass(Z_low,10/5000,1);
plot(Z_low2);
Y_low2=fft(Z_low2);
plot(real(Y_low2(2:end/2)))
n=length(Y_low2);
Y_low22=Y_low2;
Y_low22(1)=[];
pxx_low2=abs(Y_low22(1:floor(n/2))).^2 ./n .*2;
plot(pxx_low2)

X_low2= ifft(Y_low22,n,2, 'symmetric');
plot(X_low2);