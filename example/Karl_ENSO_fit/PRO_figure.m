close all
clear
w=1/(360*3.1);  % physical intrinsic freq of ENSO (fixed)  not even num
wamp=0; % ignore it (cannot make big difference)
g=-1/(360*1.6); % damping 
gamp=g*2; % + -> damping amplitude  will modulated based on this 0 value
ffreq=1/360;
n=0.015;
namp=0;
cubic=1/(360*5)*0; % cubic damping we don't care
f=1/120*0;
dt=.5;
L=360*200;
IC=[0 1];

%g=g-(3e-5)*dt;

    
%[Y,C,W,PM]=SSRMrk(w,wamp,g,gamp,ffreq,n,namp,cubic,dt,L,IC);
% [Y,C,W,PM]=PRO(w,wamp,g,gamp,ffreq,n,namp,cubic,f,dt,L,IC);
[Y,C,W,PM]=PRO(w,g,gamp,ffreq,n,namp,cubic,f,dt,L,IC);

t=1/24:1/12:L/360;

Y=reshape(Y,2,30/dt,L/30);
Y=squeeze(mean(Y,2));
% Y=Y(:,601:end);
% t=t(601:end);
Y=Y(:,1:end);
t=t(1:end);

% [Syy,f,Syylow,Syyhigh]=specprog(Y(1,:),1/12,1,3,0.25);

NFFT=256*10;
[Pxx,F] = pyulear(Y(1,:),100,NFFT,12);
[Pxx1,F1] = pyulear(Y(1,:),1,NFFT,12);


subplot(211)
plot(t,Y(1,:),'r',t,Y(2,:),'k--')


seasv=seasonalvariance(Y(1,:)',12);

subplot(234)
% plot(1:12,seasonalvariance(Y(1,:)',12)/var(Y(1,:)))
plot(1:12,seasv/var(Y(1,:)));
axis([1 12 0 2])
% monthlabel
title('Seasonal variance')
ylabel('SSTA var ($^\circ$C)')

subplot(235)
loglog(F,Pxx,'r-','Linewidth',2)
hold on
plot(f,Syy,'b')
plot(f,Syylow,'k--',f,Syyhigh,'k--')
loglog(F1,Pxx1,'k','Linewidth',1)
hold off
axis([.02 5 10^-5 10^5])
ylabel('Power')
xlabel('Frequency ($\mbox{yrs}^{-1}$)')
title('Spectrum')

subplot(236)
plot(Y(1,:),Y(2,:),Y(1,1:12:end),Y(2,1:12:end),'r.')


%kyy_
ylen=5;
function a=seasonalvariance(Y, m)
    rsp_T=reshape(Y,[m length(Y)/m]);
    a=var(rsp_T,0,2);
end