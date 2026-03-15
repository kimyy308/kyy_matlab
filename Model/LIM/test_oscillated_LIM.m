clear; close all; clc

%% =========================================================
% PARAMETERS
%% =========================================================

nyear = 1000;
nmon  = nyear*12;
t     = (1:nmon)';

lag = 21;

%% =========================================================
% Y : persistent slow oscillation
%% =========================================================

phi = 0.90;
Aosc = 0.15;
period_y = 36;

y = zeros(nmon,1);

for k = 2:nmon

    osc = Aosc*sin(2*pi*k/period_y);

    y(k) = phi*y(k-1) + osc + 0.4*randn;

end

y = y - mean(y);

%% =========================================================
% X : dynamically driven by lagged y
%% =========================================================

ax = 0.4;     % persistence in x
b  = 0.8;     % coupling strength

x = zeros(nmon,1);

for k = 2:nmon

    if k > lag
        forcing = b*y(k-lag);
    else
        forcing = 0;
    end

    x(k) = ax*x(k-1) + forcing + 0.6*randn;

end

x = x - mean(x);

%% =========================================================
% POWER SPECTRUM
%% =========================================================

[pxx,f] = pwelch(x,[],[],[],12);
[pyy,~] = pwelch(y,[],[],[],12);

period_year = 1./f;
period_year(1) = NaN;

%% =========================================================
% LEAD-LAG CORRELATION
%% =========================================================

maxlag = 60;

acf_x  = xcorr(x,maxlag,'coeff');
ccf_xy = xcorr(y,x,maxlag,'coeff');  % y leads x

lags = -maxlag:maxlag;

%% =========================================================
% LIM
%% =========================================================

x1 = x(1:end-1);
x2 = x(2:end);
A1 = (x1'*x2)/(x1'*x1);

X = [x y];

X1 = X(1:end-1,:);
X2 = X(2:end,:);

A2 = (X1'*X1)\(X1'*X2);

%% =========================================================
% FORECAST SKILL (multiple initial times)
%% =========================================================

maxlead = 48;

init_step = 1;
init_idx = 1:init_step:(nmon-maxlead);

ninit = length(init_idx);

skill1 = zeros(maxlead+1,1);
skill2 = zeros(maxlead+1,1);

for L = 0:maxlead

    pred1 = zeros(ninit,1);
    pred2 = zeros(ninit,1);
    truth = zeros(ninit,1);

    for i = 1:ninit

        t0 = init_idx(i);

        truth(i) = x(t0+L);

        pred1(i) = A1^L * x(t0);

        pred_state = A2^L * [x(t0); y(t0)];
        pred2(i) = pred_state(1);

    end

    skill1(L+1) = corr(pred1,truth);
    skill2(L+1) = corr(pred2,truth);

end

lead = 0:maxlead;

%% =========================================================
% PLOT
%% =========================================================

figure('position',[100 100 1200 900])

subplot(4,2,1)
plot(t/12,x,'k')
title('x anomaly time series')
xlabel('year')

subplot(4,2,2)
plot(t/12,y,'k')
title('y persistent oscillation')
xlabel('year')

subplot(4,2,3)
plot(period_year,pxx,'k','linewidth',1.2)
xlim([0 8])
set(gca,'xtick',1:1:8)
xlabel('period (years)')
title('Power spectrum (x)')

subplot(4,2,4)
plot(period_year,pyy,'k','linewidth',1.2)
xlim([0 8])
set(gca,'xtick',1:1:8)
xlabel('period (years)')
title('Power spectrum (y)')
xline(3,'r--')

subplot(4,2,5)
plot(lags,acf_x,'k')
title('Lead-lag corr (x vs x)')
xlabel('lag (months)')
ylim([-1 1])

subplot(4,2,6)
plot(lags,ccf_xy,'k')
title('Lead-lag corr (y leads x)')
xlabel('lag (months)')
ylim([-1 1])
xline(-21,'r--')

subplot(4,2,7)
plot(lead,skill1,'k','linewidth',1.5)
title('LIM skill (x → x)')
xlabel('lead (months)')
ylim([-0.5 1])

subplot(4,2,8)
plot(lead,skill2,'k','linewidth',1.5)
title('LIM skill (x,y → x)')
xlabel('lead (months)')
ylim([-0.5 1])