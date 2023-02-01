function [psd, psd_sig, freq, var] = Func_0027_get_PSD_siglev(data, sig_level, DOF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [psd, psd_sig, freq] = Func_0027_get_PSD_siglev(data, sig_level, DOF)
%
% Power Spectral Density using 1-d data [t] based on fft, with significant psd
%
%  input:
%  data         Time Series Data (t)
%  sig_level    Significance level (0~1), default : 0.95 (optional)
%  DOF          Degree of Freedom, default : 2 (optional)
%
%  output:
%  psd          Power Spectral Density on frequncy domain. (wave component^2 / wavenumber / cycle per period)
%  psd_sig      Significant Power spectral density level under significance level
%  freq         Frequency (cycle/samples). Period: 1/freq (ex: 50y data, 1/10 freq -> period=10y)
%  var          Variance (psd * freq)
%
%  e-mail:      kimyy308@pusan.ac.kr
%
%  Updated      29-Jan-2023 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% squeeze
data = squeeze(data)';

%% get the number of the data(n)
n=length(data);

%% determine the frequncy (if sampling rate = 1, freq: 1/n(0), 2/n, 3/n, ... n/n)
sampling_rate=1;
nyquist=sampling_rate * 1/2;
freq=(1:n/2)/(n/2)*nyquist;
period=1./freq;

%% conduct Fast Fourier Transform
pft=fft(data); % complex double
pft(1) = []; % value with 1/n freq (~ 0 freq; n period ~inifinte period)) is mean value of results
pft_real=real(pft)'; % real part of pft
pft_imag=imag(pft)'; % imaginary part of pft
psd= (1./n) .* 2 .* abs(pft(1:floor(n/2))).^2 / 0.3974;  %pft (2~n/2) only. 2/n+1 ~ n: duplicated 

[psd2, freq2]=pwelch(data,25,1,49,1);
% freq2=freq2/pi;
% pwelch(data)
% h=spectrum.welch;
% psd3=psd(h, data);

var=squeeze(freq).*squeeze(psd);

%% Statistical significance for Fourier spectrum

%% get autocorrelation coefficient(r) with lag 2
r=autocorr(data,2); 
r1=r(2);
r1=(r(2)+r(3)^(1.0/2.0))/2.0;


%% get Chi-square value
%% default Chi=5.991; %Chi-Square values 95% with DOF 2 (chi2inv(0.95, 2))
if nargin == 1
    sig_level=0.95;
    DOF=2;
end
Chi_sq=chi2inv(sig_level,DOF);  % get chi-square value
for i=1:floor(n/2)
    %% discrete power spectrum(P_k) from Markov process. If r1 == 0 : white noise
    mkv(i)=(1) ./ (1 + r1^2 - 2*r1*cos(2*pi*freq(i)));
%     mkv(i)=(1-r1^2) ./ (1 + r1^2 - 2*r1*cos(2*pi*freq(i)));
%     Red(i)=(1-r1^2) ./ (1 + r1^2 - 2*r1*cos(2*pi*freq(i)));

    %% if data(Y_t) is normally distributed, then ampltiude of Fourier component(C_k) is normally distributed.
    %% then, (C_k)^2 is chi-square distributed with 2 DOF
    %% (S_y)^2 is the variance of Y_t
    %% n* amplitude of discreted component / variance of total data
    %% n(abs(C_k)^2 / (S_y)^2 --> 1/DOF * P_k * chi_sq
end
sum1=0;
sum2=0;
for i=1:floor(n/2)
    sum1= sum1+mkv(i).*freq(i);
    sum2= sum2+var(i);
end
scale_factor = sum2 ./ sum1;
Red = mkv .* scale_factor;
psd_sig = Red ./ DOF .* Chi_sq; % Rednoise significance level


end


% Daewon's test code
% x = test; 
% h = spectrum.welch;    % Create a Welch spectral estimator.
% Hpsd = psd(h,x);             % Calculate the PSD
% f=Hpsd.Frequencies;
% power=Hpsd.Data;
% figure
% loglog(f,power)

