function [data_det] = Func_0028_detrend_linear_1d(data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [psd, psd_sig, freq] = Func_0027_get_PSD_siglev(data, sig_level, DOF)
%
% Power Spectral Density using 1-d data [t] based on fft, with significant psd
%
%  input:
%  data         Time Series Data (t)
%
%  output:
%  psd          linearly detrended data (t)
%
%  e-mail:      kimyy308@pusan.ac.kr
%
%  Updated      30-Jan-2023 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% squeeze
data = squeeze(data)';

n=length(data);
t=1:n;

p=polyfit(t,data,1); %linear fit
data_fit=p(1)*t + p(2);
data_det=data-data_fit; 

end


