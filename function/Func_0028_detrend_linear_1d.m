function [data_det, tr] = Func_0028_detrend_linear_1d(data, nanflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [psd, psd_sig, freq] = Func_0027_get_PSD_siglev(data, sig_level, DOF)
%
% Power Spectral Density using 1-d data [t] based on fft, with significant psd
%
%  input:
%  data         Time Series Data (t)
%  nanflag      'omitnan' (optional)
%
%  output:
%  data_det     linearly detrended data (t)
%  trend        trend (scalar)
%
%  e-mail:      kimyy308@pusan.ac.kr
%
%  Updated      30-Jan-2023 by Yong-Yub Kim
%  Updated      30-Aug-2023 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% squeeze
data = squeeze(data)';

if nargin < 2
    n=length(data);
    t=1:n;
    
    p=polyfit(t,data,1); %linear fit
    data_fit=p(1)*t + p(2);
    data_det=data-data_fit;

else
    data_finite=data(isfinite(data));
    n=length(data_finite);
    t=1:n;
    p=polyfit(t,data_finite,1); %linear fit
    data_fit=p(1)*t + p(2);
    data_det_finite=data_finite-data_fit;
    data_det=data;
    data_det(isfinite(data))=data_det_finite;
end

data_det=data_det';
tr=p(1);

end


