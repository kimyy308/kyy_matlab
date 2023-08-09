function [data_des_det] = Func_0029_deseasonalize_detrend_linear_3d(data, nanflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [psd, psd_sig, freq] = Func_0027_get_PSD_siglev(data, sig_level, DOF)
%
% Power Spectral Density using 1-d data [t] based on fft, with significant psd
%
%  input:
%  data         3-D Data (x,y,t)
%  nanflag      'omitnan' (optional)
%
%  output:
%  data_des_det deseasonalized, linearly detrended data (t)
%
%  e-mail:      kimyy308@pusan.ac.kr
%
%  Updated      11-Apr-2023 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% get size
    [ilen, jlen, tlen]=size(data);
    mlen=12;
    ylen=tlen/mlen;
    
    data_reshaped=reshape(data, [ilen, jlen, mlen, ylen]);

    if nargin < 2
        data_clim_mean=squeeze(mean(data_reshaped,4));
        data_des=data_reshaped-data_clim_mean;
        data_des=reshape(data_des,[ilen,jlen,tlen]);
        for i=1:ilen
            for j=1:jlen
                 if (isnan(data(i,j,1))~=1 ...
                         & nansum(data(i,j,:))~=0)
                    data_des_det(i,j,:)= Func_0028_detrend_linear_1d(data_des(i,j,:));
                 else
                    data_des_det(i,j,:)= NaN;
                 end
            end
        end
    
    
    else
        data_clim_mean=squeeze(mean(data_reshaped,4,nanflag));
        data_des=data_reshaped-data_clim_mean;
        data_des=reshape(data_des,[ilen,jlen,tlen]);
        for i=1:ilen
            for j=1:jlen
                 if (isnan(data(i,j,1))~=1 ...
                         & nansum(data(i,j,:))~=0)
                    data_des_det(i,j,1:tlen)= Func_0028_detrend_linear_1d(data_des(i,j,:), nanflag);
                 else
                    data_des_det(i,j,1:tlen)= NaN;
                 end
            end
        end
    end

end


