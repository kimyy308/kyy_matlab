function [mean_data, error_status] = Func_0021_get_vert_area_weighted_mean(data, lon, lat, thick, lonlatflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [mean_data, error_status] = Func_0021_get_vert_area_weighted_mean(data, lonlat, depth, lonlatflag)
%
% get the vertical area weighted mean
%
%  input:
%  data                   data (2-D array, [lon or lat, depth])
%  lonlat                    lon or lat (2-D array, [lon or lat, depth])
%  thick                   thickness (2-D array, [lon or lat, depth])
%  lonlatflag               flag (lon=1 lat=2)
%
%  output:
%  mean_data              area-weighted mean (scalar)
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    22-Feb-2022 by Yong-Yub Kim  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lon(:,1)
% lat(:,1)
    dist2=m_lldist(lon(:,1), lat(:,1));
    dist2=[dist2; dist2(end)];
    dist2=repmat(dist2, 1, size(lat,2),1).*1000;
    dA = dist2.*thick;

    data=squeeze(data);
    mask_rho=NaN(size(data));
    mask_rho(isfinite(data))=1;
    dA=dA.*mask_rho;
    dA_sum = sum(dA(:), 'omitnan');
    if ismatrix(data)
        mean_data=sum(data.*dA, 'all', 'omitnan')/dA_sum;
    elseif ndims(data) ==3
        mean_data=NaN(size(data,3),1);
        for i=1:size(data,3)
            mean_data(i)=sum(data(:,:,i).*dA, 'all', 'omitnan')/dA_sum;
        end
    end
    
    error_status=1;
end

