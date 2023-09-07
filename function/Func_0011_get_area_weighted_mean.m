function [mean_data, error_status] = Func_0011_get_area_weighted_mean(data, lon, lat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [mean_data, error_status] = Func_0010_get_mask_from_data(data, lon, lat)
%
% get the polygon data from nwp_polygon_point.m corresponding to regionname
%
%  input:
%  data                   data (2-D array, [lon, lat]) or 3-D array
%  lon                    lon (2-D array, [lon, lat])
%  lat                    lat (2-D array, [lon, lat])
%
%  output:
%  mean_data              area-weighted mean (scalar)
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    24-May-2021 by Yong-Yub Kim
%  Updated    07-Jun-2021 by Yong-Yub Kim   % convert 1-d lonlat to 2-d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mean_diff_lon = mean(diff(lon(:,1)));
mean_diff_lat = mean(diff(lat(1,:)));
if (mean_diff_lon == mean_diff_lat)
    dA = (m_lldist([mean_diff_lon 0], [0 0])*1e3)^2 * cosd(lat);
    error_status=1;
else
%     m_lldist(
%     dist_lon=diff(lon);
%     dist_lon(end+1,:)=dist_lon(end,:);
%     dist_lat=diff(lat')';
%     dist_lat(:,end+1)=dist_lat(:,end);
    [~,londim]=size(lon);
    [~,latdim]=size(lat);
    if londim == 1 && latdim ==1
        [temp_lat, temp_lon] = meshgrid(lat,lon);
        lat=temp_lat;
        lon=temp_lon;
    end    
    
    for j=1:size(lon,2)
        xdist(:,j)=m_lldist(lon(:,j), lat(:,j));
    end
    xdist(end+1,:)=xdist(end,:);
    for i=1:size(lon,1)
        ydist(i,:)=m_lldist(lon(i,:), lat(i,:));
    end
    ydist(:,end+1)=ydist(:,end);
    dA=xdist.*ydist;
    error_status=2;
end

mask_rho=NaN(size(data));
mask_rho(isfinite(data))=1;
if ismatrix(data)
    dA=dA.*mask_rho;
    dA_sum = sum(dA(:), 'omitnan');
    mean_data=sum(data.*dA, 'all', 'omitnan')/dA_sum;
elseif ndims(data) ==3
    mean_data=NaN(size(data,3),1);
    for i=1:size(data,3)
        dA=dA.*mask_rho(:,:,i);
        dA_sum = sum(dA(:), 'omitnan');
        mean_data(i)=sum(data(:,:,i).*dA, 'all', 'omitnan')/dA_sum;
    end
end
   

