function [region_mask] = get_mask_from_polygon(lon, lat, refpolygon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function var=get_GCMname_from_RCM(testname);
%
% get the GCM testname corresponding to RCM test name
%
%  input:
%  lon            1-d or 2-d longitude data
%  lat            1-d or 2-d latitude data
%  polygon        polygon data
% 
%  output:
%  region_mask          mask data in polygon
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    14-May-2021 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (ndims(lon)==1)
    [lat2, lon2]=meshgrid(lat, lon);
else
    lon2 = lon; lat2 = lat;
    region_mask = double(inpolygon(lon2, lat2, refpolygon(:,1), refpolygon(:,2)));
    region_mask(region_mask==0)=NaN;
end

end

