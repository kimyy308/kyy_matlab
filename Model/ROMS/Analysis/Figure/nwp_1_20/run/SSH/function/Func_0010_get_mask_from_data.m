function [mask_rho, nanmask, error_status] = Func_0010_get_mask_from_data(data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [mask, nanmask, error_status] = Func_0010_get_mask_from_data(data);
%
% get the polygon data from nwp_polygon_point.m corresponding to regionname
%
%  input:
%  data                   data (2-D array)
%
%  output:
%  mask                   ocean mask (2-D array)
%  nanmask                land mask (2-D array)
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    20-May-2021 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mask_rho=NaN(size(data));
nanmask=NaN(size(data));
mask_rho(isfinite(data))=1;
nanmask(isnan(data))=1;
error_status=1;

end

