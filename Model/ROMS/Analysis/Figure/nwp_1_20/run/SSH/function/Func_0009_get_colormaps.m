function [colormap, error_status] = Func_0009_get_colormaps(colormapname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [colormap, error_status] = Func_0009_get_colormaps(colormapname);
%
% get the polygon data from nwp_polygon_point.m corresponding to regionname
%
%  input:
%  colormap name          colormap name (string)
%
%  output:
%  colormap               colormap (2-D array, [R; G; B;])
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    20-May-2021 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch colormapname
    case('byr')
        colormap = customcolormap_preset('red-yellow-blue');
    case('yr')
        colormap = byrmap(129:256,:);
    case('gmt_ocn')
        load('C:\Users\User\Dropbox\source\matlab\Common\Figure\gmt_ocean_mod2.mat')  % % set colormap (gmt_ocean, nonwhite)
        colormap = gmt_ocean_mod2;
end

error_status=1;

end

