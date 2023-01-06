function [colormap, error_status] = Func_0009_get_colormaps(colormapname, dropboxpath)
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
%  Updated    06-Jul-2021 by Yong-Yub Kim (jet_mod)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch colormapname
    case('byr')
        colormap = customcolormap_preset('red-yellow-blue');
    case('yr')
        byrmap = customcolormap_preset('red-yellow-blue');
        colormap = byrmap(129:256,:);
    case('by')
        byrmap = customcolormap_preset('red-yellow-blue');
        colormap = byrmap(1:128,:);
    case('gmt_ocn')
        load([dropboxpath, filesep, 'source', filesep, 'matlab', filesep, 'Common', ...
            filesep, 'Figure', filesep, 'gmt_ocean_mod2.mat'])  % % set colormap (gmt_ocean, nonwhite)
        colormap = gmt_ocean_mod2;
    case('byr2')
        load([dropboxpath, filesep, 'source', filesep, 'matlab', filesep, 'Common', ...
            filesep, 'Figure', filesep, 'byrmap2.mat'])  % % set colormap (blue-yellow-red)
        colormap = byrmap2;
    case('byr3')
        load([dropboxpath, filesep, 'source', filesep, 'matlab', filesep, 'Common', ...
            filesep, 'Figure', filesep, 'byrmap3.mat'])  % % set colormap (blue-yellow-red)
        colormap = byrmap3;
    case('yr3')
        load([dropboxpath, filesep, 'source', filesep, 'matlab', filesep, 'Common', ...
            filesep, 'Figure', filesep, 'byrmap3.mat'])  % % set colormap (blue-yellow-red)
        colormap = byrmap3(501:1000,:);
    case('by3')
        load([dropboxpath, filesep, 'source', filesep, 'matlab', filesep, 'Common', ...
            filesep, 'Figure', filesep, 'byrmap3.mat'])  % % set colormap (blue-yellow-red)
        colormap = byrmap3(1:500,:);
    case('jet_mod')
        load([dropboxpath, filesep, 'source', filesep, 'matlab', filesep, 'Common', ...
            filesep, 'Figure', filesep, 'jet_mod.mat'])  % % set colormap (jet_modified)
        colormap = jet_mod;
    case('bwr')
        load([dropboxpath, filesep, 'source', filesep, 'matlab', filesep, 'Common', ...
            filesep, 'Figure', filesep, 'bwr_map.mat'])  % % set colormap (bwrmap)
        colormap = bwrmap;
    case('jet')
        colormap =jet;
end

error_status=1;

end

