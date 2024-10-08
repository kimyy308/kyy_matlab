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
%         load([dropboxpath, filesep, 'source', filesep, 'matlab', filesep, 'Common', ...
%             filesep, 'Figure', filesep, 'bwr_map.mat'])  % % set colormap (bwrmap)
        bwrmap = customcolormap_preset('red-white-blue');
        colormap = bwrmap;
    case('wr')
%         bwrmap = customcolormap_preset('red-white-blue');
        load([dropboxpath, filesep, 'source', filesep, 'matlab', filesep, 'Common', ...
            filesep, 'Figure', filesep, 'bwr_map.mat'])  % % set colormap (bwrmap)
        colormap = bwrmap(51:100,:);
    case('bwr_10')
        colormap = customcolormap_preset('red-white-blue_10');
    case('bwr_20')
        colormap = customcolormap_preset('red-white-blue_20');
    case('bwr_10_botgray')
        colormap2 = customcolormap_preset('red-white-blue_10');
        for i=1:10
            for j=1:10
                colormap((i-1)*10+j,:)=colormap2(i,:);
            end
        end
        colormap(1,:)=[0.6 0.6 0.6];
    case('wr_10_topgray')
        colormap2 = customcolormap_preset('red-white-blue_20');
        colormap2=colormap2(11:20,:);
        for i=1:10
            for j=1:10
                colormap((i-1)*10+j,:)=colormap2(i,:);
            end
        end
        colormap(end,:)=[0.6 0.6 0.6];
    case('wr_10')
        colormap = customcolormap_preset('red-white-blue_20');
        colormap=colormap(11:20,:);
    case('wr_5')
        colormap = customcolormap_preset('red-white-blue_10');
        colormap=colormap(6:10,:);
    case('bw_10')
        colormap = customcolormap_preset('red-white-blue_20');
        colormap=colormap(1:10,:);
    case('wr_08')
        colormap = customcolormap_preset('red-white-blue_16');
        colormap=colormap(9:16,:);
    case('jet')
        colormap =jet;
    case('gray')
        colormap = [0.8 0.8 0.8];
    case('NCAR_SMYLE')
        colormap = [9 48 107; ...
            7 80 155; ...
            31 111 180; ...
            66 146 197; ...
            105 173 213; ...
            158 202 225; ...
            197 219 239; ...
            220 234 247; ...
            246 251 255; ...
            255 255 255; ...
            255 255 255; ...
            225 225 225; ...
            191 191 191; ...
            160 160 160; ...
            129 129 129; ...
            252 252 202; ...
            254 219 122; ...
            253 141 61; ...
            227 30 30; ...
            128 0 38];
        colormap = colormap / 255.0;
    case('bwg_10') % brown 2 green 
%         colormap = [89 44 1; ...
%             147 75 2; ...
%             199 123 30; ...
%             227 191 119; ...
%             248 230 192; ...
%             199 238 231; ...
%             107 208 193; ...
%             5 153 144; ...
%             1 104 94; ...
%             1 61 47];
        colormap = [84 48 5; ...
            140 81 10; ...
            191 129 45; ...
            223 194 125; ...
            246 232 195; ...
            199 234 229; ...
            128 205 193; ...
            53 151 143; ...
            1 102 94; ...
            0 60 48];
        colormap = colormap / 255.0;
    case('bwg_20') % brown 2 green 
        colormap = [89 44 0; ...
            118 59 1; ...
            147 75 2; ...
            173 99 16; ...
            199 123 30; ...
            213 157 75; ...
            227 191 119; ...
            238 211 155; ...
            248 230 192; ...
            231 233 205; ...
            214 236 218; ...
            199 238 231; ...
            153 223 212; ...
            107 208 193; ...
            56 181 169; ...
            5 153 144; ...
            3 129 119; ...
            1 104 94; ...
            1 83 71; ...
            0 61 47];
        colormap = colormap / 255.0;
end

error_status=1;

end

