close all; clear all;  clc;
warning off;

testname = 'test06';
regionname = 'NWP';
inputyear = [1983:1992]; % % put year which you want to plot [year year ...]
inputmonth = [1,2]; % % put month which you want to plot [month month ...]
allyear =[1983:1992];
refyear =[1983:1987];
checktime=[15,30];
caxisval=[0 15];
addpath(genpath(['C:\Users\User\Dropbox\source\matlab\function\']))
[dropboxpath, erorr_status] = Func_0008_set_dropbox_path(computer);
addpath(genpath([dropboxpath, '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\MICT_pollack\paper\subroutine\']))
[byrmap3, error_status] = Func_0009_get_colormaps('byr3', dropboxpath);
[byrmap, error_status] = Func_0009_get_colormaps('byr2', dropboxpath);  
[refpolygon, lonlat, error_status] = Func_0007_get_polygon_data_from_regionname(regionname);

param_script =['C:\users\user/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_10/run/fig_param/fig_param_kyy_', regionname, '.m'];
figrawdir =strcat('D:\research\Ph_D_course\2021_Particle tracking of the walleye pollock\figure\paper\'); % % where figure files will be saved
filedir = strcat('D:\Data\Model\ROMS\nwp_1_10\input\test06\'); % % where data files are
obsdir=['D:\Data\Model\ROMS\nwp_1_10\input\test06\'];

LTRANS_testname='Pollock6';

tifname = [figrawdir, 'FigS01.tif'];

correction_right_fig=[-0.1600,0,0,0]; % right, up, width, height
correction_upper_fig=[0,0,0,0];
correction_large_fig=[0,0,0,0.020];


ncname = [filedir,'roms_grid_nwp_1_10_test06.nc'];
lon_rho=ncread(ncname, 'lon_rho');
lat_rho=ncread(ncname, 'lat_rho');
mask_rho=ncread(ncname, 'mask_rho');
mask_rho(mask_rho==1)=NaN;
mask_rho(mask_rho==0)=1;
RCM_land=mask_rho;

obsname = [obsdir, 'auto01_obs_20210901.nc'];
ncinfo(obsname)
rlon=ncread(obsname, 'rlon');
rlat=ncread(obsname, 'rlat');


run(param_script);
m_grid_fontsize = m_grid_fontsize -4;

% % projection
testnameind=2;
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
ax{testnameind,1}=axes;
% set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

pc{testnameind,1}=m_pcolor(lon_rho',lat_rho', RCM_land','parent',ax{testnameind,1});
colormap(ax{testnameind,1},[0.8 0.8 0.8]);
shading(gca,m_pcolor_shading_method); 
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
       'parent', ax{testnameind,1});
% % % pcolor
hold on
for i=1:length(rlon)
% h1=m_plot(rlon(i),rlat(i),'marker','o','color','r', ...
%                   'markerfacecolor',bwrmap(colind,:), 'parent', ax{testnameind, 1}); 
[indw, inde, inds, indn]=Func_0012_findind_Y(0.05, [rlon(i) rlon(i) rlat(i) rlat(i)], lon_rho, lat_rho, 1);

    if(isnan(RCM_land(indw,inds)))
        h1=m_plot(rlon(i),rlat(i),'marker','o','color','r', 'MarkerSize', 1.5, ...
                          'markerfacecolor','r', 'parent', ax{testnameind, 1}); 
    end
end

txt{testnameind,1}=m_text(118, 48, 'OBS points for DA', 'FontSize', m_grid_fontsize+4); 

% % figure save
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height*2]) 
print('-dtiff','-r500',tifname); RemoveWhiteSpace([], 'file', tifname);
hold off;
close all;

save ('D:\research\Ph_D_course\2021_Particle tracking of the walleye pollock\data_dryad\DataS01_S01.mat', ...
            'lon_rho', 'lat_rho', ...
            'rlon', 'rlat', 'RCM_land')
