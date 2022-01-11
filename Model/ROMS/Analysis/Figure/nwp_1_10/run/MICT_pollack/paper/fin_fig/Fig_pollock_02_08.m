close all; clear all;  clc;
warning off;

testname = 'test06';
regionname = 'ES';
inputyear = [1983:1987]; % % put year which you want to plot [year year ...]\
inputyear2 = [1988:1992]; % % put year which you want to plot [year year ...]
inputmonth = [1,2]; % % put month which you want to plot [month month ...]
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
filedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\'); % % where data files are
savedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\pollock_3\');
inputdir = ['/home/auto/MAMS/Data/01_NWP_1_10/Input/'];
LTRANS_testname='Pollock6';

tifname = [figrawdir, 'Fig02_08.tif'];

% correction_right_fig=[-0.1600,0,0,0]; % right, up, width, height
correction_upper_fig=[0,0,0,0];
% correction_upper_fig=[0.04,0,0,0];

% correction_large_fig=[0,0,0,0.020];
correction_large_fig=[0,0,0,0];



% % %  OISST data
run(param_script);
m_grid_fontsize = m_grid_fontsize -4
ind=1;
clear comb_data
for yearij = 1:length(inputyear)
    tempyear = inputyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
%             ncname = [filedir,testname,regionname,'_rms_clim_',num2str(1983,'%04i'),'_',num2str(2018,'%02i'),'.nc'];
        ncname = [filedir,testname,regionname,'_sq_er_clim_',num2str(1983,'%04i'),'_',num2str(2018,'%02i'),'.nc'];
% model_land=ncread(ncname, 'avhrr_land');
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
        sst_ind=(inputyear(yearij)-1983)*12+inputmonth(monthij);
        data=ncread(ncname, 'avhrr_sst',[1 1 sst_ind], [inf inf 1]);
        lon_OISST=ncread(ncname, 'lon');
        lat_OISST=ncread(ncname, 'lat');
%                         lastday_m=size(data,3);
        if (exist('comb_data')==0)
            comb_data=data;
        else
            comb_data(:,:,end+1)=data;
        end
    end
end
mean_OISST = sum(comb_data,3)/size(comb_data,3);
mean_OISST(mean_OISST==0)=NaN;
mean_OISST_early=mean_OISST;

clear comb_data
for yearij = 1:length(inputyear2)
    tempyear = inputyear2(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
%             ncname = [filedir,testname,regionname,'_rms_clim_',num2str(1983,'%04i'),'_',num2str(2018,'%02i'),'.nc'];
        ncname = [filedir,testname,regionname,'_sq_er_clim_',num2str(1983,'%04i'),'_',num2str(2018,'%02i'),'.nc'];
% model_land=ncread(ncname, 'avhrr_land');
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
        sst_ind=(inputyear2(yearij)-1983)*12+inputmonth(monthij);
        data=ncread(ncname, 'avhrr_sst',[1 1 sst_ind], [inf inf 1]);
        lon_OISST=ncread(ncname, 'lon');
        lat_OISST=ncread(ncname, 'lat');
%                         lastday_m=size(data,3);
        if (exist('comb_data')==0)
            comb_data=data;
        else
            comb_data(:,:,end+1)=data;
        end
    end
end
mean_OISST = sum(comb_data,3)/size(comb_data,3);
mean_OISST(mean_OISST==0)=NaN;
mean_OISST_late=mean_OISST;

diff_OISST=mean_OISST_late-mean_OISST_early;

% % projection
testnameind=1;
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)-0.5],'lat',[lonlat(3) lonlat(4)-0.5]);
sb1=subplot(4,5,[1 2 6 7]);  % Auto-fitted to the figure.
pos_sb{testnameind,1}=get(sb1, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [000,0,0,0];
delete(sb1); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

% % land
OISST_land=ncread(ncname, 'avhrr_land',[1 1], [inf inf]);
pc{testnameind,1}=m_pcolor(lon_OISST',lat_OISST', OISST_land','parent',ax{testnameind,1});
colormap(ax{testnameind,1},[0.8 0.8 0.8]);
shading(gca,m_pcolor_shading_method); 
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
       'xticklabels', [], 'yticklabels', [], 'parent', ax{testnameind,1});

% % % pcolor
hold on
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)-0.5],'lat',[lonlat(3) lonlat(4)-0.5]);
ax{testnameind,2}=axes;
set(ax{testnameind,2},'pos',pos_sb{testnameind});
pc{testnameind,2}=m_pcolor(lon_OISST',lat_OISST', mean_OISST_early','parent',ax{testnameind,2});
colormap(ax{testnameind,2},jet_mod);
caxis(caxisval);
shading(gca,m_pcolor_shading_method);   
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
        'xticklabels', [130, 134, 138], 'xtick',[130, 134, 138], 'XaxisLocation', 'top',...
        'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});

% % % contour
ax{testnameind,3}=axes;
set(ax{testnameind,3},'pos',pos_sb{testnameind});
[con_C,con_h]=m_contour(lon_OISST',lat_OISST', mean_OISST_early', [2, 5, 10], m_contour_color, 'linewidth', m_contour_linewidth, 'parent',ax{testnameind,3});
clabel(con_C,con_h,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
    'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
       'xticklabels', [], 'backcolor', 'none','parent', ax{testnameind,3});
caxis(caxisval);
colormap(ax{testnameind,3},jet_mod);


txt{testnameind,1}=m_text(127.5, 48, 'OISST', 'FontSize', m_grid_fontsize+4); 
txt{testnameind,1}=m_text(127.5, 46.5, '1983―1987', 'FontSize', m_grid_fontsize+4); 

% % projection 22222222222222222
testnameind=2;
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)-0.5],'lat',[lonlat(3) lonlat(4)-0.5]);
sb1=subplot(4,5,[3 4 8 9]);  % Auto-fitted to the figure.
pos_sb{testnameind,1}=get(sb1, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [000,0,0,0];
delete(sb1); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

% % land
OISST_land=ncread(ncname, 'avhrr_land',[1 1], [inf inf]);
pc{testnameind,1}=m_pcolor(lon_OISST',lat_OISST', OISST_land','parent',ax{testnameind,1});
colormap(ax{testnameind,1},[0.8 0.8 0.8]);
shading(gca,m_pcolor_shading_method); 
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
       'xticklabels', [], 'yticklabels', [], 'parent', ax{testnameind,1});

% % % pcolor
hold on
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)-0.5],'lat',[lonlat(3) lonlat(4)-0.5]);
ax{testnameind,2}=axes;
set(ax{testnameind,2},'pos',pos_sb{testnameind});
pc{testnameind,2}=m_pcolor(lon_OISST',lat_OISST', mean_OISST_late','parent',ax{testnameind,2});
colormap(ax{testnameind,2},jet_mod);
caxis(caxisval);
shading(gca,m_pcolor_shading_method);   
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
        'xticklabels', [130, 134, 138], 'xtick',[130, 134, 138], 'XaxisLocation', 'top',...
        'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});

% % % contour
ax{testnameind,3}=axes;
set(ax{testnameind,3},'pos',pos_sb{testnameind});
[con_C,con_h]=m_contour(lon_OISST',lat_OISST', mean_OISST_late', [2, 5, 10], m_contour_color, 'linewidth', m_contour_linewidth, 'parent',ax{testnameind,3});
clabel(con_C,con_h,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
    'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
       'yticklabels', [],'xticklabels', [], 'backcolor', 'none','parent', ax{testnameind,3});
caxis(caxisval);
colormap(ax{testnameind,3},jet_mod);
    
    
txt{testnameind,1}=m_text(127.5, 48, 'OISST', 'FontSize', m_grid_fontsize+4); 
txt{testnameind,1}=m_text(127.5, 46.5, '1988―1992', 'FontSize', m_grid_fontsize+4); 






% % %  RCM data
clear comb_data
for yearij = 1:length(inputyear)
    tempyear = inputyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])

        data=ncread(ncname, 'temp_surf');
        lastday_m=size(data,3);
        if (exist('comb_data')==0)
            comb_data=data;
        else
            comb_data(:,:,end+1:end+lastday_m)=data;
        end
    end
end
lon_RCM = ncread(ncname, 'lon_rho');
lat_RCM = ncread(ncname, 'lat_rho');
mean_RCM = sum(comb_data,3)/size(comb_data,3);
mean_RCM(mean_RCM==0)=NaN;
mask_model = double(inpolygon(lon_RCM,lat_RCM,refpolygon(:,1),refpolygon(:,2)));
mask_model(mask_model==0)=NaN;
mean_RCM=mean_RCM.*mask_model;
mean_RCM_early=mean_RCM;

clear comb_data
for yearij = 1:length(inputyear2)
    tempyear = inputyear2(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])

        data=ncread(ncname, 'temp_surf');
        lastday_m=size(data,3);
        if (exist('comb_data')==0)
            comb_data=data;
        else
            comb_data(:,:,end+1:end+lastday_m)=data;
        end
    end
end
mean_RCM = sum(comb_data,3)/size(comb_data,3);
mean_RCM(mean_RCM==0)=NaN;
mask_model = double(inpolygon(lon_RCM,lat_RCM,refpolygon(:,1),refpolygon(:,2)));
mask_model(mask_model==0)=NaN;
mean_RCM=mean_RCM.*mask_model;
mean_RCM_late=mean_RCM;

diff_RCM=mean_RCM_late-mean_RCM_early;


% % projection333333333333333333
testnameind=3;
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)-0.5],'lat',[lonlat(3) lonlat(4)-0.5]);
sb2=subplot(4,5,[11 12 16 17]);  % Auto-fitted to the figure.
pos_sb{testnameind,1}=get(sb2, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig+ [000,0.23,0,0];
delete(sb2); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

% % land
temp_surf=ncread(ncname, 'temp_surf', [1 1 1], [inf inf 1]);
temp_surf(isnan(temp_surf))=50000;
temp_surf(temp_surf<50000)=NaN;
RCM_land=temp_surf;
RCM_land(temp_surf==50000)=1;
pc{testnameind,1}=m_pcolor(lon_RCM',lat_RCM', RCM_land','parent',ax{testnameind,1});
colormap(ax{testnameind,1},[0.8 0.8 0.8]);
shading(gca,m_pcolor_shading_method); 
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
       'xticklabels', [], 'yticklabels', [], 'parent', ax{testnameind,1});

% % % pcolor
hold on
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)-0.5],'lat',[lonlat(3) lonlat(4)-0.5]);
ax{testnameind,2}=axes;
set(ax{testnameind,2},'pos',pos_sb{testnameind});
pc{testnameind,2}=m_pcolor(lon_RCM',lat_RCM', mean_RCM_early','parent',ax{testnameind,2});
colormap(ax{testnameind,2},jet_mod);
caxis(caxisval);
shading(gca,m_pcolor_shading_method);   

m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
        'xticklabels', [128, 132, 136, 140], 'xtick',[128, 132, 136, 140], 'XaxisLocation', 'bottom',...
        'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});

% % % contour
ax{testnameind,3}=axes;
set(ax{testnameind,3},'pos',pos_sb{testnameind});
[con_C,con_h]=m_contour(lon_RCM',lat_RCM', mean_RCM_early', [2, 5, 10], m_contour_color, 'linewidth', m_contour_linewidth, 'parent',ax{testnameind,3});
clabel(con_C,con_h,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
    'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
       'xticklabels', [], 'backcolor', 'none','parent', ax{testnameind,3});
caxis(caxisval);
colormap(ax{testnameind,3},jet_mod);  

txt{testnameind,1}=m_text(127.5, 48, 'Reanalysis', 'FontSize', m_grid_fontsize+4); 
txt{testnameind,1}=m_text(127.5, 46.5, '1983―1987', 'FontSize', m_grid_fontsize+4); 


% % projection44444444444444444
testnameind=4;
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)-0.5],'lat',[lonlat(3) lonlat(4)-0.5]);
sb2=subplot(4,5,[13 14 18 19]);  % Auto-fitted to the figure.
pos_sb{testnameind,1}=get(sb2, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig+ [000,0.23,0,0];
delete(sb2); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

% % land
temp_surf=ncread(ncname, 'temp_surf', [1 1 1], [inf inf 1]);
temp_surf(isnan(temp_surf))=50000;
temp_surf(temp_surf<50000)=NaN;
RCM_land=temp_surf;
RCM_land(temp_surf==50000)=1;
pc{testnameind,1}=m_pcolor(lon_RCM',lat_RCM', RCM_land','parent',ax{testnameind,1});
colormap(ax{testnameind,1},[0.8 0.8 0.8]);
shading(gca,m_pcolor_shading_method); 
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
       'xticklabels', [], 'yticklabels', [], 'parent', ax{testnameind,1});

% % % pcolor
hold on
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)-0.5],'lat',[lonlat(3) lonlat(4)-0.5]);
ax{testnameind,2}=axes;
set(ax{testnameind,2},'pos',pos_sb{testnameind});
pc{testnameind,2}=m_pcolor(lon_RCM',lat_RCM', mean_RCM_late','parent',ax{testnameind,2});
colormap(ax{testnameind,2},jet_mod);
caxis(caxisval);
shading(gca,m_pcolor_shading_method);   

m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
        'xticklabels', [128, 132, 136, 140], 'xtick',[128, 132, 136, 140], 'XaxisLocation', 'bottom',...
        'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});

% % % contour
ax{testnameind,3}=axes;
set(ax{testnameind,3},'pos',pos_sb{testnameind});
[con_C,con_h]=m_contour(lon_RCM',lat_RCM', mean_RCM_late', [2, 5, 10], m_contour_color, 'linewidth', m_contour_linewidth, 'parent',ax{testnameind,3});
clabel(con_C,con_h,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
    'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
       'yticklabels', [],'xticklabels', [], 'backcolor', 'none','parent', ax{testnameind,3});
caxis(caxisval);
colormap(ax{testnameind,3},jet_mod);    
    

txt{testnameind,1}=m_text(127.5, 48, 'Reanalysis', 'FontSize', m_grid_fontsize+4); 
txt{testnameind,1}=m_text(127.5, 46.5, '1988―1992', 'FontSize', m_grid_fontsize+4); 





% % colorbar
h = colorbar('eastoutside','AxisLocation','out');
caxis(caxisval);
set(h,'fontsize',m_grid_fontsize);
h_title=title(h,'^oC','fontsize',m_grid_fontsize-2);
% set(h, 'Position', [pos_sb{2}(1)+0.312, pos_sb{2}(2)+0.326, 0.0231, pos_sb{2}(3)-0.105])  % right, up, width, height
set(h, 'Position', [pos_sb{2}(1)+0.312, pos_sb{2}(2)-0.12, 0.0231, pos_sb{2}(3)+0.115])  % right, up, width, height

% % figure save
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height*2]) 
print('-dtiff','-r500',tifname); RemoveWhiteSpace([], 'file', tifname);
hold off;
close all;

save ('D:\research\Ph_D_course\2021_Particle tracking of the walleye pollock\data_dryad\Data02_7_SST.mat', ...
            'lon_OISST', 'lat_OISST', 'mean_OISST', 'OISST_land',...
            'lon_RCM', 'lat_RCM', 'mean_RCM', 'RCM_land')

