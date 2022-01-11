close all; clear all;  clc;
warning off;

testname = 'test06';
regionname = 'pollock_egg3';
earlyyear=[1983:1987];
lateyear=[1988:1992];
inputmonth = [1,2]; % % put month which you want to plot [month month ...]
allyear =[1983:1992];
refyear =[1983:1987];
checktime=[15,30];
caxisval=[0 15];
caxisval_diff=[-3 3];
addpath(genpath(['C:\Users\User\Dropbox\source\matlab\function\']))
[dropboxpath, erorr_status] = Func_0008_set_dropbox_path(computer);
addpath(genpath([dropboxpath, '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\MICT_pollack\paper\subroutine\']))
[byrmap3, error_status] = Func_0009_get_colormaps('byr3', dropboxpath);
[byrmap, error_status] = Func_0009_get_colormaps('byr2', dropboxpath);  
[refpolygon, lonlat, error_status] = Func_0007_get_polygon_data_from_regionname(regionname);

param_script =['C:\users\user/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_10/run/fig_param/fig_param_kyy_', regionname, '.m'];
figrawdir =strcat('D:\research\Ph_D_course\2021_Particle tracking of the walleye pollock\figure\paper\'); % % where figure files will be saved
filedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\'); % % where data files are
savedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\pollock_6\');
inputdir = ['/home/auto/MAMS/Data/01_NWP_1_10/Input/'];
LTRANS_testname='Pollock6';

tifname = [figrawdir, 'Fig04.tif'];

correction_right_fig=[-0.1600,0,0,0]; % right, up, width, height
correction_upper_fig=[0,0,0,0];
correction_large_fig=[0,0,0,0.020];

% % %  83-87 (early period, period I)
run(param_script);
m_grid_fontsize = m_grid_fontsize -4;
ind=1;
clear comb_data
for yearij = 1:length(earlyyear)
    tempyear = earlyyear(yearij);
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
mean_RCM_early = sum(comb_data,3)/size(comb_data,3);
mean_RCM_early(mean_RCM_early==0)=NaN;

% % projection
testnameind=1;
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
sb1=subplot(2,8,[2 3 10 11]);  % Auto-fitted to the figure.
pos_sb{testnameind,1}=get(sb1, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [-0.0800,0,0,0];
delete(sb1); % Delete the subplot axes
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
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type, ...
       'xticklabels', [128, 130, 132], 'xtick',[128, 130, 132], ...
       'backcolor', 'none','parent', ax{testnameind,1});

% % % pcolor
hold on
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
ax{testnameind,2}=axes;
set(ax{testnameind,2},'pos',pos_sb{testnameind});
pc{testnameind,2}=m_pcolor(lon_RCM',lat_RCM', mean_RCM_early','parent',ax{testnameind,2});
colormap(ax{testnameind,2},jet_mod);
caxis(caxisval);
shading(gca,m_pcolor_shading_method);   
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
        'xticklabels', [129, 131], 'xtick',[129, 131], 'XaxisLocation', 'top',...
        'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});

txt{testnameind,1}=m_text(127.15, 36.9, '(A) 83-87', 'FontSize', m_grid_fontsize-2); 

% % % contour
ax{testnameind,3}=axes;
set(ax{testnameind,3},'pos',pos_sb{testnameind});
[con_C,con_h]=m_contour(lon_RCM',lat_RCM', mean_RCM_early', [2, 5, 10], m_contour_color, 'linewidth', m_contour_linewidth-0.5, 'parent',ax{testnameind,3});
clabel(con_C,con_h,'FontSize',m_contour_label_fontsize-2,'Color',m_contour_label_color, ...
    'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type, 'yticklabels', [], ...
       'xticklabels', [], 'xtick',[], ...
       'backcolor', 'none','parent', ax{testnameind,3});
caxis(caxisval);
colormap(ax{testnameind,3},jet_mod);


% % %  88-92 (early period, period I)
run(param_script);
m_grid_fontsize = m_grid_fontsize -4;
ind=1;
clear comb_data
for yearij = 1:length(lateyear)
    tempyear = lateyear(yearij);
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
mean_RCM_late = sum(comb_data,3)/size(comb_data,3);
mean_RCM_late(mean_RCM_late==0)=NaN;

% % projection
testnameind=2;
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
sb2=subplot(2,8,[4 5 12 13]);  % Auto-fitted to the figure.
pos_sb{testnameind,1}=get(sb2, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [-0.0800,0,0,0];
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
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type, 'yticklabels', [], ...
       'xticklabels', [128, 130, 132], 'xtick',[128, 130, 132], ...
       'backcolor', 'none','parent', ax{testnameind,1});

% % % pcolor
hold on
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
ax{testnameind,2}=axes;
set(ax{testnameind,2},'pos',pos_sb{testnameind});
pc{testnameind,2}=m_pcolor(lon_RCM',lat_RCM', mean_RCM_late','parent',ax{testnameind,2});
colormap(ax{testnameind,2},jet_mod);
caxis(caxisval);
shading(gca,m_pcolor_shading_method);   
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
        'xticklabels', [129, 131], 'xtick',[129, 131], 'XaxisLocation', 'top',...
        'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});
txt{testnameind,1}=m_text(127.15, 36.9, '(B) 88-92', 'FontSize', m_grid_fontsize-2); 

% % % contour
ax{testnameind,3}=axes;
set(ax{testnameind,3},'pos',pos_sb{testnameind});
[con_C,con_h]=m_contour(lon_RCM',lat_RCM', mean_RCM_late', [2, 5, 10], m_contour_color, 'linewidth', m_contour_linewidth-0.5, 'parent',ax{testnameind,3});
clabel(con_C,con_h,'FontSize',m_contour_label_fontsize-2,'Color',m_contour_label_color, ...
    'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type, 'yticklabels', [], ...
       'xticklabels', [], 'xtick',[], ...
       'backcolor', 'none','parent', ax{testnameind,3});
caxis(caxisval);
colormap(ax{testnameind,3},jet_mod);

% % colorbar
h = colorbar('westoutside','AxisLocation','out');
caxis([caxisval]);
set(h,'fontsize',m_grid_fontsize+2);
h_title=title(h,'^oC','fontsize',m_grid_fontsize+2);
set(h, 'Position', [pos_sb{2}(1)-0.29, pos_sb{2}(2)+0.360, 0.0281, pos_sb{2}(3)-0.060])  % right, up, width, height

% % %  get difference (late-early period, period II-I)
mean_RCM_diff=mean_RCM_late-mean_RCM_early;

% % projection
testnameind=3;
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
sb2=subplot(2,8,[6 7 14 15]);  % Auto-fitted to the figure.
pos_sb{testnameind,1}=get(sb2, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [-0.0800,0,0,0];
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
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type, 'yticklabels', [], ...
       'xticklabels', [128, 130, 132], 'xtick',[128, 130, 132], ...
       'backcolor', 'none','parent', ax{testnameind,1});

% % % pcolor
hold on
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
ax{testnameind,2}=axes;
set(ax{testnameind,2},'pos',pos_sb{testnameind});
pc{testnameind,2}=m_pcolor(lon_RCM',lat_RCM', mean_RCM_diff','parent',ax{testnameind,2});
colormap(ax{testnameind,2},byrmap3);
caxis(caxisval_diff);
shading(gca,m_pcolor_shading_method);   
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
        'xticklabels', [129, 131], 'xtick',[129, 131], 'XaxisLocation', 'top',...
        'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});
% txt{testnameind,1}=m_text(127.15, 36.9, '(C) B-A', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=m_text(127.15, 36.9, '(C) (B)-(A)', 'FontSize', m_grid_fontsize-3); 


% % colorbar
h = colorbar('eastoutside','AxisLocation','out');
caxis(caxisval_diff);
set(h,'fontsize',m_grid_fontsize+2);
h_title=title(h,'^oC','fontsize',m_grid_fontsize+2);
set(h, 'Position', [pos_sb{3}(1)+0.192, pos_sb{3}(2)+0.360, 0.0231, pos_sb{2}(3)-0.060])  % right, up, width, height



% % figure save
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height*2]) 
print('-dtiff','-r500',tifname); RemoveWhiteSpace([], 'file', tifname);
hold off;
close all;

save ('D:\research\Ph_D_course\2021_particle tracking of the walleye pollock\data_dryad\Data09_SST.mat', ...
            'lon_RCM', 'lat_RCM', 'RCM_land', ...
            'mean_RCM_late', 'mean_RCM_diff', 'mean_RCM_early')

