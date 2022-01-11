close all; clear all;  clc;
warning off;

testname = 'test06';
regionname = 'pollock_egg3';
earlyyear=[1983:1987];
lateyear=[1988:1992];
inputmonth = [1,2]; % % put month which you want to plot [month month ...]
allyear =[1983:1992];
refyear =[1983:1987];
caxisval=[0 150];
caxisval_diff=[-70 70];
addpath(genpath(['C:\Users\User\Dropbox\source\matlab\function\']))
[dropboxpath, erorr_status] = Func_0008_set_dropbox_path(computer);
addpath(genpath([dropboxpath, '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\MICT_pollack\paper\subroutine\']))
[byrmap3, error_status] = Func_0009_get_colormaps('byr3', dropboxpath);
[byrmap, error_status] = Func_0009_get_colormaps('byr2', dropboxpath);  
[refpolygon, lonlat, error_status] = Func_0007_get_polygon_data_from_regionname(regionname);

param_script =['C:\users\user/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_10/run/fig_param/fig_param_kyy_', regionname, '.m'];
% figrawdir =strcat('Z:\내 드라이브\research\Ph_D_course\Particle tracking of the walleye pollock\figure\paper\'); % % where figure files will be saved
figrawdir =strcat('D:\research\Ph_D_course\2021_Particle tracking of the walleye pollock\figure\paper\'); % % where figure files will be saved
filedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\'); % % where data files are
savedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\pollock_6\');
inputdir = ['/home/auto/MAMS/Data/01_NWP_1_10/Input/'];
LTRANS_testname='Pollock6';

tifname = [figrawdir, 'Fig07_02.tif'];

correction_right_fig=[-0.1600,0,0,0]; % right, up, width, height
correction_upper_fig=[0,0,0,0];
correction_large_fig=[0,0,0,0.020];


temp_checktime=15;
% % %  83-87 (early period, period I)
run(param_script);
m_grid_fontsize = m_grid_fontsize -4;
ind=1;
clear egg_mask comb_egg_mask
for yearij = 1:length(earlyyear)
    tempyear = earlyyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
        all_checktime=ncread(ncname, 'checktime');
        ind_checktime=find(all_checktime==temp_checktime);
        egg_mask=squeeze(ncread(ncname, 'mask_par', [1 1 1 ind_checktime], [inf inf inf 1]));
        lastday_m=size(egg_mask,3);
        if (exist('comb_egg_mask')==0)
            comb_egg_mask=egg_mask;
        else
            comb_egg_mask(:,:,end+1:end+lastday_m)=egg_mask;
        end
    end
end
lon_RCM = ncread(ncname, 'lon_rho');
lat_RCM = ncread(ncname, 'lat_rho');
sum_early_egg_15d = sum(comb_egg_mask,3);

sum_early_egg_15d(sum_early_egg_15d==0)=NaN;
% [mean_early_egg, error_status] = Func_0011_get_area_weighted_mean(mean_early_egg, lon_rho, lat_rho);

% % projection
testnameind=1;
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
sb1=subplot(4,8,[2 3 10 11]);  % Auto-fitted to the figure.
% sb1=subplot(2,8,[2 3]);  % Auto-fitted to the figure.
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
       'xticklabels', [], 'xtick',[128, 130, 132], ...
       'backcolor', 'none','parent', ax{testnameind,1});

% % % pcolor
hold on
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
ax{testnameind,2}=axes;
set(ax{testnameind,2},'pos',pos_sb{testnameind});
pc{testnameind,2}=m_pcolor(lon_RCM',lat_RCM', sum_early_egg_15d','parent',ax{testnameind,2});
colormap(ax{testnameind,2},parula);
caxis(caxisval);
shading(gca,m_pcolor_shading_method);   
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
        'xticklabels', [129, 131], 'xtick',[129, 131], 'XaxisLocation', 'top',...
        'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});

txt{testnameind,1}=m_text(127.15, 36.9, '(A) 83-87', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=m_text(127.15, 40.75, '15 days', 'FontSize', m_grid_fontsize-2); 




% % %  88-92 (early period, period I)
run(param_script);
m_grid_fontsize = m_grid_fontsize -4;
ind=1;
clear egg_mask comb_egg_mask
for yearij = 1:length(lateyear)
    tempyear = lateyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
        all_checktime=ncread(ncname, 'checktime');
        ind_checktime=find(all_checktime==temp_checktime);
        egg_mask=squeeze(ncread(ncname, 'mask_par', [1 1 1 ind_checktime], [inf inf inf 1]));
        lastday_m=size(egg_mask,3);
        if (exist('comb_egg_mask')==0)
            comb_egg_mask=egg_mask;
        else
            comb_egg_mask(:,:,end+1:end+lastday_m)=egg_mask;
        end
    end
end
lon_RCM = ncread(ncname, 'lon_rho');
lat_RCM = ncread(ncname, 'lat_rho');
sum_late_egg_15d = sum(comb_egg_mask,3);

sum_late_egg_15d(sum_late_egg_15d==0)=NaN;
% [mean_early_egg, error_status] = Func_0011_get_area_weighted_mean(mean_early_egg, lon_rho, lat_rho);

% % projection
testnameind=2;
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
sb2=subplot(4,8,[4 5 12 13]);  % Auto-fitted to the figure.
% sb2=subplot(2,8,[4 5]);  % Auto-fitted to the figure.

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
       'xticklabels', [], 'xtick',[128, 130, 132], ...
       'backcolor', 'none','parent', ax{testnameind,1});

% % % pcolor
hold on
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
ax{testnameind,2}=axes;
set(ax{testnameind,2},'pos',pos_sb{testnameind});
pc{testnameind,2}=m_pcolor(lon_RCM',lat_RCM', sum_late_egg_15d','parent',ax{testnameind,2});
colormap(ax{testnameind,2},parula);
caxis(caxisval);
shading(gca,m_pcolor_shading_method);   
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
        'xticklabels', [129, 131], 'xtick',[129, 131], 'XaxisLocation', 'top',...
        'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});
txt{testnameind,1}=m_text(127.15, 36.9, '(B) 88-92', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=m_text(127.15, 40.75, '15 days', 'FontSize', m_grid_fontsize-2); 

% % %  get difference (late-early period, period II-I)
sum_diff_egg_15d=sum_late_egg_15d-sum_early_egg_15d;

% % projection
testnameind=3;
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
sb2=subplot(4,8,[6 7 14 15]);  % Auto-fitted to the figure.
% sb2=subplot(2,8,[6 7]);  % Auto-fitted to the figure.

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
       'xticklabels', [], 'xtick',[128, 130, 132], ...
       'backcolor', 'none','parent', ax{testnameind,1});

% % % pcolor
hold on
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
ax{testnameind,2}=axes;
set(ax{testnameind,2},'pos',pos_sb{testnameind});
pc{testnameind,2}=m_pcolor(lon_RCM',lat_RCM', sum_diff_egg_15d','parent',ax{testnameind,2});
colormap(ax{testnameind,2},byrmap3);
caxis(caxisval_diff);
shading(gca,m_pcolor_shading_method);   
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
        'xticklabels', [129, 131], 'xtick',[129, 131], 'XaxisLocation', 'top',...
         'YaxisLocation', 'right','backcolor', 'none','parent', ax{testnameind,2});
% txt{testnameind,1}=m_text(127.15, 36.9, '(C) B-A', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=m_text(127.15, 36.9, '(C) (B)-(A)', 'FontSize', m_grid_fontsize-3); 
txt{testnameind,1}=m_text(127.15, 40.75, '15 days', 'FontSize', m_grid_fontsize-2); 


% % % % % % % % % % 
% % % % % % % % % %  row 2
% % % % % % % % % %

correction_upper_fig = [0 0.31 0 0];

temp_checktime=30;
% % %  83-87 (early period, period I)
run(param_script);
m_grid_fontsize = m_grid_fontsize -4;
ind=1;
clear egg_mask comb_egg_mask
for yearij = 1:length(earlyyear)
    tempyear = earlyyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
        all_checktime=ncread(ncname, 'checktime');
        ind_checktime=find(all_checktime==temp_checktime);
        egg_mask=squeeze(ncread(ncname, 'mask_par', [1 1 1 ind_checktime], [inf inf inf 1]));
        lastday_m=size(egg_mask,3);
        if (exist('comb_egg_mask')==0)
            comb_egg_mask=egg_mask;
        else
            comb_egg_mask(:,:,end+1:end+lastday_m)=egg_mask;
        end
    end
end
lon_RCM = ncread(ncname, 'lon_rho');
lat_RCM = ncread(ncname, 'lat_rho');
sum_early_egg_30d = sum(comb_egg_mask,3);

sum_early_egg_30d(sum_early_egg_30d==0)=NaN;
% [mean_early_egg, error_status] = Func_0011_get_area_weighted_mean(mean_early_egg, lon_rho, lat_rho);

% % projection
testnameind=4;
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
sb1=subplot(4,8,[18 19 26 27]);  % Auto-fitted to the figure.
% sb1=subplot(2,8,[10 11]);  % Auto-fitted to the figure.

pos_sb{testnameind,1}=get(sb1, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [-0.0800,0,0,0]; % right, up, width, height
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
pc{testnameind,2}=m_pcolor(lon_RCM',lat_RCM', sum_early_egg_30d','parent',ax{testnameind,2});
colormap(ax{testnameind,2},parula);
caxis(caxisval);
shading(gca,m_pcolor_shading_method);   
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
        'xticklabels', [], 'xtick',[129, 131], 'XaxisLocation', 'top',...
        'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});

txt{testnameind,1}=m_text(127.15, 36.9, '(D) 83-87', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=m_text(127.15, 40.75, '30 days', 'FontSize', m_grid_fontsize-2); 




% % %  88-92 (early period, period I)
run(param_script);
m_grid_fontsize = m_grid_fontsize -4;
ind=1;
clear egg_mask comb_egg_mask
for yearij = 1:length(lateyear)
    tempyear = lateyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
        all_checktime=ncread(ncname, 'checktime');
        ind_checktime=find(all_checktime==temp_checktime);
        egg_mask=squeeze(ncread(ncname, 'mask_par', [1 1 1 ind_checktime], [inf inf inf 1]));
        lastday_m=size(egg_mask,3);
        if (exist('comb_egg_mask')==0)
            comb_egg_mask=egg_mask;
        else
            comb_egg_mask(:,:,end+1:end+lastday_m)=egg_mask;
        end
    end
end
lon_RCM = ncread(ncname, 'lon_rho');
lat_RCM = ncread(ncname, 'lat_rho');
sum_late_egg_30d = sum(comb_egg_mask,3);

sum_late_egg_30d(sum_late_egg_30d==0)=NaN;
% [mean_early_egg, error_status] = Func_0011_get_area_weighted_mean(mean_early_egg, lon_rho, lat_rho);

% % projection
testnameind=5;
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
sb2=subplot(4,8,[20 21 28 29]);  % Auto-fitted to the figure.
% sb2=subplot(2,8,[12 13]);  % Auto-fitted to the figure.

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
pc{testnameind,2}=m_pcolor(lon_RCM',lat_RCM', sum_late_egg_30d','parent',ax{testnameind,2});
colormap(ax{testnameind,2},parula);
caxis(caxisval);
shading(gca,m_pcolor_shading_method);   
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
        'xticklabels', [], 'xtick',[129, 131], 'XaxisLocation', 'top',...
        'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});
txt{testnameind,1}=m_text(127.15, 36.9, '(E) 88-92', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=m_text(127.15, 40.75, '30 days', 'FontSize', m_grid_fontsize-2); 

% % colorbar (# of particles)
h = colorbar('westoutside','AxisLocation','out');
caxis([caxisval]);
set(h,'fontsize',m_grid_fontsize+2);
set(h, 'Position', [pos_sb{2}(1)-0.29, pos_sb{2}(2)+0.002, 0.0281, pos_sb{2}(3)+0.080])  % right, up, width, height
h_title=title(h,'NTI','fontsize',m_grid_fontsize+2);



% % %  get difference (late-early period, period II-I)
sum_diff_egg_30d=sum_late_egg_30d-sum_early_egg_30d;

% % projection
testnameind=6;
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
sb2=subplot(4,8,[22 23 30 31]);  % Auto-fitted to the figure.
% sb2=subplot(2,8,[14 15]);  % Auto-fitted to the figure.


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
pc{testnameind,2}=m_pcolor(lon_RCM',lat_RCM', sum_diff_egg_30d','parent',ax{testnameind,2});
colormap(ax{testnameind,2},byrmap3);
caxis(caxisval_diff);
shading(gca,m_pcolor_shading_method);   
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
        'xticklabels', [], 'xtick',[129, 131], 'XaxisLocation', 'top',...
         'YaxisLocation', 'right','backcolor', 'none','parent', ax{testnameind,2});
% txt{testnameind,1}=m_text(127.15, 36.9, '(F) E-D', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=m_text(127.15, 36.9, '(F) (E)-(D)', 'FontSize', m_grid_fontsize-3); 
txt{testnameind,1}=m_text(127.15, 40.75, '30 days', 'FontSize', m_grid_fontsize-2); 





% % colorbar  (difference)
h = colorbar('eastoutside','AxisLocation','out');
caxis(caxisval_diff);
set(h,'fontsize',m_grid_fontsize+2);
set(h, 'Position', [pos_sb{3}(1)+0.242, pos_sb{3}(2)+0.002, 0.0231, pos_sb{2}(3)+0.080])  % right, up, width, height
h_title=title(h,'NTI','fontsize',m_grid_fontsize+2);



% % figure save
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y/3]);
set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height*2]) 
print('-dtiff','-r500',tifname); 
RemoveWhiteSpace([], 'file', tifname);
hold off;
close all;

save ('D:\research\Ph_D_course\2021_Particle tracking of the walleye pollock\data_dryad\Data05_numpar.mat', ...
            'lon_RCM', 'lat_RCM', 'RCM_land', ...
            'sum_early_egg_15d', 'sum_late_egg_15d', 'sum_diff_egg_15d', ...
            'sum_early_egg_30d', 'sum_late_egg_30d', 'sum_diff_egg_30d')

