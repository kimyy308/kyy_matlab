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
figrawdir =strcat('D:\research\Ph_D_course\2021_Particle tracking of the walleye pollock\figure\paper\'); % % where figure files will be saved
filedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\'); % % where data files are
savedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\pollock_6\');
inputdir = ['/home/auto/MAMS/Data/01_NWP_1_10/Input/'];
LTRANS_testname='Pollock6';

tifname = [figrawdir, 'Fig07_04.tif'];

correction_right_fig=[-0.1600,0,0,0]; % right, up, width, height
correction_upper_fig=[0,0,0,0];
correction_large_fig=[0,0,0,0.020];


% % % positioning of the figures
figs.l_scale = 1;
figs.l_unit = 'km';

figs.ylimit = [-42.4e3 inf];
figs.zlimit = [-4000 20];

figs.n_rows = 2;
figs.n_columns = 3;
figs.ax_width = 1.8633;
figs.ax_height = 1.0250;
figs.xlabel_width = 0.25;
figs.xtick_width= 0.23;
figs.right_margin = 0.07;
figs.ylabel_height = 0.24;
figs.ytick_height= 0.25;
figs.top_margin = 0.21;
figs.fig_pos  = [0,0,figs.xlabel_width+figs.n_columns*(figs.xtick_width+figs.ax_width+figs.right_margin), ...
    figs.ylabel_height+figs.n_rows*(figs.ytick_height+figs.ax_height+figs.top_margin)];
figs.fig_title = 'Figure_7';
figs.exp_title = {{'Open-Ocean'},{'Shelf &','Control','Slope'},{'No Shelf'},{'Shelf &','Gentle','Slope'},{'Shelf &','Steep','Slope'},{'Shelf &','Cliff'}};
figs.fig = figure('name',figs.fig_title,'PaperUnits','inches', ...
    'paperposition',figs.fig_pos,'position',figs.fig_pos*get(groot,'ScreenPixelsPerInch'),'visible','off');

figs.text_config.text_period.xpos=0.03;
figs.text_config.text_period.ypos=0.97;
figs.text_legend.text_period.xpos=0.02;
figs.text_legend.text_period.ypos=0.98;

for mm = 1:figs.n_rows*figs.n_columns
%     exp_num_str = num2str(exp_num(mm),'%02d');

    ii=ceil(mm/figs.n_columns);
    jj=mod(mm-1,figs.n_columns)+1;
    figs.subplot(mm).posAx = ...
        [figs.xlabel_width+figs.xtick_width+(jj-1)*(figs.ax_width+figs.right_margin+figs.xtick_width), ...
        figs.ylabel_height+figs.ytick_height+(figs.n_rows-ii)*(figs.ax_height+figs.top_margin+figs.ytick_height), ...
        figs.ax_width, ...
        figs.ax_height];
    figs.subplot(mm).ax = axes(figs.fig, 'units', 'inches', 'position',figs.subplot(mm).posAx);
    xlim(figs.subplot(mm).ax,figs.ylimit/figs.l_scale); 
    ylim(figs.subplot(mm).ax,figs.zlimit/figs.l_scale);
    figs.subplot(mm).text_period=text(figs.subplot(mm).ax, ...
        figs.text_config.text_period.xpos,figs.text_config.text_period.ypos, ...
        figs.exp_title{mm},'verticalalignment','top','horizontalalignment','left','units','normalized','interpreter','tex');
    figs.subplot(mm).text_legend=text(figs.subplot(mm).ax, ...
        figs.text_legend.text_period.xpos,figs.text_legend.text_period.ypos, ...
        ['\bf ',char(double('A')-1+mm)], ...
        'verticalalignment','bottom','horizontalalignment','left','units','normalized','fontsize',12);
    switch ii, case figs.n_rows, xlabel(figs.subplot(mm).ax,['$$y$$ (',figs.l_unit,')'],'interpreter','latex'); end
    switch jj, case 1,      ylabel(figs.subplot(mm).ax,['z (',figs.l_unit,')']); end
end

% l1 = annotation('line'); l1.Units = 'inches'; l1.X=[1,1]*0.25;	l1.Y=[0,1.5]; l1.LineWidth = 0.25; l1.Color = 'r';
% l2 = annotation('line'); l2.Units = 'inches'; l2.X=[1,1]*0.48;	l2.Y=[0,1.5]; l2.LineWidth = 0.25; l2.Color = 'r';
% l3 = annotation('line'); l3.Units = 'inches'; l3.X=[1,1]*(0.25+2.1633);	l3.Y=[0,1.5]; l3.LineWidth = 0.25; l3.Color = 'r';
% l4 = annotation('line'); l4.Units = 'inches'; l4.Y=[1,1]*ylab_h;	l4.X=[0,4.5]; l4.LineWidth = 0.25; l4.Color = 'r';
% l5 = annotation('line'); l5.Units = 'inches'; l5.Y=[1,1]*(ylab_h+ytick_h);	l5.X=[0,4.5]; l5.LineWidth = 0.25; l5.Color = 'r';
% l6 = annotation('line'); l6.Units = 'inches'; l6.Y=[1,1]*(ylab_h+ytick_h+ax_h+t_margin);	l6.X=[0,4.5]; l6.LineWidth = 0.25; l6.Color = 'r';

print(figs.fig,figs.fig_title,'-dtiff','-r300');
% print(figs.fig,figs.fig_title,'-depsc');



temp_checktime=15;
% % %  83-87 (early period, period I)
run(param_script);
m_grid_fontsize = m_grid_fontsize -8;
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
lon_sum_early_data_15d=sum(sum_early_egg_15d, 1, 'omitnan');
lat_1d=mean(lat_RCM,1,'omitnan');

% [mean_early_egg, error_status] = Func_0011_get_area_weighted_mean(mean_early_egg, lon_rho, lat_rho);

% % projection
testnameind=1;
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
sb1=subplot(6,12,[3 4 15 16]);  % Auto-fitted to the figure.
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
        'xticklabels', [128 130 132], 'xtick',[128 130 132], 'XaxisLocation', 'top',...
        'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});

txt{testnameind,1}=m_text(127.15, 36.9, '1983！1987', 'FontSize', m_grid_fontsize-4); 
txt{testnameind,1}=m_text(127.15, 40.75, '15 days', 'FontSize', m_grid_fontsize-2); 

% % barh
sb1=subplot(6,12,[5 17]);  % Auto-fitted to the figure.
pos_sb{testnameind,3}=get(sb1, 'pos'); % Get the position.
pos_sb{testnameind,3} = pos_sb{testnameind,3} + correction_upper_fig+ correction_large_fig + [-0.0800,0.2,0,-0.2];
delete(sb1); % Delete the subplot axes
ax{testnameind,3}=axes;
set(ax{testnameind,3},'pos',pos_sb{testnameind,3});

barh(lat_1d, lon_sum_early_data_15d)
%     ylabel('latitude (degree)')
xlim([0 1500])
ylim([36 41]); grid minor;
set(ax{testnameind,1}, 'fontsize', m_grid_fontsize) 
txt{testnameind,1}=text(127.15, 36.5, '1983！1987', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=text(127.15, 37.1, '15 days', 'FontSize', m_grid_fontsize-2); 
set(ax{testnameind,1}, 'XAxisLocation', 'top')
set(ax{testnameind,1}, 'XTick',[0 500 1000])
set(ax{testnameind,1}, 'XTicklabel', {})
set(ax{testnameind,1}, 'YTick',[36 38 40])
set(ax{testnameind,1}, 'YTicklabel', {'36^oN', '38^oN', '40^oN'})


% % %  88-92 (early period, period I)
run(param_script);
m_grid_fontsize = m_grid_fontsize -8;
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
sb2=subplot(6,12,[6 7 18 19]);  % Auto-fitted to the figure.
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
        'xticklabels', [128 130 132], 'xtick',[128 130 132], 'XaxisLocation', 'top',...
        'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});
txt{testnameind,1}=m_text(127.15, 36.9, '1988！1992', 'FontSize', m_grid_fontsize-4); 
txt{testnameind,1}=m_text(127.15, 40.75, '15 days', 'FontSize', m_grid_fontsize-2); 

% % %  get difference (late-early period, period II-I)
sum_diff_egg_15d=sum_late_egg_15d-sum_early_egg_15d;



% % projection
testnameind=3;
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
sb2=subplot(6,12,[9 10 21 22]);  % Auto-fitted to the figure.
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
        'xticklabels', [128 130 132], 'xtick',[128 130 132], 'XaxisLocation', 'top',...
         'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});
% txt{testnameind,1}=m_text(127.15, 36.9, '(C) B-A', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=m_text(127.15, 36.9, '(B)-(A)', 'FontSize', m_grid_fontsize-3); 
txt{testnameind,1}=m_text(127.15, 40.75, '15 days', 'FontSize', m_grid_fontsize-2); 


% % % % % % % % % % 
% % % % % % % % % %  row 2
% % % % % % % % % %

correction_upper_fig = [0 0.31 0 0];

temp_checktime=30;
% % %  83-87 (early period, period I)
run(param_script);
m_grid_fontsize = m_grid_fontsize -8;
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
sb1=subplot(6,12,[39 40 51 52]);  % Auto-fitted to the figure.
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
       'xticklabels', [], 'xtick',[], ...
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
        'xticklabels', [], 'xtick',[128 130 132], 'XaxisLocation', 'top',...
        'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});

txt{testnameind,1}=m_text(127.15, 36.9, '1983！1987', 'FontSize', m_grid_fontsize-4); 
txt{testnameind,1}=m_text(127.15, 40.75, '30 days', 'FontSize', m_grid_fontsize-2); 




% % %  88-92 (early period, period I)
run(param_script);
m_grid_fontsize = m_grid_fontsize -8;
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
sb2=subplot(6,12,[42 43 54 55]);  % Auto-fitted to the figure.
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
       'xticklabels', [], 'xtick',[], ...
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
        'xticklabels', [], 'xtick',[128 130 132], 'XaxisLocation', 'top',...
        'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});
txt{testnameind,1}=m_text(127.15, 36.9, '1988！1992', 'FontSize', m_grid_fontsize-4); 
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
sb2=subplot(6,12,[45 46 57 58]);  % Auto-fitted to the figure.
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
       'xticklabels', [], 'xtick',[], ...
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
        'xticklabels', [], 'xtick',[128 130 132], 'XaxisLocation', 'top',...
         'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});
% txt{testnameind,1}=m_text(127.15, 36.9, '(F) E-D', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=m_text(127.15, 36.9, '(E)-(D)', 'FontSize', m_grid_fontsize-3); 
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

