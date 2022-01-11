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

tifname = [figrawdir, 'Fig08.tif'];

correction_right_fig=[-0.1600,0,0,0]; % right, up, width, height
correction_upper_fig=[0,0,0,0]; % right, up, width, height
correction_large_fig=[0,0,0,-0.28]; % right, up, width, height


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
lon_rho = ncread(ncname, 'lon_rho');
lat_rho = ncread(ncname, 'lat_rho');
sum_data = sum(comb_egg_mask,3, 'omitnan');
sum_data(sum_data==0)=NaN;
testnameind=1;
lon_sum_early_data_15d=sum(sum_data, 1, 'omitnan');
lat_1d=mean(lat_rho,1,'omitnan');

% % projection
testnameind=1;
sb1=subplot(4,8,[2 3 10 11]);  % Auto-fitted to the figure.
% sb1=subplot(2,8,[2 3]);  % Auto-fitted to the figure.
pos_sb{testnameind,1}=get(sb1, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [-0.0800,0,0,0];
delete(sb1); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

barh(lat_1d, lon_sum_early_data_15d)
    ylabel('latitude (degree)')
    xlim([0 1500])
    ylim([36 41]); grid minor;
    set(ax{testnameind,1}, 'fontsize', m_grid_fontsize) 
txt{testnameind,1}=text(127.15, 36.5, '(A) 83-87', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=text(127.15, 37.1, '15 days', 'FontSize', m_grid_fontsize-2); 

set(ax{testnameind,1}, 'XAxisLocation', 'top')
set(ax{testnameind,1}, 'XTick',[0 500 1000])
set(ax{testnameind,1}, 'XTicklabel', {'0', '500', '1000'})


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
lon_rho = ncread(ncname, 'lon_rho');
lat_rho = ncread(ncname, 'lat_rho');
sum_data = sum(comb_egg_mask,3, 'omitnan');
sum_data(sum_data==0)=NaN;
testnameind=1;
lon_sum_late_data_15d=sum(sum_data, 1, 'omitnan');
lat_1d=mean(lat_rho,1,'omitnan');

% % projection
testnameind=2;
sb1=subplot(4,8,[4 5 12 13]);  % Auto-fitted to the figure.
% sb1=subplot(2,8,[2 3]);  % Auto-fitted to the figure.
pos_sb{testnameind,1}=get(sb1, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [-0.0800,0,0,0];
delete(sb1); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

barh(lat_1d, lon_sum_late_data_15d)
    xlim([0 1500])
    ylim([36 41]); grid minor;
    set(ax{testnameind,1}, 'fontsize', m_grid_fontsize) 
txt{testnameind,1}=text(127.15, 36.5, '(B) 88-92', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=text(127.15, 37.1, '15 days', 'FontSize', m_grid_fontsize-2); 

set(ax{testnameind,1}, 'XAxisLocation', 'top')
set(ax{testnameind,1}, 'XTick',[0 500 1000])
set(ax{testnameind,1}, 'XTicklabel', {'0', '500', '1000'})
set(ax{testnameind,1}, 'YTicklabel', [])

% % %  get difference (late-early period, period II-I)
lon_sum_diff_data_15d=lon_sum_late_data_15d-lon_sum_early_data_15d;

% % projection
testnameind=3;
sb2=subplot(4,8,[6 7 14 15]);  % Auto-fitted to the figure.
% sb2=subplot(2,8,[6 7]);  % Auto-fitted to the figure.

pos_sb{testnameind,1}=get(sb2, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [-0.0800,0,0,0];
delete(sb2); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

barh(lat_1d, lon_sum_diff_data_15d, 'r')
    xlim([-1000 500])
    ylim([36 41]); grid minor;
    set(ax{testnameind,1}, 'fontsize', m_grid_fontsize) 
% txt{testnameind,1}=text(-872.85, 36.5, '(C) B-A', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=text(-872.85, 36.5, '(C) (B)-(A)', 'FontSize', m_grid_fontsize-3); 
txt{testnameind,1}=text(-872.85, 37.1, '15 days', 'FontSize', m_grid_fontsize-2); 

set(ax{testnameind,1}, 'XAxisLocation', 'top')
set(ax{testnameind,1}, 'XTick',[-500 0 500])
set(ax{testnameind,1}, 'XTicklabel', {'-500', '0', '500'})
set(ax{testnameind,1}, 'YAxisLocation', 'right')


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
lon_rho = ncread(ncname, 'lon_rho');
lat_rho = ncread(ncname, 'lat_rho');
sum_data = sum(comb_egg_mask,3, 'omitnan');
sum_data(sum_data==0)=NaN;
testnameind=1;
lon_sum_early_data_30d=sum(sum_data, 1, 'omitnan');
lat_1d=mean(lat_rho,1,'omitnan');

% % projection
testnameind=4;
sb1=subplot(4,8,[18 19 26 27]);  % Auto-fitted to the figure.
% sb1=subplot(2,8,[10 11]);  % Auto-fitted to the figure.

pos_sb{testnameind,1}=get(sb1, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [-0.0800,0,0,0]; % right, up, width, height
delete(sb1); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

barh(lat_1d, lon_sum_early_data_30d)
    ylabel('latitude (degree)')
    xlabel('# of individuals')
    xlim([0 1500])
    ylim([36 41]); grid minor;
    set(ax{testnameind,1}, 'fontsize', m_grid_fontsize) 
txt{testnameind,1}=text(127.15, 36.5, '(D) 83-87', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=text(127.15, 37.1, '30 days', 'FontSize', m_grid_fontsize-2); 
set(ax{testnameind,1}, 'XTick',[0 500 1000])
set(ax{testnameind,1}, 'XTicklabel', {'0', '500', '1000'})




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
lon_rho = ncread(ncname, 'lon_rho');
lat_rho = ncread(ncname, 'lat_rho');
sum_data = sum(comb_egg_mask,3, 'omitnan');
sum_data(sum_data==0)=NaN;
testnameind=1;
lon_sum_late_data_30d=sum(sum_data, 1, 'omitnan');
lat_1d=mean(lat_rho,1,'omitnan');

% % projection
testnameind=5;
sb2=subplot(4,8,[20 21 28 29]);  % Auto-fitted to the figure.
% sb2=subplot(2,8,[12 13]);  % Auto-fitted to the figure.

pos_sb{testnameind,1}=get(sb2, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [-0.0800,0,0,0];
delete(sb2); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

barh(lat_1d, lon_sum_late_data_30d)
    xlabel('# of individuals')
    xlim([0 1500])
    ylim([36 41]); grid minor;
    set(ax{testnameind,1}, 'fontsize', m_grid_fontsize) 
txt{testnameind,1}=text(127.15, 36.5, '(E) 88-92', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=text(127.15, 37.1, '30 days', 'FontSize', m_grid_fontsize-2); 
set(ax{testnameind,1}, 'XTick',[0 500 1000])
set(ax{testnameind,1}, 'XTicklabel', {'0', '500', '1000'})
set(ax{testnameind,1}, 'YTicklabel', [])

% % %  get difference (late-early period, period II-I)
lon_sum_diff_data_30d=lon_sum_late_data_30d-lon_sum_early_data_30d;

% % projection
testnameind=6;
sb2=subplot(4,8,[22 23 30 31]);  % Auto-fitted to the figure.
% sb2=subplot(2,8,[14 15]);  % Auto-fitted to the figure.


pos_sb{testnameind,1}=get(sb2, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [-0.0800,0,0,0];
delete(sb2); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

barh(lat_1d, lon_sum_diff_data_30d, 'r')
    xlabel('# of individuals')
    xlim([-1000 500])
    ylim([36 41]); grid minor;
    set(ax{testnameind,1}, 'fontsize', m_grid_fontsize)
    set(ax{testnameind,1}, 'XTick',[-500 0 500])
    set(ax{testnameind,1}, 'XTicklabel', {'-500', '0', '500'})
    set(ax{testnameind,1}, 'YAxisLocation', 'right')

% txt{testnameind,1}=text(-872.85, 36.5, '(F) E-D', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=text(-872.85, 36.5, '(F) (E)-(D)', 'FontSize', m_grid_fontsize-3); 
txt{testnameind,1}=text(-872.85, 37.1, '30 days', 'FontSize', m_grid_fontsize-2); 



% % figure save
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y/3]);
set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height*2]) 
print('-dtiff','-r500',tifname); 
RemoveWhiteSpace([], 'file', tifname);
hold off;
close all;

save ('D:\research\Ph_D_course\2021_Particle tracking of the walleye pollock\data_dryad\Data06_numpar_lon.mat', ...
            'lat_1d', ...
            'lon_sum_early_data_15d', 'lon_sum_late_data_15d', 'lon_sum_diff_data_15d', ...
            'lon_sum_early_data_30d', 'lon_sum_late_data_30d', 'lon_sum_diff_data_30d')

