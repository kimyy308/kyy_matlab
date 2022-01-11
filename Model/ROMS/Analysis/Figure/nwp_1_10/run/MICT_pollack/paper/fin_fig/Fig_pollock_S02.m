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

tifname = [figrawdir, 'FigS02.tif'];

correction_right_fig=[-0.1600,0,0,0]; % right, up, width, height
correction_upper_fig=[0,0,0,0]; % right, up, width, height
correction_large_fig=[0,0,0,-0.15]; % right, up, width, height

dist_value = m_lldist([128 129], [38, 38]);

% % % 
% % %  row 1
% % % 
temp_checktime=15;
run(param_script);
m_grid_fontsize = m_grid_fontsize -2;
clear wind_curl comb_wind_curl u_rho comb_u_rho v_rho comb_v_rho uwind comb_uwind vwind comb_vwind ...
    egg_mask comb_egg_mask sp_ground comb_sp_ground wind_curl2 comb_wind_curl2 temp_surf comb_temp_surf curl_mask ...
    comb_mov_dist_lon_mean comb_mov_dist_lat_mean
for yearij = 1:length(allyear)
    tempyear = allyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])  

        all_checktime=ncread(ncname, 'checktime');
        ind_checktime=find(all_checktime==temp_checktime);
        ocean_time=ncread(ncname, 'time')+datenum(1900,12,31);
        mov_dist_lon_mean = ncread(ncname, 'mov_dist_lon_mean', [1 ind_checktime], [inf 1]);
        mov_dist_lat_mean = ncread(ncname, 'mov_dist_lat_mean', [1 ind_checktime], [inf 1]);
        lastday_m=length(mov_dist_lon_mean);
        if (exist('comb_mov_dist_lon_mean')==0)
            comb_ocean_time=ocean_time;
            comb_mov_dist_lon_mean= mov_dist_lon_mean;
            comb_mov_dist_lat_mean= mov_dist_lat_mean;
        else
            comb_ocean_time(end+1:end+lastday_m)=ocean_time;
            comb_mov_dist_lon_mean(end+1:end+lastday_m)=mov_dist_lon_mean;
            comb_mov_dist_lat_mean(end+1:end+lastday_m)=mov_dist_lat_mean;
        end
    end
end

half_len=round(length(comb_mov_dist_lon_mean)/2);
egg_half_1=mean(comb_mov_dist_lon_mean(1:half_len), 'omitnan');
egg_half_2=mean(comb_mov_dist_lon_mean(half_len+1:end), 'omitnan');
regime_ts_egg_mask(1:half_len)=egg_half_1;
regime_ts_egg_mask(half_len+1:length(comb_mov_dist_lon_mean))=egg_half_2;

% % projection
testnameind=1;
sb1=subplot(2,2,[1 2]);  % Auto-fitted to the figure.
pos_sb{testnameind,1}=get(sb1, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [0,0,0,0];
delete(sb1); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

mslplot{testnameind,2}=bar(comb_mov_dist_lon_mean.*dist_value, 'FaceColor', [0.6 0.6 0.6], 'parent',ax{testnameind,1});
hold on
% mslplot{testnameind,1}=plot(regime_ts_egg_mask.*dist_value,'k','parent',ax{testnameind,1});
ylabel(ax{testnameind,1},'Distance (km)')

cal_arr=1;
for cali=1:length(allyear)-1
    cal_arr=[cal_arr, cal_arr(end)+sum(eomday(allyear(cali),inputmonth))];
end
% set(ax{testnameind,1},'color','none','yaxislocation','left','xtick', cal_arr,'xticklabel',datestr(comb_ocean_time(cal_arr), 'yyyy'),'position', ax_pos+[0 0.02 -0.01 -0.02]);
set(ax{testnameind,1},'color','none','yaxislocation','left','xtick', cal_arr,'xticklabel',datestr(comb_ocean_time(cal_arr), 'yyyy'));

axis tight 
ylim(ax{testnameind,1}, [0,3].*dist_value)
set(ax{testnameind,1},'ycolor','k', 'box', 'off', 'FontSize',m_grid_fontsize);
xlabel(ax{testnameind,1}, 'Year');

set(mslplot{testnameind,1},'LineWidth',2);
grid on
% lgd=legend([mslplot{testnameind,1}], 'moved distance of individuals(lat)');
% set(lgd,'FontSize',m_grid_fontsize);
% set(lgd,'Orientation','horizontal');
% set(lgd,'Location','Northwest');

txt{testnameind,1}=text(480, 2.2.*dist_value, '(A) 15 days', 'FontSize', m_grid_fontsize); 


% % % % % % % % % % 
% % % % % % % % % %  row 2
% % % % % % % % % %
temp_checktime=30;
run(param_script);
m_grid_fontsize = m_grid_fontsize -2;
clear wind_curl comb_wind_curl u_rho comb_u_rho v_rho comb_v_rho uwind comb_uwind vwind comb_vwind ...
    egg_mask comb_egg_mask sp_ground comb_sp_ground wind_curl2 comb_wind_curl2 temp_surf comb_temp_surf curl_mask ...
    comb_mov_dist_lon_mean comb_mov_dist_lat_mean
for yearij = 1:length(allyear)
    tempyear = allyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])  

        all_checktime=ncread(ncname, 'checktime');
        ind_checktime=find(all_checktime==temp_checktime);
        ocean_time=ncread(ncname, 'time')+datenum(1900,12,31);
        mov_dist_lon_mean = ncread(ncname, 'mov_dist_lon_mean', [1 ind_checktime], [inf 1]);
        mov_dist_lat_mean = ncread(ncname, 'mov_dist_lat_mean', [1 ind_checktime], [inf 1]);
        lastday_m=length(mov_dist_lon_mean);
        if (exist('comb_mov_dist_lon_mean')==0)
            comb_ocean_time=ocean_time;
            comb_mov_dist_lon_mean= mov_dist_lon_mean;
            comb_mov_dist_lat_mean= mov_dist_lat_mean;
        else
            comb_ocean_time(end+1:end+lastday_m)=ocean_time;
            comb_mov_dist_lon_mean(end+1:end+lastday_m)=mov_dist_lon_mean;
            comb_mov_dist_lat_mean(end+1:end+lastday_m)=mov_dist_lat_mean;
        end
    end
end
half_len=round(length(comb_mov_dist_lon_mean)/2);
egg_half_1=mean(comb_mov_dist_lon_mean(1:half_len), 'omitnan');
egg_half_2=mean(comb_mov_dist_lon_mean(half_len+1:end), 'omitnan');
regime_ts_egg_mask(1:half_len)=egg_half_1;
regime_ts_egg_mask(half_len+1:length(comb_mov_dist_lon_mean))=egg_half_2;

% % projection
testnameind=2;
sb1=subplot(2,2,[3 4]);  % Auto-fitted to the figure.
pos_sb{testnameind,1}=get(sb1, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [0,0.23,0,0];
delete(sb1); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

mslplot{testnameind,2}=bar(comb_mov_dist_lon_mean.*dist_value, 'FaceColor', [0.6 0.6 0.6], 'parent',ax{testnameind,1});
hold on
% mslplot{testnameind,1}=plot(regime_ts_egg_mask.*dist_value,'k','parent',ax{testnameind,1});
ylabel(ax{testnameind,1},'Distance (km)')

cal_arr=1;
for cali=1:length(allyear)-1
    cal_arr=[cal_arr, cal_arr(end)+sum(eomday(allyear(cali),inputmonth))];
end
set(ax{testnameind,1},'color','none','yaxislocation','left','xtick', cal_arr,'xticklabel',datestr(comb_ocean_time(cal_arr), 'yyyy'));

axis tight 
ylim(ax{testnameind,1}, [0,3].*dist_value)
set(ax{testnameind,1},'ycolor','k', 'box', 'off', 'FontSize',m_grid_fontsize);
xlabel(ax{testnameind,1}, 'Year');

set(mslplot{testnameind,1},'LineWidth',2);
grid on
% lgd=legend([mslplot{testnameind,1}], 'moved distance of individuals(lat)');
% set(lgd,'FontSize',m_grid_fontsize);
% set(lgd,'Orientation','horizontal');
% set(lgd,'Location','Northwest');
txt{testnameind,1}=text(480, 2.2.*dist_value, '(B) 30 days', 'FontSize', m_grid_fontsize); 


% % figure save
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [hor_paper_size_x*2, hor_paper_size_y]);
set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height*2]) 
print('-dtiff','-r500',tifname); 
RemoveWhiteSpace([], 'file', tifname);
hold off;
close all;

% save ('Z:\내 드라이브\research\Ph_D_course\Particle tracking of the walleye pollock\data_dryad\Data06_numpar_lon.mat', ...
%             'lat_1d', ...
%             'lon_sum_early_data_15d', 'lon_sum_late_data_15d', 'lon_sum_diff_data_15d', ...
%             'lon_sum_early_data_30d', 'lon_sum_late_data_30d', 'lon_sum_diff_data_30d')

