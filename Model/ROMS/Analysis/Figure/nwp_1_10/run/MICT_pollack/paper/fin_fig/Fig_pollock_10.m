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
figrawdir =strcat('C:\Users\User\Desktop\'); % % where figure files will be saved
filedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\'); % % where data files are
savedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\pollock_6\');
inputdir = ['/home/auto/MAMS/Data/01_NWP_1_10/Input/'];
LTRANS_testname='Pollock6';

tifname = [figrawdir, 'Fig10.tif'];

correction_right_fig=[-0.1600,0,0,0]; % right, up, width, height
correction_upper_fig=[0,0,0,0]; % right, up, width, height
correction_large_fig=[0,0,0,-0.15]; % right, up, width, height


% % % 
% % %  row 1
% % % 
temp_checktime=15;
run(param_script);
m_grid_fontsize = m_grid_fontsize -2;
clear wind_curl comb_wind_curl u_rho comb_u_rho v_rho comb_v_rho uwind comb_uwind vwind comb_vwind ...
    egg_mask comb_egg_mask sp_ground comb_sp_ground wind_curl2 comb_wind_curl2 temp_surf comb_temp_surf  curl_mask
for yearij = 1:length(allyear)
    tempyear = allyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])  
        if(exist('curl_mask')==0)
            lon_rho = ncread(ncname, 'lon_rho');
            lat_rho = ncread(ncname, 'lat_rho');
            [vel_polygon, lonlat, error_status] = Func_0007_get_polygon_data_from_regionname('SK_EEZ');
            vel_mask = double(inpolygon(lon_rho,lat_rho,vel_polygon(:,1),vel_polygon(:,2)));
        end
        ocean_time=ncread(ncname, 'time')+datenum(1900,12,31);
        u_rho=ncread(ncname, 'u_rho').*vel_mask;
        u_rho(u_rho==0)=NaN;
        v_rho=ncread(ncname, 'v_rho').*vel_mask;
        v_rho(v_rho==0)=NaN;
        all_checktime=ncread(ncname, 'checktime');
        ind_checktime=find(all_checktime==temp_checktime);
        egg_mask=squeeze(ncread(ncname, 'mask_par', [1 1 1 ind_checktime], [inf inf inf 1])).*vel_mask;
%             egg_mask=ncread(ncname, 'mask_15day').*vel_mask;
        egg_mask(egg_mask==0)=NaN;

        lastday_m=size(u_rho,3);
        if (exist('comb_u_rho')==0)
            comb_u_rho=u_rho;
            comb_v_rho=v_rho;
            comb_egg_mask=egg_mask;
            comb_ocean_time=ocean_time;
        else
            comb_u_rho(:,:,end+1:end+lastday_m)=u_rho;
            comb_v_rho(:,:,end+1:end+lastday_m)=v_rho;
            comb_egg_mask(:,:,end+1:end+lastday_m)=egg_mask;
            comb_ocean_time(end+1:end+lastday_m)=ocean_time;
        end
    end
end

ts_u_rho=reshape(comb_u_rho,[size(comb_u_rho,1)*size(comb_u_rho,2), size(comb_u_rho,3)]);
mean_ts_u_rho=mean(ts_u_rho,1,'omitnan');
ts_v_rho=reshape(comb_v_rho,[size(comb_v_rho,1)*size(comb_v_rho,2), size(comb_v_rho,3)]);
mean_ts_v_rho=mean(ts_v_rho,1,'omitnan');
ts_egg_mask=reshape(comb_egg_mask,[size(comb_egg_mask,1)*size(comb_egg_mask,2), size(comb_egg_mask,3)]);
sum_ts_egg_mask=sum(ts_egg_mask,1,'omitnan');

half_len=round(length(sum_ts_egg_mask)/2);
egg_half_1=mean(sum_ts_egg_mask(1:half_len));
egg_half_2=mean(sum_ts_egg_mask(half_len+1:end));
regime_ts_egg_mask(1:half_len)=egg_half_1;
regime_ts_egg_mask(half_len+1:length(sum_ts_egg_mask))=egg_half_2;


% % projection
testnameind=1;
sb1=subplot(2,2,[1 2]);  % Auto-fitted to the figure.
pos_sb{testnameind,1}=get(sb1, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [0,0,0,0];
delete(sb1); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

mslplot{testnameind,2}=bar(sum_ts_egg_mask, 'FaceColor', [0.6 0.6 0.6], 'parent',ax{testnameind,1});
hold on
mslplot{testnameind,1}=plot(regime_ts_egg_mask,'k','parent',ax{testnameind,1});
ylabel(ax{testnameind,1},'# of individuals')

cal_arr=1;
for cali=1:length(allyear)-1
    cal_arr=[cal_arr, cal_arr(end)+sum(eomday(allyear(cali),inputmonth))];
end
% set(ax{testnameind,1},'color','none','yaxislocation','left','xtick', cal_arr,'xticklabel',datestr(comb_ocean_time(cal_arr), 'yyyy'),'position', ax_pos+[0 0.02 -0.01 -0.02]);
set(ax{testnameind,1},'color','none','yaxislocation','left','xtick', cal_arr,'xticklabel',datestr(comb_ocean_time(cal_arr), 'yyyy'));

axis tight 
% ylim(ax{testnameind,1}, [-1,1])
set(ax{testnameind,1},'ycolor','k', 'box', 'off', 'FontSize',m_grid_fontsize);
xlabel(ax{testnameind,1}, 'Year');

set(mslplot{testnameind,1},'LineWidth',2);
grid on
% lgd=legend([mslplot{testnameind,1}], 'moved distance of individuals(lat)');
% set(lgd,'FontSize',m_grid_fontsize);
% set(lgd,'Orientation','horizontal');
% set(lgd,'Location','Northwest');

txt{testnameind,1}=text(400, 30, '(A) 15 days', 'FontSize', m_grid_fontsize); 


% % % % % % % % % % 
% % % % % % % % % %  row 2
% % % % % % % % % %
temp_checktime=30;
run(param_script);
m_grid_fontsize = m_grid_fontsize -2;
clear wind_curl comb_wind_curl u_rho comb_u_rho v_rho comb_v_rho uwind comb_uwind vwind comb_vwind ...
    egg_mask comb_egg_mask sp_ground comb_sp_ground wind_curl2 comb_wind_curl2 temp_surf comb_temp_surf  curl_mask
for yearij = 1:length(allyear)
    tempyear = allyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])  
        if(exist('curl_mask')==0)
            lon_rho = ncread(ncname, 'lon_rho');
            lat_rho = ncread(ncname, 'lat_rho');
            [vel_polygon, lonlat, error_status] = Func_0007_get_polygon_data_from_regionname('SK_EEZ');
            vel_mask = double(inpolygon(lon_rho,lat_rho,vel_polygon(:,1),vel_polygon(:,2)));
        end
        ocean_time=ncread(ncname, 'time')+datenum(1900,12,31);
        u_rho=ncread(ncname, 'u_rho').*vel_mask;
        u_rho(u_rho==0)=NaN;
        v_rho=ncread(ncname, 'v_rho').*vel_mask;
        v_rho(v_rho==0)=NaN;
        all_checktime=ncread(ncname, 'checktime');
        ind_checktime=find(all_checktime==temp_checktime);
        egg_mask=squeeze(ncread(ncname, 'mask_par', [1 1 1 ind_checktime], [inf inf inf 1])).*vel_mask;
%             egg_mask=ncread(ncname, 'mask_15day').*vel_mask;
        egg_mask(egg_mask==0)=NaN;

        lastday_m=size(u_rho,3);
        if (exist('comb_u_rho')==0)
            comb_u_rho=u_rho;
            comb_v_rho=v_rho;
            comb_egg_mask=egg_mask;
            comb_ocean_time=ocean_time;
        else
            comb_u_rho(:,:,end+1:end+lastday_m)=u_rho;
            comb_v_rho(:,:,end+1:end+lastday_m)=v_rho;
            comb_egg_mask(:,:,end+1:end+lastday_m)=egg_mask;
            comb_ocean_time(end+1:end+lastday_m)=ocean_time;
        end
    end
end

ts_u_rho=reshape(comb_u_rho,[size(comb_u_rho,1)*size(comb_u_rho,2), size(comb_u_rho,3)]);
mean_ts_u_rho=mean(ts_u_rho,1,'omitnan');
ts_v_rho=reshape(comb_v_rho,[size(comb_v_rho,1)*size(comb_v_rho,2), size(comb_v_rho,3)]);
mean_ts_v_rho=mean(ts_v_rho,1,'omitnan');
ts_egg_mask=reshape(comb_egg_mask,[size(comb_egg_mask,1)*size(comb_egg_mask,2), size(comb_egg_mask,3)]);
sum_ts_egg_mask=sum(ts_egg_mask,1,'omitnan');

half_len=round(length(sum_ts_egg_mask)/2);
egg_half_1=mean(sum_ts_egg_mask(1:half_len));
egg_half_2=mean(sum_ts_egg_mask(half_len+1:end));
regime_ts_egg_mask(1:half_len)=egg_half_1;
regime_ts_egg_mask(half_len+1:length(sum_ts_egg_mask))=egg_half_2;

% % projection
testnameind=2;
sb1=subplot(2,2,[3 4]);  % Auto-fitted to the figure.
pos_sb{testnameind,1}=get(sb1, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [0,0.23,0,0];
delete(sb1); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

mslplot{testnameind,2}=bar(sum_ts_egg_mask, 'FaceColor', [0.6 0.6 0.6], 'parent',ax{testnameind,1});
hold on
mslplot{testnameind,1}=plot(regime_ts_egg_mask,'k','parent',ax{testnameind,1});
ylabel(ax{testnameind,1},'# of individuals')

cal_arr=1;
for cali=1:length(allyear)-1
    cal_arr=[cal_arr, cal_arr(end)+sum(eomday(allyear(cali),inputmonth))];
end
set(ax{testnameind,1},'color','none','yaxislocation','left','xtick', cal_arr,'xticklabel',datestr(comb_ocean_time(cal_arr), 'yyyy'));

axis tight 
% ylim(ax{testnameind,1}, [-1,1])
set(ax{testnameind,1},'ycolor','k', 'box', 'off', 'FontSize',m_grid_fontsize);
xlabel(ax{testnameind,1}, 'Year');

set(mslplot{testnameind,1},'LineWidth',2);
grid on
% lgd=legend([mslplot{testnameind,1}], 'moved distance of individuals(lat)');
% set(lgd,'FontSize',m_grid_fontsize);
% set(lgd,'Orientation','horizontal');
% set(lgd,'Location','Northwest');
txt{testnameind,1}=text(400, 45, '(B) 30 days', 'FontSize', m_grid_fontsize); 


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

