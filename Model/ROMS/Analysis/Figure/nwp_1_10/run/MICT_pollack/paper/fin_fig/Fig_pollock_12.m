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

tifname = [figrawdir, 'Fig12.tif'];

correction_right_fig=[-0.1600,0,0,0]; % right, up, width, height
correction_upper_fig=[0,0,0,0]; % right, up, width, height
correction_large_fig=[0,0,0,-0.15]; % right, up, width, height


% % % 
% % %  row 1
% % % 
run(param_script);

m_grid_fontsize = m_grid_fontsize -2;
clear wind_curl comb_wind_curl u_rho comb_u_rho v_rho comb_v_rho uwind comb_uwind vwind comb_vwind ...
    egg_mask comb_egg_mask sp_ground comb_sp_ground wind_curl2 comb_wind_curl2 temp_surf comb_temp_surf  curl_mask mean_ts_aoi
for yearij = 1:length(allyear)
    tempyear = allyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])  
        if(exist('curl_mask')==0)
            lon_rho = ncread(ncname, 'lon_rho');
            lat_rho = ncread(ncname, 'lat_rho');
            vel_polygon=[128, 37.5; 128, 39; 130, 39; 130, 37.5];
            vel_mask = double(inpolygon(lon_rho,lat_rho,vel_polygon(:,1),vel_polygon(:,2)));
        end
        ocean_time=ncread(ncname, 'time')+datenum(1900,12,31);
        u_rho=ncread(ncname, 'u_rho').*vel_mask;
        u_rho(u_rho==0)=NaN;
        lastday_m=size(u_rho,3);
        if (exist('comb_u_rho')==0)
            comb_u_rho=u_rho;
            comb_ocean_time=ocean_time;
        else
            comb_u_rho(:,:,end+1:end+lastday_m)=u_rho;
            comb_ocean_time(end+1:end+lastday_m)=ocean_time;
        end
    end
end

opts = spreadsheetImportOptions("NumVariables", 13);
opts.Sheet = "AOI";
opts.DataRange = "A2:M51";
opts.VariableNames = ["yy", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
climate2021S1 = readtable("Z:\내 드라이브\research\Ph_D_course\Particle tracking of the walleye pollock\climate_2021.xlsx", opts, "UseExcel", false);
climate2021S1 = table2array(climate2021S1);
clear opts

AOI_year = climate2021S1(:,1);
AOI_value = climate2021S1(:,2:13)';
AOI_value= AOI_value(:);
AOI_year = repmat(AOI_year, [1 12])';
AOI_year = AOI_year(:);
AOI_month = 1:12;
AOI_month = repmat(AOI_month, [length(AOI_year)/12 1])';
AOI_month = AOI_month(:);

% %                 10 year lowpass filter
AOI_value=lowpass(AOI_value, 1/120, 1, 'Steepness', 0.99);  % phase correction is applied automatically, but butterworth isn't.

AOI_year=reshape(AOI_year, [12, length(AOI_year)/12]);
AOI_month=reshape(AOI_month, [12, length(AOI_month)/12]);
AOI_value=reshape(AOI_value, [12, length(AOI_value)/12]);
AOI_year_start=find(AOI_year(1,:)==min(allyear));
AOI_year_end=find(AOI_year(1,:)==max(allyear));

cal_arr=1;
for cali=1:length(allyear)
    if cali==length(allyear)
        for monthij=2:length(inputmonth)
            cal_arr=[cal_arr, cal_arr(end)+sum(eomday(allyear(cali),inputmonth(monthij)))];
        end
    else
        for monthij=1:length(inputmonth)
            cal_arr=[cal_arr, cal_arr(end)+sum(eomday(allyear(cali),inputmonth(monthij)))];
        end
    end
end

AOI_year_yearly = mean(AOI_year(inputmonth, allyear-1969), 1);
AOI_month_yearly = mean(AOI_month(inputmonth, allyear-1969), 1);
AOI_value_yearly = mean(AOI_value(inputmonth, allyear-1969), 1);

half_len=round(length(AOI_value_yearly)/2);
egg_half_1=mean(AOI_value_yearly(1:half_len), 'omitnan');
egg_half_2=mean(AOI_value_yearly(half_len+1:end), 'omitnan');
regime_ts_egg_mask(1:half_len)=egg_half_1;
regime_ts_egg_mask(half_len+1:length(AOI_value_yearly))=egg_half_2;

% % projection
testnameind=1;
sb1=subplot(2,2,[1 2]);  % Auto-fitted to the figure.
pos_sb{testnameind,1}=get(sb1, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [0,0,0,0];
delete(sb1); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

mslplot{testnameind,2}=bar(AOI_year_yearly, AOI_value_yearly, 'FaceColor', [0.6 0.6 0.6], 'parent',ax{testnameind,1});
hold on
mslplot{testnameind,1}=plot(AOI_year_yearly, regime_ts_egg_mask, 'k','parent',ax{testnameind,1});
ylabel(ax{testnameind,1},'AOI')

axis tight 
% ylim(ax{testnameind,1}, [-1,1])
set(ax{testnameind,1},'ycolor','k', 'box', 'off', 'FontSize',m_grid_fontsize);
xlabel(ax{testnameind,1}, 'Year');

set(mslplot{testnameind,1},'LineWidth',2);
grid on

txt{testnameind,1}=text(1983, 0.5, '(A) AOI', 'FontSize', m_grid_fontsize); 


% % % % % % % % % % 
% % % % % % % % % %  row 2
% % % % % % % % % %
run(param_script);

m_grid_fontsize = m_grid_fontsize -2;
clear wind_curl comb_wind_curl u_rho comb_u_rho v_rho comb_v_rho uwind comb_uwind vwind comb_vwind ...
    egg_mask comb_egg_mask sp_ground comb_sp_ground wind_curl2 comb_wind_curl2 temp_surf comb_temp_surf  curl_mask mean_ts_EAWMI
for yearij = 1:length(allyear)
    tempyear = allyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])  
        if(exist('curl_mask')==0)
            lon_rho = ncread(ncname, 'lon_rho');
            lat_rho = ncread(ncname, 'lat_rho');
            vel_polygon=[128, 37.5; 128, 39; 130, 39; 130, 37.5];
            vel_mask = double(inpolygon(lon_rho,lat_rho,vel_polygon(:,1),vel_polygon(:,2)));
        end
        ocean_time=ncread(ncname, 'time')+datenum(1900,12,31);
        u_rho=ncread(ncname, 'u_rho').*vel_mask;
        u_rho(u_rho==0)=NaN;
        lastday_m=size(u_rho,3);
        if (exist('comb_u_rho')==0)
            comb_u_rho=u_rho;
            comb_ocean_time=ocean_time;
        else
            comb_u_rho(:,:,end+1:end+lastday_m)=u_rho;
            comb_ocean_time(end+1:end+lastday_m)=ocean_time;
        end
    end
end

opts = spreadsheetImportOptions("NumVariables", 13);
opts.Sheet = "EAWMI";
opts.DataRange = "A2:M51";
opts.VariableNames = ["yy", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
climate2021S1 = readtable("Z:\내 드라이브\research\Ph_D_course\Particle tracking of the walleye pollock\climate_2021.xlsx", opts, "UseExcel", false);
climate2021S1 = table2array(climate2021S1);
clear opts

EAWMI_year = climate2021S1(:,1);
EAWMI_value = climate2021S1(:,2:13)';
EAWMI_value= EAWMI_value(:);
EAWMI_year = repmat(EAWMI_year, [1 12])';
EAWMI_year = EAWMI_year(:);
EAWMI_month = 1:12;
EAWMI_month = repmat(EAWMI_month, [length(EAWMI_year)/12 1])';
EAWMI_month = EAWMI_month(:);

% %                 10 year lowpass filter
EAWMI_value=lowpass(EAWMI_value, 1/120, 1, 'Steepness', 0.99);  % phase correction is applied automatically, but butterworth isn't.

EAWMI_year=reshape(EAWMI_year, [12, length(EAWMI_year)/12]);
EAWMI_month=reshape(EAWMI_month, [12, length(EAWMI_month)/12]);
EAWMI_value=reshape(EAWMI_value, [12, length(EAWMI_value)/12]);
EAWMI_year_start=find(EAWMI_year(1,:)==min(allyear));
EAWMI_year_end=find(EAWMI_year(1,:)==max(allyear));

cal_arr=1;
for cali=1:length(allyear)
    if cali==length(allyear)
        for monthij=2:length(inputmonth)
            cal_arr=[cal_arr, cal_arr(end)+sum(eomday(allyear(cali),inputmonth(monthij)))];
        end
    else
        for monthij=1:length(inputmonth)
            cal_arr=[cal_arr, cal_arr(end)+sum(eomday(allyear(cali),inputmonth(monthij)))];
        end
    end
end

EAWMI_year_yearly = mean(EAWMI_year(inputmonth, allyear-1969), 1);
EAWMI_month_yearly = mean(EAWMI_month(inputmonth, allyear-1969), 1);
EAWMI_value_yearly = mean(EAWMI_value(inputmonth, allyear-1969), 1);
half_len=round(length(EAWMI_value_yearly)/2);
egg_half_1=mean(EAWMI_value_yearly(1:half_len), 'omitnan');
egg_half_2=mean(EAWMI_value_yearly(half_len+1:end), 'omitnan');
regime_ts_egg_mask(1:half_len)=egg_half_1;
regime_ts_egg_mask(half_len+1:length(EAWMI_value_yearly))=egg_half_2;
    
% % projection
testnameind=2;
sb1=subplot(2,2,[3 4]);  % Auto-fitted to the figure.
pos_sb{testnameind,1}=get(sb1, 'pos'); % Get the position.
pos_sb{testnameind,1} = pos_sb{testnameind,1} + correction_upper_fig+ correction_large_fig + [0,0.23,0,0];
delete(sb1); % Delete the subplot axes
ax{testnameind,1}=axes;
set(ax{testnameind,1},'pos',pos_sb{testnameind,1});

mslplot{testnameind,2}=bar(EAWMI_year_yearly, EAWMI_value_yearly, 'FaceColor', [0.6 0.6 0.6], 'parent',ax{testnameind,1});
hold on
mslplot{testnameind,1}=plot(EAWMI_year_yearly, regime_ts_egg_mask, 'k','parent',ax{testnameind,1});
ylabel(ax{testnameind,1},'EAWMI')

axis tight 
ylim(ax{testnameind,1}, [16,21])
set(ax{testnameind,1},'ycolor','k', 'box', 'off', 'FontSize',m_grid_fontsize);
xlabel(ax{testnameind,1}, 'Year');

set(mslplot{testnameind,1},'LineWidth',2);
grid on
% lgd=legend([mslplot{testnameind,1}], 'moved distance of individuals(lat)');
% set(lgd,'FontSize',m_grid_fontsize);
% set(lgd,'Orientation','horizontal');
% set(lgd,'Location','Northwest');
txt{testnameind,1}=text(1983, 20.5, '(B) EAWMI', 'FontSize', m_grid_fontsize); 


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

