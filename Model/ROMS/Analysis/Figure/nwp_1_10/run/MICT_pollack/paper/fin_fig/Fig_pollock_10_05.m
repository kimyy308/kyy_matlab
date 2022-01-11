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
% figrawdir =strcat('C:\Users\User\Desktop\'); % % where figure files will be saved
filedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\'); % % where data files are
savedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\pollock_6\');
inputdir = ['/home/auto/MAMS/Data/01_NWP_1_10/Input/'];
LTRANS_testname='Pollock6';

tifname = [figrawdir, 'Fig10_05.tif'];

correction_right_fig=[-0.1600,0,0,0]; % right, up, width, height
correction_upper_fig=[0,0,0,0];
correction_large_fig=[0,0,0,0.020];

% % %  83-87 (early period, period I)
run(param_script);
m_grid_fontsize = m_grid_fontsize -4;
clear uwind_early_RCM comb_uwind_early_RCM
for yearij = 1:length(earlyyear)
    tempyear = earlyyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
        uwind_early_RCM=ncread(ncname, 'uwind');
        lastday_m=size(uwind_early_RCM,3);
        if (exist('comb_uwind_early_RCM')==0)
            comb_uwind_early_RCM=uwind_early_RCM;
        else
            comb_uwind_early_RCM(:,:,end+1:end+lastday_m)=uwind_early_RCM;
        end
    end
end
mean_uwind_early_RCM=mean(comb_uwind_early_RCM,3);

clear vwind_early_RCM comb_vwind_early_RCM
for yearij = 1:length(earlyyear)
    tempyear = earlyyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
        vwind_early_RCM=ncread(ncname, 'vwind');
        lastday_m=size(vwind_early_RCM,3);
        if (exist('comb_vwind_early_RCM')==0)
            comb_vwind_early_RCM=vwind_early_RCM;
        else
            comb_vwind_early_RCM(:,:,end+1:end+lastday_m)=vwind_early_RCM;
        end
    end
end
mean_vwind_early_RCM=mean(comb_vwind_early_RCM,3);

lon_RCM = ncread(ncname, 'lon_rho');
lat_RCM = ncread(ncname, 'lat_rho');

mask_model = double(inpolygon(lon_RCM,lat_RCM,refpolygon(:,1),refpolygon(:,2)));
mask_model(mask_model==0)=NaN;
mean_uwind_early_RCM = mean_uwind_early_RCM .* mask_model;
mean_vwind_early_RCM = mean_vwind_early_RCM .* mask_model;

ref_vec_x_ind = find(abs(lon_RCM(:,1)-m_quiver_ref_text_x_location) == min(abs(lon_RCM(:,1)-m_quiver_ref_text_x_location)));
ref_vec_y_ind = find(abs(lat_RCM(1,:)-m_quiver_ref_text_y_location) == min(abs(lat_RCM(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
% ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval/2+1;
ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
mean_uwind_early_RCM(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_uwind_value;
mean_vwind_early_RCM(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_vwind_value;

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

% % % quiver
hold on
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
ax{testnameind,2}=axes;
set(ax{testnameind,2},'pos',pos_sb{testnameind});
pc{testnameind,2}=m_quiver(lon_RCM(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                    lat_RCM(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                    mean_uwind_early_RCM(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_wind_vector_size, ...
                    mean_vwind_early_RCM(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_wind_vector_size, ...
                    'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
        'xticklabels', [129, 131], 'xtick',[129, 131], 'XaxisLocation', 'top',...
        'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});

    
txt{testnameind,1}=m_text(127.15, 36.9, '(A) 83-87', 'FontSize', m_grid_fontsize-2); 
m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_wind_text, 'FontSize', m_quiver_ref_text_fontsize-8); 






% % %  88-92 (early period, period I)
run(param_script);
m_grid_fontsize = m_grid_fontsize -4;
clear u_late_RCM comb_u_late_RCM
for yearij = 1:length(lateyear)
    tempyear = lateyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
        u_late_RCM=ncread(ncname, 'uwind');
        lastday_m=size(u_late_RCM,3);
        if (exist('comb_u_late_RCM')==0)
            comb_u_late_RCM=u_late_RCM;
        else
            comb_u_late_RCM(:,:,end+1:end+lastday_m)=u_late_RCM;
        end
    end
end
mean_u_late_RCM=mean(comb_u_late_RCM,3);

clear v_late_RCM comb_v_late_RCM
for yearij = 1:length(lateyear)
    tempyear = lateyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
        v_late_RCM=ncread(ncname, 'vwind');
        lastday_m=size(v_late_RCM,3);
        if (exist('comb_v_late_RCM')==0)
            comb_v_late_RCM=v_late_RCM;
        else
            comb_v_late_RCM(:,:,end+1:end+lastday_m)=v_late_RCM;
        end
    end
end
mean_v_late_RCM=mean(comb_v_late_RCM,3);

lon_RCM = ncread(ncname, 'lon_rho');
lat_RCM = ncread(ncname, 'lat_rho');

mask_model = double(inpolygon(lon_RCM,lat_RCM,refpolygon(:,1),refpolygon(:,2)));
mask_model(mask_model==0)=NaN;
mean_u_late_RCM = mean_u_late_RCM .* mask_model;
mean_v_late_RCM = mean_v_late_RCM .* mask_model;

ref_vec_x_ind = find(abs(lon_RCM(:,1)-m_quiver_ref_text_x_location) == min(abs(lon_RCM(:,1)-m_quiver_ref_text_x_location)));
ref_vec_y_ind = find(abs(lat_RCM(1,:)-m_quiver_ref_text_y_location) == min(abs(lat_RCM(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
% ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval/2+1;
ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
mean_u_late_RCM(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_uwind_value;
mean_v_late_RCM(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_vwind_value;

% % projection
testnameind=2;
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
sb1=subplot(2,8,[4 5 12 13]);  % Auto-fitted to the figure.
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
       'xticklabels', [128, 130, 132], 'xtick',[128, 130, 132], 'yticklabels', [],...
       'backcolor', 'none','parent', ax{testnameind,1});

% % % quiver
hold on
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
ax{testnameind,2}=axes;
set(ax{testnameind,2},'pos',pos_sb{testnameind});
pc{testnameind,2}=m_quiver(lon_RCM(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                    lat_RCM(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                    mean_u_late_RCM(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_wind_vector_size, ...
                    mean_v_late_RCM(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_wind_vector_size, ...
                    'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
        'xticklabels', [129, 131], 'xtick',[129, 131], 'XaxisLocation', 'top',...
        'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});

    
txt{testnameind,1}=m_text(127.15, 36.9, '(B) 88-92', 'FontSize', m_grid_fontsize-2); 
m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_wind_text, 'FontSize', m_quiver_ref_text_fontsize-8); 



% % %  get difference (late-early period, period II-I)
mean_uwind_diff_RCM=mean_u_late_RCM-mean_uwind_early_RCM;
mean_vwind_diff_RCM=mean_v_late_RCM-mean_vwind_early_RCM;
mean_uwind_diff_RCM(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_uwind_value;
mean_vwind_diff_RCM(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_vwind_value;


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

% % % quiver
hold on
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
ax{testnameind,2}=axes;
set(ax{testnameind,2},'pos',pos_sb{testnameind});
pc{testnameind,2}=m_quiver(lon_RCM(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                    lat_RCM(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                    mean_uwind_diff_RCM(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_wind_vector_size, ...
                    mean_vwind_diff_RCM(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_wind_vector_size, ...
                    'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);
m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
        'xticklabels', [129, 131], 'xtick',[129, 131], 'XaxisLocation', 'top',...
        'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});

    
% txt{testnameind,1}=m_text(127.15, 36.9, '(C) B-A', 'FontSize', m_grid_fontsize-2); 
txt{testnameind,1}=m_text(127.15, 36.9, '(C) (B)-(A)', 'FontSize', m_grid_fontsize-3); 
m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_wind_text, 'FontSize', m_quiver_ref_text_fontsize-8); 



% % figure save
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height*2]) 
print('-dtiff','-r500',tifname); RemoveWhiteSpace([], 'file', tifname);
hold off;
close all;

save ('D:\research\Ph_D_course\2021_Particle tracking of the walleye pollock\data_dryad\Data11_wind.mat', ...
            'lon_RCM', 'lat_RCM', 'RCM_land', ...
            'mean_uwind_early_RCM', 'mean_vwind_early_RCM', ...
            'mean_u_late_RCM', 'mean_v_late_RCM', ...
            'mean_uwind_diff_RCM', 'mean_vwind_diff_RCM')
