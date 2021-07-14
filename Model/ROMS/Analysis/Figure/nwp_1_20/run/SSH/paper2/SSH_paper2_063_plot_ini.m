clear all 
close all;
dropboxpath='C:\Users\User\Dropbox';
addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\FIgure\nwp_1_20\run\SSH\function']));
system_name=computer;
er_status = Func_0008_set_dropbox_path(system_name);

scenname = 'historical'
modelnames = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR'};


% % % % CMIP5 Init plot

% for modeli= 1:length(modelnames)
% modelname = modelnames{modeli};
% % jpgname = ['Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\paper2\fin_fig_r1\ini\ini_', modelname, '_bot', '.tif'];
% jpgname = ['Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\paper2\fin_fig_r1\ini\ini_', modelname, '_surf', '.tif'];
% 
% inidir= 'D:\Data\Model\ROMS\nwp_1_20\input';
% inifilename =[inidir, '\', 'roms_nwp_ini_', modelname, '.nc']; 
% gridname= [inidir, '\test53\', 'roms_grid_nwp_1_20_test53.nc'];
% temp=ncread(inifilename,'temp', [1 1 40 1], [inf inf 1 1]);
% lon_rho=ncread(gridname,'lon_rho');
% lat_rho=ncread(gridname,'lat_rho');
% lonlat=[115 164 15 52];
% regionname = 'NWP';
% param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
% 
% run(param_script);
% 
% 
% m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
% hold on;
% m_pcolor(lon_rho',lat_rho',temp(:,:,1)'-273.15);
% shading(gca,m_pcolor_shading_method);
% m_gshhs_i('color',m_gshhs_line_color);
% m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
% % titlename = strcat('corr, ',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
% %     num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(corr_interped(:), 'omitnan'),2)));  %% + glacier contribution
% % 
% % title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
% % set colorbar 
% h = colorbar;
% colormap(jet);
% set(h,'fontsize',colorbar_fontsize);
% %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
% caxis([0 30]);
% 
% % set grid
% m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
% set(gcf, 'PaperUnits', 'points');
% set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
% saveas(gcf,jpgname,'tif');
% 
% hold off
% close all;
% end



% % % % % WOA ini plot
% jpgname = ['Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\paper2\fin_fig_r1\ini\ini_WOA98', '_surf', '.tif'];
% 
% inidir= 'D:\Data\Model\ROMS\nwp_1_20\input';
% inifilename =[inidir, '\test53\', 'roms_nwp_ini_', 'test53', '.nc']; 
% gridname= [inidir, '\test53\', 'roms_grid_nwp_1_20_test53.nc'];
% temp=ncread(inifilename,'temp', [1 1 40 1], [inf inf 1 1]);
% lon_rho=ncread(gridname,'lon_rho');
% lat_rho=ncread(gridname,'lat_rho');
% lonlat=[115 164 15 52];
% regionname = 'NWP';
% param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
% 
% run(param_script);
% 
% 
% m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
% hold on;
% m_pcolor(lon_rho',lat_rho',temp(:,:,1)');
% shading(gca,m_pcolor_shading_method);
% m_gshhs_i('color',m_gshhs_line_color);
% m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
% % titlename = strcat('corr, ',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
% %     num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(corr_interped(:), 'omitnan'),2)));  %% + glacier contribution
% % 
% % title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
% % set colorbar 
% h = colorbar;
% colormap(jet);
% set(h,'fontsize',colorbar_fontsize);
% %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
% caxis([0 30]);
% 
% % set grid
% m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
% set(gcf, 'PaperUnits', 'points');
% set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
% saveas(gcf,jpgname,'tif');
% 
% hold off
% close all;



% % % % ROMS 1976 1m plot
% 
% for modeli= 1:length(modelnames)
% modelname = modelnames{modeli};
% rcmtestname = Func_0006_get_RCMname_from_GCM(modelname, scenname);
% jpgname = ['Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\paper2\fin_fig_r1\ini\ini_', rcmtestname, '_surf', '.tif'];
% 
% inidir= 'D:\Data\Model\ROMS\nwp_1_20\input';
% outputdir = 'J:\Data\Model\ROMS\nwp_1_20';
% % J:\Data\Model\ROMS\nwp_1_20\test53\run\1976
% inifilename =[outputdir, '\', rcmtestname, '\run\1976\', rcmtestname, '_monthly_1976_01', '.nc']; 
% gridname= [inidir, '\test53\', 'roms_grid_nwp_1_20_test53.nc'];
% temp=ncread(inifilename,'temp', [1 1 40 1], [inf inf 1 1]);
% lon_rho=ncread(gridname,'lon_rho');
% lat_rho=ncread(gridname,'lat_rho');
% lonlat=[115 164 15 52];
% regionname = 'NWP';
% param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
% 
% run(param_script);
% 
% 
% m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
% hold on;
% m_pcolor(lon_rho',lat_rho',temp(:,:,1)');
% shading(gca,m_pcolor_shading_method);
% m_gshhs_i('color',m_gshhs_line_color);
% m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
% % titlename = strcat('corr, ',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
% %     num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(corr_interped(:), 'omitnan'),2)));  %% + glacier contribution
% % 
% % title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
% % set colorbar 
% h = colorbar;
% colormap(jet);
% set(h,'fontsize',colorbar_fontsize);
% %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
% caxis([0 30]);
% 
% % set grid
% m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
% set(gcf, 'PaperUnits', 'points');
% set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
% saveas(gcf,jpgname,'tif');
% 
% hold off
% close all;
% end





% % % ROMS WOA ini - GCM ini plot

for modeli= 1:length(modelnames)
modelname = modelnames{modeli};
rcmtestname = Func_0006_get_RCMname_from_GCM(modelname, scenname);
jpgname = ['Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\paper2\fin_fig_r1\ini\ini_diff_woa_', modelname, '_surf', '.tif'];

inidir= 'D:\Data\Model\ROMS\nwp_1_20\input';
outputdir = 'J:\Data\Model\ROMS\nwp_1_20';
% J:\Data\Model\ROMS\nwp_1_20\test53\run\1976

inifilename =[inidir, '\test53\', 'roms_nwp_ini_', 'test53', '.nc']; 
gridname= [inidir, '\test53\', 'roms_grid_nwp_1_20_test53.nc'];
rcm_temp=ncread(inifilename,'temp', [1 1 40 1], [inf inf 1 1]);

inifilename =[inidir, '\', 'roms_nwp_ini_', modelname, '.nc']; 
gcm_temp=ncread(inifilename,'temp', [1 1 40 1], [inf inf 1 1])-273.15;

temp= rcm_temp-gcm_temp;

lon_rho=ncread(gridname,'lon_rho');
lat_rho=ncread(gridname,'lat_rho');
lonlat=[115 164 15 52];
regionname = 'NWP';
param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];

run(param_script);


m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
hold on;
m_pcolor(lon_rho',lat_rho',temp(:,:,1)');
shading(gca,m_pcolor_shading_method);
m_gshhs_i('color',m_gshhs_line_color);
m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
% titlename = strcat('corr, ',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
%     num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(corr_interped(:), 'omitnan'),2)));  %% + glacier contribution
% 
% title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

% set colorbar 
h = colorbar;
colormap(bwrmap);
set(h,'fontsize',colorbar_fontsize);
%             title(h,'mm/y','fontsize',colorbar_title_fontsize);
caxis([-10 10]);

% set grid
m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

saveas(gcf,jpgname,'tif');

hold off
close all;
end



% % % % ROMS 1976 1m - GCM ini plot
% for modeli= 1:length(modelnames)
% modelname = modelnames{modeli};
% rcmtestname = Func_0006_get_RCMname_from_GCM(modelname, scenname);
% jpgname = ['Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\paper2\fin_fig_r1\ini\ini_diff_', rcmtestname, '_bot', '.tif'];
% 
% inidir= 'D:\Data\Model\ROMS\nwp_1_20\input';
% outputdir = 'J:\Data\Model\ROMS\nwp_1_20';
% % J:\Data\Model\ROMS\nwp_1_20\test53\run\1976
% inifilename =[outputdir, '\', rcmtestname, '\run\1976\', rcmtestname, '_monthly_1976_01', '.nc']; 
% gridname= [inidir, '\test53\', 'roms_grid_nwp_1_20_test53.nc'];
% rcm_temp=ncread(inifilename,'temp', [1 1 1 1], [inf inf 1 1]);
% 
% inifilename =[inidir, '\', 'roms_nwp_ini_', modelname, '.nc']; 
% gcm_temp=ncread(inifilename,'temp', [1 1 1 1], [inf inf 1 1])-273.15;
% 
% temp= rcm_temp-gcm_temp;
% 
% lon_rho=ncread(gridname,'lon_rho');
% lat_rho=ncread(gridname,'lat_rho');
% lonlat=[115 164 15 52];
% regionname = 'NWP';
% param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
% 
% run(param_script);
% 
% 
% m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
% hold on;
% m_pcolor(lon_rho',lat_rho',temp(:,:,1)');
% shading(gca,m_pcolor_shading_method);
% m_gshhs_i('color',m_gshhs_line_color);
% m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
% % titlename = strcat('corr, ',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
% %     num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(corr_interped(:), 'omitnan'),2)));  %% + glacier contribution
% % 
% % title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
% % set colorbar 
% h = colorbar;
% colormap(bwrmap);
% set(h,'fontsize',colorbar_fontsize);
% %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
% caxis([-10 10]);
% 
% % set grid
% m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
% set(gcf, 'PaperUnits', 'points');
% set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
% saveas(gcf,jpgname,'tif');
% 
% hold off
% close all;
% end