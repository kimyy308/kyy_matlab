% %  Created 20-Sep-2023 by Yong-Yub Kim
clc; clear all; close all;
warning off;

%% set path
[error_status, tmp.hostname] = system('hostname');
tmp.hostname=tmp.hostname(1:end-1);
switch tmp.hostname
    case 'Yong-Yubs-iMac-Pro.local'
        tmp.dropboxpath = '/Volumes/kyy_raid/kimyy/Dropbox';
        tmp.kimyypath = '/Volumes/kyy_raid/kimyy';
    case 'Yong-Yubui-MacBookPro.local'
        tmp.dropboxpath = '/Users/kimyy/Dropbox';
        tmp.kimyypath = '/Users/kimyy';
    case {'da1', 'da2', 'da3', 'da4'}
        tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
end
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

%% model configuration
% cfg.var='TS'; %SST PRECT PSL TS SSH sumChl
% cfg.vars = {'SSH', 'sumChl'};
cfg.vars = {'SST', 'PRECT', 'PSL', 'TS'};
% cfg.vars = {'SSH'};
% cfg.vars = {'SST', 'PRECT', 'PSL', 'TS', 'sumChl', 'photoC_TOT_zint'};
% cfg.vars = { 'sumChl', 'photoC_TOT_zint',  'SSH'};
% cfg.vars = {'diatChl', 'diazChl', 'spChl'};

% cfg.vars={'diatChl', 'diazChl', 'spChl', 'diatC', 'diazC', 'spC', 'SSH', 'BSF', 'HBLT', ...
%     'PAR_avg', 'NO3', 'PO4', 'SiO3', 'Fe', 'TEMP', 'UVEL', 'VVEL', 'WVEL', ...
%     'diat_agg_zint_100m', 'diat_Fe_lim_Cweight_avg_100m', ...
%     'diat_Fe_lim_surf', 'diat_light_lim_Cweight_avg_100m', 'diat_light_lim_surf', ...
%     'diat_loss_zint_100m', 'diat_N_lim_Cweight_avg_100m', 'diat_N_lim_surf', ...
%     'diat_P_lim_Cweight_avg_100m', 'diat_P_lim_surf', 'diat_SiO3_lim_Cweight_avg_100m', ...
%     'diat_SiO3_lim_surf', 'diaz_agg_zint_100m', 'diaz_Fe_lim_Cweight_avg_100m', ...
%     'diaz_Fe_lim_surf', 'diaz_light_lim_Cweight_avg_100m', 'diaz_light_lim_surf', ...
%     'diaz_loss_zint_100m', 'diaz_P_lim_Cweight_avg_100m', 'diaz_P_lim_surf', 'dustToSed', ...
%     'graze_diat_zint_100m', 'graze_diat_zoo_zint_100m', 'graze_diaz_zint_100m', ...
%     'graze_diaz_zoo_zint_100m', 'graze_sp_zint_100m', 'graze_sp_zoo_zint_100m', ...
%     'LWUP_F', 'O2_ZMIN_DEPTH', 'photoC_diat_zint_100m', 'photoC_diaz_zint_100m', ...
%     'photoC_sp_zint_100m', ...
%     'sp_agg_zint_100m', 'sp_Fe_lim_Cweight_avg_100m', 'sp_Fe_lim_surf', ...
%     'sp_light_lim_Cweight_avg_100m', 'sp_light_lim_surf', 'sp_loss_zint_100m', ...
%     'sp_N_lim_Cweight_avg_100m', 'sp_N_lim_surf', 'sp_P_lim_Cweight_avg_100m', ...
%     'sp_P_lim_surf', 'zoo_loss_zint_100m', 'photoC_TOT_zint_100m', 'sumChl'};
% cfg.vars={'graze_diaz_zint_100m', 'zoo_loss_zint_100m', 'photoC_TOT_zint_100m'};
% cfg.vars={'sp_P_lim_surf'};
% cfg.vars = {'SST', 'PRECT', 'PSL', 'TS'};
% cfg.vars = {'PAR_avg', 'NO3', 'PO4', 'SiO3', 'Fe'};
% cfg.vars = {'NO3'};
% cfg.vars = {'diatC', 'diazC', 'spC'};
% cfg.vars = {'SSH', 'BSF'};
% cfg.vars = {'photoC_diat_zint_100m', 'photoC_diaz_zint_100m', 'photoC_sp_zint_100m', 'photoC_TOT_zint_100m'};
% % 'graze_diaz_zoo_zint_100m', 
% cfg.vars = {'PAR_avg', 'NO3', 'PO4', 'SiO3', 'Fe', ...
%     'diatChl', 'diazChl', 'spChl', 'diatC', 'diazC', 'spC', 'SSH', 'BSF', ...
%     'TEMP', 'UVEL', 'VVEL', 'WVEL', 'HBLT'};
% cfg.vars = {'SSH'};
% cfg.vars = {'photoC_sp_zint_100m'};
% cfg.vars = {'diazC'};
% cfg.vars = {'TEMP', 'UVEL', 'VVEL', 'WVEL'};
% cfg.vars = {'diat_N_lim_Cweight_avg_100m'};
% cfg.vars={'photoC_TOT_zint_100m'};

% cfg.vars = {'HBLT'};
% cfg.vars={'photoC_TOT_zint_100m',  'photoC_diat_zint_100m', 'photoC_diaz_zint_100m', 'photoC_sp_zint_100m', ...
%      'diat_Fe_lim_Cweight_avg_100m', ...
%     'diat_light_lim_Cweight_avg_100m',  ...
%      'diat_N_lim_Cweight_avg_100m',  ...
%     'diat_P_lim_Cweight_avg_100m',  'diat_SiO3_lim_Cweight_avg_100m', ...
%     'sp_Fe_lim_Cweight_avg_100m',  ...
%     'sp_light_lim_Cweight_avg_100m',   ...
%     'sp_N_lim_Cweight_avg_100m', 'sp_P_lim_Cweight_avg_100m', ...
%      'diaz_Fe_lim_Cweight_avg_100m', ...
%      'diaz_light_lim_Cweight_avg_100m',  'diaz_P_lim_Cweight_avg_100m', ...
%     'zoo_loss_zint_100m', ...
%     'diaz_agg_zint_100m', 'diat_loss_zint_100m', 'diat_agg_zint_100m', 'diaz_loss_zint_100m',  ...
%     'sp_agg_zint_100m', 'sp_loss_zint_100m', ...
%     'graze_diat_zint_100m', 'graze_diat_zoo_zint_100m', 'graze_diaz_zint_100m', ...
%     'graze_diaz_zoo_zint_100m', 'graze_sp_zint_100m', 'graze_sp_zoo_zint_100m', ...
%     'diat_Fe_lim_surf', 'diat_light_lim_surf', 'diat_P_lim_surf', 'diat_SiO3_lim_surf', 'diaz_Fe_lim_surf', ...
%     'diaz_light_lim_surf', 'diaz_P_lim_surf', 'dustToSed', 'LWUP_F', 'O2_ZMIN_DEPTH', 'diat_N_lim_surf', ...
%     'sp_Fe_lim_surf', 'sp_light_lim_surf', 'sp_N_lim_surf', 'sp_P_lim_surf'  };

% cfg.vars = {'TEMPCLINE'};
% cfg.vars = {'UVEL55', 'VVEL55'};
% cfg.vars = {'IRON_FLUX', 'TEMP145', 'UVEL145', 'VVEL145', 'NO3145', 'PO4145', 'Fe145'};

% cfg.vars = {'sumChl'};

cfg.vars = {'NPP', 'GPP', 'TOTVEGC', 'TWS', 'aice', 'sithick'};
cfg.vars = {'AODDUST'};
cfg.vars = {'TS'};
cfg.vars = {'PRECT'};
% cfg.vars = {'SST'};
cfg.vars = {'GPP'};
cfg.vars = {'TWS'};
cfg.vars = {'photoC_TOT_zint_100m'};
cfg.vars = {'SST', 'PRECT', 'TS', 'PSL'};
cfg.vars = {'FIRE'};
cfg.vars = {'TEMP145'};
cfg.vars = {'photoC_TOT_zint_100m'};
cfg.vars = {'DSTFLXT'};
cfg.vars = { 'GPP', 'NPP', 'TOTVEGC', 'TWS', 'COL_FIRE_CLOSS', 'COL_FIRE_NLOSS', 'FIRE', ...
            'FPSN', 'SOILICE', 'SOILLIQ', 'TOTSOILICE', 'TOTSOILLIQ', ...
            'RAIN', 'QSOIL', 'QSOIL_ICE', 'QRUNOFF', 'QOVER', 'QRGWL', ...
            'QH2OSFC', 'NEP', 'DSTFLXT' };
cfg.vars = {'TWS'};
% cfg.vars = {'TS', 'PRECT', 'PSL', 'SST'};
% cfg.vars = {'NPP', 'GPP', 'RAIN', 'TWS'};
% cfg.vars = {'FIRE', 'TOTVEGC'};
cfg.vars = {'TS', 'SST', 'PRECT', 'PSL'};
cfg.vars = {'PSL'};

% cfg.vlayer=1:10; % 10layer. don't put more than 15
cfg.vlayer=1; % surface, vertical slice

cfg.vlayer_1st=min(cfg.vlayer);
cfg.vlayer_cnt=max(cfg.vlayer)-cfg.vlayer_1st+1;


for vari=1:length(cfg.vars)

cfg.var=cfg.vars{vari};
cfg.obs_name=f_obs_name(cfg.var);
cfg.obs_fname_mid=f_obs_name_mid(cfg.var);
cfg.obs_varname=f_obs_varname(cfg.var);
cfg.comp=Func_0025_CESM2_cmpname_var(cfg.var);
cfg.obs_fname_module=f_obs_fname_module(cfg.comp);
cfg.obs_iyears=1960:2020;

disp(cfg.var);
tic;

% dirs.hcstroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/HCST_EXP/archive/', cfg.comp, '/', cfg.var];
% dirs.obsroot=['/Volumes/kyy_raid/kimyy/Observation/', cfg.obs_name, '/monthly_reg_', cfg.obs_fname_module(2:4)];
% dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', cfg.var];
% dirs.lens2root=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/LENS2/archive_analysis/', cfg.comp, '/', cfg.var];
% dirs.assmroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/ASSM_EXP/archive_analysis/', cfg.comp, '/', cfg.var];

dirs.hcstmatroot=[tmp.kimyypath, '/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', cfg.var];
dirs.figroot=[tmp.kimyypath,'/Figure/CESM2/ESP/HCST_EXP/archive_anomaly/', cfg.comp,'/', cfg.var];
mkdir(dirs.figroot)

cfg.iyears=cfg.obs_iyears;
cfg.gnm='f09_g17';
cfg.proj_year=5;

cfg.len_t_y = length(cfg.iyears);
cfg.casename_m = ['ens_all'];

% cfg.len_t_m = length(cfg.months);
% cfg.len_t = cfg.len_t_y *cfg.len_t_m;

% %% grid set(mask from model)
% tmp.obsname=cfg.obsnames{1};
% iyear=min(cfg.iyears);
% cfg.casename_m=[cfg.gridname, '.hcst.', tmp.obsname, '-', cfg.assm_factor, 'p', cfg.ens_member];
% cfg.casename=[cfg.casename_m, '_i', num2str(iyear)];
% dirs.datadir= [dirs.hcstroot, filesep, cfg.casename_m, filesep, 'GMSV'];

% f09_g17.hcst.en4.2_ba-10p1_i2021.pop.h.once.nc

% [tmp.error_status, tmp.value]=system(['ls ', dirs.datadir, '/*once*']);  % b.e21.BHISTsmbb.f09_g17.assm.oras4_ba-10p1.pop.h.once.nc
tmp.gridname = [tmp.kimyypath, '/Model/CESM2/ESP/HCST_EXP/archive_transfer/', cfg.comp, '/grid.nc'];
tmp.maskname = [tmp.kimyypath, '/Model/CESM2/ESP/HCST_EXP/archive_transfer/ocn/RECCAP2_region_masks_all_v20210412_POP2_grid.nc'];

% grid.region_mask=ncread(tmp.gridname, 'REGION_MASK'); 
% grid.ocean_mask=NaN(size(grid.region_mask));
% grid.ocean_mask(grid.region_mask>0)=1;
% grid.tarea = ncread(tmp.gridname, 'TAREA');

switch cfg.comp
    case {'ocn', 'ice'}
        grid.tlong=ncread(tmp.gridname, 'TLONG');
        grid.tlat=ncread(tmp.gridname, 'TLAT');
        grid.mask_ocn=ncread(tmp.maskname, 'open_ocean');
        grid.mask_ocn(grid.mask_ocn<-10e10)=NaN;
        grid.mask_ocn=grid.mask_ocn./grid.mask_ocn;
        grid.tarea=ncread(tmp.gridname, 'TAREA')/1000000.0; %(m^2 -> km^2)
        grid.tarea_60=grid.tarea; grid.tarea_60(grid.tlat>60 | grid.tlat<-60)=NaN;
    case {'atm', 'lnd'}
        grid.lon=ncread(tmp.gridname, 'lon');
        grid.lat=ncread(tmp.gridname, 'lat');
        [grid.tlat grid.tlong]=meshgrid(grid.lat, grid.lon);
        grid.tarea=ncread(tmp.gridname, 'AREA');
        grid.tarea_60=grid.tarea; grid.tarea_60(grid.tlat>60 | grid.tlat<-60)=NaN;
end

grid.nlon=size(grid.tlong,1);
grid.nlat=size(grid.tlat,2);
% grid.ntime=cfg.proj_year.*12;


% % model filename example
% /mnt/lustre/proj/kimyy/Model/CESM2/ESP/HCST_EXP/archive/atm/SST/ens_all/ens_all_i2020
% SST_f09_g17.hcst.ens_all_i2020.cam.h0.2024-12.nc

fig_flags(1:100)=0;

% fig_flags(2:9)=1; %OBS ACC
% fig_flags(17)=1; %OBS trend
% fig_flags(21:26)=1; %OBS ACC
% fig_flags(31)=1; %OBS ensmean
% fig_flags(38:43)=1; %OBS ensmean bias
% fig_flags(35:37)=1;
% fig_flags(1:100)=1;
% fig_flags(48)=1; %model drift
% fig_flags(24)=1; %assm obs nRMSE
% fig_flags(24:26)=1; %assm obs RMSE
% fig_flags(41:43)=1; %assm mean bias
% fig_flags(21:30)=1; % RMSEs
% fig_flags(49)=1; %model drift ratio
% fig_flags(50:51)=1; %obs-AR1 coef, noise
% fig_flags(31:34)=1;
% fig_flags(53:57)=1;
fig_flags(58:59)=1;

%% no observation
% fig_flags(1)=1;
% fig_flags(11:16)=1;
% fig_flags(18:20)=1;
% fig_flags(27:30)=1;
% fig_flags(32:37)=1;
% fig_flags(39)=1;
% fig_flags(44:47)=1;
% fig_flags(48:49)=1;
% fig_flags(13)=1;
% fig_flags(52)=1;

S = shaperead('landareas.shp');

%% read & plot data
tmp.varname=cfg.var;

% % % %% read all data for RMSE colorbar
% % % fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
% % %             '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
% % %             '_l00', 'y.mat'];
% % % load(fig_cfg.mat_name, 'data', 'data2')
% % % data2_l0=data2;
% % % fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
% % %             '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
% % %             '_l01', 'y.mat'];
% % % load(fig_cfg.mat_name, 'data', 'data2')
% % % data2_l1=data2;
% % % fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
% % %             '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
% % %             '_l02', 'y.mat'];
% % % load(fig_cfg.mat_name, 'data', 'data2')
% % % data2_l2=data2;
% % % fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
% % %             '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
% % %             '_l03', 'y.mat'];
% % % load(fig_cfg.mat_name, 'data', 'data2')
% % % data2_l3=data2;
% % % fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
% % %             '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
% % %             '_l04', 'y.mat'];
% % % load(fig_cfg.mat_name, 'data', 'data2')
% % % data2_l4=data2;

clear tmp.ydata tmp.ydata_lens2 tmp.ydata_obs tmp.ydata_assm
for lyear=0:cfg.proj_year-1
    tmp.lyear_str=num2str(lyear, '%02i');
    fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l', tmp.lyear_str, 'y.mat'];
    load(fig_cfg.mat_name, 'data', 'data2')
%     pcolor(data2.([tmp.varname,'_corr_assm_int_l',tmp.lyear_str])'); shading flat; colorbar;

%% clim time range
    switch cfg.obs_name
       case 'GPCC'
           cfg.clim_ys=1965;
           cfg.clim_ye=2019;
       otherwise
           cfg.clim_ys=1965;
           cfg.clim_ye=2020;
   end
   cfg.clim_tlen = (cfg.clim_ys-1959)-lyear:(cfg.clim_ye-2020)+cfg.len_t_y-lyear;
   cfg.clim_tlen2=length(cfg.clim_tlen);


%% NP SVD
%% regional mca
grid.atm_mask_glo=NaN(length(grid.lon), length(grid.lat));

grid.section_tga=[100 280 10 60]; % 170W 120W 5S 5N
[CESM2_grid.obs_lon_ind, CESM2_grid.obs_lon_ind2, CESM2_grid.obs_lat_ind, CESM2_grid.obs_lat_ind2] = ...
    Func_0012_findind_Y(1, grid.section_tga, ...
    grid.lon, grid.lat, 'CESM2'); % find valid lon, lat index near station

grid.atm_mask=grid.atm_mask_glo;
grid.atm_mask(CESM2_grid.obs_lon_ind:CESM2_grid.obs_lon_ind2, CESM2_grid.obs_lat_ind:CESM2_grid.obs_lat_ind2)=1;


tmp.svd_modes=5;
[data3.([tmp.varname, '_svd_np_obs_hcst_lvmap_left', '_l', tmp.lyear_str]), ...
            data3.([tmp.varname, '_svd_np_obs_hcst_pcs_left', '_l', tmp.lyear_str]), ...
            data3.([tmp.varname, '_svd_np_obs_hcst_lvmap_right', '_l', tmp.lyear_str]), ...
            data3.([tmp.varname, '_svd_np_obs_hcst_pcs_right', '_l', tmp.lyear_str]), ...
            data3.([tmp.varname, '_svd_np_obs_hcst_lambda', '_l', tmp.lyear_str]), ...
            data3.([tmp.varname, '_svd_np_obs_hcst_scf', '_l', tmp.lyear_str])] = ...
            mca( data2.([tmp.varname, '_obs_ano', '_l', tmp.lyear_str]).*grid.atm_mask, ... 
                 data2.([tmp.varname, '_model_ano', '_l', tmp.lyear_str]).*grid.atm_mask, tmp.svd_modes);



%% SVD - NP SVD plot
fig_flags(57)=1;
if fig_flags(57)==1
    tmp.mode_want=3;

fig_cfg.name_rgn = 'Glob';
    fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
%     fig_cfg.map_proj = 'robinson';  % robinson, eqdcylin

    fig_cfg.x_lim = [-180 180];
    fig_cfg.y_lim = [-80 89];
    fig_cfg.fig_size = [0,0,6,3.5];
    fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
    fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
    fig_cfg.title_pos = [0.5,0.93];
    fig_cfg.p_lim =0.1;
    fig_cfg.c_lim = [-1 1];
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);


%% PCT (left) plot
    for SVD_mode=1:tmp.mode_want
        fig_cfg.fig_name=['l', tmp.lyear_str, ', SVD_np_obs_hcst_left_pct_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'visible','off');
    %     hold on
        plot(cfg.iyears(cfg.clim_tlen)+lyear,data3.([tmp.varname,'_svd_np_obs_hcst_pcs_left_l',tmp.lyear_str])(SVD_mode,:), 'k-', 'linewidth', 2);
    %         plot(cfg.iyears(1:end-4)+4,EOF.pc(:,SVD_mode));
        grid minor
        xlabel('year'); ylabel(['PCT-mode', num2str(SVD_mode), ', SCF:', num2str(round(data3.([tmp.varname,'_svd_np_obs_hcst_scf_l',tmp.lyear_str])(SVD_mode).*100,1)),'%']);
        set(gca, 'fontsize', 20)
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_obs'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'np_obs_hcst_left_PCT_l', tmp.lyear_str,'_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
    % hold off

%% LV (left) plot
    for SVD_mode=1:tmp.mode_want
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=squeeze(data3.([tmp.varname,'_svd_np_obs_hcst_lvmap_left_l',tmp.lyear_str])(:,:,SVD_mode));
        tmp.C=tmp.C([end, 1:end],:);
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
        
        tmp.prc95 =prctile(tmp.C(isfinite(tmp.C)), 95);
        tmp.prc05 =prctile(tmp.C(isfinite(tmp.C)), 5);

        fig_cfg.fig_name=['l',tmp.lyear_str, ', SVD_LV_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
            'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
        %% map setting
        ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
            'fontname','freeserif'); 
    
        axis off; 
        hold on;
        setm(ax_m,'origin',[0,200],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
        set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',14,'fontname','freeserif','interpreter','none')
    
        %% caxis & colorbar
        tmp.maxval=max(abs([tmp.prc05 tmp.prc95]));
        caxis(ax_m, [-tmp.maxval tmp.maxval]); 
        colormap(fig_cfg.c_map);
        cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
        set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    
        %% draw on ax_m
        h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
        shading flat;
        geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
    
    
    %% frame and label setting
        setm(ax_m,'frame','on','FLineWidth',1);
    
        label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
        label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
        mlabel; plabel;
        label_y=plabel; label_x=mlabel;
        for lxi=1:length(label_x)
            tmp.tmppos=label_x(lxi,1).Position;
            tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55; % y position correction
            label_x(lxi,1).Position=tmp.tmppos;
            label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
        end
        for lyi=1:length(label_y)
            label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
            tmp.tmppos=label_y(lyi,1).Position;
            tmp.tmppos(1)=-fig_cfg.ax_size(3)+2.6; % x position correction
            label_y(lyi,1).Position=tmp.tmppos;
        end
    
        %% save
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_obs'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'np_obs_hcst_left_LV_l', tmp.lyear_str, '_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end


%% PCT (right) plot
    for SVD_mode=1:tmp.mode_want
        fig_cfg.fig_name=['l',tmp.lyear_str, ', SVD_np_obs_hcst_right_pct_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'visible','off');
    %     hold on
        plot(cfg.iyears(cfg.clim_tlen)+lyear,data3.([tmp.varname,'_svd_np_obs_hcst_pcs_right_l',tmp.lyear_str])(SVD_mode,:), 'k-', 'linewidth', 2);
    %         plot(cfg.iyears(1:end-4)+4,EOF.pc(:,SVD_mode));
        grid minor
        xlabel('year'); ylabel(['PCT-mode', num2str(SVD_mode), ', SCF:', num2str(round(data3.([tmp.varname,'_svd_np_obs_hcst_scf_l',tmp.lyear_str])(SVD_mode).*100,1)),'%']);
        set(gca, 'fontsize', 20)
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_obs'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'np_obs_hcst_right_PCT_l', tmp.lyear_str,'_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
    % hold off

%% LV (right) plot
    for SVD_mode=1:tmp.mode_want
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=squeeze(data3.([tmp.varname,'_svd_np_obs_hcst_lvmap_right_l',tmp.lyear_str])(:,:,SVD_mode));
        tmp.C=tmp.C([end, 1:end],:);
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
        
        tmp.prc95 =prctile(tmp.C(isfinite(tmp.C)), 95);
        tmp.prc05 =prctile(tmp.C(isfinite(tmp.C)), 5);

        fig_cfg.fig_name=['l',tmp.lyear_str, ', SVD_LV_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
            'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
        %% map setting
        ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
            'fontname','freeserif'); 
    
        axis off; 
        hold on;
        setm(ax_m,'origin',[0,200],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
        set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',14,'fontname','freeserif','interpreter','none')
    
        %% caxis & colorbar
        tmp.maxval=max(abs([tmp.prc05 tmp.prc95]));
        caxis(ax_m, [-tmp.maxval tmp.maxval]); 
        colormap(fig_cfg.c_map);
        cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
        set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    
        %% draw on ax_m
        h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
        shading flat;
        geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
    
    
    %% frame and label setting
        setm(ax_m,'frame','on','FLineWidth',1);
    
        label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
        label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
        mlabel; plabel;
        label_y=plabel; label_x=mlabel;
        for lxi=1:length(label_x)
            tmp.tmppos=label_x(lxi,1).Position;
            tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55; % y position correction
            label_x(lxi,1).Position=tmp.tmppos;
            label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
        end
        for lyi=1:length(label_y)
            label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
            tmp.tmppos=label_y(lyi,1).Position;
            tmp.tmppos(1)=-fig_cfg.ax_size(3)+2.6; % x position correction
            label_y(lyi,1).Position=tmp.tmppos;
        end
    
        %% save
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_obs'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'np_obs_hcst_right_LV_l', tmp.lyear_str, '_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end

end



%% NA SVD
%% regional mca
grid.atm_mask_glo=NaN(length(grid.lon), length(grid.lat));

grid.section_tga=[280 360 20 70]; % 
[CESM2_grid.obs_lon_ind, CESM2_grid.obs_lon_ind2, CESM2_grid.obs_lat_ind, CESM2_grid.obs_lat_ind2] = ...
    Func_0012_findind_Y(1, grid.section_tga, ...
    grid.lon, grid.lat, 'CESM2'); % find valid lon, lat index near station

grid.atm_mask=grid.atm_mask_glo;
grid.atm_mask(CESM2_grid.obs_lon_ind:CESM2_grid.obs_lon_ind2, CESM2_grid.obs_lat_ind:CESM2_grid.obs_lat_ind2)=1;


tmp.svd_modes=5;
[data3.([tmp.varname, '_svd_na_obs_hcst_lvmap_left', '_l', tmp.lyear_str]), ...
            data3.([tmp.varname, '_svd_na_obs_hcst_pcs_left', '_l', tmp.lyear_str]), ...
            data3.([tmp.varname, '_svd_na_obs_hcst_lvmap_right', '_l', tmp.lyear_str]), ...
            data3.([tmp.varname, '_svd_na_obs_hcst_pcs_right', '_l', tmp.lyear_str]), ...
            data3.([tmp.varname, '_svd_na_obs_hcst_lambda', '_l', tmp.lyear_str]), ...
            data3.([tmp.varname, '_svd_na_obs_hcst_scf', '_l', tmp.lyear_str])] = ...
            mca( data2.([tmp.varname, '_obs_ano', '_l', tmp.lyear_str]).*grid.atm_mask, ... 
                 data2.([tmp.varname, '_model_ano', '_l', tmp.lyear_str]).*grid.atm_mask, tmp.svd_modes);



%% SVD - NP SVD plot
fig_flags(57)=1;
if fig_flags(57)==1
    tmp.mode_want=3;
%% PCT (left) plot
    for SVD_mode=1:tmp.mode_want
        fig_cfg.fig_name=['l', tmp.lyear_str, ', SVD_na_obs_hcst_left_pct_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'visible','off');
    %     hold on
        plot(cfg.iyears(cfg.clim_tlen)+lyear,data3.([tmp.varname,'_svd_na_obs_hcst_pcs_left_l',tmp.lyear_str])(SVD_mode,:), 'k-', 'linewidth', 2);
    %         plot(cfg.iyears(1:end-4)+4,EOF.pc(:,SVD_mode));
        grid minor
        xlabel('year'); ylabel(['PCT-mode', num2str(SVD_mode), ', SCF:', num2str(round(data3.([tmp.varname,'_svd_na_obs_hcst_scf_l',tmp.lyear_str])(SVD_mode).*100,1)),'%']);
        set(gca, 'fontsize', 20)
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_obs'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'na_obs_hcst_left_PCT_l', tmp.lyear_str,'_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
    % hold off

%% LV (left) plot
    for SVD_mode=1:tmp.mode_want
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=squeeze(data3.([tmp.varname,'_svd_na_obs_hcst_lvmap_left_l',tmp.lyear_str])(:,:,SVD_mode));
        tmp.C=tmp.C([end, 1:end],:);
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
        
        tmp.prc95 =prctile(tmp.C(isfinite(tmp.C)), 95);
        tmp.prc05 =prctile(tmp.C(isfinite(tmp.C)), 5);

        fig_cfg.fig_name=['l',tmp.lyear_str, ', SVD_LV_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
            'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
        %% map setting
        ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
            'fontname','freeserif'); 
    
        axis off; 
        hold on;
        setm(ax_m,'origin',[0,200],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
        set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',14,'fontname','freeserif','interpreter','none')
    
        %% caxis & colorbar
        tmp.maxval=max(abs([tmp.prc05 tmp.prc95]));
        caxis(ax_m, [-tmp.maxval tmp.maxval]); 
        colormap(fig_cfg.c_map);
        cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
        set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    
        %% draw on ax_m
        h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
        shading flat;
        geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
    
    
    %% frame and label setting
        setm(ax_m,'frame','on','FLineWidth',1);
    
        label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
        label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
        mlabel; plabel;
        label_y=plabel; label_x=mlabel;
        for lxi=1:length(label_x)
            tmp.tmppos=label_x(lxi,1).Position;
            tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55; % y position correction
            label_x(lxi,1).Position=tmp.tmppos;
            label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
        end
        for lyi=1:length(label_y)
            label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
            tmp.tmppos=label_y(lyi,1).Position;
            tmp.tmppos(1)=-fig_cfg.ax_size(3)+2.6; % x position correction
            label_y(lyi,1).Position=tmp.tmppos;
        end
    
        %% save
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_obs'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'na_obs_hcst_left_LV_l', tmp.lyear_str, '_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end


%% PCT (right) plot
    for SVD_mode=1:tmp.mode_want
        fig_cfg.fig_name=['l',tmp.lyear_str, ', SVD_na_obs_hcst_right_pct_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'visible','off');
    %     hold on
        plot(cfg.iyears(cfg.clim_tlen)+lyear,data3.([tmp.varname,'_svd_na_obs_hcst_pcs_right_l',tmp.lyear_str])(SVD_mode,:), 'k-', 'linewidth', 2);
    %         plot(cfg.iyears(1:end-4)+4,EOF.pc(:,SVD_mode));
        grid minor
        xlabel('year'); ylabel(['PCT-mode', num2str(SVD_mode), ', SCF:', num2str(round(data3.([tmp.varname,'_svd_na_obs_hcst_scf_l',tmp.lyear_str])(SVD_mode).*100,1)),'%']);
        set(gca, 'fontsize', 20)
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_obs'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'na_obs_hcst_right_PCT_l', tmp.lyear_str,'_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
    % hold off

%% LV (right) plot
    for SVD_mode=1:tmp.mode_want
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=squeeze(data3.([tmp.varname,'_svd_na_obs_hcst_lvmap_right_l',tmp.lyear_str])(:,:,SVD_mode));
        tmp.C=tmp.C([end, 1:end],:);
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
        
        tmp.prc95 =prctile(tmp.C(isfinite(tmp.C)), 95);
        tmp.prc05 =prctile(tmp.C(isfinite(tmp.C)), 5);

        fig_cfg.fig_name=['l',tmp.lyear_str, ', SVD_LV_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
            'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
        %% map setting
        ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
            'fontname','freeserif'); 
    
        axis off; 
        hold on;
        setm(ax_m,'origin',[0,200],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
        set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',14,'fontname','freeserif','interpreter','none')
    
        %% caxis & colorbar
        tmp.maxval=max(abs([tmp.prc05 tmp.prc95]));
        caxis(ax_m, [-tmp.maxval tmp.maxval]); 
        colormap(fig_cfg.c_map);
        cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
        set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    
        %% draw on ax_m
        h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
        shading flat;
        geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
    
    
    %% frame and label setting
        setm(ax_m,'frame','on','FLineWidth',1);
    
        label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
        label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
        mlabel; plabel;
        label_y=plabel; label_x=mlabel;
        for lxi=1:length(label_x)
            tmp.tmppos=label_x(lxi,1).Position;
            tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55; % y position correction
            label_x(lxi,1).Position=tmp.tmppos;
            label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
        end
        for lyi=1:length(label_y)
            label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
            tmp.tmppos=label_y(lyi,1).Position;
            tmp.tmppos(1)=-fig_cfg.ax_size(3)+2.6; % x position correction
            label_y(lyi,1).Position=tmp.tmppos;
        end
    
        %% save
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_obs'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'na_obs_hcst_right_LV_l', tmp.lyear_str, '_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end

end









end


%% multi lyear
tmp.lyear_str=['01_04'];
    fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l', tmp.lyear_str, 'y.mat'];
    load(fig_cfg.mat_name, 'data', 'data2')
%     pcolor(data2.([tmp.varname,'_corr_assm_int_l',tmp.lyear_str])'); shading flat; colorbar;

%% clim time range
    switch cfg.obs_name
       case 'GPCC'
           cfg.clim_ys=1965;
           cfg.clim_ye=2019;
       otherwise
           cfg.clim_ys=1965;
           cfg.clim_ye=2020;
   end
   cfg.clim_tlen = (cfg.clim_ys-1959)-lyear:(cfg.clim_ye-2020)+cfg.len_t_y-lyear;
   cfg.clim_tlen2=length(cfg.clim_tlen);


%% NP SVD
%% regional mca
grid.atm_mask_glo=NaN(length(grid.lon), length(grid.lat));

grid.section_tga=[100 280 10 60]; % 170W 120W 5S 5N
[CESM2_grid.obs_lon_ind, CESM2_grid.obs_lon_ind2, CESM2_grid.obs_lat_ind, CESM2_grid.obs_lat_ind2] = ...
    Func_0012_findind_Y(1, grid.section_tga, ...
    grid.lon, grid.lat, 'CESM2'); % find valid lon, lat index near station

grid.atm_mask=grid.atm_mask_glo;
grid.atm_mask(CESM2_grid.obs_lon_ind:CESM2_grid.obs_lon_ind2, CESM2_grid.obs_lat_ind:CESM2_grid.obs_lat_ind2)=1;


tmp.svd_modes=5;
[data3.([tmp.varname, '_svd_np_obs_hcst_lvmap_left', '_l', tmp.lyear_str]), ...
            data3.([tmp.varname, '_svd_np_obs_hcst_pcs_left', '_l', tmp.lyear_str]), ...
            data3.([tmp.varname, '_svd_np_obs_hcst_lvmap_right', '_l', tmp.lyear_str]), ...
            data3.([tmp.varname, '_svd_np_obs_hcst_pcs_right', '_l', tmp.lyear_str]), ...
            data3.([tmp.varname, '_svd_np_obs_hcst_lambda', '_l', tmp.lyear_str]), ...
            data3.([tmp.varname, '_svd_np_obs_hcst_scf', '_l', tmp.lyear_str])] = ...
            mca( data2.([tmp.varname, '_obs_ano', '_l', tmp.lyear_str]).*grid.atm_mask, ... 
                 data2.([tmp.varname, '_model_ano', '_l', tmp.lyear_str]).*grid.atm_mask, tmp.svd_modes);



%% SVD - NP SVD plot
fig_flags(57)=1;
if fig_flags(57)==1
    tmp.mode_want=3;

fig_cfg.name_rgn = 'Glob';
    fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
%     fig_cfg.map_proj = 'robinson';  % robinson, eqdcylin

    fig_cfg.x_lim = [-180 180];
    fig_cfg.y_lim = [-80 89];
    fig_cfg.fig_size = [0,0,6,3.5];
    fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
    fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
    fig_cfg.title_pos = [0.5,0.93];
    fig_cfg.p_lim =0.1;
    fig_cfg.c_lim = [-1 1];
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);


%% PCT (left) plot
    for SVD_mode=1:tmp.mode_want
        fig_cfg.fig_name=['l', tmp.lyear_str, ', SVD_np_obs_hcst_left_pct_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'visible','off');
    %     hold on
        plot(cfg.iyears(cfg.clim_tlen)+lyear,data3.([tmp.varname,'_svd_np_obs_hcst_pcs_left_l',tmp.lyear_str])(SVD_mode,:), 'k-', 'linewidth', 2);
    %         plot(cfg.iyears(1:end-4)+4,EOF.pc(:,SVD_mode));
        grid minor
        xlabel('year'); ylabel(['PCT-mode', num2str(SVD_mode), ', SCF:', num2str(round(data3.([tmp.varname,'_svd_np_obs_hcst_scf_l',tmp.lyear_str])(SVD_mode).*100,1)),'%']);
        set(gca, 'fontsize', 20)
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_obs'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'np_obs_hcst_left_PCT_l', tmp.lyear_str,'_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
    % hold off

%% LV (left) plot
    for SVD_mode=1:tmp.mode_want
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=squeeze(data3.([tmp.varname,'_svd_np_obs_hcst_lvmap_left_l',tmp.lyear_str])(:,:,SVD_mode));
        tmp.C=tmp.C([end, 1:end],:);
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
        
        tmp.prc95 =prctile(tmp.C(isfinite(tmp.C)), 95);
        tmp.prc05 =prctile(tmp.C(isfinite(tmp.C)), 5);

        fig_cfg.fig_name=['l',tmp.lyear_str, ', SVD_LV_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
            'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
        %% map setting
        ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
            'fontname','freeserif'); 
    
        axis off; 
        hold on;
        setm(ax_m,'origin',[0,200],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
        set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',14,'fontname','freeserif','interpreter','none')
    
        %% caxis & colorbar
        tmp.maxval=max(abs([tmp.prc05 tmp.prc95]));
        caxis(ax_m, [-tmp.maxval tmp.maxval]); 
        colormap(fig_cfg.c_map);
        cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
        set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    
        %% draw on ax_m
        h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
        shading flat;
        geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
    
    
    %% frame and label setting
        setm(ax_m,'frame','on','FLineWidth',1);
    
        label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
        label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
        mlabel; plabel;
        label_y=plabel; label_x=mlabel;
        for lxi=1:length(label_x)
            tmp.tmppos=label_x(lxi,1).Position;
            tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55; % y position correction
            label_x(lxi,1).Position=tmp.tmppos;
            label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
        end
        for lyi=1:length(label_y)
            label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
            tmp.tmppos=label_y(lyi,1).Position;
            tmp.tmppos(1)=-fig_cfg.ax_size(3)+2.6; % x position correction
            label_y(lyi,1).Position=tmp.tmppos;
        end
    
        %% save
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_obs'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'np_obs_hcst_left_LV_l', tmp.lyear_str, '_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end


%% PCT (right) plot
    for SVD_mode=1:tmp.mode_want
        fig_cfg.fig_name=['l',tmp.lyear_str, ', SVD_np_obs_hcst_right_pct_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'visible','off');
    %     hold on
        plot(cfg.iyears(cfg.clim_tlen)+lyear,data3.([tmp.varname,'_svd_np_obs_hcst_pcs_right_l',tmp.lyear_str])(SVD_mode,:), 'k-', 'linewidth', 2);
    %         plot(cfg.iyears(1:end-4)+4,EOF.pc(:,SVD_mode));
        grid minor
        xlabel('year'); ylabel(['PCT-mode', num2str(SVD_mode), ', SCF:', num2str(round(data3.([tmp.varname,'_svd_np_obs_hcst_scf_l',tmp.lyear_str])(SVD_mode).*100,1)),'%']);
        set(gca, 'fontsize', 20)
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_obs'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'np_obs_hcst_right_PCT_l', tmp.lyear_str,'_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
    % hold off

%% LV (right) plot
    for SVD_mode=1:tmp.mode_want
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=squeeze(data3.([tmp.varname,'_svd_np_obs_hcst_lvmap_right_l',tmp.lyear_str])(:,:,SVD_mode));
        tmp.C=tmp.C([end, 1:end],:);
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
        
        tmp.prc95 =prctile(tmp.C(isfinite(tmp.C)), 95);
        tmp.prc05 =prctile(tmp.C(isfinite(tmp.C)), 5);

        fig_cfg.fig_name=['l',tmp.lyear_str, ', SVD_LV_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
            'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
        %% map setting
        ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
            'fontname','freeserif'); 
    
        axis off; 
        hold on;
        setm(ax_m,'origin',[0,200],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
        set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',14,'fontname','freeserif','interpreter','none')
    
        %% caxis & colorbar
        tmp.maxval=max(abs([tmp.prc05 tmp.prc95]));
        caxis(ax_m, [-tmp.maxval tmp.maxval]); 
        colormap(fig_cfg.c_map);
        cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
        set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    
        %% draw on ax_m
        h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
        shading flat;
        geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
    
    
    %% frame and label setting
        setm(ax_m,'frame','on','FLineWidth',1);
    
        label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
        label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
        mlabel; plabel;
        label_y=plabel; label_x=mlabel;
        for lxi=1:length(label_x)
            tmp.tmppos=label_x(lxi,1).Position;
            tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55; % y position correction
            label_x(lxi,1).Position=tmp.tmppos;
            label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
        end
        for lyi=1:length(label_y)
            label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
            tmp.tmppos=label_y(lyi,1).Position;
            tmp.tmppos(1)=-fig_cfg.ax_size(3)+2.6; % x position correction
            label_y(lyi,1).Position=tmp.tmppos;
        end
    
        %% save
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_obs'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'np_obs_hcst_right_LV_l', tmp.lyear_str, '_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end

end



%% NA SVD
%% regional mca
grid.atm_mask_glo=NaN(length(grid.lon), length(grid.lat));

grid.section_tga=[280 360 20 70]; % 
[CESM2_grid.obs_lon_ind, CESM2_grid.obs_lon_ind2, CESM2_grid.obs_lat_ind, CESM2_grid.obs_lat_ind2] = ...
    Func_0012_findind_Y(1, grid.section_tga, ...
    grid.lon, grid.lat, 'CESM2'); % find valid lon, lat index near station

grid.atm_mask=grid.atm_mask_glo;
grid.atm_mask(CESM2_grid.obs_lon_ind:CESM2_grid.obs_lon_ind2, CESM2_grid.obs_lat_ind:CESM2_grid.obs_lat_ind2)=1;


tmp.svd_modes=5;
[data3.([tmp.varname, '_svd_na_obs_hcst_lvmap_left', '_l', tmp.lyear_str]), ...
            data3.([tmp.varname, '_svd_na_obs_hcst_pcs_left', '_l', tmp.lyear_str]), ...
            data3.([tmp.varname, '_svd_na_obs_hcst_lvmap_right', '_l', tmp.lyear_str]), ...
            data3.([tmp.varname, '_svd_na_obs_hcst_pcs_right', '_l', tmp.lyear_str]), ...
            data3.([tmp.varname, '_svd_na_obs_hcst_lambda', '_l', tmp.lyear_str]), ...
            data3.([tmp.varname, '_svd_na_obs_hcst_scf', '_l', tmp.lyear_str])] = ...
            mca( data2.([tmp.varname, '_obs_ano', '_l', tmp.lyear_str]).*grid.atm_mask, ... 
                 data2.([tmp.varname, '_model_ano', '_l', tmp.lyear_str]).*grid.atm_mask, tmp.svd_modes);



%% SVD - NP SVD plot
fig_flags(57)=1;
if fig_flags(57)==1
    tmp.mode_want=3;
%% PCT (left) plot
    for SVD_mode=1:tmp.mode_want
        fig_cfg.fig_name=['l', tmp.lyear_str, ', SVD_na_obs_hcst_left_pct_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'visible','off');
    %     hold on
        plot(cfg.iyears(cfg.clim_tlen)+lyear,data3.([tmp.varname,'_svd_na_obs_hcst_pcs_left_l',tmp.lyear_str])(SVD_mode,:), 'k-', 'linewidth', 2);
    %         plot(cfg.iyears(1:end-4)+4,EOF.pc(:,SVD_mode));
        grid minor
        xlabel('year'); ylabel(['PCT-mode', num2str(SVD_mode), ', SCF:', num2str(round(data3.([tmp.varname,'_svd_na_obs_hcst_scf_l',tmp.lyear_str])(SVD_mode).*100,1)),'%']);
        set(gca, 'fontsize', 20)
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_obs'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'na_obs_hcst_left_PCT_l', tmp.lyear_str,'_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
    % hold off

%% LV (left) plot
    for SVD_mode=1:tmp.mode_want
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=squeeze(data3.([tmp.varname,'_svd_na_obs_hcst_lvmap_left_l',tmp.lyear_str])(:,:,SVD_mode));
        tmp.C=tmp.C([end, 1:end],:);
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
        
        tmp.prc95 =prctile(tmp.C(isfinite(tmp.C)), 95);
        tmp.prc05 =prctile(tmp.C(isfinite(tmp.C)), 5);

        fig_cfg.fig_name=['l',tmp.lyear_str, ', SVD_LV_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
            'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
        %% map setting
        ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
            'fontname','freeserif'); 
    
        axis off; 
        hold on;
        setm(ax_m,'origin',[0,200],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
        set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',14,'fontname','freeserif','interpreter','none')
    
        %% caxis & colorbar
        tmp.maxval=max(abs([tmp.prc05 tmp.prc95]));
        caxis(ax_m, [-tmp.maxval tmp.maxval]); 
        colormap(fig_cfg.c_map);
        cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
        set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    
        %% draw on ax_m
        h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
        shading flat;
        geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
    
    
    %% frame and label setting
        setm(ax_m,'frame','on','FLineWidth',1);
    
        label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
        label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
        mlabel; plabel;
        label_y=plabel; label_x=mlabel;
        for lxi=1:length(label_x)
            tmp.tmppos=label_x(lxi,1).Position;
            tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55; % y position correction
            label_x(lxi,1).Position=tmp.tmppos;
            label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
        end
        for lyi=1:length(label_y)
            label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
            tmp.tmppos=label_y(lyi,1).Position;
            tmp.tmppos(1)=-fig_cfg.ax_size(3)+2.6; % x position correction
            label_y(lyi,1).Position=tmp.tmppos;
        end
    
        %% save
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_obs'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'na_obs_hcst_left_LV_l', tmp.lyear_str, '_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end


%% PCT (right) plot
    for SVD_mode=1:tmp.mode_want
        fig_cfg.fig_name=['l',tmp.lyear_str, ', SVD_na_obs_hcst_right_pct_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'visible','off');
    %     hold on
        plot(cfg.iyears(cfg.clim_tlen)+lyear,data3.([tmp.varname,'_svd_na_obs_hcst_pcs_right_l',tmp.lyear_str])(SVD_mode,:), 'k-', 'linewidth', 2);
    %         plot(cfg.iyears(1:end-4)+4,EOF.pc(:,SVD_mode));
        grid minor
        xlabel('year'); ylabel(['PCT-mode', num2str(SVD_mode), ', SCF:', num2str(round(data3.([tmp.varname,'_svd_na_obs_hcst_scf_l',tmp.lyear_str])(SVD_mode).*100,1)),'%']);
        set(gca, 'fontsize', 20)
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_obs'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'na_obs_hcst_right_PCT_l', tmp.lyear_str,'_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
    % hold off

%% LV (right) plot
    for SVD_mode=1:tmp.mode_want
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=squeeze(data3.([tmp.varname,'_svd_na_obs_hcst_lvmap_right_l',tmp.lyear_str])(:,:,SVD_mode));
        tmp.C=tmp.C([end, 1:end],:);
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
        
        tmp.prc95 =prctile(tmp.C(isfinite(tmp.C)), 95);
        tmp.prc05 =prctile(tmp.C(isfinite(tmp.C)), 5);

        fig_cfg.fig_name=['l',tmp.lyear_str, ', SVD_LV_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
            'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
        %% map setting
        ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
            'fontname','freeserif'); 
    
        axis off; 
        hold on;
        setm(ax_m,'origin',[0,200],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
        set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',14,'fontname','freeserif','interpreter','none')
    
        %% caxis & colorbar
        tmp.maxval=max(abs([tmp.prc05 tmp.prc95]));
        caxis(ax_m, [-tmp.maxval tmp.maxval]); 
        colormap(fig_cfg.c_map);
        cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
        set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    
        %% draw on ax_m
        h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
        shading flat;
        geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
    
    
    %% frame and label setting
        setm(ax_m,'frame','on','FLineWidth',1);
    
        label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
        label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
        mlabel; plabel;
        label_y=plabel; label_x=mlabel;
        for lxi=1:length(label_x)
            tmp.tmppos=label_x(lxi,1).Position;
            tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55; % y position correction
            label_x(lxi,1).Position=tmp.tmppos;
            label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
        end
        for lyi=1:length(label_y)
            label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
            tmp.tmppos=label_y(lyi,1).Position;
            tmp.tmppos(1)=-fig_cfg.ax_size(3)+2.6; % x position correction
            label_y(lyi,1).Position=tmp.tmppos;
        end
    
        %% save
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_obs'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'na_obs_hcst_right_LV_l', tmp.lyear_str, '_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end


end

end

function varname_C = f_varname_C(varname)
    switch varname
        case 'temp'
            varname_C='TEMP';
        case 'salt'
            varname_C='SALT';
    end
end



function obsname_simple = f_obs_name(varn)
    switch varn
        case 'SST'
            obsname_simple='ERSST';
        case 'PRECT'
%             obsname_simple='GPCP';
            obsname_simple='GPCC';
        case 'PSL'
            obsname_simple='ERA5';
        case 'SSH'
            obsname_simple='CMEMS';
        case 'TS'
%             obsname_simple='HadCRUT5';
            obsname_simple='ERA5';
        case 'sumChl'
            obsname_simple='OC_CCI';
        otherwise
            obsname_simple='nan';
    end
end


function obsname_simple = f_obs_name_mid(varn)
    switch varn
        case 'SST'
            obsname_simple='ersst_reg_cesm2.v5.';
        case 'PRECT'
            obsname_simple='GPCP_reg_cesm2.v5.';
        case 'PSL'
            obsname_simple='ERA5_msl_reg_cesm2.';
        case 'SSH'
            obsname_simple='CMEMS_reg_cesm2.';
        case 'TS'
            obsname_simple='HadCRUT5_reg_cesm2.';
        case 'sumChl'
            obsname_simple='OC_CCI_reg_cesm2.';
        otherwise
            obsname_simple='nan';
    end
end

function obsname_simple = f_obs_varname(varn)
    switch varn
        case 'SST'
            obsname_simple='sst';
        case 'PRECT'
            obsname_simple='precip';
        case 'PSL'
            obsname_simple='msl';
        case 'SSH'
            obsname_simple='sla';
        case 'TS'
            obsname_simple='tas_mean';
        case 'sumChl'
            obsname_simple='chlor_a';
        otherwise
            obsname_simple='nan';
    end
end

function obsname_simple = f_obs_fname_module(comp)
    switch comp
        case 'ocn'
            obsname_simple='.pop.h.';
        case 'atm'
            obsname_simple='.cam.h0.';
        case 'lnd'
            obsname_simple='.clm2.h0.';
        case 'ice'
            obsname_simple='.cice.h.';
    end
end

function obsname_simple = f_obs_iyears(varn)
    switch varn
        case 'PRECT'
            obsname_simple=1979:2020;
        case 'SSH'
            obsname_simple=1993:2020;
        case 'sumChl'
            obsname_simple=1998:2020;
        otherwise
            obsname_simple=1970:2020;
    end
end