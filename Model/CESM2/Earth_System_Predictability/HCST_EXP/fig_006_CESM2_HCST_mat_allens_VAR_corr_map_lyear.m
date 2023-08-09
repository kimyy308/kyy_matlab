% %  Created 12-Apr-2023 by Yong-Yub Kim
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
cfg.vars = {'SST'};


cfg.vlayer=1:10; % 10layer. don't put more than 15
% cfg.vlayer=1; % surface, vertical slice

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
dirs.figroot=[tmp.kimyypath,'/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', cfg.var];



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

S = shaperead('landareas.shp');

%% read & plot data
tmp.varname=cfg.var;

clear tmp.ydata tmp.ydata_lens2 tmp.ydata_obs tmp.ydata_assm
for lyear=0:cfg.proj_year-1
    tmp.lyear_str=num2str(lyear, '%02i');
    fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l', tmp.lyear_str, 'y.mat'];
    load(fig_cfg.mat_name, 'data')
%     pcolor(data.([tmp.varname,'_corr_assm_int_l',tmp.lyear_str])'); shading flat; colorbar;
    
%% obs & model corr map --------------------------------------
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

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data.([tmp.varname, '_corr_obs', '_l', tmp.lyear_str]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data.([tmp.varname, '_corr_AR1', '_l', tmp.lyear_str]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_obs', '_l', tmp.lyear_str]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=['l',tmp.lyear_str, ', corr_obs_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, fig_cfg.c_lim); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,'R','fontsize',12);

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

%         if (strcmp(tmp.season, 'INI')~=1 & strcmp(cfg.var, 'sumChl')~=1) %for Chls, it uses all available data
%             %% <AR1 area -> hatch
%             pp2 = pcolorm(tmp.Y,tmp.X,tmp.C_2, 'parent', ax_m);
%             set(pp2,'linestyle','none','Tag','HatchingRegion');
%             hp = findobj(pp2,'Tag','HatchingRegion');
%             hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);
%         end


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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 

%% obs & assm corr map --------------------------------------
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

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data.([tmp.varname, '_corr_obs_assm', '_l', tmp.lyear_str]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data.([tmp.varname, '_corr_AR1', '_l', tmp.lyear_str]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_obs_assm', '_l', tmp.lyear_str]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=['l',tmp.lyear_str, ', corr_obs_assm_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, fig_cfg.c_lim); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,'R','fontsize',12);

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

%         if (strcmp(tmp.season, 'INI')~=1 & strcmp(cfg.var, 'sumChl')~=1) %for Chls, it uses all available data
%             %% <AR1 area -> hatch
%             pp2 = pcolorm(tmp.Y,tmp.X,tmp.C_2, 'parent', ax_m);
%             set(pp2,'linestyle','none','Tag','HatchingRegion');
%             hp = findobj(pp2,'Tag','HatchingRegion');
%             hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);
%         end


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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_map', filesep, 'assm'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_assm_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 


%% obs & lens2 corr map --------------------------------------
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

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data.([tmp.varname, '_corr_obs_lens2', '_l', tmp.lyear_str]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data.([tmp.varname, '_corr_AR1', '_l', tmp.lyear_str]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_obs_lens2', '_l', tmp.lyear_str]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=['l',tmp.lyear_str, ', corr_obs_lens2_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, fig_cfg.c_lim); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,'R','fontsize',12);

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

%         if (strcmp(tmp.season, 'INI')~=1 & strcmp(cfg.var, 'sumChl')~=1) %for Chls, it uses all available data
%             %% <AR1 area -> hatch
%             pp2 = pcolorm(tmp.Y,tmp.X,tmp.C_2, 'parent', ax_m);
%             set(pp2,'linestyle','none','Tag','HatchingRegion');
%             hp = findobj(pp2,'Tag','HatchingRegion');
%             hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);
%         end


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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_map', filesep, 'lens2'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_lens2_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 


%% obs & hcst-lens2 corr map --------------------------------------
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

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data.([tmp.varname, '_corr_obs_int', '_l', tmp.lyear_str]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data.([tmp.varname, '_corr_AR1', '_l', tmp.lyear_str]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_obs_int', '_l', tmp.lyear_str]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=['l',tmp.lyear_str, ', corr_obs_int_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, fig_cfg.c_lim); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,'R','fontsize',12);

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

%         if (strcmp(tmp.season, 'INI')~=1 & strcmp(cfg.var, 'sumChl')~=1) %for Chls, it uses all available data
%             %% <AR1 area -> hatch
%             pp2 = pcolorm(tmp.Y,tmp.X,tmp.C_2, 'parent', ax_m);
%             set(pp2,'linestyle','none','Tag','HatchingRegion');
%             hp = findobj(pp2,'Tag','HatchingRegion');
%             hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);
%         end


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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_int_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;


%% model & assm corr map --------------------------------------
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

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data.([tmp.varname, '_corr_assm', '_l', tmp.lyear_str]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data.([tmp.varname, '_corr_AR1', '_l', tmp.lyear_str]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_assm', '_l', tmp.lyear_str]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=['l',tmp.lyear_str, ', corr_assm_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, fig_cfg.c_lim); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,'R','fontsize',12);

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

%         if (strcmp(tmp.season, 'INI')~=1 & strcmp(cfg.var, 'sumChl')~=1) %for Chls, it uses all available data
%             %% <AR1 area -> hatch
%             pp2 = pcolorm(tmp.Y,tmp.X,tmp.C_2, 'parent', ax_m);
%             set(pp2,'linestyle','none','Tag','HatchingRegion');
%             hp = findobj(pp2,'Tag','HatchingRegion');
%             hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);
%         end


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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_assm_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_assm_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;


 %% model-lens2 & assm corr map --------------------------------------
    fig_cfg.name_rgn = 'Glob';
    fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
    fig_cfg.x_lim = [-180 180];
    fig_cfg.y_lim = [-80 89];
    fig_cfg.fig_size = [0,0,6,3.5];
    fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
    fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
    fig_cfg.title_pos = [0.5,0.93];
    fig_cfg.p_lim =0.1;
    fig_cfg.c_lim = [-1 1];
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data.([tmp.varname, '_corr_assm_int', '_l', tmp.lyear_str]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data.([tmp.varname, '_corr_AR1', '_l', tmp.lyear_str]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_assm_int', '_l', tmp.lyear_str]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=['l',tmp.lyear_str, ', corr_assm _, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
    fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
        'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
    %% map setting
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 

    axis off; 
    hold on;
    setm(ax_m,'origin',[0,200],'MapLatLimit',fig_cfg.y_lim);
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
    text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
    'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
    'fontsize',14,'fontname','freeserif','interpreter','none')

    %% caxis & colorbar
    caxis(ax_m, fig_cfg.c_lim); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,'R','fontsize',12);

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

%         if (strcmp(tmp.season, 'INI')~=1 & strcmp(cfg.var, 'sumChl')~=1) %for Chls, it uses all available data
%             %% <AR1 area -> hatch
%             pp2 = pcolorm(tmp.Y,tmp.X,tmp.C_2, 'parent', ax_m);
%             set(pp2,'linestyle','none','Tag','HatchingRegion');
%             hp = findobj(pp2,'Tag','HatchingRegion');
%             hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);
%         end


    %% frame and label setting
    setm(ax_m,'frame','on','FLineWidth',1);

    label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
    label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
    mlabel; plabel;
    label_y=plabel; label_x=mlabel;
    for lxi=1:length(label_x)
        tmp.tmppos=label_x(lxi,1).Position;
        tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_assm_int_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_assm_int_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];

    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;


%% lens2 & assm corr map --------------------------------------
    fig_cfg.name_rgn = 'Glob';
    fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
    fig_cfg.x_lim = [-180 180];
    fig_cfg.y_lim = [-80 89];
    fig_cfg.fig_size = [0,0,6,3.5];
    fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
    fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
    fig_cfg.title_pos = [0.5,0.93];
    fig_cfg.p_lim =0.1;
    fig_cfg.c_lim = [-1 1];
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data.([tmp.varname, '_corr_assm_lens2', '_l', tmp.lyear_str]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data.([tmp.varname, '_corr_AR1', '_l', tmp.lyear_str]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_assm_lens2', '_l', tmp.lyear_str]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=['l',tmp.lyear_str, ', corr_assm_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
    fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
        'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
    %% map setting
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 

    axis off; 
    hold on;
    setm(ax_m,'origin',[0,200],'MapLatLimit',fig_cfg.y_lim);
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
    text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
    'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
    'fontsize',14,'fontname','freeserif','interpreter','none')

    %% caxis & colorbar
    caxis(ax_m, fig_cfg.c_lim); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,'R','fontsize',12);

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

%         if (strcmp(tmp.season, 'INI')~=1 & strcmp(cfg.var, 'sumChl')~=1) %for Chls, it uses all available data
%             %% <AR1 area -> hatch
%             pp2 = pcolorm(tmp.Y,tmp.X,tmp.C_2, 'parent', ax_m);
%             set(pp2,'linestyle','none','Tag','HatchingRegion');
%             hp = findobj(pp2,'Tag','HatchingRegion');
%             hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);
%         end


    %% frame and label setting
    setm(ax_m,'frame','on','FLineWidth',1);

    label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
    label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
    mlabel; plabel;
    label_y=plabel; label_x=mlabel;
    for lxi=1:length(label_x)
        tmp.tmppos=label_x(lxi,1).Position;
        tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_assm_ext_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_assm_ext_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;

    

%% assm ensemble mean map --------------------------------------
    data.assm_mean=mean(data.([tmp.varname,'_assm']),3, 'omitnan');
    if strcmp(tmp.varname, 'SST')
        data.assm_mean(data.assm_mean==0)=NaN;
        data.assm_mean=data.assm_mean-273.15;
    end

    tmp.prc95 =prctile(data.assm_mean(:), 95);
    tmp.prc05 =prctile(data.assm_mean(:), 5);
    
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

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data.assm_mean;
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data.([tmp.varname, '_corr_AR1', '_l', tmp.lyear_str]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    fig_cfg.fig_name=['l',tmp.lyear_str, ', ensmean_, ', tmp. varname];
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
    caxis(ax_m, [tmp.prc05 tmp.prc95]); 
%     colormap(fig_cfg.c_map);
    colormap(flip(autumn(10)));
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,'R','fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_ensmean_map', filesep, 'assm'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'ensmean_assm_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;




    %% model ensemble mean map --------------------------------------
    data.model_mean=mean(data.([tmp.varname,'_model_l', tmp.lyear_str]),3, 'omitnan');
    if strcmp(tmp.varname, 'SST')
        data.model_mean(data.model_mean==0)=NaN;
        data.model_mean=data.model_mean-273.15;
    end
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

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data.model_mean;
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data.([tmp.varname, '_corr_AR1', '_l', tmp.lyear_str]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    fig_cfg.fig_name=['l',tmp.lyear_str, ', ensmean_, ', tmp. varname];
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
    caxis(ax_m, [tmp.prc05 tmp.prc95]); 
%     colormap(fig_cfg.c_map);
    colormap(flip(autumn(10)));
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,'R','fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_ensmean_map', filesep, 'hcst'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'ensmean_hcst_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;



    


     %% lens2 ensemble mean map --------------------------------------
    data.lens2_mean=mean(data.([tmp.varname,'_lens2_l', tmp.lyear_str]),3, 'omitnan');
    if strcmp(tmp.varname, 'SST')
        data.model_mean(data.lens2_mean==0)=NaN;
        data.model_mean=data.lens2_mean-273.15;
    end
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

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data.lens2_mean;
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data.([tmp.varname, '_corr_AR1', '_l', tmp.lyear_str]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    fig_cfg.fig_name=['l',tmp.lyear_str, ', ensmean_, ', tmp. varname];
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
    caxis(ax_m, [tmp.prc05 tmp.prc95]); 
%     colormap(fig_cfg.c_map);
    colormap(flip(autumn(10)));
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,'R','fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_ensmean_map', filesep, 'lens2'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'ensmean_lens2_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
    


%% for histogram (assm)

    histind=1;
    for corrref=-1:0.1:-0.1
        %% hcst_int
        corr_histo.(['hcst_int_ind_l', tmp.lyear_str]){histind} = ...
            find(data.([tmp.varname, '_corr_assm_int', '_l', tmp.lyear_str])>=corrref & ...
            data.([tmp.varname, '_corr_assm_int', '_l', tmp.lyear_str])<corrref+0.1);
        corr_histo.(['hcst_int_area_l', tmp.lyear_str])(histind) = ...
            sum(grid.tarea_60(corr_histo.(['hcst_int_ind_l', tmp.lyear_str]){histind}), 'omitnan');
        %% lens2
        corr_histo.(['lens2_ind_l', tmp.lyear_str]){histind} = ...
            find(data.([tmp.varname, '_corr_assm_lens2', '_l', tmp.lyear_str])>=corrref & ...
            data.([tmp.varname, '_corr_assm_lens2', '_l', tmp.lyear_str])<corrref+0.1);
        corr_histo.(['lens2_area_l', tmp.lyear_str])(histind) = ...
            sum(grid.tarea_60(corr_histo.(['lens2_ind_l', tmp.lyear_str]){histind}), 'omitnan');
        histind=histind+1;
    end

    for corrref=0.1:0.1:1.0
        %% hcst_int
        corr_histo.(['hcst_int_ind_l', tmp.lyear_str]){histind} = ...
            find(data.([tmp.varname, '_corr_assm_int', '_l', tmp.lyear_str])>corrref-0.1 & ...
            data.([tmp.varname, '_corr_assm_int', '_l', tmp.lyear_str])<corrref);
        corr_histo.(['hcst_int_area_l', tmp.lyear_str])(histind) = ...
            sum(grid.tarea_60(corr_histo.(['hcst_int_ind_l', tmp.lyear_str]){histind}), 'omitnan');

        %% lens2
        corr_histo.(['lens2_ind_l', tmp.lyear_str]){histind} = ...
            find(data.([tmp.varname, '_corr_assm_lens2', '_l', tmp.lyear_str])>corrref-0.1 & ...
            data.([tmp.varname, '_corr_assm_lens2', '_l', tmp.lyear_str])<=corrref);
        corr_histo.(['lens2_area_l', tmp.lyear_str])(histind) = ...
            sum(grid.tarea_60(corr_histo.(['lens2_ind_l', tmp.lyear_str]){histind}), 'omitnan');
        histind=histind+1;
    end
    corr_ref_x=-0.95:0.1:0.95;

%% for histogram (obs)
    histind=1;
    for corrref=-1:0.1:-0.1
        %hcst_int
        corr_histo.(['obs_hcst_int_ind_l', tmp.lyear_str]){histind} = ...
            find(data.([tmp.varname, '_corr_obs_int', '_l', tmp.lyear_str])>=corrref & ...
            data.([tmp.varname, '_corr_obs_int', '_l', tmp.lyear_str])<corrref+0.1);
        corr_histo.(['obs_hcst_int_area_l', tmp.lyear_str])(histind) = ...
            sum(grid.tarea_60(corr_histo.(['obs_hcst_int_ind_l', tmp.lyear_str]){histind}), 'omitnan');
        %lens2
        corr_histo.(['obs_lens2_ind_l', tmp.lyear_str]){histind} = ...
            find(data.([tmp.varname, '_corr_obs_lens2', '_l', tmp.lyear_str])>=corrref & ...
            data.([tmp.varname, '_corr_obs_lens2', '_l', tmp.lyear_str])<corrref+0.1);
        corr_histo.(['obs_lens2_area_l', tmp.lyear_str])(histind) = ...
            sum(grid.tarea_60(corr_histo.(['obs_lens2_ind_l', tmp.lyear_str]){histind}), 'omitnan');
        histind=histind+1;
    end

    for corrref=0.1:0.1:1.0
        %hcst_int
        corr_histo.(['obs_hcst_int_ind_l', tmp.lyear_str]){histind} = ...
            find(data.([tmp.varname, '_corr_obs_int', '_l', tmp.lyear_str])>corrref-0.1 & ...
            data.([tmp.varname, '_corr_obs_int', '_l', tmp.lyear_str])<corrref);
        corr_histo.(['obs_hcst_int_area_l', tmp.lyear_str])(histind) = ...
            sum(grid.tarea_60(corr_histo.(['obs_hcst_int_ind_l', tmp.lyear_str]){histind}), 'omitnan');

        %lens2
        corr_histo.(['obs_lens2_ind_l', tmp.lyear_str]){histind} = ...
            find(data.([tmp.varname, '_corr_obs_lens2', '_l', tmp.lyear_str])>corrref-0.1 & ...
            data.([tmp.varname, '_corr_obs_lens2', '_l', tmp.lyear_str])<=corrref);
        corr_histo.(['obs_lens2_area_l', tmp.lyear_str])(histind) = ...
            sum(grid.tarea_60(corr_histo.(['obs_lens2_ind_l', tmp.lyear_str]){histind}), 'omitnan');
        histind=histind+1;
    end
    corr_ref_x=-0.95:0.1:0.95;

%% hcst_int histogram
    fig_h = figure('name',fig_cfg.fig_name,'visible','off');

    bar(corr_ref_x,corr_histo.(['hcst_int_area_l', tmp.lyear_str]), 'linewidth', 2);
    xlabel('R'); ylabel('Area (km^2)');
    set(gca, 'fontsize', 20)
    grid minor
%     ylim([-0.2 1])
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_histogram', filesep, 'hcst_int'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'histogram_hcst_int_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;

%% lens2 histogram
    fig_h = figure('name',fig_cfg.fig_name,'visible','off');
    bar(corr_ref_x,corr_histo.(['lens2_area_l', tmp.lyear_str]), 'linewidth', 2)
    xlabel('R'); ylabel('Area (km^2)');
    set(gca, 'fontsize', 20)
    grid minor
%     ylim([-0.2 1])
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_histogram', filesep, 'lens2'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'histogram_lens2_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
    


%% EOF
    tmp.mode_want=3;
    [EOF.lv, EOF.pc, EOF.var_exp] = ...
        Func_0024_EOF_3d((data.([tmp.varname, '_model_l', tmp.lyear_str])(:,:,1:end-lyear) - data.([tmp.varname, '_lens2_l', tmp.lyear_str])(:,:,1:end-lyear)), tmp.mode_want);


%% EOF- PCT plot
    for EOF_mode=1:tmp.mode_want
        fig_cfg.fig_name=['l',tmp.lyear_str, ', EOF_pct_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'visible','off');
    %     hold on
        plot(cfg.iyears(1:end-lyear)+lyear,EOF.pc(:,EOF_mode), 'k-', 'linewidth', 2);
    %         plot(cfg.iyears(1:end-4)+4,EOF.pc(:,EOF_mode));
        grid minor
        xlabel('year'); ylabel(['PCT-mode', num2str(EOF_mode), ', ', num2str(round(EOF.var_exp(EOF_mode),1)),'% explained']);
        set(gca, 'fontsize', 20)
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_EOF', filesep, 'hcst_int'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'PCT_hcst_l', tmp.lyear_str,'_mode',num2str(EOF_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
    % hold off

%% EOF- LV plot
    for EOF_mode=1:tmp.mode_want
        tmp.C=squeeze(EOF.lv(:,:,EOF_mode));
        tmp.C=tmp.C([end, 1:end],:);
        
        tmp.prc95 =prctile(tmp.C(isfinite(tmp.C)), 95);
        tmp.prc05 =prctile(tmp.C(isfinite(tmp.C)), 5);

        fig_cfg.fig_name=['l',tmp.lyear_str, ', ensmean_, ', tmp. varname];
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
%         tmp.maxval=max(abs(tmp.C(isfinite(tmp.C))));
%         caxis(ax_m, [-tmp.maxval tmp.maxval]); 
%         caxis(ax_m, [tmp.prc05 tmp.prc95]); 
        tmp.maxval=max(abs([tmp.prc05 tmp.prc95]));
        caxis(ax_m, [-tmp.maxval tmp.maxval]); 
        colormap(fig_cfg.c_map);
%         colormap(jet(20));
        cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
        set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    %     title(cb,'R','fontsize',12);
    
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
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_EOF', filesep, 'hcst_int'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'LV_hcst_l', tmp.lyear_str, '_mode',num2str(EOF_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end






end

%% histogram - all (assm)
fig_h = figure('name',fig_cfg.fig_name,'visible','off');

% newcolors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'}; %rainbow
newcolors = {'#F00','#F80','#0B0','#00F','#50F','#A0F'}; %rainbow (except yellow)

colororder(newcolors)

for lyear=0:cfg.proj_year-1
    tmp.lyear_str=num2str(lyear, '%02i');
    hold on 
    plot(corr_ref_x,corr_histo.(['hcst_int_area_l', tmp.lyear_str]), 'linewidth', -lyear+cfg.proj_year);
    tmp.mean_group(lyear+1)=sum(corr_ref_x.*corr_histo.(['hcst_int_area_l', tmp.lyear_str])) / ...
        sum(corr_histo.(['hcst_int_area_l', tmp.lyear_str]));
end
plot(corr_ref_x,corr_histo.lens2_area_l01, 'k--', 'linewidth', 1);
tmp.mean_group_ext = sum(corr_ref_x.*corr_histo.lens2_area_l01) / ...
        sum(corr_histo.lens2_area_l01);
fig_cfg.colororder=colororder;
xline(0)
%avg
% xline(tmp.mean_group(1), 'color', fig_cfg.colororder(1,:), 'linewidth', 4);
% xline(tmp.mean_group(2), 'color', fig_cfg.colororder(2,:), 'linewidth', 3);
% xline(tmp.mean_group(3), 'color', fig_cfg.colororder(3,:), 'linewidth', 2);
% xline(tmp.mean_group(4), 'color', fig_cfg.colororder(4,:), 'linewidth', 1);
% xline(tmp.mean_group_ext, 'k--', 'linewidth', 1);

switch cfg.comp
    case 'ocn'
        ylim([0 12*10^11])
        tmp.yref=1.*10^11;
    case 'ice'
        ylim([0 6.5*10^10])
        tmp.yref=1.*10^10;
    case 'atm'
        ylim([0 14*10^13])
        tmp.yref = 1.*10^13;
    case 'lnd'
        ylim([0 9*10^13])
        tmp.yref = 1.*10^13;
end
line([tmp.mean_group(1), tmp.mean_group(1)], [0, tmp.yref], 'color', fig_cfg.colororder(1,:), 'linewidth', 5);
line([tmp.mean_group(2), tmp.mean_group(2)], [0, tmp.yref], 'color', fig_cfg.colororder(2,:), 'linewidth', 4);
line([tmp.mean_group(3), tmp.mean_group(3)], [0, tmp.yref], 'color', fig_cfg.colororder(3,:), 'linewidth', 3);
line([tmp.mean_group(4), tmp.mean_group(4)], [0, tmp.yref], 'color', fig_cfg.colororder(4,:), 'linewidth', 2);
line([tmp.mean_group(5), tmp.mean_group(5)], [0, tmp.yref], 'color', fig_cfg.colororder(5,:), 'linewidth', 1);

line([tmp.mean_group_ext, tmp.mean_group_ext], [0, tmp.yref], 'color', 'k', 'linewidth', 1, 'LineStyle','--');

hold off
% legend('INT-L0', 'INT-L1', 'INT-L2', 'INT-L3', 'INT-L4', 'EXT', 'Location', 'northwest')
legend('INT-L1', 'INT-L2', 'INT-L3', 'INT-L4', 'INT-L5', 'LENS2', 'Location', 'northwest')

grid minor
xlim([-1 1])

xlabel('R'); ylabel('Area (km^2)');
set(gca, 'fontsize', 20)
dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_histogram', filesep, 'hcst_int'];
if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
cfg.figname=[dirs.figdir, filesep, 'plot_hcst_all_', tmp.varname, '.tif'];
print(fig_h, cfg.figname, '-dpng');
RemoveWhiteSpace([], 'file', cfg.figname);
close all;




% % % % % %% histogram - all (obs)
% % % % % fig_h = figure('name',fig_cfg.fig_name,'visible','off');
% % % % % 
% % % % % % newcolors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'}; %rainbow
% % % % % newcolors = {'#F00','#F80','#0B0','#00F','#50F','#A0F'}; %rainbow (except yellow)
% % % % % 
% % % % % colororder(newcolors)
% % % % % 
% % % % % for lyear=0:cfg.proj_year-1
% % % % %     tmp.lyear_str=num2str(lyear, '%02i');
% % % % %     hold on 
% % % % %     plot(corr_ref_x,corr_histo.(['obs_hcst_int_area_l', tmp.lyear_str]), 'linewidth', -lyear+cfg.proj_year);
% % % % %     tmp.mean_group(lyear+1)=sum(corr_ref_x.*corr_histo.(['obs_hcst_int_area_l', tmp.lyear_str])) / ...
% % % % %         sum(corr_histo.(['obs_hcst_int_area_l', tmp.lyear_str]));
% % % % % end
% % % % % plot(corr_ref_x,corr_histo.obs_lens2_area_l01, 'k--', 'linewidth', 1);
% % % % % tmp.mean_group_ext = sum(corr_ref_x.*corr_histo.obs_lens2_area_l01) / ...
% % % % %         sum(corr_histo.obs_lens2_area_l01);
% % % % % fig_cfg.colororder=colororder;
% % % % % xline(0)
% % % % % %avg
% % % % % % xline(tmp.mean_group(1), 'color', fig_cfg.colororder(1,:), 'linewidth', 4);
% % % % % % xline(tmp.mean_group(2), 'color', fig_cfg.colororder(2,:), 'linewidth', 3);
% % % % % % xline(tmp.mean_group(3), 'color', fig_cfg.colororder(3,:), 'linewidth', 2);
% % % % % % xline(tmp.mean_group(4), 'color', fig_cfg.colororder(4,:), 'linewidth', 1);
% % % % % % xline(tmp.mean_group_ext, 'k--', 'linewidth', 1);
% % % % % 
% % % % % switch cfg.comp
% % % % %     case 'ocn'
% % % % %         ylim([0 10*10^11])
% % % % %         tmp.yref=1.*10^11;
% % % % %     case 'atm'
% % % % %         ylim([0 14*10^13])
% % % % %         tmp.yref = 1.*10^13;
% % % % % end
% % % % % line([tmp.mean_group(1), tmp.mean_group(1)], [0, tmp.yref], 'color', fig_cfg.colororder(1,:), 'linewidth', 5);
% % % % % line([tmp.mean_group(2), tmp.mean_group(2)], [0, tmp.yref], 'color', fig_cfg.colororder(2,:), 'linewidth', 4);
% % % % % line([tmp.mean_group(3), tmp.mean_group(3)], [0, tmp.yref], 'color', fig_cfg.colororder(3,:), 'linewidth', 3);
% % % % % line([tmp.mean_group(4), tmp.mean_group(4)], [0, tmp.yref], 'color', fig_cfg.colororder(4,:), 'linewidth', 2);
% % % % % line([tmp.mean_group(5), tmp.mean_group(5)], [0, tmp.yref], 'color', fig_cfg.colororder(5,:), 'linewidth', 1);
% % % % % 
% % % % % line([tmp.mean_group_ext, tmp.mean_group_ext], [0, tmp.yref], 'color', 'k', 'linewidth', 1, 'LineStyle','--');
% % % % % 
% % % % % hold off
% % % % % % legend('INT-L0', 'INT-L1', 'INT-L2', 'INT-L3', 'INT-L4', 'EXT', 'Location', 'northwest')
% % % % % legend('INT-L1', 'INT-L2', 'INT-L3', 'INT-L4', 'INT-L5', 'LENS2', 'Location', 'northwest')
% % % % % 
% % % % % grid minor
% % % % % xlim([-1 1])
% % % % % 
% % % % % xlabel('R'); ylabel('Area (km^2)');
% % % % % set(gca, 'fontsize', 20)
% % % % % dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_histogram', filesep, 'obs_hcst_int'];
% % % % % if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
% % % % % cfg.figname=[dirs.figdir, filesep, 'plot_obs_hcst_all_', tmp.varname, '.tif'];
% % % % % print(fig_h, cfg.figname, '-dpng');
% % % % % RemoveWhiteSpace([], 'file', cfg.figname);
% % % % % close all;



% %% timeseries
% sta_lonlat = {[345, 75], [300, 85], [260, 85], [220, 80], [180, 80], ...
%                 [140, 85], [100, 85], [40, 75], [360, 65], [310, 55], ...
%                 [230, 50], [200, 50], [180, 50], [150, 55], [180, 35], ...
%                 [140, 33], [124, 35], [240,20], [180, 20], [150, 15], ...
%                 [180, 10], [260, 0], [160, 0], [280, -20], [240, -20], ...
%                 [200, -20], [160, -20], [160, -30], [200, -40], [200, -50], ...
%                 [110, -10], [90, 15], [60, 10], [42, -10], [40, -30], ...
%                 [40, -40], [40, -50], [340, 40], [320, 30], [300, 30], ...
%                 [285, 30], [340, 20], [300, 20], [340, 10], [320, 10], ...
%                 [270, 25], [360, 0], [320, 0], [10, -10], [330, -10], ...
%                 [330, -30], [5, 40], [60, -65], [180, -70], [240, -70, ...
%                 [320, -70]};
% 345, 75 (iceland)
% 179, 80 (north pole)
% 160, -30 (Austrailia eastern offshore

% xpoint=160;
% ypoint = -30;
% [grid.id_w, grid.id_e, grid.id_s, grid.id_n] = Func_0012_findind_Y(0.5, [xpoint, ypoint], ...
%             grid.tlong, ...
%             grid.tlat, 'CESM2'); % find valid lon, lat index near station
% 
% 
% plot(cfg.iyears+4,squeeze(data.([tmp.varname, '_model_l04'])(grid.id_w,grid.id_s,:)), 'linewidth', 2)
% hold on
% plot(cfg.iyears+4,squeeze(data.([tmp.varname, '_assm'])(grid.id_w,grid.id_s,:)), 'linewidth', 2)
% plot(cfg.iyears+4,squeeze(data.([tmp.varname, '_lens2_l04'])(grid.id_w,grid.id_s,:)), 'linewidth', 2)
% hold off
% legend ('hcst', 'assm', 'lens2')
% title(['l04, ', tmp.varname, ', ', num2str(xpoint),'E, ', num2str(ypoint), 'N'])
% grid minor
% disp('end')

% plot(squeeze(data.diatChl_model_l04(grid.id_w,grid.id_s,:)) - squeeze(data.diatChl_lens2_l04(grid.id_w,grid.id_s,:)))
% hold on
% plot(squeeze(data.diatChl_assm(grid.id_w,grid.id_s,:))- squeeze(data.diatChl_lens2_l04(grid.id_w,grid.id_s,:)))
% hold off
% legend ('hcst-lens2', 'assm-lens2')
% title('l04, diatChl, 179E, 80N')


toc;
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
            obsname_simple='GPCP';
        case 'PSL'
            obsname_simple='ERA5';
        case 'SSH'
            obsname_simple='CMEMS';
        case 'TS'
            obsname_simple='HadCRUT5';
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