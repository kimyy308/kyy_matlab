% %  Created 22-Sep-2023 by Yong-Yub Kim
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
cfg.vars = {'TS', 'SST', 'PRECT', 'PSL'};
% cfg.vars = {'TS'};
% cfg.vlayer=1:10; % 10layer. don't put more than 15
cfg.vlayer=1; % surface, vertical slice

cfg.vars = {'TS', 'SST', 'PRECT', 'PSL'};
cfg.vars={'TLAI','SOILWATER_10CM', 'COL_FIRE_CLOSS', 'FAREA_BURNED'};


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
vstr=['v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i')];
dirs.figroot=[tmp.kimyypath, '/Figure/CESM2/ESP/HCST_EXP/archive_anomaly/', cfg.comp,'/', cfg.var, filesep, vstr];

mkdir(dirs.figroot)

cfg.iyears=cfg.obs_iyears;
cfg.gnm='f09_g17';
cfg.proj_year=5;
cfg.season = {'JFM1', 'AMJ1', 'JAS1', 'OND1'};

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
fig_flags(1:100)=1;
% fig_flags(48)=1; %model drift
% fig_flags(24)=1; %assm obs nRMSE
% fig_flags(24:26)=1; %assm obs RMSE
% fig_flags(41:43)=1; %assm mean bias
% fig_flags(21:30)=1; % RMSEs
% fig_flags(49)=1; %model drift ratio
% fig_flags(50:51)=1; %obs-AR1 coef, noise
% fig_flags(31:34)=1;
% fig_flags(53:57)=1;
% fig_flags(6)=1;
% fig_flags(52)=1;

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
% fig_flags(58:59)=1;

S = shaperead('landareas.shp');

%% read & plot data
tmp.varname=cfg.var;

%% read all data for RMSE colorbar
fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_JFM1', '.mat'];
load(fig_cfg.mat_name, 'data', 'data2')
data2_JFM1=data2;
fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_AMJ1', '.mat'];
load(fig_cfg.mat_name, 'data', 'data2')
data2_AMJ1=data2;
fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_JAS1', '.mat'];
load(fig_cfg.mat_name, 'data', 'data2')
data2_JAS1=data2;
fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_OND1', '.mat'];
load(fig_cfg.mat_name, 'data', 'data2')
data2_OND1=data2;
% % % fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
% % %             '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
% % %             '_l04', 'y.mat'];
% % % load(fig_cfg.mat_name, 'data', 'data2')
% % % data2_l4=data2;

clear tmp.ydata tmp.ydata_lens2 tmp.ydata_obs tmp.ydata_assm
% for lyear=0:cfg.proj_year-1
%     tmp.lyear_str=num2str(lyear, '%02i');
%     fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
%             '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
%             '_', tmp.season, 'y.mat'];
%     load(fig_cfg.mat_name, 'data', 'data2')
%     pcolor(data2.([tmp.varname,'_corr_assm_int_', tmp.season])'); shading flat; colorbar;
for lss=1:length(cfg.season)
    tmp.season=cfg.season{lss};
    tmp.mons = f_season_mons(tmp.season);
    fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_', tmp.season, '.mat'];
    load(fig_cfg.mat_name, 'data', 'data2');

%% clim time range
    switch cfg.obs_name
       case 'GPCC'
           cfg.clim_ys=1961;
           cfg.clim_ye=2019;
       otherwise
           cfg.clim_ys=1961;
           cfg.clim_ye=2020;
    end
   lyear=0;
   cfg.clim_tlen = (cfg.clim_ys-1959)-lyear:(cfg.clim_ye-2020)+cfg.len_t_y-lyear;
   cfg.clim_tlen2=length(cfg.clim_tlen);


%% AR1_assm corr map --------------------------------------
if fig_flags(1)==1
for fake=1:1
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
    data2.([tmp.varname, '_corr_assm_AR1_ano', '_', tmp.season])(data2.([tmp.varname, '_corr_assm_AR1_ano', '_', tmp.season])==0)=NaN;
    tmp.C=data2.([tmp.varname, '_corr_assm_AR1_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 

    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_assm_AR1_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_assm_AR1_ano_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_AR1', filesep, 'assm'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_assm_AR1_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% AR1_obs corr map --------------------------------------
if fig_flags(2)==1
for fake=1:1
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
    data2.([tmp.varname, '_corr_obs_AR1', '_', tmp.season])(data2.([tmp.varname, '_corr_obs_AR1', '_', tmp.season])==0)=NaN;
    tmp.C=data2.([tmp.varname, '_corr_obs_AR1', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_obs_AR1', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_obs_AR1_ano_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_AR1', filesep, 'obs'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_AR1_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% obs & model corr map --------------------------------------
if fig_flags(3)==1
for fake=1:1
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
    tmp.C=data2.([tmp.varname, '_corr_obs_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_obs_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_obs_ano_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% obs & assm corr map --------------------------------------
if fig_flags(4)==1
for fake=1:1
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
    tmp.C=data2.([tmp.varname, '_corr_obs_assm_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_obs_assm_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_obs_assm_ano_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_assm_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% obs & lens2 corr map --------------------------------------
if fig_flags(5)==1
for fake=1:1
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
    tmp.C=data2.([tmp.varname, '_corr_obs_lens2_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_obs_lens2_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_obs_lens2_ano_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_lens2_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% obs & hcst-lens2, obs-hcst2 & hcst-lens2 corr map --------------------------------------
if fig_flags(6)==1
for fake=1:1
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
    tmp.C=data2.([tmp.varname, '_corr_obs_int_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_obs_int_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_obs_int_ano_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_int_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_int_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;

%% 2
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
    tmp.C=data2.([tmp.varname, '_corr_obs_int2_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_obs_int2_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_obs_int_ano_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_int2_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_int2_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
end
end

%% obs-det & model-det corr map --------------------------------------
if fig_flags(7)==1
for fake=1:1
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
    tmp.C=data2.([tmp.varname, '_corr_obs_det_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_obs_det_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_obs_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_det_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_det_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% obs-det & assm-det corr map --------------------------------------
if fig_flags(8)==1
for fake=1:1
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
    tmp.C=data2.([tmp.varname, '_corr_obs_assm_det_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_obs_assm_det_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_obs_assm_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_det_map', filesep, 'assm'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_assm_det_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% obs-det & lens2-det corr map --------------------------------------
if fig_flags(9)==1
for fake=1:1
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
    tmp.C=data2.([tmp.varname, '_corr_obs_lens2_det_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_obs_lens2_det_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_obs_lens2_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_det_map', filesep, 'lens2'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_lens2_det_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% model & obs (>AR1) corr map --------------------------------------
if fig_flags(10)==1
for fake=1:1
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
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10_botgray', tmp.dropboxpath);

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data2.([tmp.varname, '_corr_obs_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);

    tmp.C_H=data2.([tmp.varname, '_corr_obs_AR1', '_', tmp.season]);
    tmp.C_H=tmp.C_H([end, 1:end],:);

    tmp.C(tmp.C<tmp.C_H)=-1;

%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_obs_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_obs_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_ano_overAR1_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
end
end


%% model & assm corr map --------------------------------------
if fig_flags(11)==1
for fake=1:1
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
    tmp.C=data2.([tmp.varname, '_corr_assm_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_assm_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_assm_ano_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    cfg.figname=[dirs.figdir, filesep, 'corr_assm_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
end
end

%% model-det & assm-det corr map --------------------------------------
if fig_flags(12)==1
for fake=1:1
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
    tmp.C=data2.([tmp.varname, '_corr_assm_det_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_assm_det_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_assm_det_ano, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_assm_det_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_assm_det_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
end
end

 %% model-lens2 & assm corr, model-lens2 & assm-lens2 map --------------------------------------
if fig_flags(13)==1
for fake=1:1

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
    tmp.C=data2.([tmp.varname, '_corr_assm_int_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_assm_int_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_assm_int_ano_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    cfg.figname=[dirs.figdir, filesep, 'corr_assm_int_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];

    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;


%% 2
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
    tmp.C=data2.([tmp.varname, '_corr_assm_int2_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_assm_int2_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_assm_int2_ano_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_assm_int2_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_assm_int2_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];

    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
 end
end

%% lens2 & assm corr map --------------------------------------
if fig_flags(14)==1
for fake=1:1

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
    tmp.C=data2.([tmp.varname, '_corr_assm_lens2_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_assm_lens2_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_assm_lens2_ano_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    cfg.figname=[dirs.figdir, filesep, 'corr_assm_ext_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
end
end

%% lens2-det & assm-det corr map --------------------------------------
if fig_flags(15)==1
for fake=1:1

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
    tmp.C=data2.([tmp.varname, '_corr_assm_lens2_det_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);

    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_assm_lens2_det_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_assm_lens2_det_ano_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_assm_ext_det_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_assm_ext_det_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
end
end

%% model & assm (>AR1) corr map --------------------------------------
if fig_flags(16)==1
for fake=1:1
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
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10_botgray', tmp.dropboxpath);

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data2.([tmp.varname, '_corr_assm_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);

    tmp.C_H=data2.([tmp.varname, '_corr_assm_AR1_ano', '_', tmp.season]);
    tmp.C_H=tmp.C_H([end, 1:end],:);

    tmp.C(tmp.C<tmp.C_H)=-1;

%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_assm_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_assm_ano_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    fig_cfg.c_map(1,:)=[0.8 0.8 0.8]; % 1< : grayed
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
    cfg.figname=[dirs.figdir, filesep, 'corr_assm_ano_overAR1_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
end
end

%% obs trend map --------------------------------------
if fig_flags(17)==1
for fake=1:1
    data2.([tmp.varname, '_obs_trend2', '_', tmp.season])(data2.([tmp.varname, '_obs_trend2', '_', tmp.season])==0)=NaN;
    tmp.trend_all= [ data2.([tmp.varname, '_obs_trend2', '_', tmp.season]), ...
                 data2.([tmp.varname, '_assm_trend2', '_', tmp.season]), ...
                 data2.([tmp.varname, '_model_trend2', '_', tmp.season]), ...
                 data2.([tmp.varname, '_lens2_trend2', '_', tmp.season])];
    tmp.trend_prc99 =prctile(tmp.trend_all(:), 99);
    tmp.trend_prc01 =prctile(tmp.trend_all(:), 1);
    tmp.trend_max=max([abs(tmp.trend_prc99) abs(tmp.trend_prc01)]);

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
    tmp.C=data2.([tmp.varname, '_obs_trend2', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_obs_trend2', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', trend_obs_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, [-tmp.trend_max tmp.trend_max]); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,[data.units, '/y'],'fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_trend_map', filesep, 'obs'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'trend_obs_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% assm trend map --------------------------------------
if fig_flags(18)==1
for fake=1:1
    data2.([tmp.varname, '_obs_trend2', '_', tmp.season])(data2.([tmp.varname, '_obs_trend2', '_', tmp.season])==0)=NaN;
    tmp.trend_all= [ data2.([tmp.varname, '_obs_trend2', '_', tmp.season]), ...
                 data2.([tmp.varname, '_assm_trend2', '_', tmp.season]), ...
                 data2.([tmp.varname, '_model_trend2', '_', tmp.season]), ...
                 data2.([tmp.varname, '_lens2_trend2', '_', tmp.season])];
    tmp.trend_prc99 =prctile(tmp.trend_all(:), 99);
    tmp.trend_prc01 =prctile(tmp.trend_all(:), 1);
    tmp.trend_max=max([abs(tmp.trend_prc99) abs(tmp.trend_prc01)]);


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
    tmp.C=data2.([tmp.varname, '_assm_trend2', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_assm_trend2', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', trend_assm_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, [-tmp.trend_max tmp.trend_max]); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,[data.units, '/y'],'fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_trend_map', filesep, 'assm'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'trend_assm_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% lens2 trend map --------------------------------------
if fig_flags(19)==1
for fake=1:1

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
    tmp.C=data2.([tmp.varname, '_lens2_trend2', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_lens2_trend2', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', trend_lens2_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, [-tmp.trend_max tmp.trend_max]); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,[data.units, '/y'],'fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_trend_map', filesep, 'lens2'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'trend_lens2_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% hcst trend map --------------------------------------
if fig_flags(20)==1
for fake=1:1

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
    tmp.C=data2.([tmp.varname, '_model_trend2', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_model_trend2', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', trend_model_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, [-tmp.trend_max tmp.trend_max]); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,[data.units, '/y'],'fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_trend_map', filesep, 'hcst'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'trend_hcst_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% obs & assm rmse map --------------------------------------
if fig_flags(21)==1
for fake=1:1
    tmp.rmse_obs_all= [ data2_JFM1.([tmp.varname, '_assm_rmse2_obs', '_JFM1']), ...
                 data2_JFM1.([tmp.varname, '_model_rmse2_obs', '_JFM1']), ...
                 data2_JFM1.([tmp.varname, '_lens2_rmse2_obs', '_JFM1']), ...
                 data2_AMJ1.([tmp.varname, '_assm_rmse2_obs', '_AMJ1']), ...
                 data2_AMJ1.([tmp.varname, '_model_rmse2_obs', '_AMJ1']), ...
                 data2_AMJ1.([tmp.varname, '_lens2_rmse2_obs', '_AMJ1']), ...
                 data2_JAS1.([tmp.varname, '_assm_rmse2_obs', '_JAS1']), ...
                 data2_JAS1.([tmp.varname, '_model_rmse2_obs', '_JAS1']), ...
                 data2_JAS1.([tmp.varname, '_lens2_rmse2_obs', '_JAS1']), ...
                 data2_OND1.([tmp.varname, '_assm_rmse2_obs', '_OND1']), ...
                 data2_OND1.([tmp.varname, '_model_rmse2_obs', '_OND1']), ...
                 data2_OND1.([tmp.varname, '_lens2_rmse2_obs', '_OND1'])];
    tmp.rmse_obs_prc99 =prctile(tmp.rmse_obs_all(:), 99);
    tmp.rmse_obs_prc01 =prctile(tmp.rmse_obs_all(:), 1);

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
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr_10', tmp.dropboxpath);
    fig_cfg.c_map = flip(fig_cfg.c_map);

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data2.([tmp.varname, '_assm_rmse2_obs', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_assm_rmse2_obs', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', assm_rmse2_obs_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, [tmp.rmse_obs_prc01 tmp.rmse_obs_prc99]); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,data.units,'fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_rmse_obs_map', filesep, 'assm'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'assm_rmse_obs_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% obs & lens2 rmse map --------------------------------------
if fig_flags(22)==1
for fake=1:1

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
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr_10', tmp.dropboxpath);
    fig_cfg.c_map = flip(fig_cfg.c_map);

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data2.([tmp.varname, '_lens2_rmse2_obs', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_lens2_rmse2_obs', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', lens2_rmse2_obs_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, [tmp.rmse_obs_prc01 tmp.rmse_obs_prc99]); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,data.units,'fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_rmse_obs_map', filesep, 'lens2'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'lens2_rmse_obs_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% obs & model rmse map --------------------------------------
if fig_flags(23)==1
for fake=1:1

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
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr_10', tmp.dropboxpath);
    fig_cfg.c_map = flip(fig_cfg.c_map);

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data2.([tmp.varname, '_model_rmse2_obs', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_model_rmse2_obs', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', model_rmse2_obs_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, [tmp.rmse_obs_prc01 tmp.rmse_obs_prc99]); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,data.units,'fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_rmse_obs_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'model_rmse_obs_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% obs & assm nrmse map --------------------------------------
if fig_flags(24)==1
for fake=1:1

    tmp.nrmse_obs_all= [ data2_JFM1.([tmp.varname, '_assm_nrmse2_obs', '_JFM1']), ...
                 data2_JFM1.([tmp.varname, '_model_nrmse2_obs', '_JFM1']), ...
                 data2_JFM1.([tmp.varname, '_lens2_nrmse2_obs', '_JFM1']), ...
                 data2_AMJ1.([tmp.varname, '_assm_nrmse2_obs', '_AMJ1']), ...
                 data2_AMJ1.([tmp.varname, '_model_nrmse2_obs', '_AMJ1']), ...
                 data2_AMJ1.([tmp.varname, '_lens2_nrmse2_obs', '_AMJ1']), ...
                 data2_JAS1.([tmp.varname, '_assm_nrmse2_obs', '_JAS1']), ...
                 data2_JAS1.([tmp.varname, '_model_nrmse2_obs', '_JAS1']), ...
                 data2_JAS1.([tmp.varname, '_lens2_nrmse2_obs', '_JAS1']), ...
                 data2_OND1.([tmp.varname, '_assm_nrmse2_obs', '_OND1']), ...
                 data2_OND1.([tmp.varname, '_model_nrmse2_obs', '_OND1']), ...
                 data2_OND1.([tmp.varname, '_lens2_nrmse2_obs', '_OND1'])];



    tmp.nrmse_obs_all(isinf(tmp.nrmse_obs_all))=NaN;
    tmp.nrmse_obs_prc90 =prctile(tmp.nrmse_obs_all(:), 90);
    tmp.nrmse_obs_prc01 =prctile(tmp.nrmse_obs_all(:), 1);

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
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr_10', tmp.dropboxpath);
    fig_cfg.c_map = flip(fig_cfg.c_map);

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data2.([tmp.varname, '_assm_nrmse2_obs', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
    tmp.C(isinf(tmp.C))=NaN;
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_assm_nrmse2_obs', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', assm_nrmse2_obs_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, [tmp.nrmse_obs_prc01 tmp.nrmse_obs_prc90]); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,' ','fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_nrmse_obs_map', filesep, 'assm'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'assm_nrmse_obs_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% obs & lens2 nrmse map --------------------------------------
if fig_flags(25)==1
for fake=1:1

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
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr_10', tmp.dropboxpath);
    fig_cfg.c_map = flip(fig_cfg.c_map);
    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data2.([tmp.varname, '_lens2_nrmse2_obs', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
    tmp.C(isinf(tmp.C))=NaN;
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_lens2_nrmse2_obs', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', lens2_rmse2_obs_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, [tmp.nrmse_obs_prc01 tmp.nrmse_obs_prc90]); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,' ','fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_nrmse_obs_map', filesep, 'lens2'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'lens2_nrmse_obs_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% obs & model nrmse map --------------------------------------
if fig_flags(26)==1
for fake=1:1

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
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr_10', tmp.dropboxpath);
    fig_cfg.c_map = flip(fig_cfg.c_map);

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data2.([tmp.varname, '_model_nrmse2_obs', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
    tmp.C(isinf(tmp.C))=NaN;
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_model_nrmse2_obs', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', model_nrmse2_obs_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, [tmp.nrmse_obs_prc01 tmp.nrmse_obs_prc90]); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,' ','fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_nrmse_obs_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'model_nrmse_obs_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% assm & lens2 rmse map --------------------------------------
if fig_flags(27)==1
for fake=1:1
    tmp.rmse_all= [ data2_JFM1.([tmp.varname, '_model_rmse2', '_JFM1']), ...
                 data2_JFM1.([tmp.varname, '_lens2_rmse2', '_JFM1']), ...
                 data2_AMJ1.([tmp.varname, '_model_rmse2', '_AMJ1']), ...
                 data2_AMJ1.([tmp.varname, '_lens2_rmse2', '_AMJ1']), ...
                 data2_JAS1.([tmp.varname, '_model_rmse2', '_JAS1']), ...
                 data2_JAS1.([tmp.varname, '_lens2_rmse2', '_JAS1']), ...
                 data2_OND1.([tmp.varname, '_model_rmse2', '_OND1']), ...
                 data2_OND1.([tmp.varname, '_lens2_rmse2', '_OND1'])];
    tmp.rmse_prc99 =prctile(tmp.rmse_all(:), 99);
    tmp.rmse_prc01 =prctile(tmp.rmse_all(:), 1);

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
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr_10', tmp.dropboxpath);
    fig_cfg.c_map = flip(fig_cfg.c_map);

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data2.([tmp.varname, '_lens2_rmse2', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_lens2_rmse2', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', lens2_rmse2_assm_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, [tmp.rmse_prc01 tmp.rmse_prc99]); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,data.units,'fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_rmse_assm_map', filesep, 'lens2'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'lens2_rmse_assm_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% assm & model rmse map --------------------------------------
if fig_flags(28)==1
for fake=1:1

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
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr_10', tmp.dropboxpath);
    fig_cfg.c_map = flip(fig_cfg.c_map);

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data2.([tmp.varname, '_model_rmse2', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_model_rmse2', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', hcst_rmse2_assm_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, [tmp.rmse_prc01 tmp.rmse_prc99]); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,data.units,'fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_rmse_assm_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'model_rmse2_assm_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% assm & lens2 nrmse map --------------------------------------
if fig_flags(29)==1
for fake=1:1
    tmp.nrmse_all= [ data2_JFM1.([tmp.varname, '_model_nrmse2', '_JFM1']), ...
                 data2_JFM1.([tmp.varname, '_lens2_nrmse2', '_JFM1']), ...
                 data2_AMJ1.([tmp.varname, '_model_nrmse2', '_AMJ1']), ...
                 data2_AMJ1.([tmp.varname, '_lens2_nrmse2', '_AMJ1']), ...
                 data2_JAS1.([tmp.varname, '_model_nrmse2', '_JAS1']), ...
                 data2_JAS1.([tmp.varname, '_lens2_nrmse2', '_JAS1']), ...
                 data2_OND1.([tmp.varname, '_model_nrmse2', '_OND1']), ...
                 data2_OND1.([tmp.varname, '_lens2_nrmse2', '_OND1'])];
    tmp.nrmse_prc99 =prctile(tmp.nrmse_all(:), 99);
    tmp.nrmse_prc01 =prctile(tmp.nrmse_all(:), 1);

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
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr_10', tmp.dropboxpath);
    fig_cfg.c_map = flip(fig_cfg.c_map);

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data2.([tmp.varname, '_lens2_nrmse2', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_lens2_nrmse2', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', lens2_nrmse_assm_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, [tmp.nrmse_prc01 tmp.nrmse_prc99]); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,' ','fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_nrmse_assm_map', filesep, 'lens2'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'lens2_nrmse2_assm_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% assm & model nrmse map --------------------------------------
if fig_flags(30)==1
for fake=1:1

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
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr_10', tmp.dropboxpath);
    fig_cfg.c_map = flip(fig_cfg.c_map);

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data2.([tmp.varname, '_model_nrmse2', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_model_nrmse2', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', model_nrmse2_assm_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, [tmp.nrmse_prc01 tmp.nrmse_prc99]); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,' ','fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_nrmse_assm_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'model_nrmse2_assm_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% temporal std of obs
if fig_flags(38)==1
for fake=1:1

    data2.obs_stdt=mean(data2.([tmp.varname,'_obs_ano_det2_stdt']),3, 'omitnan');    
    data2.assm_stdt=mean(data2.([tmp.varname,'_assm_ano_det2_stdt']),3, 'omitnan');

    tmp.stdt_all= [  data2.assm_stdt, ...
                 data2.obs_stdt];
    
    tmp.stdt_prc99 =prctile(tmp.stdt_all(:), 99);
    tmp.stdt_prc01 =prctile(tmp.stdt_all(:), 1);
    
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
    tmp.C=data2.obs_stdt;
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    fig_cfg.fig_name=[tmp.season, ', ensstdt_, ', tmp. varname];
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
    caxis(ax_m, [tmp.stdt_prc01 tmp.stdt_prc99]); 
%     colormap(fig_cfg.c_map);
    colormap(flip(autumn(10)));
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,'\sigma_t','fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_ensstdt_map', filesep, 'obs'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'ensstdt_obs_ano_det2_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
end
end

%% temporal std of assm
if fig_flags(39)==1
for fake=1:1
    data2.obs_stdt=mean(data2.([tmp.varname,'_obs_ano_det2_stdt']),3, 'omitnan');    
    data2.assm_stdt=mean(data2.([tmp.varname,'_assm_ano_det2_stdt']),3, 'omitnan');

    tmp.stdt_all= [  data2.assm_stdt, ...
                 data2.obs_stdt];
    
    tmp.stdt_prc99 =prctile(tmp.stdt_all(:), 99);
    tmp.stdt_prc01 =prctile(tmp.stdt_all(:), 1);

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
    tmp.C=data2.([tmp.varname,'_assm_ano_det2_stdt']);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    fig_cfg.fig_name=[tmp.season, ', ensstdt_, ', tmp. varname];
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
    caxis(ax_m, [tmp.stdt_prc01 tmp.stdt_prc99]); 
%     colormap(fig_cfg.c_map);
    colormap(flip(autumn(10)));
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,'\sigma_t','fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_ensstdt_map', filesep, 'assm'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'ensstdt_assm_ano_det2_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
end
end

%% temporal std difference of assm (ASSM - OBS)
if fig_flags(40)==1
for fake=1:1
    
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

    tmp.C=data2.([tmp.varname,'_assm_ano_det2_stdt']) ...
        - data2.([tmp.varname,'_obs_ano_det2_stdt']);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    fig_cfg.fig_name=[tmp.season, ', diff_ensstdt_, ', tmp.varname];
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
    caxis(ax_m, [-max(abs(tmp.C(:))) max(abs(tmp.C(:)))]); 
    colormap(fig_cfg.c_map);
%     colormap(flip(autumn(10)));
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,'\sigma_t','fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_ensstdt_map', filesep, 'assm'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'ensstdt_ano_det2_err_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
end
end





%% AR1 coefficient based on assm
if fig_flags(46)==1
for fake=1:1
    if lyear==0
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
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr_10', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=-data2.([tmp.varname, '_assm_AR1_coef']);
        tmp.C=tmp.C([end, 1:end],:);
    %         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
    %         tmp.C_H=tmp.C_H([end, 1:end],:);
    %         tmp.C_2=tmp.C;
    %         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 
    
    
        fig_cfg.fig_name=[tmp.season, ', AR1_coef_, ', tmp. varname];
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
        colormap(fig_cfg.c_map);
        cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
        set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
        title(cb,' ','fontsize',12);
    
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
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_AR1_param', filesep, 'assm'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'AR1_coef_assm_map_', tmp.varname, '_', tmp.season, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
end
end

%% AR1 noise based on assm
if fig_flags(47)==1
for fake=1:1
    if lyear==0
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
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr_10', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=data2.([tmp.varname, '_assm_AR1_noise']);
        tmp.C=tmp.C([end, 1:end],:);
    %         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
    %         tmp.C_H=tmp.C_H([end, 1:end],:);
    %         tmp.C_2=tmp.C;
    %         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 
    
    
        fig_cfg.fig_name=[tmp.season, ', AR1_noise_, ', tmp. varname];
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
        colormap(fig_cfg.c_map);
        cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
        set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
        title(cb,data.units,'fontsize',12);
    
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
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_AR1_param', filesep, 'assm'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'AR1_noise_assm_map_', tmp.varname, '_', tmp.season, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
end
end

%% AR1 coefficient based on obs
if fig_flags(50)==1
for fake=1:1
    if lyear==0
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
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr_10', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=-data2.([tmp.varname, '_obs_AR1_coef']);
        tmp.C=tmp.C([end, 1:end],:);
    %         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
    %         tmp.C_H=tmp.C_H([end, 1:end],:);
    %         tmp.C_2=tmp.C;
    %         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 
    
    
        fig_cfg.fig_name=[tmp.season, ', obs_AR1_coef_, ', tmp. varname];
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
        colormap(fig_cfg.c_map);
        cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
        set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
        title(cb,' ','fontsize',12);
    
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
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_AR1_param', filesep, 'obs'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'AR1_coef_obs_map_', tmp.varname, '_', tmp.season, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
end
end

%% AR1 noise based on obs
if fig_flags(51)==1
for fake=1:1
    if lyear==0
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
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr_10', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=data2.([tmp.varname, '_obs_AR1_noise']);
        tmp.C=tmp.C([end, 1:end],:);
    %         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
    %         tmp.C_H=tmp.C_H([end, 1:end],:);
    %         tmp.C_2=tmp.C;
    %         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 
    
    
        fig_cfg.fig_name=[tmp.season, ', AR1_noise_, ', tmp. varname];
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
        colormap(fig_cfg.c_map);
        cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
        set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
        title(cb,data.units,'fontsize',12);
    
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
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_AR1_param', filesep, 'obs'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'AR1_noise_obs_map_', tmp.varname, '_', tmp.season, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
end
end


 %% model-svd_1st(model-lens2) & assm corr map --------------------------------------
if fig_flags(52)==1

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
    tmp.C=data2.([tmp.varname, '_corr_assm_int_svd', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_assm_int_svd', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_assm_int_svd_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    cfg.figname=[dirs.figdir, filesep, 'corr_assm_int_svd_map_', tmp.varname, '_', tmp.season, 'y.tif'];

    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;

%% 2
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
    tmp.C=data2.([tmp.varname, '_corr_assm_int2_svd', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_assm_int2_svd', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_assm_int2_svd_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_assm_int2_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_assm_int2_svd_map_', tmp.varname, '_', tmp.season, 'y.tif'];

    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
end

 %% model-svd_1st(model-lens2) & obs corr map --------------------------------------
if fig_flags(52)==1

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
    tmp.C=data2.([tmp.varname, '_corr_obs_int_svd', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_obs_int_svd', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_assm_int_svd_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_int_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_int_svd_map_', tmp.varname, '_', tmp.season, 'y.tif'];

    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;

%% 2
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
    tmp.C=data2.([tmp.varname, '_corr_obs_int2_svd', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_obs_int2_svd', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_assm_int2_svd_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_int2_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_int2_svd_map_', tmp.varname, '_', tmp.season, 'y.tif'];

    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
end

%% MSSS (assm base) map --------------------------------------
if fig_flags(53)==1
for fake=1:1
    
    data2.([tmp.varname, '_MSSS_assm_ano', '_', tmp.season]) = ...
        1 -  data2.([tmp.varname, '_model_rmse2', '_', tmp.season]).^2 ...
             ./ data2.([tmp.varname, '_lens2_rmse2', '_', tmp.season]).^2;
    data2.([tmp.varname, '_MSSS_obs_ano', '_', tmp.season]) = ...
        1 -  data2.([tmp.varname, '_model_rmse2_obs', '_', tmp.season]).^2 ...
             ./ data2.([tmp.varname, '_lens2_rmse2_obs', '_', tmp.season]).^2;

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
    tmp.C=data2.([tmp.varname, '_MSSS_assm_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);

%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_MSSS_assm_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', _MSSS_assm_ano_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, [-1 1]); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,' ','fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_MSSS_map'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'MSSS_assm_base_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% MSSS (obs base) map --------------------------------------
if fig_flags(54)==1
for fake=1:1
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
    tmp.C=data2.([tmp.varname, '_MSSS_obs_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);

%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_MSSS_obs_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', _MSSS_obs_ano_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    caxis(ax_m, [-1 1]); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    title(cb,' ','fontsize',12);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_MSSS_map'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'MSSS_obs_base_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end


%% SVD - lens2-model plot
if fig_flags(55)==1
    tmp.mode_want=3;
%% PCT (left) plot
    for SVD_mode=1:tmp.mode_want
        fig_cfg.fig_name=[tmp.season, ', SVD_lens2-model_left_pct_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'visible','off');
    %     hold on
        plot(cfg.iyears(cfg.clim_tlen)+lyear,data2.([tmp.varname,'_svd_lens2_model_pcs_left_', tmp.season])(SVD_mode,:), 'k-', 'linewidth', 2);
    %         plot(cfg.iyears(1:end-4)+4,EOF.pc(:,EOF_mode));
        grid minor
        xlabel('year'); ylabel(['PCT-mode', num2str(SVD_mode), ', SCF:', num2str(round(data2.([tmp.varname,'_svd_lens2_model_scf_', tmp.season])(SVD_mode).*100,1)),'%']);
        set(gca, 'fontsize', 20)
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_hcst'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'lens2-model_left_PCT_', tmp.season,'_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
    % hold off

%% LV (left) plot
    for SVD_mode=1:tmp.mode_want
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=squeeze(data2.([tmp.varname,'_svd_lens2_model_lvmap_left_', tmp.season])(:,:,SVD_mode));
        tmp.C=tmp.C([end, 1:end],:);
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
        
        tmp.prc95 =prctile(tmp.C(isfinite(tmp.C)), 95);
        tmp.prc05 =prctile(tmp.C(isfinite(tmp.C)), 5);

        fig_cfg.fig_name=[tmp.season, ', SVD_LV_, ', tmp. varname];
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
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_hcst'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'lens2-model_left_LV_', tmp.season, '_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end


%% PCT (right) plot
    for SVD_mode=1:tmp.mode_want
        fig_cfg.fig_name=[tmp.season, ', SVD_lens2-model_right_pct_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'visible','off');
    %     hold on
        plot(cfg.iyears(cfg.clim_tlen)+lyear,data2.([tmp.varname,'_svd_lens2_model_pcs_right_', tmp.season])(SVD_mode,:), 'k-', 'linewidth', 2);
    %         plot(cfg.iyears(1:end-4)+4,EOF.pc(:,SVD_mode));
        grid minor
        xlabel('year'); ylabel(['PCT-mode', num2str(SVD_mode), ', SCF:', num2str(round(data2.([tmp.varname,'_svd_lens2_model_scf_', tmp.season])(SVD_mode).*100,1)),'%']);
        set(gca, 'fontsize', 20)
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_hcst'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'lens2-model_right_PCT_', tmp.season,'_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
    % hold off

%% LV (right) plot
    for SVD_mode=1:tmp.mode_want
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=squeeze(data2.([tmp.varname,'_svd_lens2_model_lvmap_right_', tmp.season])(:,:,SVD_mode));
        tmp.C=tmp.C([end, 1:end],:);
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
        
        tmp.prc95 =prctile(tmp.C(isfinite(tmp.C)), 95);
        tmp.prc05 =prctile(tmp.C(isfinite(tmp.C)), 5);

        fig_cfg.fig_name=[tmp.season, ', SVD_LV_, ', tmp. varname];
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
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_hcst'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'lens2-model_right_LV_', tmp.season, '_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end

end

%% SVD - lens2-assm plot
if fig_flags(56)==1
    tmp.mode_want=3;
%% PCT (left) plot
    for SVD_mode=1:tmp.mode_want
        fig_cfg.fig_name=[tmp.season, ', SVD_lens2-assm_left_pct_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'visible','off');
    %     hold on
        plot(cfg.iyears(cfg.clim_tlen)+lyear,data2.([tmp.varname,'_svd_lens2_assm_pcs_left_', tmp.season])(SVD_mode,:), 'k-', 'linewidth', 2);
    %         plot(cfg.iyears(1:end-4)+4,EOF.pc(:,SVD_mode));
        grid minor
        xlabel('year'); ylabel(['PCT-mode', num2str(SVD_mode), ', SCF:', num2str(round(data2.([tmp.varname,'_svd_lens2_assm_scf_', tmp.season])(SVD_mode).*100,1)),'%']);
        set(gca, 'fontsize', 20)
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_assm'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'lens2-assm_left_PCT_', tmp.season,'_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
    % hold off

%% LV (left) plot
    for SVD_mode=1:tmp.mode_want
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=squeeze(data2.([tmp.varname,'_svd_lens2_assm_lvmap_left_', tmp.season])(:,:,SVD_mode));
        tmp.C=tmp.C([end, 1:end],:);
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
        
        tmp.prc95 =prctile(tmp.C(isfinite(tmp.C)), 95);
        tmp.prc05 =prctile(tmp.C(isfinite(tmp.C)), 5);

        fig_cfg.fig_name=[tmp.season, ', SVD_LV_, ', tmp. varname];
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
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_assm'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'lens2-assm_left_LV_', tmp.season, '_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end


%% PCT (right) plot
    for SVD_mode=1:tmp.mode_want
        fig_cfg.fig_name=[tmp.season, ', SVD_lens2-assm_right_pct_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'visible','off');
    %     hold on
        plot(cfg.iyears(cfg.clim_tlen)+lyear,data2.([tmp.varname,'_svd_lens2_assm_pcs_right_', tmp.season])(SVD_mode,:), 'k-', 'linewidth', 2);
    %         plot(cfg.iyears(1:end-4)+4,EOF.pc(:,SVD_mode));
        grid minor
        xlabel('year'); ylabel(['PCT-mode', num2str(SVD_mode), ', SCF:', num2str(round(data2.([tmp.varname,'_svd_lens2_assm_scf_', tmp.season])(SVD_mode).*100,1)),'%']);
        set(gca, 'fontsize', 20)
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_assm'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'lens2-assm_right_PCT_', tmp.season,'_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
    % hold off

%% LV (right) plot
    for SVD_mode=1:tmp.mode_want
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=squeeze(data2.([tmp.varname,'_svd_lens2_assm_lvmap_right_', tmp.season])(:,:,SVD_mode));
        tmp.C=tmp.C([end, 1:end],:);
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
        
        tmp.prc95 =prctile(tmp.C(isfinite(tmp.C)), 95);
        tmp.prc05 =prctile(tmp.C(isfinite(tmp.C)), 5);

        fig_cfg.fig_name=[tmp.season, ', SVD_LV_, ', tmp. varname];
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
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_assm'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'lens2-assm_right_LV_', tmp.season, '_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end

end

%% SVD - lens2-obs plot
if fig_flags(57)==1
    tmp.mode_want=3;
%% PCT (left) plot
    for SVD_mode=1:tmp.mode_want
        fig_cfg.fig_name=[tmp.season, ', SVD_lens2-obs_left_pct_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'visible','off');
    %     hold on
        plot(cfg.iyears(cfg.clim_tlen)+lyear,data2.([tmp.varname,'_svd_lens2_obs_pcs_left_', tmp.season])(SVD_mode,:), 'k-', 'linewidth', 2);
    %         plot(cfg.iyears(1:end-4)+4,EOF.pc(:,SVD_mode));
        grid minor
        xlabel('year'); ylabel(['PCT-mode', num2str(SVD_mode), ', SCF:', num2str(round(data2.([tmp.varname,'_svd_lens2_obs_scf_', tmp.season])(SVD_mode).*100,1)),'%']);
        set(gca, 'fontsize', 20)
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_obs'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'lens2-obs_left_PCT_', tmp.season,'_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
    % hold off

%% LV (left) plot
    for SVD_mode=1:tmp.mode_want
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=squeeze(data2.([tmp.varname,'_svd_lens2_obs_lvmap_left_', tmp.season])(:,:,SVD_mode));
        tmp.C=tmp.C([end, 1:end],:);
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
        
        tmp.prc95 =prctile(tmp.C(isfinite(tmp.C)), 95);
        tmp.prc05 =prctile(tmp.C(isfinite(tmp.C)), 5);

        fig_cfg.fig_name=[tmp.season, ', SVD_LV_, ', tmp. varname];
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
        cfg.figname=[dirs.figdir, filesep, 'lens2-obs_left_LV_', tmp.season, '_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end


%% PCT (right) plot
    for SVD_mode=1:tmp.mode_want
        fig_cfg.fig_name=[tmp.season, ', SVD_lens2-obs_right_pct_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'visible','off');
    %     hold on
        plot(cfg.iyears(cfg.clim_tlen)+lyear,data2.([tmp.varname,'_svd_lens2_obs_pcs_right_', tmp.season])(SVD_mode,:), 'k-', 'linewidth', 2);
    %         plot(cfg.iyears(1:end-4)+4,EOF.pc(:,SVD_mode));
        grid minor
        xlabel('year'); ylabel(['PCT-mode', num2str(SVD_mode), ', SCF:', num2str(round(data2.([tmp.varname,'_svd_lens2_obs_scf_', tmp.season])(SVD_mode).*100,1)),'%']);
        set(gca, 'fontsize', 20)
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_SVD', filesep, 'lens2_obs'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'lens2-obs_right_PCT_', tmp.season,'_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
    % hold off

%% LV (right) plot
    for SVD_mode=1:tmp.mode_want
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=squeeze(data2.([tmp.varname,'_svd_lens2_obs_lvmap_right_', tmp.season])(:,:,SVD_mode));
        tmp.C=tmp.C([end, 1:end],:);
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
        
        tmp.prc95 =prctile(tmp.C(isfinite(tmp.C)), 95);
        tmp.prc05 =prctile(tmp.C(isfinite(tmp.C)), 5);

        fig_cfg.fig_name=[tmp.season, ', SVD_LV_, ', tmp. varname];
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
        cfg.figname=[dirs.figdir, filesep, 'lens2-obs_right_LV_', tmp.season, '_mode',num2str(SVD_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end

end

%% obs & model - obs & lens2 corr map --------------------------------------
if fig_flags(58)==1
for fake=1:1
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
    tmp.C=data2.([tmp.varname, '_corr_obs_ano', '_', tmp.season]) - data2.([tmp.varname, '_corr_obs_lens2_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_obs_ano', '_', tmp.season]) ...
        -data2.([tmp.varname, '_corr_obs_lens2_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_obs_Rdiff_ano_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_Rdiff_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_Rdiff_model_lens2_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end

%% assm & model - assm & lens2 corr map --------------------------------------
if fig_flags(59)==1
for fake=1:1
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
    tmp.C=data2.([tmp.varname, '_corr_assm_ano', '_', tmp.season]) - data2.([tmp.varname, '_corr_assm_lens2_ano', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data2.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_assm_ano', '_', tmp.season]) ...
        -data2.([tmp.varname, '_corr_assm_lens2_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_assm_Rdiff_ano_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_assm_Rdiff_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_assm_Rdiff_model_lens2_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all; 
end
end







%% for histogram (assm)
for fake=1:1

    histind=1;
    for corrref=-1:0.1:-0.1
        %% hcst_int
        corr_histo.(['hcst_int_ind_', tmp.season]){histind} = ...
            find(data2.([tmp.varname, '_corr_assm_int2_ano', '_', tmp.season])>=corrref & ...
            data2.([tmp.varname, '_corr_assm_int2_ano', '_', tmp.season])<corrref+0.1);
        corr_histo.(['hcst_int_area_', tmp.season])(histind) = ...
            sum(grid.tarea_60(corr_histo.(['hcst_int_ind_', tmp.season]){histind}), 'omitnan');
        %% lens2
        corr_histo.(['lens2_ind_', tmp.season]){histind} = ...
            find(data2.([tmp.varname, '_corr_assm_lens2_ano', '_', tmp.season])>=corrref & ...
            data2.([tmp.varname, '_corr_assm_lens2_ano', '_', tmp.season])<corrref+0.1);
        corr_histo.(['lens2_area_', tmp.season])(histind) = ...
            sum(grid.tarea_60(corr_histo.(['lens2_ind_', tmp.season]){histind}), 'omitnan');
        histind=histind+1;
    end

    for corrref=0.1:0.1:1.0
        %% hcst_int
        corr_histo.(['hcst_int_ind_', tmp.season]){histind} = ...
            find(data2.([tmp.varname, '_corr_assm_int2_ano', '_', tmp.season])>corrref-0.1 & ...
            data2.([tmp.varname, '_corr_assm_int2_ano', '_', tmp.season])<corrref);
        corr_histo.(['hcst_int_area_', tmp.season])(histind) = ...
            sum(grid.tarea_60(corr_histo.(['hcst_int_ind_', tmp.season]){histind}), 'omitnan');

        %% lens2
        corr_histo.(['lens2_ind_', tmp.season]){histind} = ...
            find(data2.([tmp.varname, '_corr_assm_lens2_ano', '_', tmp.season])>corrref-0.1 & ...
            data2.([tmp.varname, '_corr_assm_lens2_ano', '_', tmp.season])<=corrref);
        corr_histo.(['lens2_area_', tmp.season])(histind) = ...
            sum(grid.tarea_60(corr_histo.(['lens2_ind_', tmp.season]){histind}), 'omitnan');
        histind=histind+1;
    end
    corr_ref_x=-0.95:0.1:0.95;
end

%% for histogram (obs)
for fake=1:1
    histind=1;
    for corrref=-1:0.1:-0.1
        %hcst_int
        corr_histo.(['obs_hcst_int_ind_', tmp.season]){histind} = ...
            find(data2.([tmp.varname, '_corr_obs_int_ano', '_', tmp.season])>=corrref & ...
            data2.([tmp.varname, '_corr_obs_int_ano', '_', tmp.season])<corrref+0.1);
        corr_histo.(['obs_hcst_int_area_', tmp.season])(histind) = ...
            sum(grid.tarea_60(corr_histo.(['obs_hcst_int_ind_', tmp.season]){histind}), 'omitnan');
        %lens2
        corr_histo.(['obs_lens2_ind_', tmp.season]){histind} = ...
            find(data2.([tmp.varname, '_corr_obs_lens2_ano', '_', tmp.season])>=corrref & ...
            data2.([tmp.varname, '_corr_obs_lens2_ano', '_', tmp.season])<corrref+0.1);
        corr_histo.(['obs_lens2_area_', tmp.season])(histind) = ...
            sum(grid.tarea_60(corr_histo.(['obs_lens2_ind_', tmp.season]){histind}), 'omitnan');
        histind=histind+1;
    end

    for corrref=0.1:0.1:1.0
        %hcst_int
        corr_histo.(['obs_hcst_int_ind_', tmp.season]){histind} = ...
            find(data2.([tmp.varname, '_corr_obs_int2_ano', '_', tmp.season])>corrref-0.1 & ...
            data2.([tmp.varname, '_corr_obs_int2_ano', '_', tmp.season])<corrref);
        corr_histo.(['obs_hcst_int_area_', tmp.season])(histind) = ...
            sum(grid.tarea_60(corr_histo.(['obs_hcst_int_ind_', tmp.season]){histind}), 'omitnan');

        %lens2
        corr_histo.(['obs_lens2_ind_', tmp.season]){histind} = ...
            find(data2.([tmp.varname, '_corr_obs_lens2_ano', '_', tmp.season])>corrref-0.1 & ...
            data2.([tmp.varname, '_corr_obs_lens2_ano', '_', tmp.season])<=corrref);
        corr_histo.(['obs_lens2_area_', tmp.season])(histind) = ...
            sum(grid.tarea_60(corr_histo.(['obs_lens2_ind_', tmp.season]){histind}), 'omitnan');
        histind=histind+1;
    end
    corr_ref_x=-0.95:0.1:0.95;
end

%% hcst_int histogram
    fig_cfg.fig_name='hist';
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
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr_10', tmp.dropboxpath);

    fig_h = figure('name',fig_cfg.fig_name,'visible','off');

    bar(corr_ref_x,corr_histo.(['hcst_int_area_', tmp.season]), 'linewidth', 2);
    xlabel('R'); ylabel('Area (km^2)');
    set(gca, 'fontsize', 20)
    grid minor
%     ylim([-0.2 1])
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_histogram', filesep, 'hcst_int'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'histogram_hcst_int2_ano_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;

%% lens2 histogram
    fig_h = figure('name',fig_cfg.fig_name,'visible','off');
    bar(corr_ref_x,corr_histo.(['lens2_area_', tmp.season]), 'linewidth', 2)
    xlabel('R'); ylabel('Area (km^2)');
    set(gca, 'fontsize', 20)
    grid minor
%     ylim([-0.2 1])
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_histogram', filesep, 'lens2'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'histogram_lens2_ano_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
    



%% EOF
    tmp.mode_want=3;
    [EOF.lv, EOF.pc, EOF.var_exp] = ...
        Func_0024_EOF_3d((data2.([tmp.varname, '_model_ano_', tmp.season])(:,:,1:end) - data2.([tmp.varname, '_lens2_ano_', tmp.season])(:,:,1:end)), tmp.mode_want);


%% EOF- PCT plot
    for EOF_mode=1:tmp.mode_want
        fig_cfg.fig_name=[tmp.season, ', EOF_pct_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'visible','off');
    %     hold on
        plot(cfg.iyears(cfg.clim_tlen)+lyear,EOF.pc(:,EOF_mode), 'k-', 'linewidth', 2);
    %         plot(cfg.iyears(1:end-4)+4,EOF.pc(:,EOF_mode));
        grid minor
        xlabel('year'); ylabel(['PCT-mode', num2str(EOF_mode), ', ', num2str(round(EOF.var_exp(EOF_mode),1)),'% explained']);
        set(gca, 'fontsize', 20)
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_EOF', filesep, 'hcst_int'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'PCT_hcst_', tmp.season,'_mode',num2str(EOF_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
    % hold off

%% EOF- LV plot
    for EOF_mode=1:tmp.mode_want
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=squeeze(EOF.lv(:,:,EOF_mode));
        tmp.C=tmp.C([end, 1:end],:);
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
        
        tmp.prc95 =prctile(tmp.C(isfinite(tmp.C)), 95);
        tmp.prc05 =prctile(tmp.C(isfinite(tmp.C)), 5);

        fig_cfg.fig_name=[tmp.season, ', ensmean_, ', tmp. varname];
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
        cfg.figname=[dirs.figdir, filesep, 'LV_hcst_', tmp.season, '_mode',num2str(EOF_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end


%% EOF (HCST-svd1)
    tmp.mode_want=3;
    [EOF.lv, EOF.pc, EOF.var_exp] = ...
        Func_0024_EOF_3d((data2.([tmp.varname, '_model_ano_', tmp.season])(:,:,1:end) ...
        - data2.([tmp.varname, '_svd_lens2_model_right_1mode_recon_', tmp.season])(:,:,1:end)), tmp.mode_want);



%% EOF- PCT plot (HCST-svd1)
    for EOF_mode=1:tmp.mode_want
        fig_cfg.fig_name=[tmp.season, ', EOF_pct_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'visible','off');
    %     hold on
        plot(cfg.iyears(cfg.clim_tlen)+lyear,EOF.pc(:,EOF_mode), 'k-', 'linewidth', 2);
    %         plot(cfg.iyears(1:end-4)+4,EOF.pc(:,EOF_mode));
        grid minor
        xlabel('year'); ylabel(['PCT-mode', num2str(EOF_mode), ', ', num2str(round(EOF.var_exp(EOF_mode),1)),'% explained']);
        set(gca, 'fontsize', 20)
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_EOF', filesep, 'hcst_int'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'PCT_hcst_svd_', tmp.season,'_mode',num2str(EOF_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
    % hold off

%% EOF- LV plot (HCST-svd1)
    for EOF_mode=1:tmp.mode_want
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=squeeze(EOF.lv(:,:,EOF_mode));
        tmp.C=tmp.C([end, 1:end],:);
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
        
        tmp.prc95 =prctile(tmp.C(isfinite(tmp.C)), 95);
        tmp.prc05 =prctile(tmp.C(isfinite(tmp.C)), 5);

        fig_cfg.fig_name=[tmp.season, ', ensmean_, ', tmp. varname];
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
        cfg.figname=[dirs.figdir, filesep, 'LV_hcst_svd_', tmp.season, '_mode',num2str(EOF_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end



end


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

function mons = f_season_mons(season)
    switch season
        case 'INI'
            mons = [1];
        case 'FMA'
            mons = [2,3,4];
        case 'MAM'
            mons = [3,4,5];
        case 'AMJ'
            mons = [4,5,6];
        case 'JJA'
            mons = [6,7,8];
        case 'JAS'
            mons = [7,8,9];
        case 'SON'
            mons = [9,10,11];
        case 'OND'
            mons = [10,11,12];
        case 'DJF'
            mons = [12,13,14];
        case 'JFM'
            mons = [13,14,15];
        case 'JFM1'
            mons = [1,2,3];
        case 'AMJ1'
            mons = [4,5,6]; 
        case 'JAS1'
            mons = [7,8,9];
        case 'OND1'
            mons = [10,11,12];
        case 'AMJ2'
            mons = [16,17,18];
        case 'JAS2'
            mons = [19,20,21];
        case 'OND2'
            mons = [22,23,24];
        case 'JFM2'
            mons = [25,26,27];
        case 'AMJ3'
            mons = [28,29,30];
        case 'JAS3'
            mons = [31,32,33];
        case 'OND3'
            mons = [34,35,36];
        case 'JFM3'
            mons = [37,38,39];
        case 'AMJ4'
            mons = [40,41,42];
        case 'JAS4'
            mons = [43,44,45];
        case 'OND4'
            mons = [46,47,48];
        case 'JFM4'
            mons = [49,50,51];
        otherwise
            mons = str2num(season);
    end
end