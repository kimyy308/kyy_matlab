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
cfg.vars = {'TS'};
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
cfg.obs_iyears=1970:2019;

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
cfg.season = {'MAM1', 'JJA1', 'SON1', 'DJF1', 'MAM2', 'JJA2', 'SON2'};
cfg.sstext = {'2:MAM', '5:JJA', '8:SON', '11:DJF', '14:MAM', '17:JJA', '20:SON'};
cfg.len_t_y = length(cfg.iyears);
cfg.casename_m = ['ens_all'];

tmp.gridname = [tmp.kimyypath, '/Model/CESM2/ESP/HCST_EXP/archive_transfer/', cfg.comp, '/grid.nc'];
tmp.maskname = [tmp.kimyypath, '/Model/CESM2/ESP/HCST_EXP/archive_transfer/ocn/RECCAP2_region_masks_all_v20210412_POP2_grid.nc'];

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
        
        grid.tlong(grid.tlong>=180)=grid.tlong(grid.tlong>=180)-360;

        grid.tarea=ncread(tmp.gridname, 'AREA');
        grid.tarea_60=grid.tarea; grid.tarea_60(grid.tlat>60 | grid.tlat<-60)=NaN;
end

grid.nlon=size(grid.tlong,1);
grid.nlat=size(grid.tlat,2);
% grid.ntime=cfg.proj_year.*12;


fig_flags(1:100)=0;

fig_flags(1:100)=1;


S = shaperead('landareas.shp');

%% read & plot data
tmp.varname=cfg.var;

clear tmp.ydata tmp.ydata_lens2 tmp.ydata_obs tmp.ydata_assm

for lss=1:length(cfg.season)
    tmp.season=cfg.season{lss};
    tmp.ss_text=cfg.sstext{lss};
    tmp.mons = f_season_mons(tmp.season);
    fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_y', num2str(min(cfg.iyears), '%04i'), '_y', num2str(max(cfg.iyears), '%04i'), ...
            '_', tmp.season, '.mat'];
    load(fig_cfg.mat_name, 'data', 'data2');

%% clim time range
% % %     switch cfg.obs_name
% % %        case 'GPCC'
% % %            cfg.clim_ys=1961;
% % %            cfg.clim_ye=2019;
% % %        otherwise
% % %            cfg.clim_ys=1961;
% % %            cfg.clim_ye=2020;
% % %     end
% % %    lyear=0;
% % %    cfg.clim_tlen = (cfg.clim_ys-1959)-lyear:(cfg.clim_ye-2020)+cfg.len_t_y-lyear;
% % %    cfg.clim_tlen2=length(cfg.clim_tlen);



%% model & assm corr map --------------------------------------
if fig_flags(11)==1
for fake=1:1
    fig_cfg.name_rgn = 'Glob';
    fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
%     fig_cfg.map_proj = 'robinson';  % robinson, eqdcylin

    fig_cfg.x_lim = [-180 180];
%     fig_cfg.x_lim = [0 360];    
    fig_cfg.y_lim = [-80 89];
    fig_cfg.fig_size = [0,0,6,3.5];
    fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
    fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
    fig_cfg.title_pos = [0.5,0.93];
    fig_cfg.p_lim =0.1;
    fig_cfg.c_lim = [-1 1];
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('NCAR_SMYLE', tmp.dropboxpath);

%     tmp.X=grid.tlong([end, 1:end],:);   
%     tmp.Y=grid.tlat([end, 1:end],:);
%     tmp.C=data2.([tmp.varname, '_corr_assm_ano', '_', tmp.season]);
%     tmp.C=tmp.C([end, 1:end],:);

    tmp.X=grid.tlong([end/2+1:end, 1:end/2],:);
    tmp.Y=grid.tlat([end/2+1:end, 1:end/2],:);
    tmp.C=data2.([tmp.varname, '_corr_assm_ano', '_', tmp.season]);
    tmp.D=data2.([tmp.varname, '_corr_assm_ano_p', '_', tmp.season]);
    tmp.C(tmp.D>0.1)=NaN;
    tmp.C=tmp.C([end/2+1:end, 1:end/2],:);



%     [tmp.mean_corr, tmp.err] = ...
%         Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_assm_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ',', tmp. varname];
    fig_cfg.fig_name=tmp.ss_text;
    
    fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
        'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
    %% map setting
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 

    axis off; 
    hold on;
%     setm(ax_m,'origin',[0,200],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    setm(ax_m,'origin',[0,0],'MapLonLimit',fig_cfg.x_lim, 'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
    text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
    'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
    'fontsize',14,'fontname','freeserif','interpreter','none')

%     %% caxis & colorbar
    caxis(ax_m, fig_cfg.c_lim); 
    colormap(fig_cfg.c_map);
%     cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
%     set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
%     title(cb,'R','fontsize',12);

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
% % %     setm(ax_m,'frame','on','FLineWidth',1);

% % %     label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
% % %     label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
% % %     mlabel; plabel;
% % %     label_y=plabel; label_x=mlabel;
% % %     for lxi=1:length(label_x)
% % %         tmp.tmppos=label_x(lxi,1).Position;
% % %         tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55; % y position correction
% % %         label_x(lxi,1).Position=tmp.tmppos;
% % %         label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
% % %     end
% % %     for lyi=1:length(label_y)
% % %         label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
% % %         tmp.tmppos=label_y(lyi,1).Position;
% % %         tmp.tmppos(1)=-fig_cfg.ax_size(3)+2.6; % x position correction
% % %         label_y(lyi,1).Position=tmp.tmppos;
% % %     end
% % % 
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
    fig_cfg.x_lim = [0 360];        
    fig_cfg.y_lim = [-80 89];
    fig_cfg.fig_size = [0,0,6,3.5];
    fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
    fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
    fig_cfg.title_pos = [0.5,0.93];
    fig_cfg.p_lim =0.1;
    fig_cfg.c_lim = [-1 1];
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('NCAR_SMYLE', tmp.dropboxpath);
% 
%     tmp.X=grid.tlong([end, 1:end],:);
%     tmp.Y=grid.tlat([end, 1:end],:);
%     tmp.C=data2.([tmp.varname, '_corr_assm_det_ano', '_', tmp.season]);
%     tmp.C=tmp.C([end, 1:end],:);

    tmp.X=grid.tlong([end/2+1:end, 1:end/2],:);
    tmp.Y=grid.tlat([end/2+1:end, 1:end/2],:);
    tmp.C=data2.([tmp.varname, '_corr_assm_det_ano', '_', tmp.season]);
    tmp.D=data2.([tmp.varname, '_corr_assm_det_ano_p', '_', tmp.season]);
    tmp.C(tmp.D>0.1)=NaN;
    tmp.C=tmp.C([end/2+1:end, 1:end/2],:);


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_assm_det_ano', '_', tmp.season]), grid.tlong, grid.tlat);
%     fig_cfg.fig_name=[tmp.season, ',', tmp. varname];
    fig_cfg.fig_name=tmp.ss_text;

    fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
        'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
    %% map setting
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 

    axis off; 
    hold on;
%     setm(ax_m,'origin',[0,200],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    setm(ax_m,'origin',[0,0],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)

    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
    text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
    'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
    'fontsize',14,'fontname','freeserif','interpreter','none')

%     %% caxis & colorbar
    caxis(ax_m, fig_cfg.c_lim); 
    colormap(fig_cfg.c_map);
%     cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
%     set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
%     title(cb,'R','fontsize',12);

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);


%% frame and label setting
% % %     setm(ax_m,'frame','on','FLineWidth',1);
% % % 
% % %     label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
% % %     label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
% % %     mlabel; plabel;
% % %     label_y=plabel; label_x=mlabel;
% % %     for lxi=1:length(label_x)
% % %         tmp.tmppos=label_x(lxi,1).Position;
% % %         tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55; % y position correction
% % %         label_x(lxi,1).Position=tmp.tmppos;
% % %         label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
% % %     end
% % %     for lyi=1:length(label_y)
% % %         label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
% % %         tmp.tmppos=label_y(lyi,1).Position;
% % %         tmp.tmppos(1)=-fig_cfg.ax_size(3)+2.6; % x position correction
% % %         label_y(lyi,1).Position=tmp.tmppos;
% % %     end
% % % 
    %% save
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_assm_det_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_assm_det_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
end
end



%% model & obs corr map --------------------------------------
if fig_flags(11)==1
for fake=1:1
    fig_cfg.name_rgn = 'Glob';
    fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
%     fig_cfg.map_proj = 'robinson';  % robinson, eqdcylin

    fig_cfg.x_lim = [-180 180];
%     fig_cfg.x_lim = [0 360];    
    fig_cfg.y_lim = [-80 89];
    fig_cfg.fig_size = [0,0,6,3.5];
    fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
    fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
    fig_cfg.title_pos = [0.5,0.93];
    fig_cfg.p_lim =0.1;
    fig_cfg.c_lim = [-1 1];
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('NCAR_SMYLE', tmp.dropboxpath);

%     tmp.X=grid.tlong([end, 1:end],:);   
%     tmp.Y=grid.tlat([end, 1:end],:);
%     tmp.C=data2.([tmp.varname, '_corr_assm_ano', '_', tmp.season]);
%     tmp.C=tmp.C([end, 1:end],:);

    tmp.X=grid.tlong([end/2+1:end, 1:end/2],:);
    tmp.Y=grid.tlat([end/2+1:end, 1:end/2],:);
    tmp.C=data2.([tmp.varname, '_corr_obs_ano', '_', tmp.season]);
    tmp.D=data2.([tmp.varname, '_corr_obs_ano_p', '_', tmp.season]);
    tmp.C(tmp.D>0.1)=NaN;
    tmp.C=tmp.C([end/2+1:end, 1:end/2],:);



%     [tmp.mean_corr, tmp.err] = ...
%         Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_assm_ano', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ',', tmp. varname];
    fig_cfg.fig_name=tmp.ss_text;
    
    fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
        'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
    %% map setting
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 

    axis off; 
    hold on;
%     setm(ax_m,'origin',[0,200],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    setm(ax_m,'origin',[0,0],'MapLonLimit',fig_cfg.x_lim, 'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
    text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
    'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
    'fontsize',14,'fontname','freeserif','interpreter','none')

%     %% caxis & colorbar
    caxis(ax_m, fig_cfg.c_lim); 
    colormap(fig_cfg.c_map);
%     cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
%     set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
%     title(cb,'R','fontsize',12);

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
% % %     setm(ax_m,'frame','on','FLineWidth',1);

% % %     label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
% % %     label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
% % %     mlabel; plabel;
% % %     label_y=plabel; label_x=mlabel;
% % %     for lxi=1:length(label_x)
% % %         tmp.tmppos=label_x(lxi,1).Position;
% % %         tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55; % y position correction
% % %         label_x(lxi,1).Position=tmp.tmppos;
% % %         label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
% % %     end
% % %     for lyi=1:length(label_y)
% % %         label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
% % %         tmp.tmppos=label_y(lyi,1).Position;
% % %         tmp.tmppos(1)=-fig_cfg.ax_size(3)+2.6; % x position correction
% % %         label_y(lyi,1).Position=tmp.tmppos;
% % %     end
% % % 
    %% save
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
end
end

%% model-det & obs-det corr map --------------------------------------
if fig_flags(12)==1
for fake=1:1
    fig_cfg.name_rgn = 'Glob';
    fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
%     fig_cfg.map_proj = 'robinson';  % robinson, eqdcylin

    fig_cfg.x_lim = [-180 180];
    fig_cfg.x_lim = [0 360];        
    fig_cfg.y_lim = [-80 89];
    fig_cfg.fig_size = [0,0,6,3.5];
    fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
    fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
    fig_cfg.title_pos = [0.5,0.93];
    fig_cfg.p_lim =0.1;
    fig_cfg.c_lim = [-1 1];
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('NCAR_SMYLE', tmp.dropboxpath);
% 
%     tmp.X=grid.tlong([end, 1:end],:);
%     tmp.Y=grid.tlat([end, 1:end],:);
%     tmp.C=data2.([tmp.varname, '_corr_assm_det_ano', '_', tmp.season]);
%     tmp.C=tmp.C([end, 1:end],:);

    tmp.X=grid.tlong([end/2+1:end, 1:end/2],:);
    tmp.Y=grid.tlat([end/2+1:end, 1:end/2],:);
    tmp.C=data2.([tmp.varname, '_corr_obs_det_ano', '_', tmp.season]);
    tmp.D=data2.([tmp.varname, '_corr_obs_det_ano_p', '_', tmp.season]);
    tmp.C(tmp.D>0.1)=NaN;
    tmp.C=tmp.C([end/2+1:end, 1:end/2],:);


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_assm_det_ano', '_', tmp.season]), grid.tlong, grid.tlat);
%     fig_cfg.fig_name=[tmp.season, ',', tmp. varname];
    fig_cfg.fig_name=tmp.ss_text;

    fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
        'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
    %% map setting
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 

    axis off; 
    hold on;
%     setm(ax_m,'origin',[0,200],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    setm(ax_m,'origin',[0,0],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)

    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
    text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
    'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
    'fontsize',14,'fontname','freeserif','interpreter','none')

%     %% caxis & colorbar
    caxis(ax_m, fig_cfg.c_lim); 
    colormap(fig_cfg.c_map);
%     cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
%     set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
%     title(cb,'R','fontsize',12);

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);


%% frame and label setting
% % %     setm(ax_m,'frame','on','FLineWidth',1);
% % % 
% % %     label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
% % %     label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
% % %     mlabel; plabel;
% % %     label_y=plabel; label_x=mlabel;
% % %     for lxi=1:length(label_x)
% % %         tmp.tmppos=label_x(lxi,1).Position;
% % %         tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55; % y position correction
% % %         label_x(lxi,1).Position=tmp.tmppos;
% % %         label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
% % %     end
% % %     for lyi=1:length(label_y)
% % %         label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
% % %         tmp.tmppos=label_y(lyi,1).Position;
% % %         tmp.tmppos(1)=-fig_cfg.ax_size(3)+2.6; % x position correction
% % %         label_y(lyi,1).Position=tmp.tmppos;
% % %     end
% % % 
    %% save
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_det_map', filesep, 'model'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_det_ano_map_', tmp.varname, '_', tmp.season, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
end
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