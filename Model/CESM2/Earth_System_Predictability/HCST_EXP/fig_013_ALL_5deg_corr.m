% %  Created 12-Jan-2024 by Yong-Yub Kim
clc; clear all; close all;
warning off;

%% set path
[error_status, tmp.hostname] = system('hostname');
tmp.hostname=tmp.hostname(1:end-1);
switch tmp.hostname
    case 'Yong-Yubs-iMac-Pro.local'
        tmp.dropboxpath = '/Volumes/kyy_raid/kimyy/Dropbox';
    case {'da1', 'da2', 'da3', 'da4'}
        tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
end
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Common', tmp.fs, 'mca']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Common', tmp.fs, 'order']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

cfg.ly_s = 1;
cfg.ly_e = 5;

%% model configuration

cfg.vars = {'TS', 'PSL', 'PRECT', 'SST'};
% cfg.vars = {'PRECT', 'PSL'};
% cfg.vars = {'SST'};
% cfg.vars = {'NPP', 'GPP', 'RAIN', 'TOTVEGC', 'FIRE'};
% cfg.vars = {'TWS'};
% cfg.vars = {'SOILWATER_10CM'};
cfg.vars = {'TS', 'PRECT', 'SST', 'PSL'};
% cfg.vars = {'PRECT', 'SST', 'PSL'};
% cfg.vars = {'TS'};
cfg.vars={'SOILWATER_10CM', 'TLAI', 'FAREA_BURNED', 'COL_FIRE_CLOSS', 'TWS'};
% cfg.vars={'SSH', 'photoC_TOT_zint', 'photoC_TOT_zint_100m', 'IRON_FLUX', 'HMXL', 'HBLT'};
cfg.vars={'SSH'};
% cfg.vars={'photoC_TOT_zint_100m'};
cfg.vars={'TLAI'};
% tmp.dimids= [1, 2, 4];
cfg.vars={'HMXL', 'HBLT'};
cfg.vars={'NO3', 'SALT'};
cfg.vars={'NO3', 'SALT', 'mul_VVEL_NO3', 'mul_WVEL_NO3','mul_UVEL_NO3'};
cfg.vars={'mul_VVEL_NO3', 'mul_WVEL_NO3','mul_UVEL_NO3'};

% cfg.vars={'TS', 'SST', 'PRECT', 'TWS', 'GPP', 'FAREA_BURNED', 'photoC_TOT_zint_100m'};
cfg.vars={'TS', 'SST', 'PRECT', 'PSL', 'TWS', 'GPP', ...
    'photoC_TOT_zint_100m', 'COL_FIRE_CLOSS', 'FAREA_BURNED', 'SOILWATER_10CM', 'TLAI', 'AEROD_v', 'SSH'};
cfg.vars={'SALT'};
% cfg.vars={'SSH'};
% cfg.vars={'FAREA_BURNED'};
% cfg.vars={'TS'};
% cfg.vars={'PRECT'};
% cfg.vars={'PSL'};
% cfg.vars={'AEROD_v'};
% cfg.vars={'GPP', 'FAREA_BURNED'};
% cfg.vars={'SSH'};


cfg.vlayer=1; % surf, vertical slice 

% cfg.vlayer=1:10; % 10layer. don't put more than 15

cfg.vlayer_1st=min(cfg.vlayer);
cfg.vlayer_cnt=max(cfg.vlayer)-cfg.vlayer_1st+1;


for vari=1:length(cfg.vars)
    cfg.var=cfg.vars{vari};
    cfg.obs_name=f_obs_name(cfg.var);
    cfg.obs_fname_mid=f_obs_name_mid(cfg.var);
    cfg.obs_varname=f_obs_varname(cfg.var);
    cfg.comp=Func_0025_CESM2_cmpname_var(cfg.var);
    cfg.obs_fname_module=f_obs_fname_module(cfg.comp);
    cfg.obs_iyears=1965:2020;
    cfg.obs_iyears2=f_obs_iyears(cfg.var);
    
    % dirs.hcstroot=['/mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive/', cfg.comp, '/', cfg.var];
    % dirs.obsroot=['/Volumes/kyy_raid/kimyy/Observation/', cfg.obs_name, '/monthly_reg_', cfg.obs_fname_module(2:4)];
    dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', cfg.var, '/5deg'];
    % dirs.lens2root=['/mnt/lustre/proj/earth.system.predictability/LENS2/archive_analysis/', cfg.comp, '/', cfg.var];
    % dirs.assmroot=['/mnt/lustre/proj/earth.system.predictability/ASSM_EXP/archive_analysis/', cfg.comp, '/', cfg.var];
    
    dirs.hcstroot=['/mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive_yearly_transfer/', cfg.comp, '/', cfg.var];
    dirs.obsroot=['/mnt/lustre/proj/kimyy/Observation/', cfg.obs_name];
    dirs.matroot=['/mnt/lustre/proj/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp,'/', cfg.var];
    dirs.lens2root=['/mnt/lustre/proj/earth.system.predictability/LENS2/archive_yearly_transfer/', cfg.comp, '/', cfg.var];
    dirs.assmroot=['/mnt/lustre/proj/earth.system.predictability/ASSM_EXP/archive_yearly_transfer/', cfg.comp, '/', cfg.var];    
    
    dirs.corrroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/regrid_5deg/statistics/corr'];
   
    load([dirs.corrroot, '/', 'corr_5deg_all_', cfg.var, '.mat'], 'cfg_assm', 'cfg_hcst', 'cfg_lens2', 'corrval', 'grid');



    cfg.iyears=cfg.obs_iyears;
    cfg.gnm='f09_g17';
    cfg.proj_year=5;

    cfg.max_iy=max(cfg.iyears);
    cfg.min_iy=min(cfg.iyears);
    cfg.len_t_y = length(cfg.iyears);
       
    % [tmp.error_status, tmp.value]=system(['ls ', dirs.datadir, '/*once*']);  % b.e21.BHISTsmbb.f09_g17.assm.oras4_ba-10p1.pop.h.once.nc
    tmp.gridname = [dirs.hcstroot, tmp.fs, '../grid.nc'];
    tmp.maskname = '/mnt/lustre/proj/kimyy/Model/CESM2/ESP/grids/ocn/RECCAP2_region_masks_all_v20210412_POP2_grid.nc';
    
    % plot set, S.
    S = shaperead('landareas.shp');

    %% read & plot data
    tmp.varname=cfg.var;
    
%% corr each member based, median, ASSM <-> LENS2 (5deg)
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
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=squeeze(median(corrval.assm_lens2.val, 1));
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(median(corrval.assm_lens2.val, 1)), grid.tlong, grid.tlat);
        fig_cfg.fig_name=['mem_med_as_le_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_assm_map', filesep, 'lens2'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_memb_median_assm_lens2_map_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;

%% corr each member based, median, OBS <-> LENS2 (5deg)
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
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=squeeze(median(corrval.obs_lens2.val, 1));
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(median(corrval.obs_lens2.val, 1)), grid.tlong, grid.tlat);
        fig_cfg.fig_name=['mem_me_ob_le_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_obs_map', filesep, 'lens2'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_memb_median_obs_lens2_map_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;


%% corr each member based, median, OBS <-> ASSM (5deg)
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
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=squeeze(median(corrval.obs_assm.val, 1));
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(median(corrval.obs_assm.val, 1)), grid.tlong, grid.tlat);
        fig_cfg.fig_name=['mem_me_ob_as_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_obs_map', filesep, 'assm'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_memb_median_obs_assm_map_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;


%% corr each member based, mean, ASSM <-> LENS2 (5deg)
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
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=squeeze(mean(corrval.assm_lens2.val, 1));
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(mean(corrval.assm_lens2.val, 1)), grid.tlong, grid.tlat);
        fig_cfg.fig_name=['mem_med_as_le_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_assm_map', filesep, 'lens2'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_memb_mean_assm_lens2_map_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;

%% corr each member based, mean, OBS <-> LENS2 (5deg)
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
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=squeeze(mean(corrval.obs_lens2.val, 1));
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(mean(corrval.obs_lens2.val, 1)), grid.tlong, grid.tlat);
        fig_cfg.fig_name=['mem_me_ob_le_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_obs_map', filesep, 'lens2'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_memb_mean_obs_lens2_map_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;


%% corr each member based, mean, OBS <-> ASSM (5deg)
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
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=squeeze(mean(corrval.obs_assm.val, 1));
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(mean(corrval.obs_assm.val, 1)), grid.tlong, grid.tlat);
        fig_cfg.fig_name=['mem_me_ob_as_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_obs_map', filesep, 'assm'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_memb_mean_obs_assm_map_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;


%% ensmeble mean part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% corr ensmeble mean based, ASSM <-> LENS2 (5deg)
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
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=squeeze(corrval.assm_lens2_em.val);
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(corrval.assm_lens2_em.val, grid.tlong, grid.tlat);
        fig_cfg.fig_name=['ens_me_as_le_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_assm_map', filesep, 'lens2'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_ens_mean_assm_lens2_map_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;


        
%% corr ensmeble mean based,  OBS <-> LENS2 (5deg)
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
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=corrval.obs_lens2_em.val;
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(corrval.obs_lens2_em.val, grid.tlong, grid.tlat);
        fig_cfg.fig_name=['ens_me_ob_le_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_obs_map', filesep, 'lens2'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_ens_mean_obs_lens2_map_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;


%% corr ensmeble mean based, OBS <-> ASSM (5deg)
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
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=corrval.obs_assm_em.val;
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(corrval.obs_assm_em.val, grid.tlong, grid.tlat);
        fig_cfg.fig_name=['ens_me_ob_as_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_obs_map', filesep, 'assm'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_ens_mean_obs_assm_map_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;




%% lead year loop for HCST
    for lyear=cfg.ly_s:cfg.ly_e
%     for lyear=[1,3]

        tmp.lyear_str = num2str(lyear, '%02i');
           

%% corr each member based, median, ASSM <-> HCST (5deg)
for foldi=1:1
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
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=squeeze(median(corrval.assm_hcst.(['ly',num2str(lyear)]).val, 1));
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(median(corrval.assm_hcst.(['ly',num2str(lyear)]).val, 1)), grid.tlong, grid.tlat);
        fig_cfg.fig_name=['l',tmp.lyear_str, ', mem_med_as_hc_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_assm_map', filesep, 'model'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_memb_median_assm_hcst_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
end

%% corr each member based, median, OBS <-> HCST (5deg)
for foldi=1:1
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
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=squeeze(median(corrval.obs_hcst.(['ly',num2str(lyear)]).val, 1));
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(median(corrval.obs_hcst.(['ly',num2str(lyear)]).val, 1)), grid.tlong, grid.tlat);
        fig_cfg.fig_name=['l',tmp.lyear_str, ', mem_med_ob_hc_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_obs_map', filesep, 'model'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_memb_median_obs_hcst_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
end

%% corr each member based, mean, ASSM <-> HCST (5deg)
for foldi=1:1

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
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=squeeze(mean(corrval.assm_hcst.(['ly',num2str(lyear)]).val, 1));
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(median(corrval.assm_hcst.(['ly',num2str(lyear)]).val, 1)), grid.tlong, grid.tlat);
        fig_cfg.fig_name=['l',tmp.lyear_str, ', mem_me_as_hc_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_assm_map', filesep, 'model'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_memb_mean_assm_hcst_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
end

%% corr each member based, mean, OBS <-> HCST (5deg)
for foldi=1:1

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
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=squeeze(mean(corrval.obs_hcst.(['ly',num2str(lyear)]).val, 1));
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(mean(corrval.obs_hcst.(['ly',num2str(lyear)]).val, 1)), grid.tlong, grid.tlat);
        fig_cfg.fig_name=['l',tmp.lyear_str, ', mem_me_ob_hc_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_obs_map', filesep, 'model'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_memb_mean_obs_hcst_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
end



%% corr each member based, mean, ASSM <-> HCST-LENS2(internal signal) (5deg)
for foldi=1:1

        fig_cfg.name_rgn = 'Glob';
        fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
    
        fig_cfg.x_lim = [-180 180];
        fig_cfg.y_lim = [-80 89];
        fig_cfg.fig_size = [0,0,6,3.5];
        fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
        fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
        fig_cfg.title_pos = [0.5,0.93];
        fig_cfg.p_lim =0.1;
        fig_cfg.c_lim = [-0.5 0.5];
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=squeeze(mean(corrval.assm_hcst_int.(['ly',num2str(lyear)]).val, 1));
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(median(corrval.assm_hcst_int.(['ly',num2str(lyear)]).val, 1)), grid.tlong, grid.tlat);
        fig_cfg.fig_name=['l',tmp.lyear_str, ', mem_me_as_hc_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_assm_map', filesep, 'model'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_memb_mean_assm_hcst_int_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
end

%% corr each member based, median, ASSM <-> HCST - ASSM <-> LENS2(5deg)
for foldi=1:1
    fig_cfg.name_rgn = 'Glob';
        fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
    
        fig_cfg.x_lim = [-180 180];
        fig_cfg.y_lim = [-80 89];
        fig_cfg.fig_size = [0,0,6,3.5];
        fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
        fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
        fig_cfg.title_pos = [0.5,0.93];
        fig_cfg.p_lim =0.1;
        fig_cfg.c_lim = [-0.5 0.5];
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=squeeze(median(corrval.assm_hcst.(['ly',num2str(lyear)]).val, 1)) - squeeze(median(corrval.assm_lens2.val, 1));
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(...
            squeeze(median(corrval.assm_hcst.(['ly',num2str(lyear)]).val, 1)) - squeeze(median(corrval.assm_lens2.val, 1)), ...
            grid.tlong, grid.tlat);
        fig_cfg.fig_name=['l',tmp.lyear_str, ', mem_me_as_hc_int_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_assm_map', filesep, 'model'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_memb_median_assm_hcst_diff_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
end

%% corr each member based, mean, ASSM <-> HCST - ASSM <-> LENS2(5deg)
for foldi=1:1
    fig_cfg.name_rgn = 'Glob';
        fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
    
        fig_cfg.x_lim = [-180 180];
        fig_cfg.y_lim = [-80 89];
        fig_cfg.fig_size = [0,0,6,3.5];
        fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
        fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
        fig_cfg.title_pos = [0.5,0.93];
        fig_cfg.p_lim =0.1;
        fig_cfg.c_lim = [-0.5 0.5];
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=squeeze(mean(corrval.assm_hcst.(['ly',num2str(lyear)]).val, 1)) - squeeze(mean(corrval.assm_lens2.val, 1));
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(...
            squeeze(mean(corrval.assm_hcst.(['ly',num2str(lyear)]).val, 1)) - squeeze(mean(corrval.assm_lens2.val, 1)), ...
            grid.tlong, grid.tlat);
        fig_cfg.fig_name=['l',tmp.lyear_str, ', mem_me_as_hc_int_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_assm_map', filesep, 'model'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_memb_mean_assm_hcst_diff_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
end

%% corr each member based, median, OBS <-> HCST - OBS <-> LENS2(5deg)
for foldi=1:1
        fig_cfg.name_rgn = 'Glob';
        fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
    
        fig_cfg.x_lim = [-180 180];
        fig_cfg.y_lim = [-80 89];
        fig_cfg.fig_size = [0,0,6,3.5];
        fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
        fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
        fig_cfg.title_pos = [0.5,0.93];
        fig_cfg.p_lim =0.1;
        fig_cfg.c_lim = [-0.5 0.5];
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=squeeze(median(corrval.obs_hcst.(['ly',num2str(lyear)]).val, 1)) - squeeze(median(corrval.obs_lens2.val, 1));
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(...
            squeeze(median(corrval.obs_hcst.(['ly',num2str(lyear)]).val, 1)) - squeeze(median(corrval.obs_lens2.val, 1)), ...
            grid.tlong, grid.tlat);
        fig_cfg.fig_name=['l',tmp.lyear_str, ', mem_me_ob_hc_int_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_obs_map', filesep, 'model'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_memb_median_obs_hcst_diff_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
end

%% corr each member based, mean, OBS <-> HCST - OBS <-> LENS2(5deg)
for foldi=1:1
        fig_cfg.name_rgn = 'Glob';
        fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
    
        fig_cfg.x_lim = [-180 180];
        fig_cfg.y_lim = [-80 89];
        fig_cfg.fig_size = [0,0,6,3.5];
        fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
        fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
        fig_cfg.title_pos = [0.5,0.93];
        fig_cfg.p_lim =0.1;
        fig_cfg.c_lim = [-0.5 0.5];
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=squeeze(mean(corrval.obs_hcst.(['ly',num2str(lyear)]).val, 1)) - squeeze(mean(corrval.obs_lens2.val, 1));
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(...
            squeeze(mean(corrval.obs_hcst.(['ly',num2str(lyear)]).val, 1)) - squeeze(mean(corrval.obs_lens2.val, 1)), ...
            grid.tlong, grid.tlat);
        fig_cfg.fig_name=['l',tmp.lyear_str, ', mem_me_ob_hc_int_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_obs_map', filesep, 'model'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_memb_mean_obs_hcst_diff_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
end




%% ensemble mean based %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% corr ensemble mean based,  ASSM <-> HCST (5deg)
for foldi=1:1
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
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=corrval.assm_hcst_em.(['ly',num2str(lyear)]).val;
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(corrval.assm_hcst_em.(['ly',num2str(lyear)]).val, grid.tlong, grid.tlat);
        fig_cfg.fig_name=['l',tmp.lyear_str, ', ens_me_as_hc_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_assm_map', filesep, 'model'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_ens_mean_assm_hcst_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
end

%% corr ensemble mean based,  OBS <-> HCST (5deg)
for foldi=1:1
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
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=corrval.obs_hcst_em.(['ly',num2str(lyear)]).val;
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(corrval.obs_hcst_em.(['ly',num2str(lyear)]).val, grid.tlong, grid.tlat);
        fig_cfg.fig_name=['l',tmp.lyear_str, ', ens_me_ob_hc_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_obs_map', filesep, 'model'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_ens_mean_obs_hcst_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
end

%% corr ensemble mean based,  ASSM <-> HCST-LENS2(internal signal) (5deg)
for foldi=1:1
        fig_cfg.name_rgn = 'Glob';
        fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
    
        fig_cfg.x_lim = [-180 180];
        fig_cfg.y_lim = [-80 89];
        fig_cfg.fig_size = [0,0,6,3.5];
        fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
        fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
        fig_cfg.title_pos = [0.5,0.93];
        fig_cfg.p_lim =0.1;
        fig_cfg.c_lim = [-0.5 0.5];
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=corrval.assm_hcst_em_int.(['ly',num2str(lyear)]).val;
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(corrval.assm_hcst_em_int.(['ly',num2str(lyear)]).val, grid.tlong, grid.tlat);
        fig_cfg.fig_name=['l',tmp.lyear_str, ', ens_me_as_hc_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_assm_map', filesep, 'model'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_ens_mean_assm_hcst_int_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
end


%% corr ensemble mean based,  ASSM <-> HCST - ASSM <-> LENS2(5deg)
for foldi=1:1
        fig_cfg.name_rgn = 'Glob';
        fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
    
        fig_cfg.x_lim = [-180 180];
        fig_cfg.y_lim = [-80 89];
        fig_cfg.fig_size = [0,0,6,3.5];
        fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
        fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
        fig_cfg.title_pos = [0.5,0.93];
        fig_cfg.p_lim =0.1;
        fig_cfg.c_lim = [-0.5 0.5];
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=corrval.assm_hcst_em.(['ly',num2str(lyear)]).val - corrval.assm_lens2_em.val;
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(...
            corrval.assm_hcst_em.(['ly',num2str(lyear)]).val - corrval.assm_lens2_em.val, ...
            grid.tlong, grid.tlat);
        fig_cfg.fig_name=['l',tmp.lyear_str, ', ens_me_as_hc_int_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_assm_map', filesep, 'model'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_ens_mean_assm_hcst_diff_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
end

%% corr ensemble mean based,  OBS <-> HCST - OBS <-> LENS2(5deg)
for foldi=1:1
        fig_cfg.name_rgn = 'Glob';
        fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
    
        fig_cfg.x_lim = [-180 180];
        fig_cfg.y_lim = [-80 89];
        fig_cfg.fig_size = [0,0,6,3.5];
        fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
        fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
        fig_cfg.title_pos = [0.5,0.93];
        fig_cfg.p_lim =0.1;
        fig_cfg.c_lim = [-0.5 0.5];
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        
        tmp.C=corrval.obs_hcst_em.(['ly',num2str(lyear)]).val - corrval.obs_lens2_em.val;
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(...
            corrval.obs_hcst_em.(['ly',num2str(lyear)]).val - corrval.obs_lens2_em.val, ...
            grid.tlong, grid.tlat);
        fig_cfg.fig_name=['l',tmp.lyear_str, ', ens_me_ob_hc_int_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_obs_map', filesep, 'model'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_ens_mean_obs_hcst_diff_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
end


    end
end




function scenname = f_scen(year)
    if year<=2014
        scenname='BHISTsmbb';
    else
        scenname='BSSP370smbb';
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
            obsname_simple='GPCC';
        case 'RAIN'
            obsname_simple='GPCC';
        case 'PSL'
            obsname_simple='ERA5';
        case 'SOILWATER_10CM'
%             obsname_simple='CMEMS';
            obsname_simple='GLEAM';
        case 'TWS'
            obsname_simple='NOAA';
        case 'SSH'
            obsname_simple='CMEMS';
        case 'TS'
%             obsname_simple='HadCRUT5';
            obsname_simple='ERA5';
        case 'sumChl'
            obsname_simple='OC_CCI';
        case 'TLAI'
            obsname_simple='NOAA'; % NOAA LAI
        case 'FAREA_BURNED'
            obsname_simple='MODIS'; % MODIS Fire_cci v5.1
            obsname_simple='AVHRR'; % AVHRR-LTDR
        case 'COL_FIRE_CLOSS'
            obsname_simple='GFED'; % GFED
        case 'photoC_TOT_zint'
%             obsname_simple='VGPM'; % VGPM
            obsname_simple='CMEMS'; %Globcolour            
        case 'photoC_TOT_zint_100m'
%             obsname_simple='VGPM'; % VGPM
            obsname_simple='CMEMS'; %Globcolour
        case 'GPP'
%             obsname_simple='ORNL_DAAC';
            obsname_simple='VODCA2GPP';
        otherwise
            obsname_simple='nan';
    end
end


function obsname_simple = f_obs_name_mid(varn)
    switch varn
        case 'SST'
            obsname_simple='ersst_reg_5deg.v5.';
        case 'PRECT'
            obsname_simple='GPCC_reg_5deg.v5.';
        case 'RAIN'
            obsname_simple='GPCC_reg_5deg.v5.';
        case 'PSL'
            obsname_simple='ERA5_msl_reg_5deg.';
        case 'SOILWATER_10CM'
%             obsname_simple='SM_reg_5deg.';
            obsname_simple='GLEAM_reg_5deg.v5.';
        case 'TWS'
            obsname_simple='TSW_reg_5deg.';
        case 'SSH'
            obsname_simple='CMEMS_reg_5deg.';
        case 'TS'
%             obsname_simple='HadCRUT5_reg_5deg.';
            obsname_simple='ERA5_t2m_reg_5deg.';
        case 'sumChl'
            obsname_simple='OC_CCI_reg_5deg.';
        case 'TLAI'
            obsname_simple='LAI_reg_5deg.'; % NOAA LAI
        case 'FAREA_BURNED'
            obsname_simple='burned_area_reg_5deg.'; % MODIS Fire_cci v5.1
            obsname_simple='AVHRR-LTDR_reg_5deg.'; % AVHRR-LTDR
        case 'COL_FIRE_CLOSS'
            obsname_simple='FIRE_CLOSS_reg_5deg.'; % GFED
        case 'photoC_TOT_zint'
%             obsname_simple='s_vgpm_reg_5deg.'; % VGPM
%             obsname_simple='ensmean_reg_5deg.'; % VGPM   
            obsname_simple='CMEMS_reg_5deg.'; %Globcolour;
        case 'photoC_TOT_zint_100m'
%             obsname_simple='s_vgpm_reg_5deg.'; % VGPM
%             obsname_simple='ensmean_reg_5deg.'; % VGPM
            obsname_simple='CMEMS_reg_5deg.'; %Globcolour;
        case 'GPP'
%             obsname_simple='ORNL_DAAC_reg_5deg.'; % ORNL_DAAC
            obsname_simple='VODCA2GPP_reg_5deg.'; % VODCA2GPP
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
        case 'RAIN'
            obsname_simple='precip';  % GPCC
        case 'PSL'
            obsname_simple='msl';
        case 'SOILWATER_10CM'
%             obsname_simple='sm';
            obsname_simple='SMsurf'; %GLEAM
        case 'TWS'
            obsname_simple='w';
        case 'SSH'
            obsname_simple='sla';
        case 'TS'
%             obsname_simple='tas_mean'; %HadCRUT5
            obsname_simple='t2m';  %ERA5
        case 'sumChl'
            obsname_simple='chlor_a';
        case 'TLAI'
            obsname_simple='LAI'; % NOAA LAI
        case 'FAREA_BURNED'
            obsname_simple='burned_area'; % MODIS Fire_cci v5.1, AVHRR-LTDR
        case 'COL_FIRE_CLOSS'
            obsname_simple='C'; % GFED
        case 'photoC_TOT_zint'
%             obsname_simple='npp'; % VGPM
            obsname_simple='PP'; %Globcolour
        case 'photoC_TOT_zint_100m'
%             obsname_simple='npp'; % VGPM
            obsname_simple='PP'; %Globcolour
        case 'GPP'
            obsname_simple='GPP';
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
            obsname_simple=1979:2020; % GPCP
            obsname_simple=1960:2019; % GPCC
        case 'RAIN'
            obsname_simple=1960:2019;
        case 'SSH'
            obsname_simple=1993:2020;
        case 'sumChl'
            obsname_simple=1998:2020;
        case 'SOILWATER_10CM'
%             obsname_simple=1979:2020; % C3S Surface Soil Moisture
            obsname_simple=1980:2020; % GLEAM SMsurf
        case 'TLAI'
            obsname_simple=1982:2018; % NOAA LAI
        case 'FAREA_BURNED'
%             obsname_simple=2001:2020; % MODIS Fire_cci v5.1
            obsname_simple=1982:2018; % AVHRR-LTDR
        case 'COL_FIRE_CLOSS'
            obsname_simple=1997:2020; % GFED
        case 'photoC_TOT_zint'
            obsname_simple=2003:2020; % VGPM (ens)
            obsname_simple=1998:2020; % GlobColour            
        case 'photoC_TOT_zint_100m'
%             obsname_simple=2003:2020; % VGPM (ens)
            obsname_simple=1998:2020; % GlobColour
        case 'GPP'
%             obsname_simple=1982:2016; % ORNL_DAAC
            obsname_simple=1989:2019; % VODCA2GPP
        otherwise
            obsname_simple=1960:2020;
    end
end
