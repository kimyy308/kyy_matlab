% %  Created 22-Jan-2024 by Yong-Yub Kim
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


cfg.vars={'SST'};

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

dirs.hcstmatroot=[tmp.kimyypath, '/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', cfg.var];
dirs.figroot=[tmp.kimyypath,'/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', cfg.var];



cfg.iyears=cfg.obs_iyears;
cfg.gnm='f09_g17';
cfg.proj_year=5;

cfg.len_t_y = length(cfg.iyears);
cfg.casename_m = ['ens_all'];


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


S = shaperead('landareas.shp');

%% read & plot data
tmp.varname=cfg.var;

clear tmp.ydata tmp.ydata_lens2 tmp.ydata_obs tmp.ydata_assm
    fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'pred_ratio_hcst_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l60m', '.mat'];
    load(fig_cfg.mat_name, 'data_lm')
    
    
    lyt=0.5:1/12:4.5;
%% for ratio ------------------------------------------------------------------------------------------------------------
    tmp.mode_want=3;
    tmp.ddd=movmean(data_lm.ratio_for_EOF(:,:,:),12,3,'Endpoints', 'discard');
    tmp.ddd(tmp.ddd==0)=NaN;
    [grid.id_w, grid.id_e, grid.id_s, grid.id_n] = Func_0012_findind_Y(0.5, [-360 360, -60 60], ...
            grid.tlong, ...
            grid.tlat, 'CESM2'); % find valid lon, lat index near station
    tmp.ddd(:,1:grid.id_s-1,:)=NaN;
    tmp.ddd(:,grid.id_n+1:end,:)=NaN;

    [EOF.lv, EOF.pc, EOF.var_exp] = ...
        Func_0024_EOF_3d(tmp.ddd, tmp.mode_want);
%     [EOF.lv, EOF.pc, EOF.var_exp] = ...
%         Func_0024_EOF_3d(movmean(data_lm.ratio_for_EOF(:,100:170,:),12,3,'Endpoints', 'discard'), tmp.mode_want);
    for mi=1:tmp.mode_want
        if mean(EOF.pc(:,mi),'omitnan') < 0
            EOF.pc(:,mi)=EOF.pc(:,mi) .* -1;
            EOF.lv(:,:,mi)=EOF.lv(:,:,mi) .* -1;
        end
    end

% % %     figure(1); plot(lyt, EOF.pc(:,1));
% % %     figure(2); pcolor(EOF.lv(:,:,1)'); shading flat; colorbar; caxis([0 0.01])
% % %     EOF.var_exp(1)

    %% EOF- PCT plot
    for EOF_mode=1:tmp.mode_want
        fig_cfg.fig_name=[ 'EOF_ratio_pct_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'visible','off');
        
        if mean(EOF.pc(:,EOF_mode),'omitnan')<0
            EOF.pc(:,EOF_mode)=-EOF.pc(:,EOF_mode);
            EOF.lv(:,:,EOF_mode)=-EOF.lv(:,:,EOF_mode);
        end

        plot(lyt,EOF.pc(:,EOF_mode), 'k-', 'linewidth', 2);
        grid minor
        xlabel('year'); ylabel(['PCT-mode', num2str(EOF_mode), ', ', num2str(round(EOF.var_exp(EOF_mode),1)),'% explained']);
        set(gca, 'fontsize', 20)
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_pred_ratio_EOF', filesep, 'ratio'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'PCT_hcst_pred_ratio','_mode',num2str(EOF_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end

    %% EOF- LV plot
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

    for EOF_mode=1:tmp.mode_want
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=squeeze(EOF.lv(:,:,EOF_mode));
        tmp.C=tmp.C([end, 1:end],:);
        
        tmp.prc95 =prctile(tmp.C(isfinite(tmp.C)), 95);
        tmp.prc05 =prctile(tmp.C(isfinite(tmp.C)), 5);

        fig_cfg.fig_name=['LV_ratio_, ', tmp. varname];
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
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_pred_ratio_EOF', filesep, 'ratio'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'LV_hcst_pred_ratio', '_mode',num2str(EOF_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
%% --------------------------------------------------------------------------------------------------------------------



%% for inc ------------------------------------------------------------------------------------------------------------
    tmp.mode_want=3;
    tmp.ddd=movmean(data_lm.inc_for_EOF(:,:,:),12,3,'Endpoints', 'discard');
    tmp.ddd(tmp.ddd==0)=NaN;
    [grid.id_w, grid.id_e, grid.id_s, grid.id_n] = Func_0012_findind_Y(0.5, [-360 360, -60 60], ...
            grid.tlong, ...
            grid.tlat, 'CESM2'); % find valid lon, lat index near station
    tmp.ddd(:,1:grid.id_s-1,:)=NaN;
    tmp.ddd(:,grid.id_n+1:end,:)=NaN;

    [EOF.lv, EOF.pc, EOF.var_exp] = ...
        Func_0024_EOF_3d(tmp.ddd, tmp.mode_want);
%     [EOF.lv, EOF.pc, EOF.var_exp] = ...
%         Func_0024_EOF_3d(movmean(data_lm.inc_for_EOF(:,100:170,:),12,3,'Endpoints', 'discard'), tmp.mode_want);
    for mi=1:tmp.mode_want
        if mean(EOF.pc(:,mi),'omitnan') < 0
            EOF.pc(:,mi)=EOF.pc(:,mi) .* -1;
            EOF.lv(:,:,mi)=EOF.lv(:,:,mi) .* -1;
        end
    end

% % %     figure(1); plot(lyt, EOF.pc(:,1));
% % %     figure(2); pcolor(EOF.lv(:,:,1)'); shading flat; colorbar; caxis([0 0.01])
% % %     EOF.var_exp(1)

    %% EOF- PCT plot
    for EOF_mode=1:tmp.mode_want
        fig_cfg.fig_name=[ 'EOF_inc_pct_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'visible','off');
        
        if mean(EOF.pc(:,EOF_mode),'omitnan')<0
            EOF.pc(:,EOF_mode)=-EOF.pc(:,EOF_mode);
            EOF.lv(:,:,EOF_mode)=-EOF.lv(:,:,EOF_mode);
        end

        plot(lyt,EOF.pc(:,EOF_mode), 'k-', 'linewidth', 2);
        grid minor
        xlabel('year'); ylabel(['PCT-mode', num2str(EOF_mode), ', ', num2str(round(EOF.var_exp(EOF_mode),1)),'% explained']);
        set(gca, 'fontsize', 20)
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_pred_ratio_EOF', filesep, 'inc'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'PCT_hcst_pred_inc','_mode',num2str(EOF_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end

    %% EOF- LV plot
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

    for EOF_mode=1:tmp.mode_want
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=squeeze(EOF.lv(:,:,EOF_mode));
        tmp.C=tmp.C([end, 1:end],:);
        
        tmp.prc95 =prctile(tmp.C(isfinite(tmp.C)), 95);
        tmp.prc05 =prctile(tmp.C(isfinite(tmp.C)), 5);

        fig_cfg.fig_name=['LV_inc_, ', tmp. varname];
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
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_pred_ratio_EOF', filesep, 'inc'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'LV_hcst_pred_inc', '_mode',num2str(EOF_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
%% --------------------------------------------------------------------------------------------------------------------

%% for spread ------------------------------------------------------------------------------------------------------------
    tmp.mode_want=3;
    tmp.ddd=movmean(data_lm.stde_for_EOF(:,:,:),12,3,'Endpoints', 'discard');
    tmp.ddd(tmp.ddd==0)=NaN;
    [grid.id_w, grid.id_e, grid.id_s, grid.id_n] = Func_0012_findind_Y(0.5, [-360 360, -60 60], ...
            grid.tlong, ...
            grid.tlat, 'CESM2'); % find valid lon, lat index near station
    tmp.ddd(:,1:grid.id_s-1,:)=NaN;
    tmp.ddd(:,grid.id_n+1:end,:)=NaN;

    [EOF.lv, EOF.pc, EOF.var_exp] = ...
        Func_0024_EOF_3d(tmp.ddd, tmp.mode_want);
%     [EOF.lv, EOF.pc, EOF.var_exp] = ...
%         Func_0024_EOF_3d(movmean(data_lm.stde_for_EOF(:,100:170,:),12,3,'Endpoints', 'discard'), tmp.mode_want);
    for mi=1:tmp.mode_want
        if mean(EOF.pc(:,mi),'omitnan') < 0
            EOF.pc(:,mi)=EOF.pc(:,mi) .* -1;
            EOF.lv(:,:,mi)=EOF.lv(:,:,mi) .* -1;
        end
    end

% % %     figure(1); plot(lyt, EOF.pc(:,1));
% % %     figure(2); pcolor(EOF.lv(:,:,1)'); shading flat; colorbar; caxis([0 0.01])
% % %     EOF.var_exp(1)

    %% EOF- PCT plot
    for EOF_mode=1:tmp.mode_want
        fig_cfg.fig_name=[ 'EOF_stde_pct_, ', tmp. varname];
        fig_h = figure('name',fig_cfg.fig_name,'visible','off');
        
        if mean(EOF.pc(:,EOF_mode),'omitnan')<0
            EOF.pc(:,EOF_mode)=-EOF.pc(:,EOF_mode);
            EOF.lv(:,:,EOF_mode)=-EOF.lv(:,:,EOF_mode);
        end

        plot(lyt,EOF.pc(:,EOF_mode), 'k-', 'linewidth', 2);
        grid minor
        xlabel('year'); ylabel(['PCT-mode', num2str(EOF_mode), ', ', num2str(round(EOF.var_exp(EOF_mode),1)),'% explained']);
        set(gca, 'fontsize', 20)
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_pred_ratio_EOF', filesep, 'stde'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'PCT_hcst_pred_stde','_mode',num2str(EOF_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end

    %% EOF- LV plot
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

    for EOF_mode=1:tmp.mode_want
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=squeeze(EOF.lv(:,:,EOF_mode));
        tmp.C=tmp.C([end, 1:end],:);
        
        tmp.prc95 =prctile(tmp.C(isfinite(tmp.C)), 95);
        tmp.prc05 =prctile(tmp.C(isfinite(tmp.C)), 5);

        fig_cfg.fig_name=['LV_stde_, ', tmp. varname];
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
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_pred_ratio_EOF', filesep, 'stde'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'LV_hcst_pred_stde', '_mode',num2str(EOF_mode),'_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end
%% --------------------------------------------------------------------------------------------------------------------




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