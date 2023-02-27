% %  Created 30-Jan-2023 by Yong-Yub Kim
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

cfg.levs=[1, 5, 10, 15];

%% model configuration
% cfg.vars={'ALK', 'BSF', 'DIC', 'HMXL', 'NO3', 'PO4', 'SALT', 'SSH', 'SiO3', 'TEMP', ...
%     'WVEL', 'diatChl', 'diazChl', 'photoC_TOT_zint', 'photoC_TOT_zint_100m', 'spChl', ...
% 'UVEL', 'VVEL', 'PD'};
% cfg.vars={'UVEL', 'VVEL', 'PD'};
cfg.vars={'NO3', 'SiO3', 'PO4', 'WVEL', 'diatChl', 'diazChl', 'spChl', ...
    'DIC', 'ALK', 'TEMP', 'SALT', 'BSF', 'HMXL', 'photoC_TOT_zint_100m', 'photoC_TOT_zint', ...
    'SSH', 'UVEL', 'VVEL', 'UISOP', 'VISOP', 'PD','Fe', 'HBLT', 'HMXL_DR', 'TBLT', 'TMXL', 'XBLT', 'XMXL', ...
    'TMXL_DR', 'XMXL_DR', 'IRON_FLUX', 'diatC_zint_100m_2', 'diazC_zint_100m_2', 'spC_zint_100m_2', ...
    'UE_PO4', 'VN_PO4', 'WT_PO4'};
% abc=dir('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/ASSM_EXP/archive_yearly_analysis/ocn')
% cfg.vars=abc(3:end).name;

cfg.obs={'en4.2_ba', 'projdv7.3_ba'};

for obsi= 1:length(cfg.obs)
    name.obs=cfg.obs{obsi};
for vari= 1:length(cfg.vars)
    name.var=cfg.vars{vari};
    dirs.assmroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/ASSM_EXP/archive_yearly_analysis/ocn/', name.var];
    dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/ASSM_EXP/archive_yearly_analysis/ocn/', name.var];
    dirs.grid='/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/ASSM_EXP/archive';
    dirs.ens=[dirs.assmroot, tmp.fs, 'ens_all'];
    dirs.corr=[dirs.assmroot, tmp.fs, 'corr', tmp.fs, 'ens'];
    
    fname.grid=[dirs.grid, tmp.fs, 'pop.h.once.nc'];    
    S = shaperead('landareas.shp');

    for levi=1:length(cfg.levs)
        tmp.lev=cfg.levs(levi);
        str.lev=num2str(tmp.lev, '%02i');
        
%% Ensemble set  --------------------------------------
%% Signal to Noise Ratio (SNR) map  --------------------------------------
        name.fign = 'SNR';
        cfg.dim_var=Func_0026_CESM2_dim_var(name.var);
        dirs.figdir= [dirs.figroot, filesep, 'ens_set', filesep, name.obs];
        if (strcmp(cfg.dim_var, '4d'))
            cfg.figname=[dirs.figdir, filesep, name.fign, '_', name.var, '_', name.obs, '_lev', str.lev, '.tif'];
            tmp.title=[name.var, ', ', name.fign, ', ', str.lev];
        elseif (strcmp(cfg.dim_var, '3d'))
            cfg.figname=[dirs.figdir, filesep, name.fign, '_', name.var, '_', name.obs, '.tif'];
            tmp.title=[name.var, ', ', name.fign];
        end
        if ~exist(cfg.figname)
            
            if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
            
            fname.(name.fign)=[dirs.ens, tmp.fs, name.var, '_', name.obs, '_y_det_', name.fign, '.nc'];
            grid.tlong=ncread(fname.grid, 'TLONG');
            grid.tlat=ncread(fname.grid, 'TLAT');
            if (strcmp(cfg.dim_var, '4d'))
                tmp.C=ncread(fname.(name.fign), name.var, [1 1 tmp.lev 1], [inf inf 1 1]);
            elseif (strcmp(cfg.dim_var, '3d'))
                tmp.C=ncread(fname.(name.fign), name.var, [1 1 1], [inf inf 1]);
            end
            
            tmp.X=grid.tlong([end, 1:end],:);
            tmp.Y=grid.tlat([end, 1:end],:);
            tmp.C=tmp.C([end, 1:end],:);
            
            tmp.C(tmp.C<1)=1;
            tmp.maxC=max(max(tmp.C));
            tmp.prc95= prctile(tmp.C(:), 95);
            fig_cfg=fig_param_set('yr', [1 tmp.prc95], tmp.dropboxpath);  % cmap, clim, droppath
             fig_cfg.c_map(1,:)=[0.8 0.8 0.8]; % 1< : grayed
            fig_cfg.fig_name=tmp.title;
            fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
                'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
            %% map setting
            ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
                'fontname','freeserif'); 
            axis off; 
            hold on;
            setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
            set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
            text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
            'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
            'fontsize',14,'fontname','freeserif','interpreter','none')
        
            %% caxis & colorbar
            if diff(fig_cfg.c_lim)~=0, caxis(ax_m, fig_cfg.c_lim); end 
            colormap(fig_cfg.c_map);
            cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
            set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
            % title(cb,'R','fontsize',12);
        
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
            end
        
            %% save
            print(fig_h, cfg.figname, '-r200','-dpng');
            RemoveWhiteSpace([], 'file', cfg.figname);
    %         fprintf('%7.1f sec\n', toc(lap_time) );
            close all;
        end
        

%% ens_mean_magnitude map  -----------------------------------------------
        name.fign = 'ENS_mean';
        cfg.dim_var=Func_0026_CESM2_dim_var(name.var);
        dirs.figdir= [dirs.figroot, filesep, 'ens_set', filesep, name.obs];
        if (strcmp(cfg.dim_var, '4d'))
            cfg.figname=[dirs.figdir, filesep, name.fign, '_', name.var, '_', name.obs, '_lev', str.lev, '.tif'];
            tmp.title=[name.var, ', ', name.fign, ', ', str.lev];
        elseif (strcmp(cfg.dim_var, '3d'))
            cfg.figname=[dirs.figdir, filesep, name.fign, '_', name.var, '_', name.obs, '.tif'];
            tmp.title=[name.var, ', ', name.fign];
        end
        if ~exist(cfg.figname)
            if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
            
            fname.(name.fign)=[dirs.ens, tmp.fs, name.var, '_', name.obs, '_y_det_', 'time_mean_ensmean', '.nc'];
            grid.tlong=ncread(fname.grid, 'TLONG');
            grid.tlat=ncread(fname.grid, 'TLAT');
            if (strcmp(cfg.dim_var, '4d'))
                tmp.C=ncread(fname.(name.fign), name.var, [1 1 tmp.lev 1], [inf inf 1 1]);
                if strcmp(name.var(end-1:end), 'hl')
                    tmp.C(tmp.C==-1)=0;
                end
            elseif (strcmp(cfg.dim_var, '3d'))
                tmp.C=ncread(fname.(name.fign), name.var, [1 1 1], [inf inf 1]);
            end      

            tmp.X=grid.tlong([end, 1:end],:);
            tmp.Y=grid.tlat([end, 1:end],:);
            tmp.C=tmp.C([end, 1:end],:);
            
%             tmp.C(tmp.C<1)=1;
            tmp.maxC=max(max(tmp.C));
            tmp.prc01= prctile(tmp.C(:), 1);
%             tmp.prc01= prctile(tmp.C(:), 1);
            tmp.prc99= prctile(tmp.C(:), 99);
            fig_cfg=fig_param_set('jet', [tmp.prc01 tmp.prc99], tmp.dropboxpath);  % cmap, clim, droppath
%              fig_cfg.c_map(1,:)=[0.8 0.8 0.8]; % 1< : grayed
            fig_cfg.fig_name=tmp.title;
            fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
                'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
            %% map setting
            ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
                'fontname','freeserif'); 
            axis off; 
            hold on;
            setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
            set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
            text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
            'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
            'fontsize',14,'fontname','freeserif','interpreter','none')
        
            %% caxis & colorbar
            if diff(fig_cfg.c_lim)~=0, caxis(ax_m, fig_cfg.c_lim); end 

            colormap(fig_cfg.c_map);
            cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
            set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
%             % title(cb,'R','fontsize',12);
        
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
            end
        
            %% save
            print(fig_h, cfg.figname, '-r200','-dpng');
            RemoveWhiteSpace([], 'file', cfg.figname);
    %         fprintf('%7.1f sec\n', toc(lap_time) );
            close all;
        end




%% ens_spread map  -----------------------------------------------
        name.fign = 'ENS_spread';
        cfg.dim_var=Func_0026_CESM2_dim_var(name.var);
        dirs.figdir= [dirs.figroot, filesep, 'ens_set', filesep, name.obs];
        if (strcmp(cfg.dim_var, '4d'))
            cfg.figname=[dirs.figdir, filesep, name.fign, '_', name.var, '_', name.obs, '_lev', str.lev, '.tif'];
            tmp.title=[name.var, ', ', name.fign, ', ', str.lev];
        elseif (strcmp(cfg.dim_var, '3d'))
            cfg.figname=[dirs.figdir, filesep, name.fign, '_', name.var, '_', name.obs, '.tif'];
            tmp.title=[name.var, ', ', name.fign];
        end
        if ~exist(cfg.figname)
            if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
            
            fname.(name.fign)=[dirs.ens, tmp.fs, name.var, '_', name.obs, '_y_det_', 'time_mean_spread', '.nc'];
            grid.tlong=ncread(fname.grid, 'TLONG');
            grid.tlat=ncread(fname.grid, 'TLAT');
            if (strcmp(cfg.dim_var, '4d'))
                tmp.C=ncread(fname.(name.fign), name.var, [1 1 tmp.lev 1], [inf inf 1 1]);
            elseif (strcmp(cfg.dim_var, '3d'))
                tmp.C=ncread(fname.(name.fign), name.var, [1 1 1], [inf inf 1]);
            end
            
            tmp.X=grid.tlong([end, 1:end],:);
            tmp.Y=grid.tlat([end, 1:end],:);
            tmp.C=tmp.C([end, 1:end],:);
            
%             tmp.C(tmp.C<1)=1;
            tmp.maxC=max(max(tmp.C));
            tmp.prc01= prctile(tmp.C(:), 1);
            tmp.prc99= prctile(tmp.C(:), 99);
            fig_cfg=fig_param_set('jet', [tmp.prc01 tmp.prc99], tmp.dropboxpath);  % cmap, clim, droppath
%              fig_cfg.c_map(1,:)=[0.8 0.8 0.8]; % 1< : grayed
            fig_cfg.fig_name=tmp.title;
            fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
                'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
            %% map setting
            ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
                'fontname','freeserif'); 
            axis off; 
            hold on;
            setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
            set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
            text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
            'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
            'fontsize',14,'fontname','freeserif','interpreter','none')
        
            %% caxis & colorbar
            if diff(fig_cfg.c_lim)~=0, caxis(ax_m, fig_cfg.c_lim); end 
            colormap(fig_cfg.c_map);
            cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
            set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
            % title(cb,'R','fontsize',12);
        
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
            end
        
            %% save
            print(fig_h, cfg.figname, '-r200','-dpng');
            RemoveWhiteSpace([], 'file', cfg.figname);
    %         fprintf('%7.1f sec\n', toc(lap_time) );
            close all;
        end




%% ens_amplitude map  -----------------------------------------------
        name.fign = 'ENS_variability';
        cfg.dim_var=Func_0026_CESM2_dim_var(name.var);
        dirs.figdir= [dirs.figroot, filesep, 'ens_set', filesep, name.obs];
        if (strcmp(cfg.dim_var, '4d'))
            cfg.figname=[dirs.figdir, filesep, name.fign, '_', name.var, '_', name.obs, '_lev', str.lev, '.tif'];
            tmp.title=[name.var, ', ', name.fign, ', ', str.lev];
        elseif (strcmp(cfg.dim_var, '3d'))
            cfg.figname=[dirs.figdir, filesep, name.fign, '_', name.var, '_', name.obs, '.tif'];
            tmp.title=[name.var, ', ', name.fign];
        end
        if ~exist(cfg.figname)
            if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
            
            fname.(name.fign)=[dirs.ens, tmp.fs, name.var, '_', name.obs, '_y_det_', 'time_std_ensmean', '.nc'];
            grid.tlong=ncread(fname.grid, 'TLONG');
            grid.tlat=ncread(fname.grid, 'TLAT');
            if (strcmp(cfg.dim_var, '4d'))
                tmp.C=ncread(fname.(name.fign), name.var, [1 1 tmp.lev 1], [inf inf 1 1]);
            elseif (strcmp(cfg.dim_var, '3d'))
                tmp.C=ncread(fname.(name.fign), name.var, [1 1 1], [inf inf 1]);
            end
            
            tmp.X=grid.tlong([end, 1:end],:);
            tmp.Y=grid.tlat([end, 1:end],:);
            tmp.C=tmp.C([end, 1:end],:);
            
%             tmp.C(tmp.C<1)=1;
            tmp.maxC=max(max(tmp.C));
            tmp.prc01= prctile(tmp.C(:), 1);
            tmp.prc99= prctile(tmp.C(:), 99);
            fig_cfg=fig_param_set('jet', [tmp.prc01 tmp.prc99], tmp.dropboxpath);  % cmap, clim, droppath
%              fig_cfg.c_map(1,:)=[0.8 0.8 0.8]; % 1< : grayed
            fig_cfg.fig_name=tmp.title;
            fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
                'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
            %% map setting
            ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
                'fontname','freeserif'); 
            axis off; 
            hold on;
            setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
            set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
            text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
            'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
            'fontsize',14,'fontname','freeserif','interpreter','none')
        
            %% caxis & colorbar
            if diff(fig_cfg.c_lim)~=0, caxis(ax_m, fig_cfg.c_lim); end 
            colormap(fig_cfg.c_map);
            cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
            set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
            % title(cb,'R','fontsize',12);
        
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
            end
        
            %% save
            print(fig_h, cfg.figname, '-r200','-dpng');
            RemoveWhiteSpace([], 'file', cfg.figname);
    %         fprintf('%7.1f sec\n', toc(lap_time) );
            close all;
        end



%% Correlation set  --------------------------------------
%% correlation (mean) map  --------------------------------------
%ttest  p= 1-tcdf(r*sqrt((mu)/(1-r^2)), mu)

        name.fign = 'corr_mean';
        cfg.dim_var=Func_0026_CESM2_dim_var(name.var);
        dirs.figdir= [dirs.figroot, filesep, 'corr_set', filesep, name.obs];
        if (strcmp(cfg.dim_var, '4d'))
            cfg.figname=[dirs.figdir, filesep, name.fign, '_', name.var, '_', name.obs, '_lev', str.lev, '.tif'];
            tmp.title=[name.var, ', ', name.fign, ', ', str.lev];
        elseif (strcmp(cfg.dim_var, '3d'))
            cfg.figname=[dirs.figdir, filesep, name.fign, '_', name.var, '_', name.obs, '.tif'];
            tmp.title=[name.var, ', ', name.fign];
        end

        if ~exist(cfg.figname)
            if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
            
            fname.(name.fign)=[dirs.corr, tmp.fs, name.var, '_', name.obs, '_y_det_', 'corr_ensmean', '.nc'];
            grid.tlong=ncread(fname.grid, 'TLONG');
            grid.tlat=ncread(fname.grid, 'TLAT');
            if (strcmp(cfg.dim_var, '4d'))
                tmp.C=ncread(fname.(name.fign), name.var, [1 1 tmp.lev 1], [inf inf 1 1]);
            elseif (strcmp(cfg.dim_var, '3d'))
                tmp.C=ncread(fname.(name.fign), name.var, [1 1 1], [inf inf 1]);
            end
            
            tmp.X=grid.tlong([end, 1:end],:);
            tmp.Y=grid.tlat([end, 1:end],:);
            tmp.C=tmp.C([end, 1:end],:);
            
%             tmp.C(tmp.C<1)=1;
            tmp.maxC=max(max(tmp.C));
            tmp.prc95= prctile(tmp.C(:), 95);
            fig_cfg=fig_param_set('byr', [-1 1], tmp.dropboxpath);  % cmap, clim, droppath
%              fig_cfg.c_map(1,:)=[0.8 0.8 0.8]; % 1< : grayed
            fig_cfg.fig_name=tmp.title;
            fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
                'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
            %% map setting
            ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
                'fontname','freeserif'); 
            axis off; 
            hold on;
            setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
            set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
            text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
            'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
            'fontsize',14,'fontname','freeserif','interpreter','none')
        
            %% caxis & colorbar
            if diff(fig_cfg.c_lim)~=0, caxis(ax_m, fig_cfg.c_lim); end 
            colormap(fig_cfg.c_map);
            cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
            set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
            % title(cb,'R','fontsize',12);
        
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
            end
        
            %% save
            print(fig_h, cfg.figname, '-r200','-dpng');
            RemoveWhiteSpace([], 'file', cfg.figname);
            close all;
        end



%% correlation (max) map  --------------------------------------

        name.fign = 'corr_max';
        cfg.dim_var=Func_0026_CESM2_dim_var(name.var);
        dirs.figdir= [dirs.figroot, filesep, 'corr_set', filesep, name.obs];
        if (strcmp(cfg.dim_var, '4d'))
            cfg.figname=[dirs.figdir, filesep, name.fign, '_', name.var, '_', name.obs, '_lev', str.lev, '.tif'];
            tmp.title=[name.var, ', ', name.fign, ', ', str.lev];
        elseif (strcmp(cfg.dim_var, '3d'))
            cfg.figname=[dirs.figdir, filesep, name.fign, '_', name.var, '_', name.obs, '.tif'];
            tmp.title=[name.var, ', ', name.fign];
        end

        if ~exist(cfg.figname)
            if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
            
            fname.(name.fign)=[dirs.corr, tmp.fs, name.var, '_', name.obs, '_y_det_', 'corr_ensmax', '.nc'];
            grid.tlong=ncread(fname.grid, 'TLONG');
            grid.tlat=ncread(fname.grid, 'TLAT');
            if (strcmp(cfg.dim_var, '4d'))
                tmp.C=ncread(fname.(name.fign), name.var, [1 1 tmp.lev 1], [inf inf 1 1]);
            elseif (strcmp(cfg.dim_var, '3d'))
                tmp.C=ncread(fname.(name.fign), name.var, [1 1 1], [inf inf 1]);
            end
            
            tmp.X=grid.tlong([end, 1:end],:);
            tmp.Y=grid.tlat([end, 1:end],:);
            tmp.C=tmp.C([end, 1:end],:);
            
%             tmp.C(tmp.C<1)=1;
            tmp.maxC=max(max(tmp.C));
            tmp.prc95= prctile(tmp.C(:), 95);
            fig_cfg=fig_param_set('byr', [-1 1], tmp.dropboxpath);  % cmap, clim, droppath
%              fig_cfg.c_map(1,:)=[0.8 0.8 0.8]; % 1< : grayed
            fig_cfg.fig_name=tmp.title;
            fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
                'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
            %% map setting
            ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
                'fontname','freeserif'); 
            axis off; 
            hold on;
            setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
            set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
            text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
            'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
            'fontsize',14,'fontname','freeserif','interpreter','none')
        
            %% caxis & colorbar
            if diff(fig_cfg.c_lim)~=0, caxis(ax_m, fig_cfg.c_lim); end 
            colormap(fig_cfg.c_map);
            cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
            set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
            % title(cb,'R','fontsize',12);
        
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
            end
        
            %% save
            print(fig_h, cfg.figname, '-r200','-dpng');
            RemoveWhiteSpace([], 'file', cfg.figname);
            close all;
        end



%% correlation (min) map  --------------------------------------

        name.fign = 'corr_min';
        cfg.dim_var=Func_0026_CESM2_dim_var(name.var);
        dirs.figdir= [dirs.figroot, filesep, 'corr_set', filesep, name.obs];
        if (strcmp(cfg.dim_var, '4d'))
            cfg.figname=[dirs.figdir, filesep, name.fign, '_', name.var, '_', name.obs, '_lev', str.lev, '.tif'];
            tmp.title=[name.var, ', ', name.fign, ', ', str.lev];
        elseif (strcmp(cfg.dim_var, '3d'))
            cfg.figname=[dirs.figdir, filesep, name.fign, '_', name.var, '_', name.obs, '.tif'];
            tmp.title=[name.var, ', ', name.fign];
        end

        if ~exist(cfg.figname)
            if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
            
            fname.(name.fign)=[dirs.corr, tmp.fs, name.var, '_', name.obs, '_y_det_', 'corr_ensmin', '.nc'];
            grid.tlong=ncread(fname.grid, 'TLONG');
            grid.tlat=ncread(fname.grid, 'TLAT');
            if (strcmp(cfg.dim_var, '4d'))
                tmp.C=ncread(fname.(name.fign), name.var, [1 1 tmp.lev 1], [inf inf 1 1]);
            elseif (strcmp(cfg.dim_var, '3d'))
                tmp.C=ncread(fname.(name.fign), name.var, [1 1 1], [inf inf 1]);
            end
            
            tmp.X=grid.tlong([end, 1:end],:);
            tmp.Y=grid.tlat([end, 1:end],:);
            tmp.C=tmp.C([end, 1:end],:);
            
%             tmp.C(tmp.C<1)=1;
            tmp.maxC=max(max(tmp.C));
            tmp.prc95= prctile(tmp.C(:), 95);
            fig_cfg=fig_param_set('byr', [-1 1], tmp.dropboxpath);  % cmap, clim, droppath
%              fig_cfg.c_map(1,:)=[0.8 0.8 0.8]; % 1< : grayed
            fig_cfg.fig_name=tmp.title;
            fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
                'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
            %% map setting
            ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
                'fontname','freeserif'); 
            axis off; 
            hold on;
            setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
            set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
            text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
            'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
            'fontsize',14,'fontname','freeserif','interpreter','none')
        
            %% caxis & colorbar
            if diff(fig_cfg.c_lim)~=0, caxis(ax_m, fig_cfg.c_lim); end 
            colormap(fig_cfg.c_map);
            cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
            set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
            % title(cb,'R','fontsize',12);
        
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
            end
        
            %% save
            print(fig_h, cfg.figname, '-r200','-dpng');
            RemoveWhiteSpace([], 'file', cfg.figname);
            close all;
        end



%% correlation (std) map  --------------------------------------

        name.fign = 'corr_std';
        cfg.dim_var=Func_0026_CESM2_dim_var(name.var);
        dirs.figdir= [dirs.figroot, filesep, 'corr_set', filesep, name.obs];
        if (strcmp(cfg.dim_var, '4d'))
            cfg.figname=[dirs.figdir, filesep, name.fign, '_', name.var, '_', name.obs, '_lev', str.lev, '.tif'];
            tmp.title=[name.var, ', ', name.fign, ', ', str.lev];
        elseif (strcmp(cfg.dim_var, '3d'))
            cfg.figname=[dirs.figdir, filesep, name.fign, '_', name.var, '_', name.obs, '.tif'];
            tmp.title=[name.var, ', ', name.fign];
        end

        if ~exist(cfg.figname)
            if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
            
            fname.(name.fign)=[dirs.corr, tmp.fs, name.var, '_', name.obs, '_y_det_', 'corr_ensstd', '.nc'];
            grid.tlong=ncread(fname.grid, 'TLONG');
            grid.tlat=ncread(fname.grid, 'TLAT');
            if (strcmp(cfg.dim_var, '4d'))
                tmp.C=ncread(fname.(name.fign), name.var, [1 1 tmp.lev 1], [inf inf 1 1]);
            elseif (strcmp(cfg.dim_var, '3d'))
                tmp.C=ncread(fname.(name.fign), name.var, [1 1 1], [inf inf 1]);
            end
            
            tmp.X=grid.tlong([end, 1:end],:);
            tmp.Y=grid.tlat([end, 1:end],:);
            tmp.C=tmp.C([end, 1:end],:);
            
%             tmp.C(tmp.C<1)=1;
            tmp.maxC=max(max(tmp.C));
            tmp.prc01= prctile(tmp.C(:), 1);
            tmp.prc99= prctile(tmp.C(:), 99);
            fig_cfg=fig_param_set('jet', [tmp.prc01 tmp.prc99], tmp.dropboxpath);  % cmap, clim, droppath
%              fig_cfg.c_map(1,:)=[0.8 0.8 0.8]; % 1< : grayed
            fig_cfg.fig_name=tmp.title;
            fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
                'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
            %% map setting
            ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
                'fontname','freeserif'); 
            axis off; 
            hold on;
            setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
            set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
            text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
            'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
            'fontsize',14,'fontname','freeserif','interpreter','none')
        
            %% caxis & colorbar
            if diff(fig_cfg.c_lim)~=0, caxis(ax_m, fig_cfg.c_lim); end 
            colormap(fig_cfg.c_map);
            cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
            set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
            % title(cb,'R','fontsize',12);
        
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
            end
        
            %% save
            print(fig_h, cfg.figname, '-r200','-dpng');
            RemoveWhiteSpace([], 'file', cfg.figname);
            close all;
        end


    end
end
end

%-------------------------------------------------------- function

function gmsst = f_gm_var(var_2d, area)
    mask_var=NaN(size(var_2d));
    mask_var(isfinite(var_2d))=1;
    area=area.*mask_var;
    var_2d=squeeze(var_2d);
    var_2d_sum=var_2d.*area;
    var_2d_sum=sum(var_2d_sum(:), 'omitnan');
    gmsst=var_2d_sum ./ sum(area(:), 'omitnan');
end

function fig_cfg = fig_param_set(cmap, c_lim, droppath)
    % cmap : 'byr', 'jet', ...
    fig_cfg.name_rgn = 'Glob';
    fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
    fig_cfg.x_lim = [-180 180];
    fig_cfg.y_lim = [-80 89];
    fig_cfg.fig_size = [0,0,6,3.5];
    fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
    fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
    fig_cfg.title_pos = [0.5,0.93];
    fig_cfg.p_lim =0.1;
    fig_cfg.c_lim = c_lim;
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps(cmap, droppath);
end

function component = f_cmpname_var(var_model)
    switch var_model
        case {'TEMP', 'SALT', 'SSH', 'UVEL', 'VVEL', 'diatChl', 'diazChl', 'spChl'...
               'dust_FLUX_IN', 'dust_REMIN', 'dustToSed', 'pH_3D', 'Fe', 'O2', ...
               'P_iron_FLUX_100m', 'POC_FLUX_100m', 'POP_FLUX_100m', ...
               'SiO2_FLUX_100m', 'TMXL_DR', 'XMXL_DR', ...
               'photoC_TOT', 'photoC_TOT_zint', 'photoC_TOT_zint_100m'}
            component='ocn'; %t.component
        case {'AODDUST', 'AODVIS', 'dst_a1SF', 'dst_a2SF', 'dst_a3SF', ...
                'PRECT', 'PSL', 'SST', 'dst_a1', 'dst_a2', 'dst_a3', ...
                'U', 'V', 'ATM_COARSE_DUST_FLUX_CPL', 'ATM_FINE_DUST_FLUX_CPL', 'SEAICE_DUST_FLUX_CPL'}
            component='atm';
        case {'DSTFLXT', 'DSTDEP', 'DSL', 'QSOIL', 'TWS'} %TOTVEGC, FIRE
            component='lnd';
        case {'TOTAL_DISCHARGE_TO_OCEAN_LIQ', ...
                'TOTAL_DISCHARGE_TO_OCEAN_ICE'}
            component='rof';
    end
end

