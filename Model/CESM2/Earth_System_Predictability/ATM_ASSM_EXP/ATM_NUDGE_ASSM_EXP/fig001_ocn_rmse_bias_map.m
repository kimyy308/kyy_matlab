% %  Created 27-Feb-2023 by Yong-Yub Kim
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
% cfg.levs=[32];

%% model configuration

%% ocn
% cfg.vars={'WVEL', 'diatChl', 'diazChl', 'spChl', 'DIC', 'ALK', 'TEMP', 'SALT', ...
%     'NO3', 'SiO3', 'PO4', 'UVEL', 'VVEL', 'PD', 'UE_PO4', 'VN_PO4', 'WT_PO4', ...
%     'Fe', 'diatC', 'diazC', 'spC', 'BSF', 'HMXL', 'photoC_TOT_zint_100m', ...
%     'photoC_TOT_zint', 'SSH', 'IRON_FLUX', 'HBLT', 'HMXL_DR', 'TBLT', 'TMXL', 'XBLT', 'XMXL', ...
%     'TMXL_DR', 'XMXL_DR', 'diatC_zint_100m_2', 'diazC_zint_100m_2', 'spC_zint_100m_2' };

%% ocn_bgc
% cfg.vars={'diat_agg', 'diatFe', 'diatP', 'diatSi', 'diaz_agg', 'diaz_Nfix', 'diazFe', 'diazP', ...
%     'dust_FLUX_IN', 'dust_REMIN', 'O2' ...
%     'photoFe_diat', 'photoFe_diaz', 'photoFe_sp', 'photoNO3_diat', 'photoNO3_diaz', ...
%     'photoNO3_sp', 'PO4_diat_uptake', 'PO4_diaz_uptake', 'PO4_sp_uptake', 'sp_agg', ...
%     'spFe', 'spP', 'zooC', 'diat_agg_zint_100m', 'diat_Fe_lim_Cweight_avg_100m', ...
%     'diat_Fe_lim_surf', 'diat_light_lim_Cweight_avg_100m', 'diat_light_lim_surf', ...
%     'diat_loss_zint_100m', 'diat_N_lim_Cweight_avg_100m', 'diat_N_lim_surf', ...
%     'diat_P_lim_Cweight_avg_100m', 'diat_P_lim_surf', 'diat_SiO3_lim_Cweight_avg_100m', ...
%     'diat_SiO3_lim_surf', 'diaz_agg_zint_100m', 'diaz_Fe_lim_Cweight_avg_100m', ...
%     'diaz_Fe_lim_surf', 'diaz_light_lim_Cweight_avg_100m', 'diaz_light_lim_surf', ...
%     'diaz_loss_zint', 'diaz_P_lim_Cweight_avg_100m', 'diaz_P_lim_surf', 'dustToSed', ...
%     'graze_diat_zint_100m', 'graze_diat_zoo_zint_100m', 'graze_diaz_zint_100m', ...
%     'graze_diaz_zoo_zint_100m', 'graze_sp_zint_100m', 'graze_sp_zoo_zint_100m', ...
%     'LWUP_F', 'O2_ZMIN_DEPTH', 'photoC_diat_zint_100m', 'photoC_diaz_zint_100m', ...
%     'photoC_NO3_TOT', 'photoC_NO3_TOT_zint_100m', 'photoC_sp_zint_100m', ...
%     'photoC_TOT', 'POC_FLUX_100m', 'sp_agg_zint_100m', 'sp_Fe_lim_Cweight_avg_100m', 'sp_Fe_lim_surf', ...
%     'sp_light_lim_Cweight_avg_100m', 'sp_light_lim_surf', 'sp_loss_zint_100m', ...
%     'sp_N_lim_Cweight_avg_100m', 'sp_N_lim_surf', 'sp_P_lim_Cweight_avg_100m', ...
%     'sp_P_lim_surf', 'zoo_loss_zint_100m'};
cfg.vars = {'photoC_TOT_zint_100m'};

%% atm
% cfg.vars={'AODDUST', 'CLOUD', 'RELHUM', 'SFdst_a2', 'SST', 'dst_a2SF', 'dst_a2_SRF', ...
%     'FSDS', 'FSNS', 'TAUX', 'TAUY', ...
%     'IVT', 'PRECL', 'PS', 'Q','U', 'UBOT', 'V850', 'uIVT', ...
%     'PRECC', 'PRECT', 'PSL', 'T', 'U850', 'V', 'VBOT', 'vIVT', 'U10'};


% cfg.vars={'WVEL'};
% abc=dir('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/ASSM_EXP/archive_yearly_analysis/ocn')
% cfg.vars=abc(3:end).name;

% cfg.obs={'f09_g17.assm.projdv7.3_ba-20p1'};
% cfg.exps={'ATM_NUDGE_ASSM_EXP', 'ATM_NUDGE2_ASSM_EXP','ATM_ASSM_EXP'};

cfg.exps={'ATM_NOASSM_ASSM_EXP'};


% for obsi= 1:length(cfg.obs)
%     name.obs=cfg.obs{obsi};
for expi= 1:length(cfg.exps)
    name.exp=cfg.exps{expi};
    name.exp_ref='ATM_ASSM_EXP';

    name.obs_ref = 'f09_g17.assm.projdv7.3_ba-20p1';
switch name.exp
    case {'ATM_NUDGE_ASSM_EXP', 'ATM_NUDGE2_ASSM_EXP','ATM_ASSM_EXP'}
        name.obs = 'f09_g17.assm.projdv7.3_ba-20p1';
    case 'ATM_NOASSM_ASSM_EXP'
        name.obs = 'LE2-1231.011';
end
for vari= 1:length(cfg.vars)
    name.var=cfg.vars{vari};
    name.comp= Func_0025_CESM2_cmpname_var(name.var);
    [name.freq, name.outf]=func_var_freq(name.var);
    dirs.assmroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/ATM_TEST/preliminary/', ...
        name.exp, tmp.fs, 'archive_analysis/', name.comp, tmp.fs, name.freq, tmp.fs, name.var, tmp.fs, name.obs];
    dirs.assm_refroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/ATM_TEST/preliminary/', ...
        name.exp_ref, tmp.fs, 'archive_analysis/', name.comp, tmp.fs, name.freq, tmp.fs, name.var, tmp.fs, name.obs_ref];
    dirs.figroot = strrep(dirs.assmroot, 'kimyy/Model', 'kimyy/Figure');
%     dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/ATM_ASSM_EXP/', ...
%         name.exp, tmp.fs, 'archive_analysis/', name.comp, tmp.fs, name.var];
    dirs.grid='/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/ASSM_EXP/archive';
    
    dirs.rmse=[dirs.assmroot, tmp.fs, 'RMSE'];
    dirs.rmse_ref=[dirs.assm_refroot, tmp.fs, 'RMSE'];
    
    fname.grid=[dirs.grid, tmp.fs, 'pop.h.once.nc'];    
    S = shaperead('landareas.shp');

    for levi=1:length(cfg.levs)
        tmp.lev=cfg.levs(levi);
        str.lev=num2str(tmp.lev, '%02i');
        
%% RMSE map  --------------------------------------
        %% name setting
        name.fign = 'RMSE';
        cfg.dim_var=Func_0026_CESM2_dim_var(name.var);
        dirs.figdir= [dirs.figroot, filesep, 'RMSE'];
        if (strcmp(cfg.dim_var, '4d'))
            cfg.figname=[dirs.figdir, filesep, name.fign, '_', name.var, '_', name.obs, '_lev', str.lev, '.tif'];
            tmp.title=[name.var, ', ', name.fign, ', ', str.lev];
        elseif (strcmp(cfg.dim_var, '3d'))
            cfg.figname=[dirs.figdir, filesep, name.fign, '_', name.var, '_', name.obs, '.tif'];
            tmp.title=[name.var, ', ', name.fign];
        end
        %% figure
        if ~exist(cfg.figname)
            if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end

            fname.(name.fign)=[dirs.rmse, tmp.fs, ...
                'RMSE_', name.var, '_', name.obs, '.', name.outf, '.nc'];
            fname_ref.(name.fign)=[dirs.rmse_ref, tmp.fs, ...
                'RMSE_', name.var, '_', name.obs_ref, '.', name.outf, '.nc'];
            [grid.tlong, grid.tlat] = func_get_grid(fname.grid, fname.(name.fign), name.comp);
            
            if (strcmp(cfg.dim_var, '4d'))
                tmp.C=ncread(fname.(name.fign), name.var, [1 1 tmp.lev 1], [inf inf 1 1]);
                tmp.C_ref=ncread(fname_ref.(name.fign), name.var, [1 1 tmp.lev 1], [inf inf 1 1]);
            elseif (strcmp(cfg.dim_var, '3d'))
                tmp.C=ncread(fname.(name.fign), name.var, [1 1 1], [inf inf 1]);
                tmp.C_ref=ncread(fname_ref.(name.fign), name.var, [1 1 1], [inf inf 1]);
            end
            
            tmp.X=grid.tlong([end, 1:end],:);
            tmp.Y=grid.tlat([end, 1:end],:);
            tmp.C=tmp.C([end, 1:end],:);
            
%             tmp.C(tmp.C<1)=1;
            tmp.maxC=max(max(tmp.C));
            tmp.minC=min(min(tmp.C));
            tmp.maxC_ref=max(max(tmp.C_ref));
            tmp.minC_ref=min(min(tmp.C_ref));
            tmp.prc05= prctile(tmp.C(:), 5);
            tmp.prc95= prctile(tmp.C(:), 95);
            tmp.prc05_ref= prctile(tmp.C_ref(:), 5);
            tmp.prc95_ref= prctile(tmp.C_ref(:), 95);
%             fig_cfg=fig_param_set('yr', [tmp.minC tmp.maxC], tmp.dropboxpath);  % cmap, clim, droppath
%             fig_cfg=fig_param_set('yr', [0e-4 2e-4], tmp.dropboxpath);  % cmap, clim, droppath
            fig_cfg=fig_param_set('jet', [tmp.prc05_ref tmp.prc95_ref], tmp.dropboxpath);  % cmap, clim, droppath

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


%% Bias map  --------------------------------------
        %% name setting
        name.fign = 'BIAS';
        cfg.dim_var=Func_0026_CESM2_dim_var(name.var);
        dirs.figdir= [dirs.figroot, filesep, 'BIAS'];
        if (strcmp(cfg.dim_var, '4d'))
            cfg.figname=[dirs.figdir, filesep, name.fign, '_', name.var, '_', name.obs, '_lev', str.lev, '.tif'];
            tmp.title=[name.var, ', ', name.fign, ', ', str.lev];
        elseif (strcmp(cfg.dim_var, '3d'))
            cfg.figname=[dirs.figdir, filesep, name.fign, '_', name.var, '_', name.obs, '.tif'];
            tmp.title=[name.var, ', ', name.fign];
        end
        %% figure
        if ~exist(cfg.figname)
            if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end

            fname.(name.fign)=[dirs.rmse, tmp.fs, ...
                'mean_bias', '_', name.var, '_', name.obs, '.', name.outf, '.nc'];
            fname_ref.(name.fign)=[dirs.rmse_ref, tmp.fs, ...
                'mean_bias', '_', name.var, '_', name.obs_ref, '.', name.outf, '.nc'];
            [grid.tlong, grid.tlat] = func_get_grid(fname.grid, fname.(name.fign), name.comp);
            
            if (strcmp(cfg.dim_var, '4d'))
                tmp.C=ncread(fname.(name.fign), name.var, [1 1 tmp.lev 1], [inf inf 1 1]);
                tmp.C_ref=ncread(fname_ref.(name.fign), name.var, [1 1 tmp.lev 1], [inf inf 1 1]);
            elseif (strcmp(cfg.dim_var, '3d'))
                tmp.C=ncread(fname.(name.fign), name.var, [1 1 1], [inf inf 1]);
                tmp.C_ref=ncread(fname_ref.(name.fign), name.var, [1 1 1], [inf inf 1]);
            end
            
            tmp.X=grid.tlong([end, 1:end],:);
            tmp.Y=grid.tlat([end, 1:end],:);
            tmp.C=tmp.C([end, 1:end],:);
            
%             tmp.C(tmp.C<1)=1;
            tmp.maxC=max(max(tmp.C));
            tmp.minC=min(min(tmp.C));
            tmp.maxC_ref=max(max(tmp.C_ref));
            tmp.minC_ref=min(min(tmp.C_ref));
            tmp.prc05= prctile(tmp.C(:), 5);
            tmp.prc95= prctile(tmp.C(:), 95);
            tmp.prc05_ref= prctile(tmp.C_ref(:), 5);
            tmp.prc95_ref= prctile(tmp.C_ref(:), 95);
            tmp.abs_prc95_ref= prctile(abs(tmp.C_ref(:)), 95);

            fig_cfg=fig_param_set('byr', [-tmp.abs_prc95_ref tmp.abs_prc95_ref], tmp.dropboxpath);  % cmap, clim, droppath

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








    end







end
end
% end


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



function [freq, outf]= func_var_freq(var_model)

switch var_model
    case{'U', 'V', 'IVT', 'PRECL', 'PRECC', 'PRECT', 'PS', 'PSL', 'Q', 'T', 'U850', 'UBOT', 'V850', ...
            'VBOT', 'uIVT', 'vIVT'}
        freq='6hr';
        outf='cam.h2';
    case {'FSDS', 'FSNS', 'TAUX', 'TAUY' }
        freq='daily'; 
        outf='cam.h1';
    case { 'CLOUD','OMEGA', 'RELHUM', ...
           'AODDUST', 'AODVIS', 'dst_a1SF', 'dst_a2SF', 'dst_a3SF', ...
            'SST', 'DSTFLXT', 'DSTDEP', 'DSL', 'QSOIL', 'TWS',...
           'TOTAL_DISCHARGE_TO_OCEAN_LIQ', ...
            'TOTAL_DISCHARGE_TO_OCEAN_ICE', ...
            'ATM_COARSE_DUST_FLUX_CPL', 'ATM_FINE_DUST_FLUX_CPL', 'SEAICE_DUST_FLUX_CPL', ...
             'dry_deposition_NOy_as_N', 'dry_deposition_NHx_as_N', ...
            'dst_a1_SRF', 'dst_a2DDF', 'dst_a2SFWET', 'dst_a2_SRF',  ...
            'dst_a3_SRF', 'SFdst_a1', 'SFdst_a2', 'SFdst_a3', ...
            'TS', 'U10', 'wet_deposition_NHx_as_N', 'wet_deposition_NOy_as_N'}
        freq='mon'; 
        outf='cam.h0';
    case{'TEMP', 'SALT', 'UVEL', 'VVEL', 'diatChl', 'diazChl', 'spChl', 'PD', ...
           'dustToSed', 'SSH', 'dust_FLUX_IN', 'dust_REMIN',  'pH_3D', 'Fe', 'O2', ...
           'dst_a1', 'dst_a2', 'dst_a3',  'ALK',  'DIC',  'NO3', 'PO4', 'SiO3', 'WVEL', ...
           'UISOP', 'VISOP', 'UE_PO4', 'VN_PO4', 'WT_PO4', ...
           'BSF', 'HMXL', 'photoC_TOT_zint', 'photoC_TOT_zint_100m', ...
           'P_iron_FLUX_100m', 'POC_FLUX_100m', 'POP_FLUX_100m', 'IRON_FLUX', ...
           'SiO2_FLUX_100m', 'TMXL_DR', 'XMXL_DR', ...
           'diatC_zint_100m_2', 'diazC_zint_100m_2', 'spC_zint_100m_2', ...
           'NOx_FLUX', 'HBLT', 'HMXL_DR', 'TBLT', 'TMXL', 'XBLT', 'XMXL', 'diat_agg', 'diatC', 'diazC', 'spC', ...
           'diatFe', 'diatP', 'diatSi', 'diaz_agg', 'diaz_Nfix', 'diazFe', 'diazP', ...
            'photoFe_diat', 'photoFe_diaz', 'photoFe_sp', 'photoNO3_diat', 'photoNO3_diaz', ...
            'photoNO3_sp', 'PO4_diat_uptake', 'PO4_diaz_uptake', 'PO4_sp_uptake', 'sp_agg', ...
            'spFe', 'spP', 'zooC', 'diat_agg_zint_100m', 'diat_Fe_lim_Cweight_avg_100m', ...
            'diat_Fe_lim_surf', 'diat_light_lim_Cweight_avg_100m', 'diat_light_lim_surf', ...
            'diat_loss_zint_100m', 'diat_N_lim_Cweight_avg_100m', 'diat_N_lim_surf', ...
            'diat_P_lim_Cweight_avg_100m', 'diat_P_lim_surf', 'diat_SiO3_lim_Cweight_avg_100m', ...
            'diat_SiO3_lim_surf', 'diaz_agg_zint_100m', 'diaz_Fe_lim_Cweight_avg_100m', ...
            'diaz_Fe_lim_surf', 'diaz_light_lim_Cweight_avg_100m', 'diaz_light_lim_surf', ...
            'diaz_loss_zint', 'diaz_P_lim_Cweight_avg_100m', 'diaz_P_lim_surf', ...
            'graze_diat_zint_100m', 'graze_diat_zoo_zint_100m', 'graze_diaz_zint_100m', ...
            'graze_diaz_zoo_zint_100m', 'graze_sp_zint_100m', 'graze_sp_zoo_zint_100m', ...
            'LWUP_F', 'O2_ZMIN_DEPTH', 'photoC_diat_zint_100m', 'photoC_diaz_zint_100m', ...
            'photoC_NO3_TOT', 'photoC_NO3_TOT_zint_100m', 'photoC_sp_zint_100m', ...
            'photoC_TOT', 'sp_agg_zint_100m', 'sp_Fe_lim_Cweight_avg_100m', 'sp_Fe_lim_surf', ...
            'sp_light_lim_Cweight_avg_100m', 'sp_light_lim_surf', 'sp_loss_zint_100m', ...
            'sp_N_lim_Cweight_avg_100m', 'sp_N_lim_surf', 'sp_P_lim_Cweight_avg_100m', ...
            'sp_P_lim_surf', 'zoo_loss_zint_100m'}
        freq='mon'; 
        outf='pop.h';

end

end

function [tlong, tlat]= func_get_grid(gridname, dataname, comp)
switch comp
    case 'ocn'
        tlong=ncread(gridname, 'TLONG');
        tlat=ncread(gridname, 'TLAT');
    case 'atm'
        lon=ncread(dataname, 'lon');
        lat=ncread(dataname, 'lat');
        [tlat, tlong] = meshgrid(lat, lon);
end
end