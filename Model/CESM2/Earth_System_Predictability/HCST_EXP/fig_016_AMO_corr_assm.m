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
cfg.ly_e = 1; %5

%% model configuration

cfg.vars={'DpCO2_ALT_CO2'};
cfg.vars={'DpCO2', 'DpCO2_ALT_CO2', 'FG_CO2', 'SSH'};

cfg.vars={'DIC_ALT_CO2'};
% cfg.vars={'DIC'};
cfg.vars={'SALT'};

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
    cfg.obs_iyears=1960:2020;
    cfg.obs_iyears2=f_obs_iyears(cfg.var);
    
    % dirs.hcstroot=['/mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive/', cfg.comp, '/', cfg.var];
    % dirs.obsroot=['/Volumes/kyy_raid/kimyy/Observation/', cfg.obs_name, '/monthly_reg_', cfg.obs_fname_module(2:4)];
    dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', cfg.var];
    % dirs.lens2root=['/mnt/lustre/proj/earth.system.predictability/LENS2/archive_analysis/', cfg.comp, '/', cfg.var];
    % dirs.assmroot=['/mnt/lustre/proj/earth.system.predictability/ASSM_EXP/archive_analysis/', cfg.comp, '/', cfg.var];
    
    dirs.hcstroot=['/mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive_yearly_transfer/', cfg.comp, '/', cfg.var];
    dirs.obsroot=['/mnt/lustre/proj/kimyy/Observation/', cfg.obs_name];
    dirs.matroot=['/mnt/lustre/proj/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp,'/', cfg.var];
    dirs.lens2root=['/mnt/lustre/proj/earth.system.predictability/LENS2/archive_yearly_transfer/', cfg.comp, '/', cfg.var];
    dirs.assmroot=['/mnt/lustre/proj/earth.system.predictability/ASSM_EXP/archive_yearly_transfer/', cfg.comp, '/', cfg.var];    
    
    dirs.corrroot='/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/corr';
   
    load([dirs.corrroot, '/', 'corr_all_', cfg.var, '_v', num2str(cfg.vlayer_1st), '_v', num2str(cfg.vlayer_cnt), '.mat'], 'cfg_assm', 'corrval', 'data_assm_em', 'grid', ...
        'clim_index_ym', 'clim_index_jfm', 'clim_index_djf');


% pcolor(corrval.obs_assm_em_ym.val'); shading flat; colorbar;

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
  

% %     %% for time series
% %     for fold=1
% %     sta_lonlat = {[125 160 30 40], [180 240 70 90]};
% %     
% %     for stai=1:length(sta_lonlat)
% %         grids=grid;
% %         if length(sta_lonlat{stai})==2
% %             xpoint=sta_lonlat{stai}(1);
% %             ypoint = sta_lonlat{stai}(2);
% %             
% %             [grids.id_w, grids.id_e, grids.id_s, grids.id_n] = Func_0012_findind_Y(3, [xpoint, ypoint], ...
% %                     grids.tlong, ...
% %                     grids.tlat, 'CESM2'); % find valid lon, lat index near station
% %         elseif length(sta_lonlat{stai})==4
% %             xpoint1=sta_lonlat{stai}(1);
% %             xpoint2=sta_lonlat{stai}(2);
% %             ypoint1 = sta_lonlat{stai}(3);
% %             ypoint2 = sta_lonlat{stai}(4);
% %             [grids.id_w, grids.id_e, grids.id_s, grids.id_n] = Func_0012_findind_Y(2, [xpoint1,xpoint2, ypoint1,ypoint2], ...
% %                     grids.tlong, ...
% %                     grids.tlat, 'CESM2'); % find valid lon, lat index near station
% %         end
% %         
% %         grids.tlong_cut=grids.tlong(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:);
% %         grids.tlat_cut=grids.tlat(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:);
% % 
% %         tmp.ASSM_mean = Func_0011_get_area_weighted_mean( ...
% %             data_assm_em.([tmp.varname, '_ym'])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:), ...
% %             grids.tlong_cut, grids.tlat_cut);
% %         
% % %     plot(cfg.obs_iyears, (tmp.ASSM_mean-mean(tmp.ASSM_mean))/std(tmp.ASSM_mean));
% % %     hold on
% % %     plot(cfg.obs_iyears, (clim_index_ym-mean(clim_index_ym,'omitnan'))/std(clim_index_ym,'omitnan'));
% % %     hold off
% % %     legend('ASSM', 'AMO')
% %     
% % %     fig_h=figure;
% % 
% %     %% time series -------
% %     fig_weig=2;
% %     set(gcf,'PaperPosition',[0 0 30 20]);
% % 
% %     %% left plot -------
% %     yyaxis left
% %     plot_legend(1)=plot(cfg.obs_iyears, tmp.ASSM_mean, ...
% %     'r-', 'linewidth', 2.*fig_weig);
% %     
% %     set(gca, 'Ycolor', 'r', 'box', 'off', ...
% %         'FontSize', 22.*fig_weig, ...
% %         'TickLabelInterpreter', 'latex');
% % %     axLH.TickLabelInterpreter = 'latex';
% %     
% %     ylabel(gca, ['$$ ', tmp.varname, ' $$'], 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
% %     xlabel(gca,'$$ Year $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
% % 
% %     %% right plot -------
% %     yyaxis right
% %     plot_legend(2)=plot(cfg.obs_iyears, clim_index_ym, ...
% %     'k-', 'linewidth', 2.*fig_weig);
% %     set(gca, 'color', 'none', ...
% %     'ycolor', 'k', 'box', 'off', ...
% %     'FontSize', 22.*fig_weig, ...
% %     'TickLabelInterpreter', 'latex');
% %    
% %     tmp.corr=corrcoef(tmp.ASSM_mean, clim_index_ym, 'Rows', 'complete');
% %     text(2000, -0.25, ['R: ', num2str(round(tmp.corr(1,2),2))], 'fontsize', 25.*fig_weig)
% % 
% %     % ylabel(axRH, '\alpha', 'fontsize', 25.*fig_weig)
% %     ylabel(gca, 'AMO Index', 'fontsize', 25.*fig_weig);
% %     xlabel(gca,'$$ Year $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig);
% %     grid on
% % 
% %     [hLg]=legend(plot_legend, ...
% %     { '$$ ASSM $$', ...
% %     '$$ AMO $$',}, ...
% %     'location', 'northoutside', 'interpreter', 'latex', ...
% %     'Orientation', 'horizontal', 'Fontsize', 35.*fig_weig, ...
% %     'NumColumns', 2);  %% for median
% %     
% %     %% save
% %     dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_TS_AMO_map', filesep, 'assm'];
% %     if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
% %     cfg.figname=[dirs.figdir, filesep, 'TS_ym_AMO_assm_em_', '_', num2str(xpoint1), 'E_', num2str(xpoint2), 'E_', ...
% %                 num2str(ypoint1), 'N_', num2str(ypoint2), 'N_', tmp.varname, '.tif'];
% %     print(gcf, cfg.figname, '-dpng');
% %     RemoveWhiteSpace([], 'file', cfg.figname);
% %     close all;
% % 
% %     end
% %     end
% % 
%% corr ensmeble mean based, OBS <-> ASSM (5deg) (yearly mean)
for fold=1
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
        
        tmp.C=corrval.obs_assm_em_ym.val;
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(corrval.obs_assm_em_ym.val, grid.tlong, grid.tlat);
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_AMO_map', filesep, 'assm'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_em_AMO_assm_ym_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
end

%% corr ensmeble mean based, OBS <-> ASSM (5deg) (jfm mean)
for fold=1
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
        
        tmp.C=corrval.obs_assm_em_jfm.val;
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(corrval.obs_assm_em_jfm.val, grid.tlong, grid.tlat);
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_AMO_map', filesep, 'assm'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_em_AMO_assm_jfm_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
end


%% corr ensmeble mean based, OBS <-> ASSM (5deg) (djf mean)
for fold=1
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
        
        tmp.C=corrval.obs_assm_em_djf.val;
        tmp.C=tmp.C([end, 1:end],:);
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(corrval.obs_assm_em_djf.val, grid.tlong, grid.tlat);
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
        dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_AMO_map', filesep, 'assm'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_em_AMO_assm_djf_', tmp.varname, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
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
