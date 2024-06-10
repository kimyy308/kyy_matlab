% %  Created 01-May-2024 by Yong-Yub Kim
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


hatchflag=1;


%% model configuration
cfg.vlayer=1:31; % surf, vertical slice 

cfg.vlayer_1st=min(cfg.vlayer);
cfg.vlayer_cnt=max(cfg.vlayer)-cfg.vlayer_1st+1;

str.vl_1=num2str(cfg.vlayer_1st, '%02i');
str.vl_max=num2str(max(cfg.vlayer), '%02i');


cfg.regions_set = { [130 130 15 50], ...
                    [143 143 15 50], ...
                    [130 190 20 20], ...
                    [130 190 30 30], ...
                    [130 190 35 35], ...
                    [130 190 42 42]};
for ri=1:length(cfg.regions_set)

    cfg.regions=cfg.regions_set{ri};
%     cfg.regions = [137 137 15 50];

cfg.vars= {'TEMP', 'SALT', 'DIC', 'DIC_ALT_CO2'};
% cfg.vars= {'FG_CO2'};



for vari=1:length(cfg.vars)

cfg.var=cfg.vars{vari};

cfg.comp=Func_0025_CESM2_cmpname_var(cfg.var);

tmp.fdir=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/tr_sysong/mat/', cfg.comp, filesep, cfg.var];
tmp.flist=dir( [tmp.fdir, '/', '*', cfg.var, '*', ...
    num2str(cfg.regions(1)),'*',  num2str(cfg.regions(2)),'*', ...
    num2str(cfg.regions(3)),'*', num2str(cfg.regions(4)),'*', ...
    '_v', str.vl_1, '_v', str.vl_max, '*'] );
    for kk = 1: length (tmp.flist)
        tmp.fname_in = tmp.flist(kk).name;
    end
    output=load([tmp.fdir, filesep, tmp.fname_in]);


    grid=output.grid;
       
    cfg.comp=Func_0025_CESM2_cmpname_var(cfg.var);
    dirs.hcstroot=['/mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive_yearly_transfer/', cfg.comp, '/', cfg.var];
    tmp.gridname = [dirs.hcstroot, tmp.fs, '../grid.nc'];
    tmp.maskname = '/mnt/lustre/proj/kimyy/Model/CESM2/ESP/grids/ocn/RECCAP2_region_masks_all_v20210412_POP2_grid.nc';
    

    %% read & plot data
    tmp.varname=cfg.var;
    
    fig_cfg.name_rgn = 'Glob';
    fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin

    fig_cfg.x_lim = [130 180];
    fig_cfg.y_lim = [15 50];
%     fig_cfg.c_lim = [-1 1];
%     fig_cfg.c_lim2 = [-0.5 0.5];
%     [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
%     [fig_cfg.c_map2, tmp.err_stat] = Func_0009_get_colormaps('bwg_10', tmp.dropboxpath);

loc_column_first=1;
loc_row_first=11;


%% 2-d grid set
if cfg.regions(1)==cfg.regions(2)
    [cut_hor, z_t_2d]=meshgrid(output.grid.cut_tlat, output.grid.z_t(cfg.vlayer_1st:max(cfg.vlayer)));
    str.xlabel='Lat';
elseif cfg.regions(3)==cfg.regions(4)
    [cut_hor, z_t_2d]=meshgrid(output.grid.cut_tlong, output.grid.z_t(cfg.vlayer_1st:max(cfg.vlayer)));
    str.xlabel='Lon';
end
% pcolor(cut_hor,-z_t_2d,squeeze(output.data.TEMP_obs_clim)'); 
% shading flat; 
% colorbar;

%% SUBPLOT(3,2,1); OBS mean
for subi=1:1
    fig_cfg.fig_size = [0,0,13,15]; %% paper size (original)
    fig_cfg.ax_size = [loc_column_first, loc_row_first, 5.4, 2.7];
    fig_cfg.cb_size = [2, 1, 9, 0.3];
    fig_cfg.title_pos = [0.5,1.07];
    fig_cfg.c_map=jet;
    [fig_cfg.c_map, tmp.err]=Func_0009_get_colormaps('bwg_20', tmp.dropboxpath);

    tmp.X=cut_hor;
    tmp.Y=-z_t_2d;
    d_all=[output.data.([cfg.var, '_obs_clim'])(:), output.data.([cfg.var, '_assm_clim'])(:), output.data.([cfg.var, '_lens2_clim'])(:)];
    d_all(abs(d_all)>10e10)=NaN;
    d_all(d_all==0)=NaN;
    fig_caxis1(1)=min(d_all(:), [], 'omitnan');
    fig_caxis1(2)=max(d_all(:), [], 'omitnan');

    tmp.C=squeeze(output.data.([cfg.var, '_obs_clim']))';

    fig_cfg.fig_name='$$ (a) \hspace{1mm}  M(OBS) $$';
    fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
        'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','on');

    %% map setting
    ax_m=subplot(3,2,1);
    set(ax_m, 'Parent', fig_h , 'fontname', 'freeserif');
    hold on;
    grid on;

    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolor(tmp.X,tmp.Y,tmp.C,'parent',ax_m); 
    shading flat;
    axis tight;

    %% frame and label setting
    xlabel(str.xlabel);
    ylabel('Depth');

    %% color set
    caxis(ax_m, fig_caxis1);
    colormap(ax_m,fig_cfg.c_map);
    set(ax_m, 'Fontsize', 15)
%     disp('a')
end


%% SUBPLOT(3,2,2); corr, OBS std
for subi=1:1
    fig_cfg.ax_size = [loc_column_first+5.7, loc_row_first, 5.4, 2.7];

    tmp.X=cut_hor;
    tmp.Y=-z_t_2d;
    d_all=[output.data.([cfg.var, '_obs_stdt'])(:), output.data.([cfg.var, '_assm_stdt'])(:), output.data.([cfg.var, '_lens2_stdt'])(:)];
    d_all(abs(d_all)>10e10)=NaN;
    d_all(d_all==0)=NaN;
    fig_caxis2(1)=min(d_all(:), [], 'omitnan');
    fig_caxis2(2)=max(d_all(:), [], 'omitnan');

%     fig_caxis2(1)=min([min(output.data.([cfg.var, '_obs_stdt'])(:)), min(output.data.([cfg.var, '_assm_stdt'])(:)), min(output.data.([cfg.var, '_lens2_stdt'])(:))]);
%     fig_caxis2(2)=max([max(output.data.([cfg.var, '_obs_stdt'])(:)), max(output.data.([cfg.var, '_assm_stdt'])(:)), max(output.data.([cfg.var, '_lens2_stdt'])(:))]);

    tmp.C=squeeze(output.data.([cfg.var, '_obs_stdt']))';

    fig_cfg.fig_name='$$ (b) \hspace{1mm}  \sigma (OBS) $$';

    %% map setting
    ax_m=subplot(3,2,2);
    set(ax_m, 'Parent', fig_h , 'fontname', 'freeserif');
    hold on;
    grid on;

    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolor(tmp.X,tmp.Y,tmp.C,'parent',ax_m); 
    shading flat;
    axis tight;

    %% frame and label setting
    xlabel(str.xlabel);
    ylabel('Depth');
    ax_m.YAxisLocation='right';

    %% color set
    caxis(ax_m, fig_caxis2);
    colormap(ax_m,'jet');
    set(ax_m, 'Fontsize', 15)
end

%% SUBPLOT(3,2,3); corr, ASSM mean
for subi=1:1
    fig_cfg.ax_size = [loc_column_first, loc_row_first-3.8, 5.4, 2.7];

    tmp.X=cut_hor;
    tmp.Y=-z_t_2d;

    tmp.C=squeeze(output.data.([cfg.var, '_assm_clim']))';

    fig_cfg.fig_name='$$ (c) \hspace{1mm}  M(ASSM) $$';
   
    %% map setting
    ax_m=subplot(3,2,3);
    set(ax_m, 'Parent', fig_h , 'fontname', 'freeserif');
    hold on;
    grid on;

    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolor(tmp.X,tmp.Y,tmp.C,'parent',ax_m); 
    shading flat;
    axis tight;

    %% frame and label setting
    xlabel(str.xlabel);
    ylabel('Depth');

    %% color set
    caxis(ax_m, fig_caxis1);
    colormap(ax_m,fig_cfg.c_map);
    set(ax_m, 'Fontsize', 15)
end

%% caxis & colorbar (ACC)
    caxis(ax_m, fig_caxis1); 
    colormap(ax_m,fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','Location', 'southoutside', 'position',fig_cfg.cb_size + [0, 1, 0, 0]);
    set(cb,'fontsize',15,'fontname','freeserif','TickDir','both');
%     cb.TickLabelPosition = 'top';
%     cb_title=title(cb,'$$ r $$','fontsize', 22, 'Position', [880, 0, 0]); % hor, ver, ?
    cb_title=title(cb,['$$ ', output.data.units, ' $$'],'fontsize', 22, 'Position', [680, 0, 0]); % hor, ver, ?
%     cb_title=title(cb,['$$ ', output.data.units, ' $$'],'fontsize', 22, 'Position', [660, 0, 0]); % hor, ver, ?

    set(cb_title, 'interpreter', 'latex');


%% SUBPLOT(3,2,4); ASSM std
for subi=1:1
    fig_cfg.ax_size = [loc_column_first+5.7, loc_row_first-3.8, 5.4, 2.7];

    tmp.X=cut_hor;
    tmp.Y=-z_t_2d;

    tmp.C=squeeze(output.data.([cfg.var, '_assm_stdt']))';

    fig_cfg.fig_name='$$ (d) \hspace{1mm}  \sigma (ASSM) $$';

    %% map setting
    ax_m=subplot(3,2,4);
    set(ax_m, 'Parent', fig_h , 'fontname', 'freeserif');
    hold on;
    grid on;

    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolor(tmp.X,tmp.Y,tmp.C,'parent',ax_m); 
    shading flat;
    axis tight;

    %% frame and label setting
    xlabel(str.xlabel);
    ylabel('Depth');
    ax_m.YAxisLocation='right';

    %% color set
    caxis(ax_m, fig_caxis2);
    colormap(ax_m,'jet');
    set(ax_m, 'Fontsize', 15)
end



%% SUBPLOT(3,2,5); LENS2 mean
for subi=1:1
    fig_cfg.ax_size = [loc_column_first, loc_row_first-7.6, 5.4, 2.7];

    tmp.X=cut_hor;
    tmp.Y=-z_t_2d;

    tmp.C=squeeze(output.data.([cfg.var, '_lens2_clim']))';

    fig_cfg.fig_name='$$ (e) \hspace{1mm}  M(LENS2) $$';
   
    %% map setting
    ax_m=subplot(3,2,5);
    set(ax_m, 'Parent', fig_h , 'fontname', 'freeserif');
    hold on;
    grid on;

    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolor(tmp.X,tmp.Y,tmp.C,'parent',ax_m); 
    shading flat;
    axis tight;

    %% frame and label setting
    xlabel(str.xlabel);
    ylabel('Depth');

    %% color set
    caxis(ax_m, fig_caxis1);
    colormap(ax_m,fig_cfg.c_map);
    set(ax_m, 'Fontsize', 15)
end

%% SUBPLOT(3,2,6); corr, ASSM<->HCST - ASSM <-> LENS2, median(individual) (LY1)
for subi=1:1
    fig_cfg.ax_size = [loc_column_first+5.7, loc_row_first-7.6, 5.4, 2.7];

    tmp.X=cut_hor;
    tmp.Y=-z_t_2d;

    tmp.C=squeeze(output.data.([cfg.var, '_lens2_stdt']))';

    fig_cfg.fig_name='$$ (f) \hspace{1mm}  \sigma (LENS2) $$';

    %% map setting
    ax_m=subplot(3,2,6);
    set(ax_m, 'Parent', fig_h , 'fontname', 'freeserif');
    hold on;
    grid on;

    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolor(tmp.X,tmp.Y,tmp.C,'parent',ax_m); 
    shading flat;
    axis tight;

    %% frame and label setting
    xlabel(str.xlabel);
    ylabel('Depth');
    ax_m.YAxisLocation='right';

    %% color set
    caxis(ax_m, fig_caxis2);
    colormap(ax_m,'jet');
    set(ax_m, 'Fontsize', 15)
end

%% caxis & colorbar (dACC)
    caxis(ax_m, fig_caxis2); 
    colormap(ax_m, jet);
    cb = colorbar(ax_m,'units','inches','Location', 'southoutside', 'position',fig_cfg.cb_size + [0 0 0 0]);
    set(cb,'fontsize',15,'fontname','freeserif','TickDir','both');
%     cb.TickLabelPosition = 'top';
%     cb_title=title(cb,'$$ \Delta r $$','fontsize', 22, 'Position', [680, -4, 0]); % hor, ver, ?
    cb_title=title(cb,['$$ ', output.data.units, ' $$'],'fontsize', 22, 'Position', [680, -4, 0]); % hor, ver, ?

    set(cb_title, 'interpreter', 'latex');    

    %% save
    dirs.figdir= ['/Users/kimyy/Desktop/backup/Research/Postdoc/03_IBS/2023_Kuroshio_Extension_Carbon_Cycle/Figure', filesep, cfg.var, filesep, 'vert_section'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, ...
        filesep, cfg.var,'_clim_std_', ...
        'h_', num2str(cfg.regions(1)), '_', num2str(cfg.regions(2)), '_', ...
        num2str(cfg.regions(3)), '_', num2str(cfg.regions(4)), '_', ...
        'v_', num2str(cfg.vlayer_1st,'%02i'), '_', ...
        num2str(max(cfg.vlayer),'%02i'),'.tif'];
    print(fig_h, cfg.figname, '-dpng');
%     RemoveWhiteSpace([], 'file', cfg.figname);
    close all;

% end

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
