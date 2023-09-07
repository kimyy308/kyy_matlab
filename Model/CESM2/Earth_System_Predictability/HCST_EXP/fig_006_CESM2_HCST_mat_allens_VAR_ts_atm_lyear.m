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

cfg.vars = {'SST', 'TS', 'PSL', 'PRECT'};

cfg.vlayer=1; % surface, vertical slice
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

disp(cfg.var);
tic;

% dirs.hcstroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/HCST_EXP/archive/', cfg.comp, '/', cfg.var];
% dirs.obsroot=['/Volumes/kyy_raid/kimyy/Observation/', cfg.obs_name, '/monthly_reg_', cfg.obs_fname_module(2:4)];
% dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', cfg.var];
% dirs.lens2root=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/LENS2/archive_analysis/', cfg.comp, '/', cfg.var];
% dirs.assmroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/ASSM_EXP/archive_analysis/', cfg.comp, '/', cfg.var];

dirs.hcstmatroot=[tmp.kimyypath, '/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', cfg.var];
dirs.figroot=[tmp.kimyypath, '/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', cfg.var];



cfg.iyears=cfg.obs_iyears;
cfg.gnm='f09_g17';
% cfg.proj_year=1;
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

% grids.region_mask=ncread(tmp.gridname, 'REGION_MASK'); 
% grids.ocean_mask=NaN(size(grids.region_mask));
% grids.ocean_mask(grids.region_mask>0)=1;
% grids.tarea = ncread(tmp.gridname, 'TAREA');

switch cfg.comp
    case {'ocn', 'ice'}
        grids.tlong=ncread(tmp.gridname, 'TLONG');
        grids.tlat=ncread(tmp.gridname, 'TLAT');
        grids.mask_ocn=ncread(tmp.maskname, 'open_ocean');
        grids.mask_ocn(grids.mask_ocn<-10e10)=NaN;
        grids.mask_ocn=grids.mask_ocn./grids.mask_ocn;
        grids.tarea=ncread(tmp.gridname, 'TAREA')/1000000.0; %(m^2 -> km^2)
        grids.tarea_60=grids.tarea; grids.tarea_60(grids.tlat>60 | grids.tlat<-60)=NaN;
    case {'atm', 'lnd'}
        grids.lon=ncread(tmp.gridname, 'lon');
        grids.lat=ncread(tmp.gridname, 'lat');
        [grids.tlat grids.tlong]=meshgrid(grids.lat, grids.lon);
        grids.tarea=ncread(tmp.gridname, 'AREA')/1000000.0; %(m^2 -> km^2)
        grids.tarea_60=grids.tarea; grids.tarea_60(grids.tlat>60 | grids.tlat<-60)=NaN;
end

grids.nlon=size(grids.tlong,1);
grids.nlat=size(grids.tlat,2);
% grids.ntime=cfg.proj_year.*12;


sta_lonlat = {[20, 15], [20, 80], [40, 3], [50, 9], [62, 32], ...
              [50, 65], [70, 42], [65, 60], [82, 48], [82, 57], ...
              [105, 57], [110, 41], [120, 50], [122, 67], [130, 60], ...
              [120, -30], [120, -21], [131, -30], [131, -20], [140,-30], ...
              [147, -22], [140, -8], [200, 65], [220, 61], [230, 55], ...
              [248, 56], [248, 30], [252, 36], [260, 38], [270, 75], ...
              [270, 52], [285, -5], [287, 5], [292, 5], [293, -42], ...
              [296, -35], [308, -12], [320, -9], [320, 65], [338, 75], ...
              [340, 65], [350, 16], [360, 15], [360, 28], [10, 15], [8, 30], ...
              [70, -75], [125, -80], [240, -80], [340, -80], [360, -80], [18, 78], ...
              [345, 75], [300, 85], [260, 85], [220, 80], [180, 80], ...
                [140, 85], [100, 85], [40, 75], [360, 65], [310, 55], ...
                [230, 50], [200, 50], [180, 50], [150, 55], [180, 35], ...
                [140, 33], [124, 35], [240,20], [180, 20], [150, 15], ...
                [180, 10], [260, 0], [160, 0], [280, -20], [240, -20], ...
                [200, -20], [160, -20], [160, -30], [200, -40], [200, -50], ...
                [110, -10], [90, 15], [60, 10], [42, -10], [40, -30], ...
                [40, -40], [40, -50], [340, 40], [320, 30], [300, 30], ...
                [285, 30], [340, 20], [300, 20], [340, 10], [320, 10], ...
                [270, 25], [360, 0], [320, 0], [10, -10], [330, -10], ...
                [330, -30], [5, 40], [60, -65], [60, -25], [180, -70], ...
                [240, -70], [340, -50], [150, 30], [5, 35], ...
                [320, -70], [200, 21], [238, 32], [200, 21], ...
                [5, 60], [151, 34], [190,0], [160, 47], ...
                [137, 3], [137, 5], [137,10], [137, 15], [137, 20], ...
                [137,25], [137,30], [137, 34],[297, 24], [297, 30], ...
                 [297, 35], [300, 40], [60,-30], [90, 15], [140, 15], ...
                 [200, 20], [200,50], [250, 20], [160, -10], [180,-10],...
                 [210,-25], [300, 25], [330, -10], [360, -10], [330, -30]};

%% caution!!! if you use regional mean, spread (regional mean of spread of each grids) is wrong, spread should not be used
sta_lonlat = {[0 360 -90 90], [190 240 -5 5]};


sta_lonlat = {[205 0]};


% CalCOFI : 238, 32
% K2 : 160, 47`
% JMA : 137, 3 ~ 34
% BATS : 294 ~ 300, 24 ~ 35
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
    disp(fig_cfg.mat_name)


    for stai=1:length(sta_lonlat)
        fig_h = figure('name','ts','visible','off');
        fig_h.Position= [0,0,1000,500];
        if length(sta_lonlat{stai})==2
            xpoint=sta_lonlat{stai}(1);
            ypoint = sta_lonlat{stai}(2);
            [grids.id_w, grids.id_e, grids.id_s, grids.id_n] = Func_0012_findind_Y(2, [xpoint, ypoint], ...
                    grids.tlong, ...
                    grids.tlat, 'CESM2'); % find valid lon, lat index near station
        elseif length(sta_lonlat{stai})==4
            xpoint1=sta_lonlat{stai}(1);
            xpoint2=sta_lonlat{stai}(2);
            ypoint1 = sta_lonlat{stai}(3);
            ypoint2 = sta_lonlat{stai}(4);
            [grids.id_w, grids.id_e, grids.id_s, grids.id_n] = Func_0012_findind_Y(2, [xpoint1,xpoint2, ypoint1,ypoint2], ...
                    grids.tlong, ...
                    grids.tlat, 'CESM2'); % find valid lon, lat index near station
        end
        
        %% colormap set
        val_transparent = 0.15;

        cmap_HCST = [0,0,1]; % blue
        cmap_HCST_b = rgb2hsv(cmap_HCST);
        cmap_HCST_b(:,2) =  val_transparent;
        cmap_HCST_b = hsv2rgb(cmap_HCST_b);

        cmap_ASSM = [1, 0, 0]; % Red [1, 0, 0], Yellow [1, 1, 0]
        cmap_ASSM_b = rgb2hsv(cmap_ASSM);
        cmap_ASSM_b(:,2) =  val_transparent;
        cmap_ASSM_b = hsv2rgb(cmap_ASSM_b);
        
        cmap_LENS2 = [0.5, 1.0, 0.5]; % Red [1, 0, 0], Yellow [1, 1, 0], Green [0.5, 1.0, 0.5]
        cmap_LENS2_b = rgb2hsv(cmap_LENS2);
        cmap_LENS2_b(:,2) =  val_transparent;
        cmap_LENS2_b = hsv2rgb(cmap_LENS2_b);

        cmap_OBS = [0 0 0]; % black [0 0 0];
        
        grids.tlong_cut=grids.tlong(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:);
        grids.tlat_cut=grids.tlat(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:);
        
        tmp.HCST_mean = Func_0011_get_area_weighted_mean(data.([tmp.varname, '_model_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:), grids.tlong_cut, grids.tlat_cut);
        tmp.HCST_lower= Func_0011_get_area_weighted_mean(data.([tmp.varname, '_model_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:) - ...
            squeeze(data.([tmp.varname, '_model_stde_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:))/2.0, grids.tlong_cut, grids.tlat_cut);
        tmp.HCST_upper= Func_0011_get_area_weighted_mean(data.([tmp.varname, '_model_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:) + ...
            squeeze(data.([tmp.varname, '_model_stde_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:))/2.0, grids.tlong_cut, grids.tlat_cut);

        tmp.ASSM_mean = Func_0011_get_area_weighted_mean(data.([tmp.varname, '_assm'])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:), grids.tlong_cut, grids.tlat_cut);
        tmp.ASSM_lower= Func_0011_get_area_weighted_mean(data.([tmp.varname, '_assm'])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:) - ...
            squeeze(data.([tmp.varname, '_assm_stde'])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:))/2.0, grids.tlong_cut, grids.tlat_cut);
        tmp.ASSM_upper= Func_0011_get_area_weighted_mean(data.([tmp.varname, '_assm'])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:) + ...
            squeeze(data.([tmp.varname, '_assm_stde'])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:))/2.0, grids.tlong_cut, grids.tlat_cut);

        tmp.LENS2_mean = Func_0011_get_area_weighted_mean(data.([tmp.varname, '_lens2_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:), grids.tlong_cut, grids.tlat_cut);
        tmp.LENS2_lower= Func_0011_get_area_weighted_mean(data.([tmp.varname, '_lens2_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:) - ...
            squeeze(data.([tmp.varname, '_lens2_stde_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:))/2.0, grids.tlong_cut, grids.tlat_cut);
        tmp.LENS2_upper= Func_0011_get_area_weighted_mean(data.([tmp.varname, '_lens2_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:) + ...
            squeeze(data.([tmp.varname, '_lens2_stde_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:))/2.0, grids.tlong_cut, grids.tlat_cut);
        
        tmp.OBS = Func_0011_get_area_weighted_mean(data.([tmp.varname, '_obs'])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:), grids.tlong_cut, grids.tlat_cut);
        
        tmp.time= cfg.iyears+lyear;
        
        hold on
        if length(sta_lonlat{stai})==2
            fig_ts.LENS2_range=fill([tmp.time(1:end-lyear), flip(tmp.time(1:end-lyear))], ...
                [tmp.LENS2_lower(1:end-lyear); flip(tmp.LENS2_upper(1:end-lyear))], cmap_LENS2_b);
            fig_ts.LENS2_range.EdgeColor = 'none';
            set(get(get(fig_ts.LENS2_range,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            
            fig_ts.HCST_range=fill([tmp.time, flip(tmp.time)], ...
                [tmp.HCST_lower; flip(tmp.HCST_upper)], cmap_HCST_b);
            fig_ts.HCST_range.EdgeColor = 'none';
            set(get(get(fig_ts.HCST_range,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
            fig_ts.ASSM_range=fill([tmp.time(1:end-lyear), flip(tmp.time(1:end-lyear))], ...
                [tmp.ASSM_lower(1:end-lyear); flip(tmp.ASSM_upper(1:end-lyear))], cmap_ASSM_b);
            fig_ts.ASSM_range.EdgeColor = 'none';
            set(get(get(fig_ts.ASSM_range,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end


        plot(tmp.time(1:end-lyear), tmp.LENS2_mean(1:end-lyear), 'linewidth', 2, 'color', cmap_LENS2)
        plot(tmp.time(1:end-lyear), tmp.ASSM_mean(1:end-lyear), 'linewidth', 2, 'color', cmap_ASSM)
        plot(tmp.time, tmp.HCST_mean, 'linewidth', 2, 'color', cmap_HCST)
        plot(tmp.time(1:end-lyear), tmp.OBS(1:end-lyear), 'linewidth', 2, 'color', cmap_OBS)

        hold off
        
%         plot(cfg.iyears+lyear,squeeze(data.([tmp.varname, '_model_l', tmp.lyear_str])(grids.id_w,grids.id_s,:)), 'linewidth', 2)
%         hold on
%         plot(cfg.iyears+lyear,squeeze(data.([tmp.varname, '_assm'])(grids.id_w,grids.id_s,:)), 'linewidth', 2)
%         plot(cfg.iyears+lyear,squeeze(data.([tmp.varname, '_lens2_l', tmp.lyear_str])(grids.id_w,grids.id_s,:)), 'linewidth', 2)
%         hold off
        legend ('LENS2', 'ASSM', 'HCST', 'OBS', 'Location', 'Northwest', 'Orientation', 'Horizontal')
        if length(sta_lonlat{stai})==2
            title(['l', tmp.lyear_str, ', ', tmp.varname, ', ', num2str(xpoint),'E, ', num2str(ypoint), 'N'])
        elseif length(sta_lonlat{stai})==4
            title(['l', tmp.lyear_str, ', ', tmp.varname, ', ', num2str(xpoint1),'E~', num2str(xpoint2),'E, ', ...
                num2str(ypoint1), 'N~',num2str(ypoint2), 'N'])
        end
        grid minor
        xlim([1960 2025])
        set(gca, 'fontsize', 20)
        
        %% corr skills
        [tmp.pot_skill_hcst, tmp.pot_skill_hcst_p]=corrcoef(tmp.ASSM_mean(1:end-lyear), tmp.HCST_mean(1:end-lyear));
        if (tmp.pot_skill_hcst_p>0.05) tmp.pot_skill_hcst=NaN(2,2); end
        [tmp.pot_skill_lens2, tmp.pot_skill_lens2_p]=corrcoef(tmp.ASSM_mean(1:end-lyear), tmp.LENS2_mean(1:end-lyear));
        if (tmp.pot_skill_lens2_p>0.05) tmp.pot_skill_lens2=NaN(2,2); end

        [tmp.skill_assm, tmp.skill_assm_p]=corrcoef(tmp.OBS(1:end-lyear), tmp.ASSM_mean(1:end-lyear));
        if (tmp.skill_assm_p>0.05) tmp.skill_assm=NaN(2,2); end
        [tmp.skill_hcst, tmp.skill_hcst_p]=corrcoef(tmp.OBS(1:end-lyear), tmp.HCST_mean(1:end-lyear));
        if (tmp.skill_hcst_p>0.05) tmp.skill_hcst=NaN(2,2); end
        [tmp.skill_lens2, tmp.skill_lens2_p]=corrcoef(tmp.OBS(1:end-lyear), tmp.LENS2_mean(1:end-lyear));
        if (tmp.skill_lens2_p>0.05) tmp.skill_lens2=NaN(2,2); end
        
        yl=ylim;
        text(1960, min(yl)+diff(yl)/30, ['AS-HC:', num2str(round(tmp.pot_skill_hcst(1,2),2))])
        text(1967, min(yl)+diff(yl)/30, ['AS-LE:', num2str(round(tmp.pot_skill_lens2(1,2),2))])
        
        text(1980, min(yl)+diff(yl)/30, ['OB-AS:', num2str(round(tmp.skill_assm(1,2),2))])
        text(1987, min(yl)+diff(yl)/30, ['OB-HC:', num2str(round(tmp.skill_hcst(1,2),2))])
        text(1994, min(yl)+diff(yl)/30, ['OB-LE:', num2str(round(tmp.skill_lens2(1,2),2))])


        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_time_series', filesep, 'l',tmp.lyear_str];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        if length(sta_lonlat{stai})==2
            cfg.figname=[dirs.figdir, filesep, 'ts_all_l',tmp.lyear_str, '_', num2str(xpoint), 'E_', num2str(ypoint), 'N_', tmp.varname, '.tif'];
        elseif length(sta_lonlat{stai})==4
            cfg.figname=[dirs.figdir, filesep, 'ts_all_l',tmp.lyear_str, '_', num2str(xpoint1), 'E_', num2str(xpoint2), 'E_', ...
                num2str(ypoint1), 'N_', num2str(ypoint2), 'N_', tmp.varname, '.tif'];
        end
            print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;


        %% anomaly plot
        fig_h = figure('name','ts','visible','off');
        fig_h.Position= [0,0,1000,500];

        tmp.HCST_tmean=mean(tmp.HCST_mean, 'omitnan');
        tmp.HCST_mean_ano=tmp.HCST_mean-tmp.HCST_tmean;
        tmp.HCST_lower_ano= tmp.HCST_lower-tmp.HCST_tmean;
        tmp.HCST_upper_ano= tmp.HCST_upper-tmp.HCST_tmean;

        tmp.ASSM_tmean=mean(tmp.ASSM_mean, 'omitnan');
        tmp.ASSM_mean_ano=tmp.ASSM_mean-tmp.ASSM_tmean;
        tmp.ASSM_lower_ano= tmp.ASSM_lower-tmp.ASSM_tmean;
        tmp.ASSM_upper_ano= tmp.ASSM_upper-tmp.ASSM_tmean;

        tmp.LENS2_tmean=mean(tmp.LENS2_mean, 'omitnan');
        tmp.LENS2_mean_ano=tmp.LENS2_mean-tmp.LENS2_tmean;
        tmp.LENS2_lower_ano= tmp.LENS2_lower-tmp.LENS2_tmean;
        tmp.LENS2_upper_ano= tmp.LENS2_upper-tmp.LENS2_tmean;

        tmp.OBS_tmean=mean(tmp.OBS, 'omitnan');
        tmp.OBS_ano=tmp.OBS-tmp.OBS_tmean;

        hold on

        if length(sta_lonlat{stai})==2
            fig_ts.LENS2_range=fill([tmp.time(1:end-lyear), flip(tmp.time(1:end-lyear))], ...
                [tmp.LENS2_lower_ano(1:end-lyear); flip(tmp.LENS2_upper_ano(1:end-lyear))], cmap_LENS2_b);
            fig_ts.LENS2_range.EdgeColor = 'none';
            set(get(get(fig_ts.LENS2_range,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
            fig_ts.HCST_range=fill([tmp.time, flip(tmp.time)], ...
                [tmp.HCST_lower_ano; flip(tmp.HCST_upper_ano)], cmap_HCST_b);
            fig_ts.HCST_range.EdgeColor = 'none';
            set(get(get(fig_ts.HCST_range,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
            fig_ts.ASSM_range=fill([tmp.time(1:end-lyear), flip(tmp.time(1:end-lyear))], ...
                [tmp.ASSM_lower_ano(1:end-lyear); flip(tmp.ASSM_upper_ano(1:end-lyear))], cmap_ASSM_b);
            fig_ts.ASSM_range.EdgeColor = 'none';
            set(get(get(fig_ts.ASSM_range,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end

        plot(tmp.time(1:end-lyear), tmp.LENS2_mean_ano(1:end-lyear), 'linewidth', 2, 'color', cmap_LENS2)
        plot(tmp.time(1:end-lyear), tmp.ASSM_mean_ano(1:end-lyear), 'linewidth', 2, 'color', cmap_ASSM)
        plot(tmp.time, tmp.HCST_mean_ano, 'linewidth', 2, 'color', cmap_HCST)
        plot(tmp.time(1:end-lyear), tmp.OBS_ano(1:end-lyear), 'linewidth', 2, 'color', cmap_OBS)

        hold off

        legend ('LENS2', 'ASSM', 'HCST', 'OBS', 'Location', 'Northwest')
        if length(sta_lonlat{stai})==2
            title(['anom, l', tmp.lyear_str, ', ', tmp.varname, ', ', num2str(xpoint),'E, ', num2str(ypoint), 'N'])
        elseif length(sta_lonlat{stai})==4
            title(['anom, l', tmp.lyear_str, ', ', tmp.varname, ', ', num2str(xpoint1),'E~', num2str(xpoint2),'E, ', ...
                num2str(ypoint1), 'N~',num2str(ypoint2), 'N'])
        end
        grid minor
        xlim([1960 2025])
        set(gca, 'fontsize', 20)
        
        %% corr skills
        [tmp.pot_skill_hcst, tmp.pot_skill_hcst_p]=corrcoef(tmp.ASSM_mean_ano(1:end-lyear), tmp.HCST_mean_ano(1:end-lyear));
        if (tmp.pot_skill_hcst_p>0.05) tmp.pot_skill_hcst=NaN(2,2); end
        [tmp.pot_skill_lens2, tmp.pot_skill_lens2_p]=corrcoef(tmp.ASSM_mean_ano(1:end-lyear), tmp.LENS2_mean_ano(1:end-lyear));
        if (tmp.pot_skill_lens2_p>0.05) tmp.pot_skill_lens2=NaN(2,2); end

        [tmp.skill_assm, tmp.skill_assm_p]=corrcoef(tmp.OBS_ano(1:end-lyear), tmp.ASSM_mean_ano(1:end-lyear));
        if (tmp.skill_assm_p>0.05) tmp.skill_assm=NaN(2,2); end
        [tmp.skill_hcst, tmp.skill_hcst_p]=corrcoef(tmp.OBS_ano(1:end-lyear), tmp.HCST_mean_ano(1:end-lyear));
        if (tmp.skill_hcst_p>0.05) tmp.skill_hcst=NaN(2,2); end
        [tmp.skill_lens2, tmp.skill_lens2_p]=corrcoef(tmp.OBS_ano(1:end-lyear), tmp.LENS2_mean_ano(1:end-lyear));
        if (tmp.skill_lens2_p>0.05) tmp.skill_lens2=NaN(2,2); end
        
        yl=ylim;
        text(1960, min(yl)+diff(yl)/30, ['AS-HC:', num2str(round(tmp.pot_skill_hcst(1,2),2))])
        text(1967, min(yl)+diff(yl)/30, ['AS-LE:', num2str(round(tmp.pot_skill_lens2(1,2),2))])
        
        text(1980, min(yl)+diff(yl)/30, ['OB-AS:', num2str(round(tmp.skill_assm(1,2),2))])
        text(1987, min(yl)+diff(yl)/30, ['OB-HC:', num2str(round(tmp.skill_hcst(1,2),2))])
        text(1994, min(yl)+diff(yl)/30, ['OB-LE:', num2str(round(tmp.skill_lens2(1,2),2))])


        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_time_series_ano', filesep, 'l',tmp.lyear_str];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        if length(sta_lonlat{stai})==2
            cfg.figname=[dirs.figdir, filesep, 'ts_ano_all_l',tmp.lyear_str, '_', num2str(xpoint), 'E_', num2str(ypoint), 'N_', tmp.varname, '.tif'];
        elseif length(sta_lonlat{stai})==4
            cfg.figname=[dirs.figdir, filesep, 'ts_ano_all_l',tmp.lyear_str, '_', num2str(xpoint1), 'E_', num2str(xpoint2), 'E_', ...
                num2str(ypoint1), 'N_', num2str(ypoint2), 'N_', tmp.varname, '.tif'];
        end
            print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;

        
        %% detrended plot
        fig_h = figure('name','ts','visible','off');
        fig_h.Position= [0,0,1000,500];

        [tmp.HCST_mean_det, tmp.HCST_mean_trend] = Func_0028_detrend_linear_1d(tmp.HCST_mean, 'omitnan');
        tmp.HCST_mean_tr_ano=tmp.HCST_mean-tmp.HCST_mean_det;
        tmp.HCST_lower_det = tmp.HCST_lower-tmp.HCST_mean_tr_ano;
        tmp.HCST_upper_det = tmp.HCST_upper-tmp.HCST_mean_tr_ano;

        [tmp.ASSM_mean_det, tmp.ASSM_mean_trend] = Func_0028_detrend_linear_1d(tmp.ASSM_mean, 'omitnan');
        tmp.ASSM_mean_tr_ano=tmp.ASSM_mean-tmp.ASSM_mean_det;
        tmp.ASSM_lower_det = tmp.ASSM_lower-tmp.ASSM_mean_tr_ano;
        tmp.ASSM_upper_det = tmp.ASSM_upper-tmp.ASSM_mean_tr_ano;

        [tmp.LENS2_mean_det, tmp.LENS2_mean_trend] = Func_0028_detrend_linear_1d(tmp.LENS2_mean, 'omitnan');
        tmp.LENS2_mean_tr_ano=tmp.LENS2_mean-tmp.LENS2_mean_det;
        tmp.LENS2_lower_det = tmp.LENS2_lower-tmp.LENS2_mean_tr_ano;
        tmp.LENS2_upper_det = tmp.LENS2_upper-tmp.LENS2_mean_tr_ano;

        [tmp.OBS_det, tmp.OBS_trend] =Func_0028_detrend_linear_1d(tmp.OBS, 'omitnan');

        hold on
        if length(sta_lonlat{stai})==2
            fig_ts.LENS2_range=fill([tmp.time(1:end-lyear), flip(tmp.time(1:end-lyear))], ...
                [tmp.LENS2_lower_det(1:end-lyear); flip(tmp.LENS2_upper_det(1:end-lyear))], cmap_LENS2_b);
            fig_ts.LENS2_range.EdgeColor = 'none';
            set(get(get(fig_ts.LENS2_range,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
            fig_ts.HCST_range=fill([tmp.time, flip(tmp.time)], ...
                [tmp.HCST_lower_det; flip(tmp.HCST_upper_det)], cmap_HCST_b);
            fig_ts.HCST_range.EdgeColor = 'none';
            set(get(get(fig_ts.HCST_range,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
            fig_ts.ASSM_range=fill([tmp.time(1:end-lyear), flip(tmp.time(1:end-lyear))], ...
                [tmp.ASSM_lower_det(1:end-lyear); flip(tmp.ASSM_upper_det(1:end-lyear))], cmap_ASSM_b);
            fig_ts.ASSM_range.EdgeColor = 'none';
            set(get(get(fig_ts.ASSM_range,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end


        plot(tmp.time(1:end-lyear), tmp.LENS2_mean_det(1:end-lyear), 'linewidth', 2, 'color', cmap_LENS2)
        plot(tmp.time(1:end-lyear), tmp.ASSM_mean_det(1:end-lyear), 'linewidth', 2, 'color', cmap_ASSM)
        plot(tmp.time, tmp.HCST_mean_det, 'linewidth', 2, 'color', cmap_HCST)
        plot(tmp.time(1:end-lyear), tmp.OBS_det(1:end-lyear), 'linewidth', 2, 'color', cmap_OBS)

        hold off

        legend ('LENS2', 'ASSM', 'HCST', 'OBS', 'Location', 'Northwest')
        if length(sta_lonlat{stai})==2
            title(['det, l', tmp.lyear_str, ', ', tmp.varname, ', ', num2str(xpoint),'E, ', num2str(ypoint), 'N'])
        elseif length(sta_lonlat{stai})==4
            title(['det, l', tmp.lyear_str, ', ', tmp.varname, ', ', num2str(xpoint1),'E~', num2str(xpoint2),'E, ', ...
                num2str(ypoint1), 'N~',num2str(ypoint2), 'N'])
        end
        grid minor
        xlim([1960 2025])
        set(gca, 'fontsize', 20)
        
        %% corr skills
        [tmp.pot_skill_hcst, tmp.pot_skill_hcst_p]=corrcoef(tmp.ASSM_mean_det(1:end-lyear), tmp.HCST_mean_det(1:end-lyear));
        if (tmp.pot_skill_hcst_p>0.05) tmp.pot_skill_hcst=NaN(2,2); end
        [tmp.pot_skill_lens2, tmp.pot_skill_lens2_p]=corrcoef(tmp.ASSM_mean_det(1:end-lyear), tmp.LENS2_mean_det(1:end-lyear));
        if (tmp.pot_skill_lens2_p>0.05) tmp.pot_skill_lens2=NaN(2,2); end

        [tmp.skill_assm, tmp.skill_assm_p]=corrcoef(tmp.OBS_det(1:end-lyear), tmp.ASSM_mean_det(1:end-lyear));
        if (tmp.skill_assm_p>0.05) tmp.skill_assm=NaN(2,2); end
        [tmp.skill_hcst, tmp.skill_hcst_p]=corrcoef(tmp.OBS_det(1:end-lyear), tmp.HCST_mean_det(1:end-lyear));
        if (tmp.skill_hcst_p>0.05) tmp.skill_hcst=NaN(2,2); end
        [tmp.skill_lens2, tmp.skill_lens2_p]=corrcoef(tmp.OBS_det(1:end-lyear), tmp.LENS2_mean_det(1:end-lyear));
        if (tmp.skill_lens2_p>0.05) tmp.skill_lens2=NaN(2,2); end
        
        yl=ylim;
        text(1960, min(yl)+diff(yl)/30, ['AS-HC:', num2str(round(tmp.pot_skill_hcst(1,2),2))])
        text(1967, min(yl)+diff(yl)/30, ['AS-LE:', num2str(round(tmp.pot_skill_lens2(1,2),2))])
        
        text(1980, min(yl)+diff(yl)/30, ['OB-AS:', num2str(round(tmp.skill_assm(1,2),2))])
        text(1987, min(yl)+diff(yl)/30, ['OB-HC:', num2str(round(tmp.skill_hcst(1,2),2))])
        text(1994, min(yl)+diff(yl)/30, ['OB-LE:', num2str(round(tmp.skill_lens2(1,2),2))])

        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_time_series_det', filesep, 'l',tmp.lyear_str];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        if length(sta_lonlat{stai})==2
            cfg.figname=[dirs.figdir, filesep, 'ts_det_all_l',tmp.lyear_str, '_', num2str(xpoint), 'E_', num2str(ypoint), 'N_', tmp.varname, '.tif'];
        elseif length(sta_lonlat{stai})==4
            cfg.figname=[dirs.figdir, filesep, 'ts_det_all_l',tmp.lyear_str, '_', num2str(xpoint1), 'E_', num2str(xpoint2), 'E_', ...
                num2str(ypoint1), 'N_', num2str(ypoint2), 'N_', tmp.varname, '.tif'];
        end
            print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    end


end





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

% plot(squeeze(data.diatChl_model_l04(grids.id_w,grids.id_s,:)) - squeeze(data.diatChl_lens2_l04(grids.id_w,grids.id_s,:)))
% hold on
% plot(squeeze(data.diatChl_assm(grids.id_w,grids.id_s,:))- squeeze(data.diatChl_lens2_l04(grids.id_w,grids.id_s,:)))
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


