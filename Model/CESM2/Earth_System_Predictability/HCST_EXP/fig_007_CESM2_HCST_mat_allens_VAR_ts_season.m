% %  Created 12-Apr-2023 by Yong-Yub Kim
clc; clear all; close all;
warning off;

%% set path
[error_status, tmp.hostname] = system('hostname');
tmp.hostname=tmp.hostname(1:end-1);
switch tmp.hostname
    case 'Yong-Yubs-iMac-Pro.local'
        tmp.dropboxpath = '/Volumes/kyy_raid/kimyy/Dropbox';
    case 'Yong-Yubui-MacBookPro.local'
        tmp.dropboxpath = '/Users/kimyy/Dropbox';
    case {'da1', 'da2', 'da3', 'da4'}
        tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
end
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

%% model configuration
% cfg.var='TS'; %SST PRECT PSL TS SSH sumChl
% cfg.vars = {'SSH', 'sumChl'};
% % cfg.vars = {'SST', 'PRECT', 'PSL', 'TS'};
% cfg.vars = {'SSH'};
% cfg.vars = {'SST', 'PRECT', 'PSL', 'TS', 'sumChl', 'photoC_TOT_zint'};
% cfg.vars = { 'sumChl', 'photoC_TOT_zint',  'SSH'};
% cfg.vars = {'diatChl', 'diazChl', 'spChl'};

cfg.vars={'NO3', 'PO4', 'SiO3', 'Fe', ...
    'diat_light_lim_Cweight_avg_100m', 'diat_P_lim_Cweight_avg_100m', ...
    'diat_SiO3_lim_Cweight_avg_100m', 'diat_N_lim_Cweight_avg_100m', ...
    'diatChl', 'diazChl', 'spChl', 'diatC', 'diazC', 'spC', 'SSH', 'BSF', ...
    'PAR_avg',  ...
    'diat_agg_zint_100m', 'diat_Fe_lim_Cweight_avg_100m', ...
    'diat_Fe_lim_surf',  'diat_light_lim_surf', ...
    'diat_loss_zint_100m',  'diat_N_lim_surf', ...
     'diat_P_lim_surf',  ...
    'diat_SiO3_lim_surf', 'diaz_agg_zint_100m', 'diaz_Fe_lim_Cweight_avg_100m', ...
    'diaz_Fe_lim_surf', 'diaz_light_lim_Cweight_avg_100m', 'diaz_light_lim_surf', ...
    'diaz_loss_zint_100m', 'diaz_P_lim_Cweight_avg_100m', 'diaz_P_lim_surf', 'dustToSed', ...
    'graze_diat_zint_100m', 'graze_diat_zoo_zint_100m', 'graze_diaz_zint_100m', ...
    'graze_diaz_zoo_zint_100m', 'graze_sp_zint_100m', 'graze_sp_zoo_zint_100m', ...
    'LWUP_F', 'O2_ZMIN_DEPTH', 'photoC_diat_zint_100m', 'photoC_diaz_zint_100m', ...
    'photoC_sp_zint_100m', ...
    'sp_agg_zint_100m', 'sp_Fe_lim_Cweight_avg_100m', 'sp_Fe_lim_surf', ...
    'sp_light_lim_Cweight_avg_100m', 'sp_light_lim_surf', 'sp_loss_zint_100m', ...
    'sp_N_lim_Cweight_avg_100m', 'sp_N_lim_surf', 'sp_P_lim_Cweight_avg_100m', ...
    'sp_P_lim_surf', 'zoo_loss_zint_100m', 'photoC_TOT_zint_100m'};



% cfg.vars={'graze_diaz_zint_100m', 'zoo_loss_zint_100m', 'photoC_TOT_zint_100m'};
% cfg.vars={'photoC_TOT_zint_100m'};
% cfg.vars = {'SST', 'PRECT', 'PSL', 'TS'};
% cfg.vars = {'PAR_avg', 'NO3', 'PO4', 'SiO3', 'Fe'};
% cfg.vars = {'NO3'};
% cfg.vars = {'diatC', 'diazC', 'spC'};
% cfg.vars = {'SSH', 'BSF'};
% cfg.vars = {'photoC_diat_zint_100m', 'photoC_diaz_zint_100m', 'photoC_sp_zint_100m', 'photoC_TOT_zint_100m'};
% % 'graze_diaz_zoo_zint_100m', 
% cfg.vars = {'diatChl', 'diazChl', 'spChl', 'diatC', 'diazC', 'spC', 'SSH', 'BSF'};
% cfg.vars = {'SSH'};
% cfg.vars = {'photoC_sp_zint_100m'};
% cfg.vars = {'diazC'};
cfg.vars = {'TEMP', 'UVEL', 'VVEL', 'WVEL', 'SSH', 'PAR_avg', 'NO3', 'PO4', 'SiO3', 'Fe', ...
    'SSH', 'BSF', 'diatC', 'diazC', 'spC', ...
    'diatChl', 'diazChl', 'spChl', 'HBLT'};
% cfg.vars = {'sp_P_lim_Cweight_avg_100m'};
% cfg.vars = {'graze_sp_zint_100m'};
% cfg.vars = {'diat_N_lim_Cweight_avg_100m'};
% cfg.vars = {'diat_N_lim_Cweight_avg_100m'};

% cfg.vars = {'photoC_diat_zint_100m', ...
%     'photoC_diaz_zint_100m', 'photoC_sp_zint_100m', 'photoC_TOT_zint_100m', 'diatChl', 'diazChl', 'spChl', 'diat_agg_zint_100m', 'diat_Fe_lim_Cweight_avg_100m', 'diat_light_lim_Cweight_avg_100m', ...
%     'diat_loss_zint_100m', 'diat_N_lim_Cweight_avg_100m', 'diat_P_lim_Cweight_avg_100m', ...
%     'diat_SiO3_lim_Cweight_avg_100m', 'diaz_agg_zint_100m', 'diaz_Fe_lim_Cweight_avg_100m', ...
%     'diaz_light_lim_Cweight_avg_100m', 'diaz_loss_zint_100m', 'diaz_P_lim_Cweight_avg_100m', ...
%     'sp_agg_zint_100m', 'sp_Fe_lim_Cweight_avg_100m', 'sp_light_lim_Cweight_avg_100m', ...
%     'sp_loss_zint_100m', 'sp_N_lim_Cweight_avg_100m', 'SSH'};
% cfg.vars = { 'diatC', 'Fe', 'spC', 'NO3', 'PO4', 'SiO3', 'TEMP', 'UVEL', 'VVEL', 'WVEL', 'PAR_avg', 'BSF', 'diazC', 'sp_P_lim_Cweight_avg_100m' };
% cfg.vars = {'TS', 'SST', 'PRECT', 'PSL'};
% cfg.vars = {'HBLT'};
% cfg.vars={'TEMP', 'UVEL', 'VVEL', 'WVEL'};
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

cfg.vars = {'sumChl'};

cfg.vlayer=1:10; % 10layer. don't put more than 15
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

dirs.hcstmatroot=['/Users/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', cfg.var];
vstr=['v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i')];
dirs.figroot=[tmp.kimyypath, '/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', cfg.var, filesep, vstr];




cfg.iyears=cfg.obs_iyears;
cfg.gnm='f09_g17';
% cfg.proj_year=1;
cfg.season = {'AMJ', 'JAS', 'OND', 'JFM'};


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
tmp.gridname = ['/Users/kimyy/Model/CESM2/ESP/HCST_EXP/archive_transfer/', cfg.comp, '/grid.nc'];
tmp.maskname = '/Users/kimyy/Model/CESM2/ESP/HCST_EXP/archive_transfer/ocn/RECCAP2_region_masks_all_v20210412_POP2_grid.nc';

% grid.region_mask=ncread(tmp.gridname, 'REGION_MASK'); 
% grid.ocean_mask=NaN(size(grid.region_mask));
% grid.ocean_mask(grid.region_mask>0)=1;
% grid.tarea = ncread(tmp.gridname, 'TAREA');

switch cfg.comp
    case 'ocn'
        grid.tlong=ncread(tmp.gridname, 'TLONG');
        grid.tlat=ncread(tmp.gridname, 'TLAT');
        grid.mask_ocn=ncread(tmp.maskname, 'open_ocean');
        grid.mask_ocn(grid.mask_ocn<-10e10)=NaN;
        grid.mask_ocn=grid.mask_ocn./grid.mask_ocn;
        grid.tarea=ncread(tmp.gridname, 'TAREA')/1000000.0; %(m^2 -> km^2)
        grid.tarea_60=grid.tarea; grid.tarea_60(grid.tlat>60 | grid.tlat<-60)=NaN;
    case 'atm'
        grid.lon=ncread(tmp.gridname, 'lon');
        grid.lat=ncread(tmp.gridname, 'lat');
        [grid.tlat grid.tlong]=meshgrid(grid.lat, grid.lon);
        grid.tarea=ncread(tmp.gridname, 'AREA')/1000000.0; %(m^2 -> km^2)
        grid.tarea_60=grid.tarea; grid.tarea_60(grid.tlat>60 | grid.tlat<-60)=NaN;
end

grid.nlon=size(grid.tlong,1);
grid.nlat=size(grid.tlat,2);
% grid.ntime=cfg.proj_year.*12;


sta_lonlat = {[345, 75], [300, 85], [260, 85], [220, 80], [180, 80], ...
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

% sta_lonlat = { [238,32], [160, 47], [137, 3], [137, 5], [137,10], ...
%                 [137, 15], [137, 20], [137,25], [137,30], [137, 34], ...
%                 [297, 24], [297, 30], [297, 35]};
% sta_lonlat = { [60,-30], [90, 15], [140, 15], [200, 20], [200,50], ...
%                 [250, 20], [160, -10], [180,-10], [210,-25], [300, 25], ...
%                 [330, -10], [360, -10], [330, -30]};
% sta_lonlat = { [300, 40]};
% CalCOFI : 238, 32
% K2 : 160, 47`
% JMA : 137, 3 ~ 34
% BATS : 294 ~ 300, 24 ~ 35
S = shaperead('landareas.shp');

%% read & plot data
tmp.varname=cfg.var;

clear tmp.ydata tmp.ydata_lens2 tmp.ydata_obs tmp.ydata_assm
for lss=1:length(cfg.season)
    tmp.season=cfg.season{lss};
    tmp.mons = f_season_mons(tmp.season);
    fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_', tmp.season, '.mat'];
    load(fig_cfg.mat_name, 'data')
    
    tmp.mon=tmp.mons(2);
    if tmp.mon>=12 && mod(tmp.mon,12)~=0
        tmp.mon_str=num2str(tmp.mon-12*floor(tmp.mon/12), '%02i');
        tmp.fy = cfg.iyears+floor(tmp.mon/12);
%         tmp.fy_str=num2str(iyear+floor(tmp.mon/12), '%04i');
    elseif tmp.mon>=12 && mod(tmp.mon,12) ==0
        tmp.mon_str='12';
        tmp.fy=cfg.iyears+floor(tmp.mon/12)-1;
%         tmp.fy_str=num2str(iyear+floor(tmp.mon/12)-1, '%04i');
    else
        tmp.mon_str=num2str(tmp.mon, '%02i');
        tmp.fy=cfg.iyears;
%         tmp.fy_str=num2str(iyear, '%04i');
    end


    for stai=1:length(sta_lonlat)
        fig_h = figure('name','ts','visible','off');

        xpoint=sta_lonlat{stai}(1);
        ypoint = sta_lonlat{stai}(2);
        [grid.id_w, grid.id_e, grid.id_s, grid.id_n] = Func_0012_findind_Y(2, [xpoint, ypoint], ...
                    grid.tlong, ...
                    grid.tlat, 'CESM2'); % find valid lon, lat index near station
        
        
        plot(tmp.fy,squeeze(data.([tmp.varname, '_model_', tmp.season])(grid.id_w,grid.id_s,:)), 'linewidth', 2)
        hold on
        plot(tmp.fy,squeeze(data.([tmp.varname, '_assm_', tmp.season])(grid.id_w,grid.id_s,:)), 'linewidth', 2)
        plot(tmp.fy,squeeze(data.([tmp.varname, '_lens2_', tmp.season])(grid.id_w,grid.id_s,:)), 'linewidth', 2)
        hold off
        legend ('HCST', 'ASSM', 'LENS2', 'Location', 'Northwest')
        title([tmp.season, ', ', tmp.varname, ', ', num2str(xpoint),'E, ', num2str(ypoint), 'N'])
        grid minor

        set(gca, 'fontsize', 20)

        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_time_series', filesep, tmp.season];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'ts_all_',tmp.season, '_', num2str(xpoint), 'E_', num2str(ypoint), 'N_', tmp.varname, '.tif'];
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

function [mons] = f_season_mons(season)
    switch season
        case 'INI'
            mons = [1];
        case 'FMA'
            mons = [2,3,4];
            lyear = 0;
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
