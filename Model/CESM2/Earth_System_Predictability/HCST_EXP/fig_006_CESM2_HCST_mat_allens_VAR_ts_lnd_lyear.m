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

cfg.vars = {'NPP', 'GPP', 'TWS', 'TOTVEGC'};

% cfg.vlayer=1; % surface, vertical slice
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
        grid.tarea=ncread(tmp.gridname, 'AREA')/1000000.0; %(m^2 -> km^2)
        grid.tarea_60=grid.tarea; grid.tarea_60(grid.tlat>60 | grid.tlat<-60)=NaN;
end

grid.nlon=size(grid.tlong,1);
grid.nlat=size(grid.tlat,2);
% grid.ntime=cfg.proj_year.*12;


sta_lonlat = {[20, 15], [20, 80], [40, 3], [50, 9], [62, 32], ...
              [50, 65], [70, 42], [65, 60], [82, 48], [82, 57], ...
              [105, 57], [110, 41], [120, 50], [122, 67], [130, 60], ...
              [120, -30], [120, -21], [131, -30], [131, -20], [140,-30], ...
              [147, -22], [140, -8], [200, 65], [220, 61], [230, 55], ...
              [248, 56], [248, 30], [252, 36], [260, 38], [270, 75], ...
              [270, 52], [285, -5], [287, 5], [292, 5], [293, -42], ...
              [296, -35], [308, -12], [320, -9], [320, 65], [338, 75], ...
              [340, 65], [350, 16], [360, 15], [360, 28], [10, 15], [8, 30], ...
              [70, -75], [125, -80], [240, -80], [340, -80], [360, -80], [18, 78]};
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

        xpoint=sta_lonlat{stai}(1);
        ypoint = sta_lonlat{stai}(2);
        [grid.id_w, grid.id_e, grid.id_s, grid.id_n] = Func_0012_findind_Y(2, [xpoint, ypoint], ...
                    grid.tlong, ...
                    grid.tlat, 'CESM2'); % find valid lon, lat index near station
        
        
        plot(cfg.iyears+lyear,squeeze(data.([tmp.varname, '_model_l', tmp.lyear_str])(grid.id_w,grid.id_s,:)), 'linewidth', 2)
        hold on
        plot(cfg.iyears+lyear,squeeze(data.([tmp.varname, '_assm'])(grid.id_w,grid.id_s,:)), 'linewidth', 2)
        plot(cfg.iyears+lyear,squeeze(data.([tmp.varname, '_lens2_l', tmp.lyear_str])(grid.id_w,grid.id_s,:)), 'linewidth', 2)
        hold off
        legend ('HCST', 'ASSM', 'LENS2', 'Location', 'Northwest')
        title(['l', tmp.lyear_str, ', ', tmp.varname, ', ', num2str(xpoint),'E, ', num2str(ypoint), 'N'])
        grid minor

        set(gca, 'fontsize', 20)

        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_time_series', filesep, 'l',tmp.lyear_str];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'ts_all_l',tmp.lyear_str, '_', num2str(xpoint), 'E_', num2str(ypoint), 'N_', tmp.varname, '.tif'];
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


