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
% cfg.vars = {'TWS', 'RAIN'};
% cfg.vars = {'RAIN'};
% cfg.vars = {'TOTVEGC'};
% cfg.vars = {'TWS'};
% cfg.vars = {'RAIN'};
% cfg.vars = {'TS', 'TWS', 'SOILWATER_10CM'};
% cfg.vars = {'RAIN', 'TWS'};
% cfg.vars = {'TS', 'TWS', 'FIRE', 'SOILWATER_10CM'};
% cfg.vars = {'SST'};
% cfg.vars = {'photoC_TOT_zint_100m'};
cfg.vars = {'PSL'};

cfg.vars={'PRECT', 'SST', 'TS', 'PSL', 'TWS', 'SOILWATER_10CM', 'COL_FIRE_CLOSS', 'TLAI','FAREA_BURNED'};
cfg.vars={'SOILWATER_10CM'};
cfg.vars={'PRECT', 'SST', 'TS', 'PSL', 'TWS', 'SOILWATER_10CM', 'COL_FIRE_CLOSS', 'TLAI','FAREA_BURNED'};
cfg.vars={'photoC_TOT_zint_100m'};
% cfg.vars={'SSH'};
% cfg.vars = {'SST', 'TS', 'PSL', 'PRECT'};
% cfg.vars={'HBLT', 'HMXL'}
% cfg.vars={'TS'};
cfg.vars={'SOILWATER_10CM', 'TLAI', 'FAREA_BURNED', 'COL_FIRE_CLOSS', 'TWS', 'GPP'};
cfg.vars={'GPP'};
cfg.vlayer=1; % surface, vertical slice
% cfg.vlayer=1:10; % 10layer. don't put more than 15
cfg.vlayer_1st=min(cfg.vlayer);
cfg.vlayer_cnt=max(cfg.vlayer)-cfg.vlayer_1st+1;

cfg.ly_s=1;
cfg.ly_e=3;


for vari=1:length(cfg.vars)

cfg.var=cfg.vars{vari};
cfg.obs_name=f_obs_name(cfg.var);
cfg.obs_fname_mid=f_obs_name_mid(cfg.var);
cfg.obs_varname=f_obs_varname(cfg.var);
cfg.comp=Func_0025_CESM2_cmpname_var(cfg.var);
cfg.obs_fname_module=f_obs_fname_module(cfg.comp);
cfg.obs_iyears=1960:2020-cfg.ly_e+1;

disp(cfg.var);
tic;

if strcmp(cfg.var, 'PRECT')==1
    tlag=-1;
else
    tlag=0;
end

% dirs.hcstroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/HCST_EXP/archive/', cfg.comp, '/', cfg.var];
% dirs.obsroot=['/Volumes/kyy_raid/kimyy/Observation/', cfg.obs_name, '/monthly_reg_', cfg.obs_fname_module(2:4)];
% dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', cfg.var];
% dirs.lens2root=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/LENS2/archive_analysis/', cfg.comp, '/', cfg.var];
% dirs.assmroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/ASSM_EXP/archive_analysis/', cfg.comp, '/', cfg.var];

dirs.hcstmatroot=[tmp.kimyypath, '/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', cfg.var];
vstr=['v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i')];
dirs.figroot=[tmp.kimyypath, '/Figure/CESM2/ESP/HCST_EXP/archive_anomaly/', cfg.comp,'/', cfg.var, filesep, vstr];



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


sta_lonlat = {[2, -2], [20, 15], [22, -30], [22, 0], [20, 80], [30, -10], [40, 3], [50, 9], [62, 32], ...
              [50, 65], [70, 42], [65, 60], [78, 15], [82, -18], [82, 48], [82, 57], [90, 25], ...
              [105, 57], [110, 41], [112, 30], [120, 50], [122, 67], [130, 60], [140, 10], ...
              [120, -30], [120, -21], [131, -30], [131, -20], [140,-30], ...
              [147, -22], [140, -8], [200, 65], [220, 61], [230, 55], ...
              [248, 56], [248, 30], [252, 36], [260, 38], [270, 75], ...
              [270, 52], [285, -5], [287, 5], [292, 5], [293, -42], ...
              [296, -35], [302, -12], [308, -12], [320, -9], [320, 65], [338, 75], ...
              [340, 65], [350, 16], [359, 15], [359, 28], [10, 15], [8, 30], ...
              [70, -75], [125, -80], [240, -80], [340, -80], [359, -80], [18, 78], ...
              [345, 75], [300, 85], [260, 85], [220, 80], [160,20], ...
              [160, 35], [160, 40], [180, 35], [180, 40], [180, 80], ...
                [140, 85], [100, 85], [40, 75], [359, 65], [310, 55], ...
                [230, 50], [200, 50], [180, 50], [150, 55], [180, 35], ...
                [120, -40], [140, 33], [124, 35], [240,20], [180, 20], [150, 15], ...
                [180, 10], [200, 0], [260, 0], [160, 0], [280, -20], [240, -20], ...
                [200, -20], [160, -20], [160, -30], [200, -40], [200, -50], ...
                [110, -10], [90, 15], [60, 10], [42, -10], [40, -30], ...
                [40, -40], [40, -50], [340, 40], [320, 30], [300, 30], ...
                [240,35], [285, 30], [340, 20], [300, 20], [340, 10], [320, 10], ...
                [270, 25], [359, 0], [320, 0], [10, -10], [330, -10], ...
                [330, -30], [5, 40], [60, -65], [60, -25], [180, -70], ...
                [240, -70], [340, -50], [150, 30], [5, 35], ...
                [320, -70], [200, 21], [238, 32], [200, 21], ...
                [5, 60], [151, 34], [190,0], [160, 47], ...
                [137, 3], [137, 5], [137,10], [137, 15], [137, 20], ...
                [137,25], [137,30], [137, 34],[297, 24], [297, 30], ...
                 [297, 35], [300, 40], [60,-30], [90, 15], [140, 15], ...
                 [200, 20], [200,50], [250, 20], [160, -10], [180,-10],...
                 [210,-25], [300, 25], [320, -70], [330, -10], [359, -10], [330, -30], [340, -30], ...
                 [240, -70], [250, -50], [310, 50], [330, 60], [303, -70]};

%% caution!!! if you use regional mean, spread (regional mean of spread of each grids) is wrong, spread should not be used
%%% lnd
sta_lonlat = {[10 40 -20 0], [20 40 -20 5], [20 40 -40 -20], [20 40 -40 -10], [20 40 5 20], [20 50 -40 -10], [30 50 -10 5], ...
     [20 40 35 48], [20 60 20 60], [20 180 55 70], [20 90 50 55], [30 40 -10 0], [32 50 -10 10], [20 65 40 60], ...
     [60 100 8 30], [60 100 40 60], [70 130 50 70], ...
     [80 100 30 45], [80 120 20 55], [90 120 10 30], [90 160 -10 20], [100 160 -10 15], [100 160 60 80], [100 130 -40 -10], ...
     [100 140 20 40], [100 140 -10 10], [110 160 55 70], [235 250 20 45], ...
     [100 160 -40 -10], [100 140 -40 -10], [130 160 60 75], [130 160 -40 -10], [140 160 -42 -20], [140 160, -40 -15], ...
     [100 140 50 65], [190 300 50 70], [220 250 30 40], [220 260 50 70], [230 268 15 45], ...
     [240 270 10 30], [240 270 10 35], [240 270 20 35], [240 280 20 35], [238 260 20 45], [238 280 10 45], [250 270 25 35], ...
     [265 285 30 40],[280 300 0 10], [280 310 -50 -30], [280 305 40 50], [290 320 -20 -5], ...
     [280 325 -10 12], [300 325 -10 0], [300 325 -10 10], [280 300 -60 -30], [278 320 -60 -20], [300 340 15 30], ...
     [300 340 60 85], [310 320 -20 0], [340 350 0 35], [340 359 0 15], [340 359 10 20], [340 359 20 30], ...
     [1 20 -40 10], [1 20 20 35], [1 20 40 60], [1 40 -40 -20], [1 40 -40 -15], ...
     [100 130 -30 -10], [240 270 30 40], [280 300, -10 5], [305 330, -15 0], [340 359, 5 15], ...
    [20 60 -10 10], [35 50 -10 5], [20 40 -40, -22], [100 145 30 70], [90 110 10 25], [285 305 -10 10]};
% %      ... % ~ocn
% %      sta_lonlat = { ...
% %      [1 40 15 25], [1 80 -70 55], [1 120 -80 -60], [1 359 -40 40],   ...
% %      [20 40 30 40], [160 280 -70 -50], [30 60 -25 -10], [30 120 -25 -15], ...
% %      [80 120 -35 -15], [60 100 -10 5], [120 160 5 20], [120 250 15 20], ...
% %      [140 170 30 50], [140 220 25 30], [140 220 25 35], [140 180 -20 -5], [140 220 -35 -10], ...
% %      [160 220 50 60], [140 240 -10 50], [190 240 15 60], [190 280 -30 20], ... 
% %      [200 240 30 60], [235 260 10 30], [240 300 -15 0], [240 280 -10 5], [290 340 20 30], ...
% %      [290 310 5 15], [320 340 10 20], [280 359 50 80], [50 90 -40 -15], [40 120 -40 -20], ...
% %      [40 80 -35 -25], [40 120 -40 0], [40 120 -30 20], [80 100 0 20], ...
% %      [100 280 -30 20], [100 120 -40 -20],  [140 200 15 22],  [130 162 5 13], ...
% %      [180 260 -10 30], [220 280 -70 -40], [280 320 35 50], [280 340 10 30], ...
% %      [280 359 -20 30], [300 359 60 90], [305 345 15 25], [310 359 -60 -30], [320 359 -20 10], ...
% %      [340 360 -5 5], [340 360 -10 0], [300 360 -40 -20], [330 350 55 65], [350 20 60 70], [1 20 -25 5]};

%  sta_lonlat = {[140 200 15 22]};
% sta_lonlat = {[130 162 5 13]};
% sta_lonlat = {[235 250 20 45]};
% 
% sta_lonlat = {[190 0], [300 30],[160 20]};
% sta_lonlat = {[200,80]};
% sta_lonlat = {[320,60]};
% sta_lonlat = {[300 340 50 65]};
% sta_lonlat = {[140 220 15 30]};
% sta_lonlat = {[160 190 -20 -10]};
% sta_lonlat = {[310 330 40 50]};
% sta_lonlat = {[325 50]};

% sta_lonlat = {[0 360 -90 90], [190 240 -5 5]};
% 
% 
% % sta_lonlat = {[205 0]};
% % sta_lonlat = {[50 9]};
% % sta_lonlat = {[120 40]};
% sta_lonlat = {[50 9], [30, -15], [70 30], [95 40], [130, -30], [150, -30], [240 40], [250 35], [295, -10], [310, -10], [358, 15], [310, 75]}; %TWS
% %  sta_lonlat = {[310 72]};
%  sta_lonlat = {[110 -3]};
% %  sta_lonlat = {[242 35]};
% %  sta_lonlat = {[20 40 -40 -10], [60 100 40 60], [100 160 -10 15], [100 160 -40 -10], [220 250 30 40], [300 325 -10 10], [300 340 60 85]};
%  sta_lonlat = {[40 50 0 12]};
%  sta_lonlat = {[0 45 10 20]};
%  sta_lonlat = {[90 140 50 75]};
%  sta_lonlat = {[130 150 55 75]};
%   sta_lonlat = {[300 85]};
% sta_lonlat = {[240,35]};
% sta_lonlat = {[160, 20], [160,35], [240,35]};
% sta_lonlat = {[320,-70]};

sta_lonlat= {[145 -25]};
% CalCOFI : 238, 32
% K2 : 160, 47`
% JMA : 137, 3 ~ 34
% BATS : 294 ~ 300, 24 ~ 35
S = shaperead('landareas.shp');

%% read & plot data
tmp.varname=cfg.var;

clear tmp.ydata tmp.ydata_lens2 tmp.ydata_obs tmp.ydata_assm

lyear=0;
tmp.lyear_str=[num2str(cfg.ly_s-1,'%02i'),'_',num2str(cfg.ly_e-1,'%02i')];

% for lyear=0:cfg.proj_year-1
%     tmp.lyear_str=num2str(lyear, '%02i');
    fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l', tmp.lyear_str, 'y.mat'];
    load(fig_cfg.mat_name, 'data', 'data2')
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
            [grids.id_w, grids.id_e, grids.id_s, grids.id_n] = Func_0012_findind_Y(1, [xpoint1,xpoint2, ypoint1,ypoint2], ...
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
        

        if (grids.id_w+grids.id_e+grids.id_s+grids.id_n)==4
            tmp.HCST_mean=NaN(1,size(data2.([tmp.varname, '_model_ano_l', tmp.lyear_str]),3));
            tmp.ASSM_mean=tmp.HCST_mean;
            tmp.LENS2_mean=tmp.HCST_mean;
            tmp.HCST_lower=tmp.HCST_mean;
            tmp.HCST_upper=tmp.HCST_mean;
            tmp.ASSM_lower=tmp.HCST_mean;
            tmp.ASSM_upper=tmp.HCST_mean;
            tmp.LENS2_lower=tmp.HCST_mean;
            tmp.LENS2_upper=tmp.HCST_mean;
            tmp.OBS = tmp.HCST_mean;
        else
            tmp.HCST_mean = Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_model_ano_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:), grids.tlong_cut, grids.tlat_cut);
            tmp.ASSM_mean = Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_assm_ano_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:), grids.tlong_cut, grids.tlat_cut);
            tmp.LENS2_mean = Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_lens2_ano_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:), grids.tlong_cut, grids.tlat_cut);
                    
            if length(sta_lonlat{stai})==2
                tmp.HCST_lower= Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_model_ano_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:) - ...
                    squeeze(data.([tmp.varname, '_model_stde_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,1+4-cfg.ly_s:end))/2.0, grids.tlong_cut, grids.tlat_cut);
                tmp.HCST_upper= Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_model_ano_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:) + ...
                    squeeze(data.([tmp.varname, '_model_stde_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,1+4-cfg.ly_s:end))/2.0, grids.tlong_cut, grids.tlat_cut);
                tmp.ASSM_lower= Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_assm_ano_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:) - ...
                    squeeze(data.([tmp.varname, '_assm_stde'])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,1+4-cfg.ly_s:end))/2.0, grids.tlong_cut, grids.tlat_cut);
                tmp.ASSM_upper= Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_assm_ano_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:) + ...
                    squeeze(data.([tmp.varname, '_assm_stde'])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,1+4-cfg.ly_s:end))/2.0, grids.tlong_cut, grids.tlat_cut);
                tmp.LENS2_lower= Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_lens2_ano_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:) - ...
                    squeeze(data.([tmp.varname, '_lens2_stde_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,1+4-cfg.ly_s:end))/2.0, grids.tlong_cut, grids.tlat_cut);
                tmp.LENS2_upper= Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_lens2_ano_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:) + ...
                    squeeze(data.([tmp.varname, '_lens2_stde_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,1+4-cfg.ly_s:end))/2.0, grids.tlong_cut, grids.tlat_cut);
            end

            tmp.OBS = Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_obs_ano_l', tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:), grids.tlong_cut, grids.tlat_cut);
            
%             tmp.OBS = tmp.OBS/86400.0; % GPP

% % %             for loni=1:size(grids.tlong,1)
% % %                 for lati=1:size(grids.tlat,2)
% % %                     tmp.corr=corrcoef(tmp.HCST_mean, data2.([tmp.varname, '_obs_ano_l', tmp.lyear_str])(loni,lati,:), 'Rows', 'complete');
% % %                     tmp.corrmap(loni,lati)=tmp.corr(1,2);
% % %                 end
% % %             end
% % % %             pcolor(tmp.corrmap'); shading flat; colorbar;
% % %             fig_cfg.name_rgn = 'Glob';
% % %             fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
% % %         
% % %             fig_cfg.x_lim = [-180 180];
% % %             fig_cfg.y_lim = [-80 89];
% % %             fig_cfg.fig_size = [0,0,6,3.5].*2;
% % %             fig_cfg.ax_size = [0.3,0.7,5.4,2.7].*2;
% % %             fig_cfg.cb_size = [5.15,0.8,0.15,2.3].*2;
% % %             fig_cfg.title_pos = [0.5,0.93];
% % %             fig_cfg.p_lim =0.1;
% % %             fig_cfg.c_lim = [-1 1];
% % %             [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
% % %         
% % %             tmp.X=grids.tlong([end, 1:end],:);
% % %             tmp.Y=grids.tlat([end, 1:end],:);
% % %             tmp.C=tmp.corrmap;
% % %             tmp.C=tmp.C([end, 1:end],:);
% % %              [tmp.mean_corr, tmp.err] = ...
% % %             Func_0011_get_area_weighted_mean(tmp.corrmap, grids.tlong, grids.tlat);
% % % 
% % %             fig_cfg.fig_name=['l', tmp.lyear_str, ', corr_assm_Rdiff_ano_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
% % %             fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
% % %                 'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','on');
% % %             ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
% % %                 'fontname','freeserif'); 
% % % 
% % %             axis off; 
% % %             hold on;
% % %             setm(ax_m,'origin',[0,200],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
% % %             set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
% % %             text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
% % %             'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
% % %             'fontsize',14,'fontname','freeserif','interpreter','none')
% % %         
% % %             %% caxis & colorbar
% % %             caxis(ax_m, [-1 1]); 
% % %             colormap(fig_cfg.c_map);
% % %             cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
% % %             set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
% % %             title(cb,'R','fontsize',12);
% % %         
% % %             %% draw on ax_m
% % %             h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
% % %             shading flat;
% % %             geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
            


%             if strcmp(tmp.varname, 'photoC_TOT_zint_100m')==1 || strcmp(tmp.varname, 'photoC_TOT_zint')==1
                
                %% mean
%                 tmp.ASSM_order=Func_0011_get_area_weighted_mean(data.([tmp.varname, '_assm_clim'])(grids.id_w:grids.id_e,grids.id_s:grids.id_n), grids.tlong_cut, grids.tlat_cut);
%                 tmp.LENS2_order=Func_0011_get_area_weighted_mean(data.([tmp.varname, '_lens2_clim'])(grids.id_w:grids.id_e,grids.id_s:grids.id_n), grids.tlong_cut, grids.tlat_cut);
%                 tmp.OBS_order=Func_0011_get_area_weighted_mean(data.([tmp.varname, '_obs_clim'])(grids.id_w:grids.id_e,grids.id_s:grids.id_n), grids.tlong_cut, grids.tlat_cut);
%                 tmp.HCST_order=Func_0011_get_area_weighted_mean(data.([tmp.varname, '_model_clim_l',tmp.lyear_str])(grids.id_w:grids.id_e,grids.id_s:grids.id_n), grids.tlong_cut, grids.tlat_cut);
                
                %% std
                tmp.ASSM_order=std(tmp.ASSM_mean, 'omitnan');
                tmp.LENS2_order=std(tmp.LENS2_mean, 'omitnan');
                tmp.OBS_order=std(tmp.OBS, 'omitnan');
                tmp.HCST_order=std(tmp.HCST_mean, 'omitnan');

                tmp.OBS_norm=tmp.OBS./tmp.OBS_order;
                tmp.HCST_mean_norm=tmp.HCST_mean./tmp.HCST_order;
                tmp.LENS2_mean_norm=tmp.LENS2_mean./tmp.LENS2_order;
                tmp.ASSM_mean_norm=tmp.ASSM_mean./tmp.ASSM_order;

%                 tmp.OBS=tmp.OBS./std(tmp.OBS,'omitnan');
%                 tmp.HCST_mean = tmp.HCST_mean ./ std(tmp.HCST_mean,'omitnan');
%                 tmp.LENS2_mean = tmp.LENS2_mean ./ std(tmp.LENS2_mean,'omitnan');
%                 tmp.ASSM_mean = tmp.ASSM_mean ./ std(tmp.ASSM_mean,'omitnan');
%                 if length(sta_lonlat{stai})==2
%                     tmp.HCST_upper = tmp.HCST_upper ./ std(tmp.HCST_mean,'omitnan');
%                     tmp.HCST_lower = tmp.HCST_lower ./ std(tmp.HCST_mean,'omitnan');
%                     tmp.ASSM_upper = tmp.ASSM_upper ./ std(tmp.ASSM_mean,'omitnan');
%                     tmp.ASSM_lower = tmp.ASSM_lower ./ std(tmp.ASSM_mean,'omitnan');
%                     tmp.LENS2_upper = tmp.LENS2_upper ./ std(tmp.LENS2_mean,'omitnan');
%                     tmp.LENS2_lower = tmp.LENS2_lower ./ std(tmp.LENS2_mean,'omitnan');
%                 end
%             end
        end


%         tmp.OBS2 = Func_0011_get_area_weighted_mean(data.([tmp.varname, '_obs'])(grids.id_w:grids.id_e,grids.id_s:grids.id_n,:), grids.tlong_cut, grids.tlat_cut);

        switch tmp.varname
           case {'PRECT', 'RAIN', 'SNOW'}
               cfg.clim_ys=1965-cfg.ly_e+1;
               cfg.clim_ye=2019-cfg.ly_e+1;
           otherwise
               cfg.clim_ys=1965-cfg.ly_e+1;
               cfg.clim_ye=2020-cfg.ly_e+1;
        end
        cfg.len_t_y=length(cfg.iyears);
        cfg.clim_tlen = (cfg.clim_ys-1959):(cfg.clim_ye-(2020-cfg.ly_e+1))+cfg.len_t_y;
       
%         tmp.time= cfg.iyears+(cfg.ly_s+cfg.ly_e)/2-1;
        tmp.time= (cfg.clim_ys:cfg.clim_ye)+(cfg.ly_s+cfg.ly_e)/2-1;


        %% normal plot
        hold on
        if length(sta_lonlat{stai})==2
            fig_ts.LENS2_range=fill([tmp.time, flip(tmp.time)], ...
                [tmp.LENS2_lower(1:end-lyear); flip(tmp.LENS2_upper(1:end-lyear))], cmap_LENS2_b);
            fig_ts.LENS2_range.EdgeColor = 'none';
            set(get(get(fig_ts.LENS2_range,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            
            fig_ts.HCST_range=fill([tmp.time, flip(tmp.time)], ...
                [tmp.HCST_lower; flip(tmp.HCST_upper)], cmap_HCST_b);
            fig_ts.HCST_range.EdgeColor = 'none';
            set(get(get(fig_ts.HCST_range,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
            fig_ts.ASSM_range=fill([tmp.time, flip(tmp.time)], ...
                [tmp.ASSM_lower(1:end-lyear); flip(tmp.ASSM_upper(1:end-lyear))], cmap_ASSM_b);
            fig_ts.ASSM_range.EdgeColor = 'none';
            set(get(get(fig_ts.ASSM_range,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end

        plot(tmp.time, tmp.LENS2_mean(1:end-lyear), 'linewidth', 2, 'color', cmap_LENS2)
        plot(tmp.time, tmp.ASSM_mean(1:end-lyear), 'linewidth', 2, 'color', cmap_ASSM)
        plot(tmp.time, tmp.HCST_mean, 'linewidth', 2, 'color', cmap_HCST)
        plot(tmp.time, tmp.OBS(1:end-lyear), 'linewidth', 2, 'color', cmap_OBS)

        hold off
        
%         plot(cfg.iyears+lyear,squeeze(data.([tmp.varname, '_model_l', tmp.lyear_str])(grids.id_w,grids.id_s,:)), 'linewidth', 2)
%         hold on
%         plot(cfg.iyears+lyear,squeeze(data.([tmp.varname, '_assm'])(grids.id_w,grids.id_s,:)), 'linewidth', 2)
%         plot(cfg.iyears+lyear,squeeze(data.([tmp.varname, '_lens2_l', tmp.lyear_str])(grids.id_w,grids.id_s,:)), 'linewidth', 2)
%         hold off
%         legend ('LENS2', 'ASSM', 'HCST', 'OBS', 'Location', 'Northwest', 'Orientation', 'Horizontal')
        legend ('LENS2', 'ASSM', 'HCST', 'OBS', 'Location', 'Southoutside', 'Orientation', 'Horizontal')
        
        if length(sta_lonlat{stai})==2
            title(['l', tmp.lyear_str, ', ', tmp.varname, ', ', num2str(xpoint),'E, ', num2str(ypoint), 'N'])
        elseif length(sta_lonlat{stai})==4
            title(['l', tmp.lyear_str, ', ', tmp.varname, ', ', num2str(xpoint1),'E~', num2str(xpoint2),'E, ', ...
                num2str(ypoint1), 'N~',num2str(ypoint2), 'N'])
        end
        grid minor
        xlim([1960 2025])
        set(gca, 'fontsize', 20)
        
% % % % %         %% NPGO correlations;
% % % % % 
% % % % %         [NPGO_year, NPGO_month, NPGO]=Func_0032_get_NPGO(3);
% % % % % %         NPGO(NPGO_year)
% % % % %         NPGO_year(NPGO_year<min(tmp.time))=NaN;
% % % % %         NPGO_year(NPGO_year>max(tmp.time))=NaN;
% % % % %         [tmp.skill_NPGO_assm, tmp.skill_NPGO_assm_p]=corrcoef(NPGO(isfinite(NPGO_year)), tmp.ASSM_mean(1:end-lyear), 'Rows', 'complete');
% % % % %         if (tmp.skill_NPGO_assm_p>0.1) tmp.skill_NPGO_assm=NaN(2,2); end
% % % % %         [tmp.skill_NPGO_hcst, tmp.skill_NPGO_hcst_p]=corrcoef(NPGO(isfinite(NPGO_year)), tmp.HCST_mean(1:end-lyear), 'Rows', 'complete');
% % % % %         if (tmp.skill_NPGO_hcst_p>0.1) tmp.skill_NPGO_hcst=NaN(2,2); end
% % % % % 
% % % % %         plot(tmp.time, NPGO(isfinite(NPGO_year)))
% % % % %         hold on
% % % % %         plot(tmp.time, tmp.ASSM_mean(1:end-lyear));
% % % % %         hold off
% % % % %         
% % % % %         for loni=1:size(grids.tlong,1)
% % % % %             for lati=1:size(grids.tlat,2)
% % % % %                 tmp.corr=corrcoef(NPGO(isfinite(NPGO_year)), data2.([tmp.varname, '_assm_ano_l', tmp.lyear_str])(loni,lati,:), 'Rows', 'complete');
% % % % %                 tmp.corrmap(loni,lati)=tmp.corr(1,2);
% % % % %             end
% % % % %         end
% % % % %         pcolor(tmp.corrmap'); shading flat; colorbar;
% % % % % 
% % % % % 
% % % % %         %% AMO correlations
% % % % %         [AMO_year, AMO_month, AMO] = Func_0033_get_AMO(2);
% % % % %         AMO_year(AMO_year<min(tmp.time))=NaN;
% % % % %         AMO_year(AMO_year>max(tmp.time))=NaN;
% % % % %         [tmp.skill_AMO_assm, tmp.skill_AMO_assm_p]=corrcoef(AMO(isfinite(AMO_year)), tmp.ASSM_mean(1:end-lyear), 'Rows', 'complete');
% % % % %         if (tmp.skill_AMO_assm_p>0.1) tmp.skill_AMO_assm=NaN(2,2); end
% % % % %         [tmp.skill_AMO_hcst, tmp.skill_AMO_hcst_p]=corrcoef(AMO(isfinite(AMO_year)), tmp.HCST_mean(1:end-lyear), 'Rows', 'complete');
% % % % %         if (tmp.skill_AMO_hcst_p>0.1) tmp.skill_AMO_hcst=NaN(2,2); end
% % % % % 
% % % % %         plot(tmp.time, AMO(isfinite(AMO_year)))
% % % % %         hold on
% % % % %         plot(tmp.time, tmp.ASSM_mean(1:end-lyear));
% % % % %         hold off
% % % % %         
% % % % %         for loni=1:size(grids.tlong,1)
% % % % %             for lati=1:size(grids.tlat,2)
% % % % %                 tmp.corr=corrcoef(AMO(isfinite(AMO_year)), data2.([tmp.varname, '_assm_ano_l', tmp.lyear_str])(loni,lati,:), 'Rows', 'complete');
% % % % %                 tmp.corrmap(loni,lati)=tmp.corr(1,2);
% % % % %             end
% % % % %         end
% % % % %         pcolor(tmp.corrmap'); shading flat; colorbar;
% % % % % 
% % % % % 
% % % % %         %% PDO correlations
% % % % %         [PDO_year, PDO_month, PDO] = Func_0034_get_PDO(3);
% % % % %         PDO_year(PDO_year<min(tmp.time))=NaN;
% % % % %         PDO_year(PDO_year>max(tmp.time))=NaN;
% % % % %         [tmp.skill_PDO_assm, tmp.skill_PDO_assm_p]=corrcoef(PDO(isfinite(PDO_year)), tmp.ASSM_mean(1:end-lyear), 'Rows', 'complete');
% % % % %         if (tmp.skill_PDO_assm_p>0.1) tmp.skill_PDO_assm=NaN(2,2); end
% % % % %         [tmp.skill_PDO_hcst, tmp.skill_PDO_hcst_p]=corrcoef(PDO(isfinite(PDO_year)), tmp.HCST_mean(1:end-lyear), 'Rows', 'complete');
% % % % %         if (tmp.skill_PDO_hcst_p>0.1) tmp.skill_PDO_hcst=NaN(2,2); end
% % % % % 
% % % % %         plot(tmp.time, PDO(isfinite(PDO_year)))
% % % % %         hold on
% % % % %         plot(tmp.time, tmp.ASSM_mean(1:end-lyear));
% % % % %         hold off
% % % % %         
% % % % %         for loni=1:size(grids.tlong,1)
% % % % %             for lati=1:size(grids.tlat,2)
% % % % %                 tmp.corr=corrcoef(PDO(isfinite(PDO_year)), data2.([tmp.varname, '_assm_ano_l', tmp.lyear_str])(loni,lati,:), 'Rows', 'complete');
% % % % %                 tmp.corrmap(loni,lati)=tmp.corr(1,2);
% % % % %             end
% % % % %         end
% % % % %         pcolor(tmp.corrmap'); shading flat; colorbar;

% % %         %% Clim ind correlations ENSO, IPO, ....
% % %         [IPO_year, IPO_month, IPO] =Func_0032_get_clim_indices('SAM', 3);
% % %         IPO_year(IPO_year<min(tmp.time))=NaN;
% % %         IPO_year(IPO_year>max(tmp.time))=NaN;
% % %         [tmp.skill_IPO_assm, tmp.skill_IPO_assm_p]=corrcoef(IPO(isfinite(IPO_year)), tmp.ASSM_mean(1:end-lyear), 'Rows', 'complete');
% % %         if (tmp.skill_IPO_assm_p>0.1) tmp.skill_IPO_assm=NaN(2,2); end
% % %         [tmp.skill_IPO_hcst, tmp.skill_IPO_hcst_p]=corrcoef(IPO(isfinite(IPO_year)), tmp.HCST_mean(1:end-lyear), 'Rows', 'complete');
% % %         if (tmp.skill_IPO_hcst_p>0.1) tmp.skill_IPO_hcst=NaN(2,2); end
% % % 
% % %         plot(tmp.time, IPO(isfinite(IPO_year)))
% % %         hold on
% % %         plot(tmp.time, tmp.ASSM_mean(1:end-lyear));
% % %         hold off
% % %         
% % %         for loni=1:size(grids.tlong,1)
% % %             for lati=1:size(grids.tlat,2)
% % %                 tmp.corr=corrcoef(IPO(isfinite(IPO_year)), data2.([tmp.varname, '_assm_ano_l', tmp.lyear_str])(loni,lati,:), 'Rows', 'complete');
% % %                 tmp.corrmap(loni,lati)=tmp.corr(1,2);
% % %             end
% % %         end
% % %         pcolor(tmp.corrmap'); shading flat; colorbar;

        %% corr skills
        [tmp.pot_skill_hcst, tmp.pot_skill_hcst_p]=corrcoef(tmp.ASSM_mean(1:end-lyear), tmp.HCST_mean(1:end-lyear));
        if (tmp.pot_skill_hcst_p>0.1) tmp.pot_skill_hcst=NaN(2,2); end
        [tmp.pot_skill_lens2, tmp.pot_skill_lens2_p]=corrcoef(tmp.ASSM_mean(1:end-lyear), tmp.LENS2_mean(1:end-lyear));
        if (tmp.pot_skill_lens2_p>0.1) tmp.pot_skill_lens2=NaN(2,2); end

        [tmp.skill_assm, tmp.skill_assm_p]=corrcoef(tmp.OBS(1:end-lyear), tmp.ASSM_mean(1:end-lyear), 'Rows', 'complete');
        if (tmp.skill_assm_p>0.1) tmp.skill_assm=NaN(2,2); end
        [tmp.skill_hcst, tmp.skill_hcst_p]=corrcoef(tmp.OBS(1:end-lyear), tmp.HCST_mean(1:end-lyear), 'Rows', 'complete');
        if (tmp.skill_hcst_p>0.1) tmp.skill_hcst=NaN(2,2); end
        [tmp.skill_lens2, tmp.skill_lens2_p]=corrcoef(tmp.OBS(1:end-lyear), tmp.LENS2_mean(1:end-lyear), 'Rows', 'complete');
        if (tmp.skill_lens2_p>0.1) tmp.skill_lens2=NaN(2,2); end
        
        yl=ylim;
        text(1960, min(yl)+diff(yl)/30, ['AS-HC:', num2str(round(tmp.pot_skill_hcst(1,2),2))])
        text(1967, min(yl)+diff(yl)/30, ['AS-LE:', num2str(round(tmp.pot_skill_lens2(1,2),2))])
        
        text(1980, min(yl)+diff(yl)/30, ['OB-AS:', num2str(round(tmp.skill_assm(1,2),2))])
        text(1987, min(yl)+diff(yl)/30, ['OB-HC:', num2str(round(tmp.skill_hcst(1,2),2))])
        text(1994, min(yl)+diff(yl)/30, ['OB-LE:', num2str(round(tmp.skill_lens2(1,2),2))])

        if length(sta_lonlat{stai})==4
            dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_time_series', filesep, 'regional_mean',  filesep, 'l',tmp.lyear_str];
        else
            dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_time_series', filesep, 'l',tmp.lyear_str];
        end

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
        
        

        %% detrended plot
        fig_h = figure('name','ts','visible','off');
        fig_h.Position= [0,0,1000,500];

        [tmp.HCST_mean_det, tmp.HCST_mean_trend] = Func_0028_detrend_linear_1d(tmp.HCST_mean, 'omitnan');
        tmp.HCST_mean_tr_ano=tmp.HCST_mean-tmp.HCST_mean_det;


        [tmp.ASSM_mean_det, tmp.ASSM_mean_trend] = Func_0028_detrend_linear_1d(tmp.ASSM_mean, 'omitnan');
        tmp.ASSM_mean_tr_ano=tmp.ASSM_mean-tmp.ASSM_mean_det;
        
        [tmp.LENS2_mean_det, tmp.LENS2_mean_trend] = Func_0028_detrend_linear_1d(tmp.LENS2_mean, 'omitnan');
        tmp.LENS2_mean_tr_ano=tmp.LENS2_mean-tmp.LENS2_mean_det;
                
        if length(sta_lonlat{stai})==2
            tmp.HCST_lower_det = tmp.HCST_lower-tmp.HCST_mean_tr_ano;
            tmp.HCST_upper_det = tmp.HCST_upper-tmp.HCST_mean_tr_ano;
            tmp.ASSM_lower_det = tmp.ASSM_lower-tmp.ASSM_mean_tr_ano;
            tmp.ASSM_upper_det = tmp.ASSM_upper-tmp.ASSM_mean_tr_ano;
            tmp.LENS2_lower_det = tmp.LENS2_lower-tmp.LENS2_mean_tr_ano;
            tmp.LENS2_upper_det = tmp.LENS2_upper-tmp.LENS2_mean_tr_ano;
        end

        [tmp.OBS_det, tmp.OBS_trend] =Func_0028_detrend_linear_1d(tmp.OBS, 'omitnan');

        hold on
        if length(sta_lonlat{stai})==2
            fig_ts.LENS2_range=fill([tmp.time, flip(tmp.time)], ...
                [tmp.LENS2_lower_det(1:end-lyear); flip(tmp.LENS2_upper_det(1:end-lyear))], cmap_LENS2_b);
            fig_ts.LENS2_range.EdgeColor = 'none';
            set(get(get(fig_ts.LENS2_range,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
            fig_ts.HCST_range=fill([tmp.time, flip(tmp.time)], ...
                [tmp.HCST_lower_det; flip(tmp.HCST_upper_det)], cmap_HCST_b);
            fig_ts.HCST_range.EdgeColor = 'none';
            set(get(get(fig_ts.HCST_range,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
            fig_ts.ASSM_range=fill([tmp.time, flip(tmp.time)], ...
                [tmp.ASSM_lower_det(1:end-lyear); flip(tmp.ASSM_upper_det(1:end-lyear))], cmap_ASSM_b);
            fig_ts.ASSM_range.EdgeColor = 'none';
            set(get(get(fig_ts.ASSM_range,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end


        plot(tmp.time, tmp.LENS2_mean_det(1:end-lyear), 'linewidth', 2, 'color', cmap_LENS2)
        plot(tmp.time, tmp.ASSM_mean_det(1:end-lyear), 'linewidth', 2, 'color', cmap_ASSM)
        plot(tmp.time, tmp.HCST_mean_det, 'linewidth', 2, 'color', cmap_HCST)
        plot(tmp.time, tmp.OBS_det(1:end-lyear), 'linewidth', 2, 'color', cmap_OBS)

        hold off

%         legend ('LENS2', 'ASSM', 'HCST', 'OBS', 'Location', 'Northwest')
        legend ('LENS2', 'ASSM', 'HCST', 'OBS', 'Location', 'Southoutside', 'Orientation', 'Horizontal')

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
        if (tmp.pot_skill_hcst_p>0.1) tmp.pot_skill_hcst=NaN(2,2); end
        [tmp.pot_skill_lens2, tmp.pot_skill_lens2_p]=corrcoef(tmp.ASSM_mean_det(1:end-lyear), tmp.LENS2_mean_det(1:end-lyear));
        if (tmp.pot_skill_lens2_p>0.1) tmp.pot_skill_lens2=NaN(2,2); end

        [tmp.skill_assm, tmp.skill_assm_p]=corrcoef(tmp.OBS_det(1:end-lyear), tmp.ASSM_mean_det(1:end-lyear), 'Rows', 'complete');
        if (tmp.skill_assm_p>0.1) tmp.skill_assm=NaN(2,2); end
        [tmp.skill_hcst, tmp.skill_hcst_p]=corrcoef(tmp.OBS_det(1:end-lyear), tmp.HCST_mean_det(1:end-lyear), 'Rows', 'complete');
        if (tmp.skill_hcst_p>0.1) tmp.skill_hcst=NaN(2,2); end
        [tmp.skill_lens2, tmp.skill_lens2_p]=corrcoef(tmp.OBS_det(1:end-lyear), tmp.LENS2_mean_det(1:end-lyear), 'Rows', 'complete');
        if (tmp.skill_lens2_p>0.1) tmp.skill_lens2=NaN(2,2); end
        
        yl=ylim;
        text(1960, min(yl)+diff(yl)/30, ['AS-HC:', num2str(round(tmp.pot_skill_hcst(1,2),2))])
        text(1967, min(yl)+diff(yl)/30, ['AS-LE:', num2str(round(tmp.pot_skill_lens2(1,2),2))])
        
        text(1980, min(yl)+diff(yl)/30, ['OB-AS:', num2str(round(tmp.skill_assm(1,2),2))])
        text(1987, min(yl)+diff(yl)/30, ['OB-HC:', num2str(round(tmp.skill_hcst(1,2),2))])
        text(1994, min(yl)+diff(yl)/30, ['OB-LE:', num2str(round(tmp.skill_lens2(1,2),2))])
        

        if length(sta_lonlat{stai})==4
            dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_time_series_det', filesep, 'regional_mean',  filesep, 'l',tmp.lyear_str];
        else
            dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_time_series_det', filesep, 'l',tmp.lyear_str];            
        end
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




        %% normalized plot
        fig_h = figure('name','ts','visible','off');
        fig_h.Position= [0,0,1000,500];

        tmp.HCST_mean_tr_ano=tmp.HCST_mean_norm;
        tmp.ASSM_mean_tr_ano=tmp.ASSM_mean_norm;
        tmp.LENS2_mean_tr_ano=tmp.LENS2_mean_norm;
                
        if length(sta_lonlat{stai})==2
            tmp.HCST_lower_norm = tmp.HCST_lower-tmp.HCST_mean_tr_ano;
            tmp.HCST_upper_norm = tmp.HCST_upper-tmp.HCST_mean_tr_ano;
            tmp.ASSM_lower_norm = tmp.ASSM_lower-tmp.ASSM_mean_tr_ano;
            tmp.ASSM_upper_norm = tmp.ASSM_upper-tmp.ASSM_mean_tr_ano;
            tmp.LENS2_lower_norm = tmp.LENS2_lower-tmp.LENS2_mean_tr_ano;
            tmp.LENS2_upper_norm = tmp.LENS2_upper-tmp.LENS2_mean_tr_ano;
        end

        

        hold on
        if length(sta_lonlat{stai})==2
            fig_ts.LENS2_range=fill([tmp.time, flip(tmp.time)], ...
                [tmp.LENS2_lower_norm(1:end-lyear); flip(tmp.LENS2_upper_norm(1:end-lyear))], cmap_LENS2_b);
            fig_ts.LENS2_range.EdgeColor = 'none';
            set(get(get(fig_ts.LENS2_range,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
            fig_ts.HCST_range=fill([tmp.time, flip(tmp.time)], ...
                [tmp.HCST_lower_norm; flip(tmp.HCST_upper_norm)], cmap_HCST_b);
            fig_ts.HCST_range.EdgeColor = 'none';
            set(get(get(fig_ts.HCST_range,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
            fig_ts.ASSM_range=fill([tmp.time, flip(tmp.time)], ...
                [tmp.ASSM_lower_norm(1:end-lyear); flip(tmp.ASSM_upper_norm(1:end-lyear))], cmap_ASSM_b);
            fig_ts.ASSM_range.EdgeColor = 'none';
            set(get(get(fig_ts.ASSM_range,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end


        plot(tmp.time, tmp.LENS2_mean_norm(1:end-lyear), 'linewidth', 2, 'color', cmap_LENS2)
        plot(tmp.time, tmp.ASSM_mean_norm(1:end-lyear), 'linewidth', 2, 'color', cmap_ASSM)
        plot(tmp.time, tmp.HCST_mean_norm, 'linewidth', 2, 'color', cmap_HCST)
        plot(tmp.time, tmp.OBS_norm(1:end-lyear), 'linewidth', 2, 'color', cmap_OBS)

        hold off

%         legend ('LENS2', 'ASSM', 'HCST', 'OBS', 'Location', 'Northwest')
        legend ('LENS2', 'ASSM', 'HCST', 'OBS', 'Location', 'Southoutside', 'Orientation', 'Horizontal')

        if length(sta_lonlat{stai})==2
            title(['norm, l', tmp.lyear_str, ', ', tmp.varname, ', ', num2str(xpoint),'E, ', num2str(ypoint), 'N'])
        elseif length(sta_lonlat{stai})==4
            title(['norm, l', tmp.lyear_str, ', ', tmp.varname, ', ', num2str(xpoint1),'E~', num2str(xpoint2),'E, ', ...
                num2str(ypoint1), 'N~',num2str(ypoint2), 'N'])
        end
        grid minor
        xlim([1960 2025])
        set(gca, 'fontsize', 20)
        
        %% corr skills
        [tmp.pot_skill_hcst, tmp.pot_skill_hcst_p]=corrcoef(tmp.ASSM_mean_norm(1:end-lyear), tmp.HCST_mean_norm(1:end-lyear));
        if (tmp.pot_skill_hcst_p>0.1) tmp.pot_skill_hcst=NaN(2,2); end
        [tmp.pot_skill_lens2, tmp.pot_skill_lens2_p]=corrcoef(tmp.ASSM_mean_norm(1:end-lyear), tmp.LENS2_mean_norm(1:end-lyear));
        if (tmp.pot_skill_lens2_p>0.1) tmp.pot_skill_lens2=NaN(2,2); end

        [tmp.skill_assm, tmp.skill_assm_p]=corrcoef(tmp.OBS_norm(1:end-lyear), tmp.ASSM_mean_norm(1:end-lyear), 'Rows', 'complete');
        if (tmp.skill_assm_p>0.1) tmp.skill_assm=NaN(2,2); end
        [tmp.skill_hcst, tmp.skill_hcst_p]=corrcoef(tmp.OBS_norm(1:end-lyear), tmp.HCST_mean_norm(1:end-lyear), 'Rows', 'complete');
        if (tmp.skill_hcst_p>0.1) tmp.skill_hcst=NaN(2,2); end
        [tmp.skill_lens2, tmp.skill_lens2_p]=corrcoef(tmp.OBS_norm(1:end-lyear), tmp.LENS2_mean_norm(1:end-lyear), 'Rows', 'complete');
        if (tmp.skill_lens2_p>0.1) tmp.skill_lens2=NaN(2,2); end
        
        yl=ylim;
        text(1960, min(yl)+diff(yl)/30, ['AS-HC:', num2str(round(tmp.pot_skill_hcst(1,2),2))])
        text(1967, min(yl)+diff(yl)/30, ['AS-LE:', num2str(round(tmp.pot_skill_lens2(1,2),2))])
        
        text(1980, min(yl)+diff(yl)/30, ['OB-AS:', num2str(round(tmp.skill_assm(1,2),2))])
        text(1987, min(yl)+diff(yl)/30, ['OB-HC:', num2str(round(tmp.skill_hcst(1,2),2))])
        text(1994, min(yl)+diff(yl)/30, ['OB-LE:', num2str(round(tmp.skill_lens2(1,2),2))])
        

        if length(sta_lonlat{stai})==4
            dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_time_series_norm', filesep, 'regional_mean',  filesep, 'l',tmp.lyear_str];
        else
            dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_time_series_norm', filesep, 'l',tmp.lyear_str];            
        end
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        if length(sta_lonlat{stai})==2
            cfg.figname=[dirs.figdir, filesep, 'ts_norm_all_l',tmp.lyear_str, '_', num2str(xpoint), 'E_', num2str(ypoint), 'N_', tmp.varname, '.tif'];
        elseif length(sta_lonlat{stai})==4
            cfg.figname=[dirs.figdir, filesep, 'ts_norm_all_l',tmp.lyear_str, '_', num2str(xpoint1), 'E_', num2str(xpoint2), 'E_', ...
                num2str(ypoint1), 'N_', num2str(ypoint2), 'N_', tmp.varname, '.tif'];
        end
            print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;




    end


% end





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


