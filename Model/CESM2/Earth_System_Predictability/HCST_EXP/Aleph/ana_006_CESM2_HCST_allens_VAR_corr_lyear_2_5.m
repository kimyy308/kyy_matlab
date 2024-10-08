% %  Created 12-Apr-2023 by Yong-Yub Kim
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

cfg.ly_s = 2;
cfg.ly_e = 5;

%% model configuration
% cfg.var='TS'; %SST PRECT PSL TS SSH sumChl
% cfg.vars = {'SSH', 'sumChl'};
% % cfg.vars = {'SST', 'PRECT', 'PSL', 'TS'};
% cfg.vars = {'SSH'};
% cfg.vars = {'SST', 'PRECT', 'PSL', 'TS', 'sumChl', 'photoC_TOT_zint'};
% cfg.vars = { 'sumChl', 'photoC_TOT_zint',  'SSH'};
% cfg.vars = {  'SST', 'PRECT', 'PSL'};
% % 
% cfg.vars = { 'PAR_avg', 'SiO3', 'Fe', 'PO4', 'NO3'};
% cfg.vars = {'Fe', 'PO4', 'NO3'};

cfg.vars = {'diatChl', 'diazChl', 'spChl'};

% cfg.vars = {'diatChl'};
% cfg.vars = {'spChl'};
cfg.vars = {'SST', 'PRECT', 'PSL', 'TS'};
% cfg.vars = {'WVEL2'};
% cfg.vars = {'diatC', 'diazC', 'spC'};
% cfg.vars = {'SSH', 'BSF', 'diatChl', 'diazChl', 'spChl', 'diatC',
% 'diazC', 'spC', 'SST', 'PRECT', 'PSL', 'TS'};
cfg.vars = {'graze_diaz_zint_100m'};
cfg.vars = { 'TEMP', 'UVEL', 'VVEL', 'WVEL' };
cfg.vars = { 'diazC' };
% cfg.vars = { 'zooC' };
cfg.vars = {'sp_P_lim_Cweight_avg_100m'};
% cfg.vars = {'diat_light_lim_Cweight_avg_100m'};
cfg.vars = { 'HBLT' };
% cfg.vars = { 'diatC', 'Fe', 'spC', 'NO3', 'PO4', 'SiO3', 'TEMP', 'UVEL', 'VVEL', 'WVEL', 'PAR_avg', 'BSF', 'diazC' };
cfg.vars={ 'SST', 'PRECT', 'TS',  ...
    'GPP', 'NPP', 'TOTVEGC', 'TWS', 'PSL', 'SST'}; 
    
    cfg.vars = {'aice', 'sithick', ...
    'photoC_TOT_zint_100m', 'sumChl', ...
    'diatC', 'Fe', 'spC', 'NO3', 'PO4', 'SiO3', 'TEMP', 'UVEL', 'VVEL', 'WVEL', ...
    'PAR_avg', 'BSF', 'diazC', 'zooC', 'PD', ...
    'SSH', 'HBLT', 'diatChl', 'spChl', 'diazChl', ...
    'diat_agg_zint_100m', 'diat_Fe_lim_Cweight_avg_100m', ...
    'diat_Fe_lim_surf', 'diat_light_lim_Cweight_avg_100m', 'diat_light_lim_surf', ...
    'diat_loss_zint_100m', 'diat_N_lim_Cweight_avg_100m', 'diat_N_lim_surf', ...
    'diat_P_lim_Cweight_avg_100m', 'diat_P_lim_surf', 'diat_SiO3_lim_Cweight_avg_100m', ...
    'diat_SiO3_lim_surf', 'diaz_agg_zint_100m', 'diaz_Fe_lim_Cweight_avg_100m', ...
    'diaz_Fe_lim_surf', 'diaz_light_lim_Cweight_avg_100m', 'diaz_light_lim_surf', ...
    'diaz_loss_zint_100m', 'diaz_P_lim_Cweight_avg_100m', 'diaz_P_lim_surf',  ...
    'graze_diat_zint_100m', 'graze_diat_zoo_zint_100m', 'graze_diaz_zint_100m', ...
    'graze_diaz_zoo_zint_100m', 'graze_sp_zint_100m', 'graze_sp_zoo_zint_100m', ...
    'LWUP_F', 'photoC_diat_zint_100m', 'photoC_diaz_zint_100m', ...
    'photoC_sp_zint_100m', ...
    'sp_agg_zint_100m', 'sp_Fe_lim_Cweight_avg_100m', 'sp_Fe_lim_surf', ...
    'sp_light_lim_Cweight_avg_100m', 'sp_light_lim_surf', 'sp_loss_zint_100m', ...
    'sp_N_lim_Cweight_avg_100m', 'sp_N_lim_surf', 'sp_P_lim_Cweight_avg_100m', ...
    'sp_P_lim_surf', 'zoo_loss_zint_100m', 'IRON_FLUX',  'TEMPCLINE', 'NO3CLINE', ...
    'UVEL145', 'VVEL145', 'NO3145', 'PO4145', 'Fe145', 'TEMP145'};
% cfg.vars = { 'RAIN', 'QSOIL', 'QSOIL_ICE', 'QRUNOFF', 'QOVER', 'QRGWL', 'QH2OSFC', 'NEP', 'DSTFLXT'} %TOTVEGC, FIRE
cfg.vars = { 'GPP', 'NPP', 'TOTVEGC', 'TWS', 'COL_FIRE_CLOSS', 'COL_FIRE_NLOSS', 'FIRE', ...
            'FPSN', 'SOILICE', 'SOILLIQ', 'TOTSOILICE', 'TOTSOILLIQ', ...
            'RAIN', 'QSOIL', 'QSOIL_ICE', 'QRUNOFF', 'QOVER', 'QRGWL', ...
            'QH2OSFC', 'NEP', 'DSTFLXT'}; %TOTVEGC, FIRE
% cfg.vars = { 'TEMP', 'PD'} %TOTVEGC, FIRE
% cfg.vars = { 'TEMP145'} %TOTVEGC, FIRE
% cfg.vars = { 'PD145'} %TOTVEGC, FIRE
cfg.vars = {'SST', 'PRECT', 'TS', 'PSL', 'AEROD_v', 'FSDS', 'FSNS', ...
    'SFdst_a1', 'SFdst_a2', 'SFdst_a3', 'U10', 'SFCO2', 'CLDTOT'};
% cfg.vars = { 'zooC'}; %TOTVEGC, FIRE

% % % cfg.vars={'photoC_TOT_zint_100m',  'photoC_diat_zint_100m', 'photoC_diaz_zint_100m', 'photoC_sp_zint_100m', ...
% % %      'diat_Fe_lim_Cweight_avg_100m', ...
% % %     'diat_light_lim_Cweight_avg_100m',  ...
% % %      'diat_N_lim_Cweight_avg_100m',  ...
% % %     'diat_P_lim_Cweight_avg_100m',  'diat_SiO3_lim_Cweight_avg_100m', ...
% % %     'sp_Fe_lim_Cweight_avg_100m',  ...
% % %     'sp_light_lim_Cweight_avg_100m',   ...
% % %     'sp_N_lim_Cweight_avg_100m', 'sp_P_lim_Cweight_avg_100m', ...
% % %      'diaz_Fe_lim_Cweight_avg_100m', ...
% % %      'diaz_light_lim_Cweight_avg_100m',  'diaz_P_lim_Cweight_avg_100m', ...
% % %     'zoo_loss_zint_100m', ...
% % %     'diaz_agg_zint_100m', 'diat_loss_zint_100m', 'diat_agg_zint_100m', 'diaz_loss_zint_100m',  ...
% % %     'sp_agg_zint_100m', 'sp_loss_zint_100m', ...
% % %     'graze_diat_zint_100m', 'graze_diat_zoo_zint_100m', 'graze_diaz_zint_100m', ...
% % %     'graze_diaz_zoo_zint_100m', 'graze_sp_zint_100m', 'graze_sp_zoo_zint_100m', ...
% % %     'diat_Fe_lim_surf', 'diat_light_lim_surf', 'diat_P_lim_surf', 'diat_SiO3_lim_surf', 'diaz_Fe_lim_surf', ...
% % %     'diaz_light_lim_surf', 'diaz_P_lim_surf', 'dustToSed', 'LWUP_F', 'O2_ZMIN_DEPTH', 'diat_N_lim_surf', ...
% % %     'sp_Fe_lim_surf', 'sp_light_lim_surf', 'sp_N_lim_surf', 'sp_P_lim_surf'  };
% % % 
% % % % cfg.vars = {'sumChl'};
% % % % cfg.vars = {'TEMPCLINE', 'NO3CLINE', 'UVEL55', 'VVEL55'};
% % % cfg.vars = {'UVEL55', 'VVEL55'};
% % % cfg.vars = {'IRON_FLUX'};
% % % % cfg.vars = {'UVEL145', 'VVEL145', 'NO3145', 'PO4145', 'Fe145', 'TEMP145'};
% % % % cfg.vars = {'NO3145', 'PO4145', 'Fe145', 'TEMP145'};
% % % cfg.vars = {'TS'};
% % % cfg.vars = {'GPP', 'NPP', 'TOTVEGC', 'TWS', 'aice', 'sithick'};
% % % cfg.vars = {'PRECT'};


% tmp.dimids= [1, 2, 4];
% cfg.vlayer=1; % surf, vertical slice 

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
    % cfg.obs_iyears=f_obs_iyears(cfg.var);
    
    % dirs.hcstroot=['/mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive/', cfg.comp, '/', cfg.var];
    % dirs.obsroot=['/Volumes/kyy_raid/kimyy/Observation/', cfg.obs_name, '/monthly_reg_', cfg.obs_fname_module(2:4)];
    % dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', cfg.var];
    % dirs.lens2root=['/mnt/lustre/proj/earth.system.predictability/LENS2/archive_analysis/', cfg.comp, '/', cfg.var];
    % dirs.assmroot=['/mnt/lustre/proj/earth.system.predictability/ASSM_EXP/archive_analysis/', cfg.comp, '/', cfg.var];
    
    dirs.hcstroot=['/mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive_transfer/', cfg.comp, '/', cfg.var];
    dirs.obsroot=['/mnt/lustre/proj/kimyy/Observation/', cfg.obs_name, '/monthly_reg_', cfg.obs_fname_module(2:4)];
    dirs.matroot=['/mnt/lustre/proj/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp,'/', cfg.var];
    dirs.lens2root=['/mnt/lustre/proj/earth.system.predictability/LENS2/archive_transfer/', cfg.comp, '/', cfg.var];
    dirs.assmroot=['/mnt/lustre/proj/earth.system.predictability/ASSM_EXP/archive_transfer/', cfg.comp, '/', cfg.var];
    
    
    
    cfg.iyears=cfg.obs_iyears;
    cfg.gnm='f09_g17';
%     cfg.proj_year=1;
    cfg.proj_year=5;

    cfg.max_iy=max(cfg.iyears);
    cfg.min_iy=min(cfg.iyears);
    cfg.len_t_y = length(cfg.iyears);
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
    tmp.gridname = [dirs.hcstroot, tmp.fs, '../grid.nc'];
    tmp.maskname = '/mnt/lustre/proj/kimyy/Model/CESM2/ESP/grids/ocn/RECCAP2_region_masks_all_v20210412_POP2_grid.nc';
    
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
            grid.z_t=ncread(tmp.gridname, 'z_t')./100; % meter
            grid.dz=ncread(tmp.gridname, 'dz')./100; % meter
            grid.dz_res=reshape(grid.dz(cfg.vlayer), [1 1 length(grid.dz(cfg.vlayer))]);
        case {'atm', 'lnd'}
            grid.lon=ncread(tmp.gridname, 'lon');
            grid.lat=ncread(tmp.gridname, 'lat');
            [grid.tlat grid.tlong]=meshgrid(grid.lat, grid.lon);
    end
    
    grid.nlon=size(grid.tlong,1);
    grid.nlat=size(grid.tlat,2);
    % grid.ntime=cfg.proj_year.*12;
    
    % % model filename example
    % /mnt/lustre/proj/kimyy/Model/CESM2/ESP/HCST_EXP/archive/atm/SST/ens_all/ens_all_i2020
    % SST_f09_g17.hcst.ens_all_i2020.cam.h0.2024-12.nc
        
    %% read & plot data
    tmp.varname=cfg.var;
    if length(tmp.varname)>=3
     tmp.varn3=tmp.varname(end-2:end);
     switch tmp.varn3
         case '145'
             tmp.fvarname=tmp.varname(1:end-3);
         otherwise
            tmp.fvarname=cfg.var;
     end
    else
        tmp.fvarname=cfg.var;
    end

    clear tmp.ydata tmp.ydata_lens2 tmp.ydata_obs tmp.ydata_assm

    for lyear_f=cfg.ly_s:cfg.ly_e
        lyear=lyear_f-1;
        tmp.lyear_str=num2str(lyear, '%02i');
        cfg.casename_m=['all'];
    
        dirs.datadir= [dirs.hcstroot, filesep, 'ens_all', filesep];
        dirs.assmdir= [dirs.assmroot, filesep, 'ens_all', filesep];
        dirs.lens2dir= [dirs.lens2root, filesep, 'ens_smbb', filesep];
        fprintf('%02d_%s_%s  ',lyear, ',', cfg.casename_m,'_', tmp.varname); lap_time = tic;
    
        
        %% variables initialization
        data.([tmp.varname, '_model', '_l', tmp.lyear_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data.([tmp.varname, '_lens2', '_l', tmp.lyear_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data.([tmp.varname, '_bias', '_l', tmp.lyear_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data.([tmp.varname, '_obs', '_l', tmp.lyear_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data.([tmp.varname, '_assm', '_l', tmp.lyear_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data.([tmp.varname, '_AR1', '_l', tmp.lyear_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);  %% AR1 initialization
    
    %     data.([tmp.varname, '_assm'])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        
        %% read variables
        for iyear=min(cfg.iyears):max(cfg.iyears)
            tmp.iyear_str=num2str(iyear, '%04i');
            cfg.casename=['ensmean_', cfg.casename_m, '_i', tmp.iyear_str];
            tmp.fy=iyear+lyear;
            tmp.fy_str=num2str(tmp.fy, '%04i');
            %% monthly filename
            for mon=1:12
                tmp.mon_str=num2str(mon, '%02i');
    
                %% HCST
                cfg.mod_fnm=[dirs.datadir, tmp.fs, 'ens_all_i',tmp.iyear_str, tmp.fs, ...
                tmp.varname, '_', cfg.gnm, '.hcst.', cfg.casename, cfg.obs_fname_module, tmp.fy_str, '-', tmp.mon_str, '.nc'];
                try
                    ncid=netcdf.open(cfg.mod_fnm, 'NOWRITE');
                catch
                    disp(['csh ', '/mnt/lustre/proj/kimyy/Scripts/Model/CESM2/HCST/', 'hcst_extract_fix_matlab_var.csh ', ...
                        tmp.varname, ' ', num2str(length(tmp.dimids)-1), ' ', tmp.iyear_str, ' ', tmp.fy_str, ' ', tmp.mon_str, ' ', cfg.comp]);
                     disp(['csh ', '/mnt/lustre/proj/kimyy/Scripts/Model/CESM2/HCST/', 'hcst_ens_fix_matlab.csh ', ...
                        tmp.varname, ' ', tmp.iyear_str, ' ', tmp.fy_str, ' ', tmp.mon_str, ' ', cfg.comp]);
                      system(['csh ', '/mnt/lustre/proj/kimyy/Scripts/Model/CESM2/HCST/', 'hcst_extract_fix_matlab_var.csh ', ...
                        tmp.varname, ' ', num2str(length(tmp.dimids)-1), ' ', tmp.iyear_str, ' ', tmp.fy_str, ' ', tmp.mon_str, ' ', cfg.comp]);
                     system(['csh ', '/mnt/lustre/proj/kimyy/Scripts/Model/CESM2/HCST/', 'hcst_ens_fix_matlab.csh ', ...
                        tmp.varname, ' ', tmp.iyear_str, ' ', tmp.fy_str, ' ', tmp.mon_str, ' ', cfg.comp]);
                     ncid=netcdf.open(cfg.mod_fnm, 'NOWRITE');

                end
                tmpvarid=netcdf.inqVarID(ncid,tmp.fvarname);
                [tmp.fvarname,tmp.xtype,tmp.dimids,tmp.natts]= netcdf.inqVar(ncid,tmpvarid);
                if length(tmp.dimids)>3
                     tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 0], [grid.nlon grid.nlat cfg.vlayer_cnt 1]));
%                      tmp.dd=tmp.dd.*grid.mask_ocn;
                     tmp.dd(abs(tmp.dd)>1e30)=NaN;
                     ddd=mean(tmp.dd,3,'omitnan'); % depth mean
                     tmp.ydata(1:grid.nlon,1:grid.nlat,mon) = ddd; %depth averaged value
                else
                    tmp.dd=  netcdf.getVar(ncid,tmpvarid);
                    tmp.dd(abs(tmp.dd)>1e30)=NaN;
                    tmp.ydata(1:grid.nlon,1:grid.nlat,mon) = tmp.dd;
                end
                netcdf.close(ncid);
                
                tmp.ydd=tmp.ydata(1:grid.nlon,1:grid.nlat,mon);
                tmp.stdydd=max(abs(tmp.ydd(isfinite(tmp.ydd))));

                if tmp.stdydd > 10e5 % if extracted files were crashed
                      system(['csh ', '/mnt/lustre/proj/kimyy/Scripts/Model/CESM2/HCST/', 'hcst_extract_fix_matlab_var.csh ', ...
                        tmp.varname, ' ', num2str(length(tmp.dimids)-1), ' ', tmp.iyear_str, ' ', tmp.fy_str, ' ', tmp.mon_str, ' ', cfg.comp]);
                     system(['csh ', '/mnt/lustre/proj/kimyy/Scripts/Model/CESM2/HCST/', 'hcst_ens_fix_matlab.csh ', ...
                        tmp.varname, ' ', tmp.iyear_str, ' ', tmp.fy_str, ' ', tmp.mon_str, ' ', cfg.comp]);
                     ncid=netcdf.open(cfg.mod_fnm, 'NOWRITE');
                    tmpvarid=netcdf.inqVarID(ncid,tmp.fvarname);
                    [tmp.fvarname,tmp.xtype,tmp.dimids,tmp.natts]= netcdf.inqVar(ncid,tmpvarid);
                    if length(tmp.dimids)>3
                         tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 0], [grid.nlon grid.nlat cfg.vlayer_cnt 1]));
    %                      tmp.dd=tmp.dd.*grid.mask_ocn;
                         tmp.dd(abs(tmp.dd)>1e30)=NaN;
                         ddd=mean(tmp.dd,3,'omitnan'); % depth mean
                         tmp.ydata(1:grid.nlon,1:grid.nlat,mon) = ddd; %depth averaged value
                    else
                        tmp.dd=  netcdf.getVar(ncid,tmpvarid);
                        tmp.dd(abs(tmp.dd)>1e30)=NaN;
                        tmp.ydata(1:grid.nlon,1:grid.nlat,mon) = tmp.dd;
                    end
                    netcdf.close(ncid);
                end

    %             tmp.ydata(:,:,mon)=ncread(cfg.mod_fnm, tmp.varname);
                
                %% LENS2
                if tmp.fy <= cfg.max_iy
                    cfg.lens2_fnm=[dirs.lens2dir, tmp.fs, ...
                        tmp.varname, '_ensmean_', tmp.fy_str, '-', tmp.mon_str, '.nc'];
                    try
                        ncid=netcdf.open(cfg.lens2_fnm, 'NOWRITE');
                    catch
                         disp(['csh ', '/mnt/lustre/proj/kimyy/Scripts/Model/CESM2/LENS2/', 'lens2_extract_ens_fix_var_parallel.csh ', ...
                            tmp.varname, ' ', tmp.fy_str, ' ', tmp.mon_str, ' ', cfg.comp]);                    
                         system(['csh ', '/mnt/lustre/proj/kimyy/Scripts/Model/CESM2/LENS2/', 'lens2_extract_ens_fix_var_parallel.csh ', ...
                            tmp.varname, ' ', tmp.fy_str,  ' ', tmp.mon_str, ' ', cfg.comp]);
                         ncid=netcdf.open(cfg.lens2_fnm, 'NOWRITE');
                    end
                    tmpvarid=netcdf.inqVarID(ncid,tmp.fvarname);
                    if length(tmp.dimids)>3
                         tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 0], [grid.nlon grid.nlat cfg.vlayer_cnt 1]));
%                          ddd=tmp.dd.*grid.dz_res; %% weight depth
%                          ddd=sum(ddd,3); % depth sum
%                          ddd=ddd./sum(grid.dz_res(:)); % divide total depth for avg
                         tmp.dd(abs(tmp.dd)>1e30)=NaN;
                         ddd=mean(tmp.dd,3,'omitnan'); % depth mean
                         tmp.ydata_lens2(1:grid.nlon,1:grid.nlat,mon) = ddd; %depth averaged value
                    else
                        tmp.dd=netcdf.getVar(ncid,tmpvarid);
                         tmp.dd(abs(tmp.dd)>1e30)=NaN;                        
                        tmp.ydata_lens2(1:grid.nlon,1:grid.nlat,mon) =tmp.dd;
                    end
                    netcdf.close(ncid);
                else
                    tmp.ydata_lens2(1:grid.nlon,1:grid.nlat,mon) = NaN;                    
                end
    
                %% OBS
                if tmp.fy <= cfg.max_iy
                    cfg.obs_fnm=[dirs.obsroot, tmp.fs, cfg.obs_fname_mid,tmp.fy_str,tmp.mon_str, '.nc'];
                    if exist(cfg.obs_fnm)~=0
                        ncid=netcdf.open(cfg.obs_fnm, 'NOWRITE');
                        tmpvarid=netcdf.inqVarID(ncid,cfg.obs_varname);
                        tmp.dd =  netcdf.getVar(ncid,tmpvarid);
                         tmp.dd(abs(tmp.dd)>1e30)=NaN;                    
                        tmp.ydata_obs(1:grid.nlon,1:grid.nlat,mon) = tmp.dd;
                        netcdf.close(ncid);
                    else
                        tmp.ydata_obs(1:grid.nlon,1:grid.nlat,mon) = NaN;                    
                    end
                else
                    tmp.ydata_obs(1:grid.nlon,1:grid.nlat,mon) = NaN;          
                end
    %             tmp.ydata_obs(:,:,mon)=ncread(cfg.obs_fnm, cfg.obs_varname);
    
                %% ASSM
                if tmp.fy <= cfg.max_iy
                    cfg.assm_fnm=[dirs.assmdir, tmp.fs, tmp.varname, '_ensmean_', tmp.fy_str,'-', tmp.mon_str, '.nc'];
                    ncid=netcdf.open(cfg.assm_fnm, 'NOWRITE');
                    tmpvarid=netcdf.inqVarID(ncid,tmp.fvarname);
                    if length(tmp.dimids)>3
                        tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 0], [grid.nlon grid.nlat cfg.vlayer_cnt 1]));
%                          ddd=tmp.dd.*grid.dz_res; %% weight depth
%                          ddd=sum(ddd,3); % depth sum
%                          ddd=ddd./sum(grid.dz_res(:)); % divide total depth for avg
                         tmp.dd(abs(tmp.dd)>1e30)=NaN;
                         ddd=mean(tmp.dd,3,'omitnan'); % depth mean
                         tmp.ydata_assm(1:grid.nlon,1:grid.nlat,mon) = ddd; %depth averaged value             
                    else
                        tmp.ydata_assm(1:grid.nlon,1:grid.nlat,mon) = netcdf.getVar(ncid,tmpvarid);
                    end
                    netcdf.close(ncid);
                else
                    tmp.ydata_assm(1:grid.nlon,1:grid.nlat,mon) = NaN;                    
                end
    
                
                if (strcmp(cfg.var, 'sumChl') || strcmp(cfg.var, 'TS')) %for Chls, it uses all available data
                    tmp.obs_mask(:,:)=tmp.ydata_obs(:,:,mon)./tmp.ydata_obs(:,:,mon);
                    tmp.ydata_mod_obs_masked(:,:,mon) = tmp.ydata(:,:,mon) .* tmp.obs_mask(:,:); % masking with available observation
                    tmp.ydata_assm_obs_masked(:,:,mon) = tmp.ydata_assm(:,:,mon) .* tmp.obs_mask(:,:); % masking with available observation
                    tmp.ydata_lens2_obs_masked(:,:,mon) = tmp.ydata_lens2(:,:,mon) .* tmp.obs_mask(:,:); % masking with available observation
                end
    
    %             %% AR1 (integration)
    %             tmp.ydata_AR1(:,:,mon) = ...
    %                 Func_0030_AR1_prog(tmp.all_data_obs(:,:,(tmp.iind-1)*12+1), data.([tmp.varname, '_AC_lag1']), lyear*12+mon-1, 0);
            end
            
    %         data.time=ncread(cfg.datafilename, 'time');
    %         tmp.ymean= mean(tmp.ydata,3);
    %         data.([tmp.varname, '_model', '_l', tmp.lyear_str])(:,:,iyear-min(cfg.iyears)+1+lyear)= tmp.ymean;
    %         tmp.ymean_obs= mean(tmp.ydata_obs,3);
    %         data.([tmp.varname, '_obs'])(:,:,iyear-min(cfg.iyears)+1)= tmp.ymean_obs;
    
    %         if ( strcmp(cfg.var, 'sumChl') || strcmp(cfg.var, 'TS')) %for Chls, it uses all available data
    %             tmp.ymean= mean(tmp.ydata,3, 'omitnan');
    %             tmp.ymean_obs= mean(tmp.ydata_obs,3, 'omitnan');
    %         else
                tmp.ymean= mean(tmp.ydata,3);
                if isfield(tmp, 'ydata_mod_obs_masked')
                    tmp.ymean_mod_obs_masked= mean(tmp.ydata_mod_obs_masked,3, 'omitnan');
                    tmp.ymean_assm_obs_masked= mean(tmp.ydata_assm_obs_masked,3, 'omitnan');
                    tmp.ymean_lens2_obs_masked= mean(tmp.ydata_lens2_obs_masked,3, 'omitnan');
                end
                tmp.ymean_obs= mean(tmp.ydata_obs,3, 'omitnan');
                tmp.ymean_assm= mean(tmp.ydata_assm,3, 'omitnan');
                tmp.ymean_assm(tmp.ymean_assm>10e30)=NaN;
                tmp.ymean_lens2= mean(tmp.ydata_lens2,3, 'omitnan');   
    
                switch cfg.comp
                    case 'ocn'
                        tmp.ymean = tmp.ymean .* grid.mask_ocn;
                        if isfield(tmp, 'ydata_mod_obs_masked')
                            tmp.ymean_mod_obs_masked = tmp.ymean_mod_obs_masked .* grid.mask_ocn;
                            tmp.ymean_assm_obs_masked = tmp.ymean_assm_obs_masked .* grid.mask_ocn;
                            tmp.ymean_lens2_obs_masked = tmp.ymean_lens2_obs_masked .* grid.mask_ocn;
                        end
                        tmp.ymean_obs = tmp.ymean_obs .* grid.mask_ocn;
                        tmp.ymean_assm = tmp.ymean_assm .* grid.mask_ocn;
                        tmp.ymean_lens2 = tmp.ymean_lens2 .* grid.mask_ocn;                        
                 end
    
    %         end
            data.([tmp.varname, '_model', '_l', tmp.lyear_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean;
            if isfield(tmp, 'ydata_mod_obs_masked')
                data.([tmp.varname, '_mod_obs_masked', '_l', tmp.lyear_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_mod_obs_masked;
                data.([tmp.varname, '_assm_obs_masked', '_l', tmp.lyear_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_assm_obs_masked;
                data.([tmp.varname, '_lens2_obs_masked', '_l', tmp.lyear_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_lens2_obs_masked;
            end
            data.([tmp.varname, '_obs', '_l', tmp.lyear_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_obs;
            data.([tmp.varname, '_lens2', '_l', tmp.lyear_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_lens2;
            data.([tmp.varname, '_assm', '_l', tmp.lyear_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_assm;
    
    %         %% AR1 (assign)
    %         tmp.ymean_AR1= mean(tmp.ydata_AR1, 3);
    %         data.([tmp.varname, '_AR1', '_l', tmp.lyear_str])(:,:,iyear-min(cfg.iyears)+1+lyear)= tmp.ymean_AR1;
    %         tmp.ymean= mean(ncread(cfg.datafilename, ['assm_', tmp.varname], [1, 1, 1], [inf, inf, 12]), 3);
    %         data.([tmp.varname, '_assm'])(:,:,iyear-min(cfg.iyears)+1)= tmp.ymean;
        end
    end
    
    data_all.([tmp.varname, '_model'])=zeros(size(data.([tmp.varname, '_assm', '_l', tmp.lyear_str])));
    data_all.([tmp.varname, '_obs'])=zeros(size(data.([tmp.varname, '_assm', '_l', tmp.lyear_str])));
    data_all.([tmp.varname, '_assm'])=zeros(size(data.([tmp.varname, '_assm', '_l', tmp.lyear_str])));
    data_all.([tmp.varname, '_lens2'])=zeros(size(data.([tmp.varname, '_assm', '_l', tmp.lyear_str])));
    if isfield(tmp, 'ydata_mod_obs_masked')
        data_all.([tmp.varname, '_mod_obs_masked'])=data_all.([tmp.varname, '_model']);
        data_all.([tmp.varname, '_assm_obs_masked'])=data_all.([tmp.varname, '_model']);
        data_all.([tmp.varname, '_lens2_obs_masked'])=data_all.([tmp.varname, '_model']); 
    end

     for lyear_f=cfg.ly_s:cfg.ly_e
        lyear=lyear_f-1;
        tmp.lyear_str=num2str(lyear, '%02i');
        data_all.([tmp.varname, '_model'])=data_all.([tmp.varname, '_model'])+data.([tmp.varname, '_model', '_l', tmp.lyear_str]);
        data_all.([tmp.varname, '_obs'])=data_all.([tmp.varname, '_obs'])+data.([tmp.varname, '_obs', '_l', tmp.lyear_str]);
        data_all.([tmp.varname, '_assm'])=data_all.([tmp.varname, '_assm'])+data.([tmp.varname, '_assm', '_l', tmp.lyear_str]);
        data_all.([tmp.varname, '_lens2'])=data_all.([tmp.varname, '_lens2'])+data.([tmp.varname, '_lens2', '_l', tmp.lyear_str]);
        if isfield(tmp, 'ydata_mod_obs_masked')
            data_all.([tmp.varname, '_mod_obs_masked'])=data_all.([tmp.varname, '_mod_obs_masked']) + ...
                data.([tmp.varname, '_mod_obs_masked', '_l', tmp.lyear_str]);
            data_all.([tmp.varname, '_assm_obs_masked'])=data_all.([tmp.varname, '_assm_obs_masked']) + ...
                data.([tmp.varname, '_assm_obs_masked', '_l', tmp.lyear_str]);
            data_all.([tmp.varname, '_lens2_obs_masked'])=data_all.([tmp.varname, '_lens2_obs_masked']) + ...
                data.([tmp.varname, '_lens2_obs_masked', '_l', tmp.lyear_str]);
        end
     end
    tmp.ly_size=length(cfg.ly_s:cfg.ly_e);
    data_all.([tmp.varname, '_model'])=data_all.([tmp.varname, '_model'])/tmp.ly_size;
    data_all.([tmp.varname, '_obs'])=data_all.([tmp.varname, '_obs'])/tmp.ly_size;
    data_all.([tmp.varname, '_assm'])=data_all.([tmp.varname, '_assm'])/tmp.ly_size;
    data_all.([tmp.varname, '_lens2'])=data_all.([tmp.varname, '_lens2'])/tmp.ly_size;

    if isfield(tmp, 'ydata_mod_obs_masked')
        data_all.([tmp.varname, '_mod_obs_masked']) = data_all.([tmp.varname, '_mod_obs_masked']) / tmp.ly_size;
        data_all.([tmp.varname, '_assm_obs_masked']) = data_all.([tmp.varname, '_assm_obs_masked']) / tmp.ly_size;
        data_all.([tmp.varname, '_lens2_obs_masked']) = data_all.([tmp.varname, '_lens2_obs_masked']) / tmp.ly_size;        
    end

    %% get correlation coefficient
    for loni=1:grid.nlon
        for lati=1:grid.nlat
             if (isnan(data_all.([tmp.varname, '_assm'])(loni,lati,1))~=1 ...
                     & nansum(data_all.([tmp.varname, '_model'])(loni,lati,:))~=0)

             %% corr assm
                 tmp.data_assm = data_all.([tmp.varname, '_assm'])(loni,lati,:);
                 tmp.data = squeeze(data_all.([tmp.varname, '_model'])(loni,lati,:));
                 [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_assm)), tmp.data_assm(isfinite(tmp.data_assm)));
                 tmp.data_det = Func_0028_detrend_linear_1d(squeeze(data_all.([tmp.varname, '_model'])(loni,lati,:)), 'omitnan');                         
                 tmp.data_assm_det = Func_0028_detrend_linear_1d(squeeze(data_all.([tmp.varname, '_assm'])(loni,lati,:)), 'omitnan');
                 [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data_det(isfinite(tmp.data_assm_det)), tmp.data_assm_det(isfinite(tmp.data_assm_det)));

                 data_all.([tmp.varname, '_corr_assm'])(loni,lati)=tmp.corr(1,2);
                 data_all.([tmp.varname, '_corr_assm_p'])(loni,lati)=tmp.corr_p(1,2);
                 data_all.([tmp.varname, '_corr_assm_det'])(loni,lati)=tmp.corr_det(1,2);
                 data_all.([tmp.varname, '_corr_assm_det_p'])(loni,lati)=tmp.corr_det_p(1,2);

             %% hcst-lens2 ~ assm
                 tmp.data_hcst_lens2 = squeeze(data_all.([tmp.varname, '_model'])(loni,lati,:)) - ...
                                squeeze(data_all.([tmp.varname, '_lens2'])(loni,lati,:));
                 [tmp.corr_int, tmp.corr_int_p]=corrcoef(tmp.data_hcst_lens2(isfinite(tmp.data_assm)), tmp.data_assm(isfinite(tmp.data_assm)));

                 
                 data_all.([tmp.varname, '_corr_assm_int'])(loni,lati)=tmp.corr_int(1,2);
                 data_all.([tmp.varname, '_corr_assm_int_p'])(loni,lati)=tmp.corr_int_p(1,2);

             %% corr lens2 ~ assm
                 tmp.lens2 = squeeze(data_all.([tmp.varname, '_lens2'])(loni,lati,:));
                 [tmp.corr, tmp.corr_p]=corrcoef(tmp.lens2(isfinite(tmp.data_assm)), tmp.data_assm(isfinite(tmp.data_assm)));
                 tmp.lens2_det = Func_0028_detrend_linear_1d(squeeze(data_all.([tmp.varname, '_lens2'])(loni,lati,:)), 'omitnan');                         
                 [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.lens2_det(isfinite(tmp.data_assm_det)), tmp.data_assm_det(isfinite(tmp.data_assm_det)));

                 data_all.([tmp.varname, '_corr_assm_lens2'])(loni,lati)=tmp.corr(1,2);
                 data_all.([tmp.varname, '_corr_assm_lens2_p'])(loni,lati)=tmp.corr_p(1,2);
                 data_all.([tmp.varname, '_corr_assm_lens2_det'])(loni,lati)=tmp.corr_det(1,2);
                 data_all.([tmp.varname, '_corr_assm_lens2_det_p'])(loni,lati)=tmp.corr_det_p(1,2);

            
             %% corr obs ~ hcst
                 if isfield(tmp, 'ydata_mod_obs_masked')
                    tmp.data = squeeze(data_all.([tmp.varname, '_mod_obs_masked'])(loni,lati,:));
                 else
                    tmp.data = squeeze(data_all.([tmp.varname, '_model'])(loni,lati,:));
                 end
                 tmp.data_obs = squeeze(data_all.([tmp.varname, '_obs'])(loni,lati,:));
                 [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
                 if isfield(tmp, 'ydata_mod_obs_masked')
                    tmp.data_det = Func_0028_detrend_linear_1d(squeeze(data_all.([tmp.varname, '_mod_obs_masked'])(loni,lati,:)), 'omitnan');                                              
                 else
                    tmp.data_det = Func_0028_detrend_linear_1d(squeeze(data_all.([tmp.varname, '_model'])(loni,lati,:)), 'omitnan');                         
                 end
                 tmp.data_obs_det = Func_0028_detrend_linear_1d(data_all.([tmp.varname, '_obs'])(loni,lati,:));
                 [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data_det(isfinite(tmp.data_obs_det)), tmp.data_obs_det(isfinite(tmp.data_obs_det)));
%                  [tmp.corr, tmp.corr_p]=corrcoef(squeeze(data_all.([tmp.varname, '_model'])(loni,lati,1+lyear:end-cfg.proj_year+1)), ...
%                      squeeze(data_all.([tmp.varname, '_obs'])(loni,lati,1+lyear:end)));

                 if sum(squeeze(isfinite(tmp.data_obs)))<10
                    tmp.corr=NaN(2,2);
                    tmp.corr_p=NaN(2,2);
                    tmp.corr_det=NaN(2,2);
                    tmp.corr_det_p=NaN(2,2);
                 end

                 data_all.([tmp.varname, '_corr_obs'])(loni,lati)=tmp.corr(1,2);
                 data_all.([tmp.varname, '_corr_obs_p'])(loni,lati)=tmp.corr_p(1,2);
                 data_all.([tmp.varname, '_corr_obs_det'])(loni,lati)=tmp.corr_det(1,2);
                 data_all.([tmp.varname, '_corr_obs_det_p'])(loni,lati)=tmp.corr_det_p(1,2);


             %% corr obs ~ assm
                if isfield(tmp, 'ydata_mod_obs_masked')
                    tmp.data_assm = squeeze(data_all.([tmp.varname, '_assm_obs_masked'])(loni,lati,:));
                 else
                    tmp.data_assm = squeeze(data_all.([tmp.varname, '_assm'])(loni,lati,:));
                 end
                 [tmp.corr, tmp.corr_p]=corrcoef(tmp.data_assm(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
                 
                 if sum(squeeze(isfinite(tmp.data_obs)))<10
                     tmp.data_assm_det = NaN(size(tmp.data));
                 else
                     tmp.data_assm_det = Func_0028_detrend_linear_1d(tmp.data_assm(isfinite(tmp.data_obs)));
                 end
                 [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data_assm_det(isfinite(tmp.data_obs_det)), tmp.data_obs_det(isfinite(tmp.data_obs_det)));
%                  [tmp.corr, tmp.corr_p]=corrcoef(squeeze(data_all.([tmp.varname, '_model', '_', tmp.season])(loni,lati,1+lyear:end-cfg.proj_year+1)), ...
%                      squeeze(data_all.([tmp.varname, '_obs'])(loni,lati,1+lyear:end)));

                 if sum(squeeze(isfinite(tmp.data_obs)))<10
                    tmp.corr=NaN(2,2);
                    tmp.corr_p=NaN(2,2);
                    tmp.corr_det=NaN(2,2);
                    tmp.corr_det_p=NaN(2,2);
                 end

                 data_all.([tmp.varname, '_corr_obs_assm'])(loni,lati)=tmp.corr(1,2);
                 data_all.([tmp.varname, '_corr_obs_assm_p'])(loni,lati)=tmp.corr_p(1,2);
                 data_all.([tmp.varname, '_corr_obs_assm_det'])(loni,lati)=tmp.corr_det(1,2);
                 data_all.([tmp.varname, '_corr_obs_assm_det_p'])(loni,lati)=tmp.corr_det_p(1,2);

             %% corr obs ~ lens2
                if isfield(tmp, 'ydata_mod_obs_masked')
                    tmp.data_lens2 = squeeze(data_all.([tmp.varname, '_lens2_obs_masked'])(loni,lati,:));
                 else
                    tmp.data_lens2 = squeeze(data_all.([tmp.varname, '_lens2'])(loni,lati,:));
                 end
                 [tmp.corr, tmp.corr_p]=corrcoef(tmp.data_lens2(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));

                 if sum(squeeze(isfinite(tmp.data_obs)))<10
                    tmp.corr=NaN(2,2);
                    tmp.corr_p=NaN(2,2);
                 end

                 data_all.([tmp.varname, '_corr_obs_lens2'])(loni,lati)=tmp.corr(1,2);
                 data_all.([tmp.varname, '_corr_obs_lens2_p'])(loni,lati)=tmp.corr_p(1,2);


             %% corr obs ~ hcst-lens2
                tmp.data_hcst_lens2 = squeeze(tmp.data) - squeeze(tmp.data_lens2);
                 [tmp.corr, tmp.corr_p]=corrcoef(tmp.data_hcst_lens2(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
                 if sum(squeeze(isfinite(tmp.data_obs)))<10
                    tmp.corr=NaN(2,2);
                    tmp.corr_p=NaN(2,2);
                 end
                 data_all.([tmp.varname, '_corr_obs_int'])(loni,lati)=tmp.corr(1,2);
                 data_all.([tmp.varname, '_corr_obs_int_p'])(loni,lati)=tmp.corr_p(1,2);


%                 %% AR1 (corr)
%                  tmp.data_AR1 = squeeze(data_all.([tmp.varname, '_AR1'])(loni,lati,1+lyear:end-cfg.proj_year+1));
%                  [tmp.corr_AR1, tmp.corr_AR1_p]=corrcoef(tmp.data_AR1, tmp.data_obs);
%                  if (isnan(tmp.corr_AR1(1,2)))
%                     tmp.corr_AR1(1,2)=0;
%                  end
%                  data_all.([tmp.varname, '_corr_AR1'])(loni,lati)=tmp.corr_AR1(1,2);
%                  data_all.([tmp.varname, '_corr_AR1_p'])(loni,lati)=tmp.corr_AR1_p(1,2);
%                  
%                  tmp.data_AR1_det = Func_0028_detrend_linear_1d(data_all.([tmp.varname, '_AR1'])(loni,lati,1+lyear:end-cfg.proj_year+1));
%                  [tmp.corr_AR1_det, tmp.corr_AR1_det_p]=corrcoef(tmp.data_AR1_det, tmp.data_obs_det);
%                  if (isnan(tmp.corr_AR1_det(1,2)))
%                     tmp.corr_AR1_det(1,2)=0;
%                  end
%                  data_all.([tmp.varname, '_corr_AR1_det'])(loni,lati)=tmp.corr_AR1_det(1,2);
%                  data_all.([tmp.varname, '_corr_AR1_det_p'])(loni,lati)=tmp.corr_AR1_det_p(1,2);

             else
                 data_all.([tmp.varname, '_corr_obs'])(loni,lati)=NaN;
                 data_all.([tmp.varname, '_corr_obs_p'])(loni,lati)=NaN;
                 data_all.([tmp.varname, '_corr_obs_det'])(loni,lati)=NaN;
                 data_all.([tmp.varname, '_corr_obs_det_p'])(loni,lati)=NaN;
%                  data_all.([tmp.varname, '_corr_assm'])(loni,lati)=NaN;
%                  data_all.([tmp.varname, '_corr_assm_p'])(loni,lati)=NaN;

%                 %% AR1 (NaN)
%                  data_all.([tmp.varname, '_corr_AR1'])(loni,lati)=NaN;
%                  data_all.([tmp.varname, '_corr_AR1_p'])(loni,lati)=NaN;
%                  data_all.([tmp.varname, '_corr_AR1_det'])(loni,lati)=NaN;
%                  data_all.([tmp.varname, '_corr_AR1_det_p'])(loni,lati)=NaN;

                %% ASSM (NaN)
                 data_all.([tmp.varname, '_corr_assm'])(loni,lati)=NaN;
                 data_all.([tmp.varname, '_corr_assm_p'])(loni,lati)=NaN;
                 data_all.([tmp.varname, '_corr_assm_det'])(loni,lati)=NaN;
                 data_all.([tmp.varname, '_corr_assm_det_p'])(loni,lati)=NaN;
                 data_all.([tmp.varname, '_corr_assm_int'])(loni,lati)=NaN;
                 data_all.([tmp.varname, '_corr_assm_int_p'])(loni,lati)=NaN;

                %% LENS2 (NaN)
                 data_all.([tmp.varname, '_corr_assm_lens2'])(loni,lati)=NaN;
                 data_all.([tmp.varname, '_corr_assm_lens2_p'])(loni,lati)=NaN;
                 data_all.([tmp.varname, '_corr_assm_lens2_det'])(loni,lati)=NaN;
                 data_all.([tmp.varname, '_corr_assm_lens2_det_p'])(loni,lati)=NaN;
                 %% obs ~ ASSM (NaN)
                 data_all.([tmp.varname, '_corr_obs_assm'])(loni,lati)=NaN;
                 data_all.([tmp.varname, '_corr_obs_assm_p'])(loni,lati)=NaN;
                 data_all.([tmp.varname, '_corr_obs_assm_det'])(loni,lati)=NaN;
                 data_all.([tmp.varname, '_corr_obs_assm_det_p'])(loni,lati)=NaN;
                 
                 %% obs ~ LENS2 (NaN)
                 data_all.([tmp.varname, '_corr_obs_lens2'])(loni,lati)=NaN;
                 data_all.([tmp.varname, '_corr_obs_lens2_p'])(loni,lati)=NaN;

                 %% obs ~ HCST-LENS2 (NaN)
                 data_all.([tmp.varname, '_corr_obs_int'])(loni,lati)=NaN;
                 data_all.([tmp.varname, '_corr_obs_int_p'])(loni,lati)=NaN;
             end
        end
    end
    disp('abc')
    if isfield(tmp,'ydata_mod_obs_masked')
        tmp=rmfield(tmp, 'ydata_mod_obs_masked');
    end

    tmp.lyear_str=[num2str(cfg.ly_s-1, '%02i'), '_', num2str(cfg.ly_e-1, '%02i')];
    
    if ~exist(dirs.matroot,'dir'), mkdir(dirs.matroot); end
    fig_cfg.mat_name=[dirs.matroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
        '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
        '_l', tmp.lyear_str, 'y.mat'];
    save(fig_cfg.mat_name, 'data_all')
    clear data
    clear data_all

    fprintf('%7.1f sec\n', toc(lap_time) );




    


% data.([tmp.varname, '_corr_obs', '_l3_4'])= ( data.([tmp.varname, '_corr_obs', '_l03']) + ...
%     data.([tmp.varname, '_corr_obs', '_l04']) ) / 2;
% data.([tmp.varname, '_corr_obs', '_l5_9'])= ( data.([tmp.varname, '_corr_obs', '_l05']) + ...
%     data.([tmp.varname, '_corr_obs', '_l06']) + ..._
%     data.([tmp.varname, '_corr_obs', '_l07']) + ...
%     data.([tmp.varname, '_corr_obs', '_l08']) + ...
%     data.([tmp.varname, '_corr_obs', '_l09']) ) / 5;
% 
% data.([tmp.varname, '_corr_assm', '_l3_4'])= ( data.([tmp.varname, '_corr_assm', '_l03']) + ...
%     data.([tmp.varname, '_corr_assm', '_l04']) ) / 2;
% data.([tmp.varname, '_corr_assm', '_l5_9'])= ( data.([tmp.varname, '_corr_assm', '_l05']) + ...
%     data.([tmp.varname, '_corr_assm', '_l06']) + ...
%     data.([tmp.varname, '_corr_assm', '_l07']) + ...
%     data.([tmp.varname, '_corr_assm', '_l08']) + ...
%     data.([tmp.varname, '_corr_assm', '_l09']) ) / 5;


    
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
