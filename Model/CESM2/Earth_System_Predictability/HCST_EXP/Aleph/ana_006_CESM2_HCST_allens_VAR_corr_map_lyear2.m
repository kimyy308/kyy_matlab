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
    otherwise
        tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
end
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Common', tmp.fs, 'mca']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Common', tmp.fs, 'order']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

%% model configuration

cfg.vars = {'TS', 'PSL', 'PRECT', 'SST'};
% cfg.vars = {'PRECT', 'PSL'};
% cfg.vars = {'SST'};
% cfg.vars = {'NPP', 'GPP', 'RAIN', 'TOTVEGC', 'FIRE'};
% cfg.vars = {'TWS'};
% cfg.vars = {'SOILWATER_10CM'};
cfg.vars = {'TS', 'PRECT', 'SST', 'PSL'};
cfg.vars = {'PRECT', 'SST', 'PSL'};

cfg.vars={'COL_FIRE_CLOSS', 'COL_FIRE_NLOSS','DSTFLXT','FAREA_BURNED','FIRE','FPSN','GPP','NEP', ...
'NFIRE','NPP','Q2M','QH2OSFC','QOVER','QRGWL','QRUNOFF','QSOIL','QSOIL_ICE','RAIN','RH2M', ...
'SNOW','SOILICE','SOILLIQ','SOILWATER_10CM','TLAI','TOTSOILICE','TOTSOILLIQ','TOTVEGC','TWS'};

cfg.vars={'FAREA_BURNED','FIRE','FPSN','GPP','NEP', ...
'NFIRE','NPP','Q2M','QH2OSFC','QOVER','QRGWL','QRUNOFF','QSOIL','QSOIL_ICE','RAIN','RH2M', ...
'SNOW','SOILICE','SOILLIQ','SOILWATER_10CM','TLAI','TOTSOILICE','TOTSOILLIQ','TOTVEGC','TWS'};

cfg.vars={'SSH', 'photoC_TOT_zint', 'photoC_TOT_zint_100m', 'IRON_FLUX', 'HMXL', 'HBLT', ...
    'NO3', 'PO4'}; 
    
cfg.vars={'SiO3', 'Fe', 'PD', 'SALT', 'TEMP', 'UVEL', 'VVEL', 'WVEL', 'zooC'};

cfg.vars={'SOILWATER_10CM', 'TLAI', 'FAREA_BURNED', 'COL_FIRE_CLOSS', 'TWS', 'SSH', 'photoC_TOT_zint', 'photoC_TOT_zint_100m'};
% cfg.vars={'SSH'};
% cfg.vars={'FAREA_BURNED'};
cfg.vars={'SOILWATER_10CM'};
cfg.vars={'NO3', 'SALT', 'mul_VVEL_NO3', 'mul_WVEL_NO3'};
cfg.vars={'mul_VVEL_NO3', 'mul_WVEL_NO3','mul_UVEL_NO3'};
cfg.vars={'pCO2SURF'};
cfg.vars={'NO3', 'SALT', 'PD'};

cfg.vars={'GPP'};
cfg.vars={'photoC_TOT_zint_100m'};
cfg.vars ={'ALK', 'DpCO2', 'PH', 'DIC'};
cfg.vars={'FAREA_BURNED'};
cfg.vars={'Jint_100m_NO3', 'tend_zint_100m_NO3'};
cfg.vars={'subt_tend_zint_100m_NO3_Jint_100m_NO3'};
cfg.vars={'NO3','TEMP'};
cfg.vars={'zooC'};
cfg.vars={'FG_CO2'};
cfg.vars={'DIC'};
cfg.vars={'DIC_ALT_CO2', 'DpCO2_ALT_CO2'};
% tmp.dimids= [1, 2, 4];
cfg.vlayer=1; % surf, vertical slice 
% cfg.vlayer=1:10; % 10layer. don't put more than 15
% cfg.vlayer=10; % 100m, vertical slice 
% cfg.vlayer=15; %150m
% cfg.vlayer=27; %305m
cfg.vlayer_1st=min(cfg.vlayer);
cfg.vlayer_cnt=max(cfg.vlayer)-cfg.vlayer_1st+1;


for vari=1:length(cfg.vars)
    cfg.var=cfg.vars{vari};
    cfg.obs_name=f_obs_name(cfg.var);
    cfg.obs_fname_mid=f_obs_name_mid(cfg.var);
    cfg.obs_varname=f_obs_varname(cfg.var);
    cfg.comp=Func_0025_CESM2_cmpname_var(cfg.var);
    cfg.obs_fname_module=f_obs_fname_module(cfg.comp);

%     cfg.obs_iyears=1960:2020;
    
    cfg.obs_iyears=1965:2020;

    cfg.obs_iyears2=f_obs_iyears(cfg.var);
  
    % dirs.hcstroot=['/mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive/', cfg.comp, '/', cfg.var];
    % dirs.obsroot=['/Volumes/kyy_raid/kimyy/Observation/', cfg.obs_name, '/monthly_reg_', cfg.obs_fname_module(2:4)];
    % dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', cfg.var];
    % dirs.lens2root=['/mnt/lustre/proj/earth.system.predictability/LENS2/archive_analysis/', cfg.comp, '/', cfg.var];
    % dirs.assmroot=['/mnt/lustre/proj/earth.system.predictability/ASSM_EXP/archive_analysis/', cfg.comp, '/', cfg.var];
    
    dirs.hcstroot=['/mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive_yearly_transfer/', cfg.comp, '/', cfg.var];
    dirs.obsroot=['/mnt/lustre/proj/kimyy/Observation/', cfg.obs_name];
    dirs.matroot=['/mnt/lustre/proj/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp,'/', cfg.var];
    dirs.lens2root=['/mnt/lustre/proj/earth.system.predictability/LENS2/archive_yearly_transfer/', cfg.comp, '/', cfg.var];
    dirs.assmroot=['/mnt/lustre/proj/earth.system.predictability/ASSM_EXP/archive_yearly_transfer/', cfg.comp, '/', cfg.var];    
    
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
            [grid.tlat, grid.tlong]=meshgrid(grid.lat, grid.lon);
            grid.area=ncread(tmp.gridname, 'AREA');
            grid.lfrac=ncread(tmp.gridname, 'LANDFRAC');
            grid.area=grid.area.*grid.lfrac;
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
    %% lead year loop
    for lyear=0:cfg.proj_year-1
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
        data.([tmp.varname, '_obs'])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data.([tmp.varname, '_assm'])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        

    %     data.([tmp.varname, '_assm'])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        
%% read variables
        for iyear=min(cfg.iyears):max(cfg.iyears)
            tmp.iyear_str=num2str(iyear, '%04i');
            cfg.casename=['ensmean_', cfg.casename_m, '_i', tmp.iyear_str];
            tmp.fy=iyear+lyear;
            tmp.fy_str=num2str(tmp.fy, '%04i');
%% yearly filename    
%% HCST mean
            cfg.mod_fnm=[dirs.datadir, tmp.fs, 'ens_all_i',tmp.iyear_str, tmp.fs, ...
            tmp.varname, '_', cfg.gnm, '.hcst.', cfg.casename, cfg.obs_fname_module, tmp.fy_str, '.nc'];
            ncid=netcdf.open(cfg.mod_fnm, 'NOWRITE');
            tmpvarid=netcdf.inqVarID(ncid,tmp.fvarname);
            [tmp.fvarname,tmp.xtype,tmp.dimids,tmp.natts]= netcdf.inqVar(ncid,tmpvarid);
            if length(tmp.dimids)>3
                 tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 0], [grid.nlon grid.nlat cfg.vlayer_cnt 1]));
%                      tmp.dd=tmp.dd.*grid.mask_ocn;
                 tmp.dd(abs(tmp.dd)>1e30)=NaN;
                 ddd=mean(tmp.dd,3,'omitnan'); % depth mean
                 tmp.ymean(1:grid.nlon,1:grid.nlat) = ddd; %depth averaged value
            else
                tmp.dd=  netcdf.getVar(ncid,tmpvarid);
                tmp.dd(abs(tmp.dd)>1e30)=NaN;
                tmp.ymean(1:grid.nlon,1:grid.nlat) = tmp.dd;
            end
            netcdf.close(ncid);            

%% HCST ensstd (spread)
            cfg.casename_std=['ensstd_', cfg.casename_m, '_i', tmp.iyear_str];
            cfg.mod_fnm=[dirs.datadir, tmp.fs, 'ens_all_i',tmp.iyear_str, tmp.fs, ...
                tmp.varname, '_', cfg.gnm, '.hcst.', cfg.casename_std, cfg.obs_fname_module, tmp.fy_str '.nc'];
            if exist(cfg.mod_fnm)~=0
                ncid=netcdf.open(cfg.mod_fnm, 'NOWRITE');
                 tmpvarid=netcdf.inqVarID(ncid,tmp.fvarname);
                [tmp.fvarname,tmp.xtype,tmp.dimids,tmp.natts]= netcdf.inqVar(ncid,tmpvarid);
                if length(tmp.dimids)>3
                     tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 0], [grid.nlon grid.nlat cfg.vlayer_cnt 1]));
    %                      tmp.dd=tmp.dd.*grid.mask_ocn;
                     tmp.dd(abs(tmp.dd)>1e30)=NaN;
                     ddd=mean(tmp.dd,3,'omitnan'); % depth mean
                     tmp.ystde(1:grid.nlon,1:grid.nlat) = ddd; %depth averaged value
                else
                    tmp.dd=  netcdf.getVar(ncid,tmpvarid);
                    tmp.dd(abs(tmp.dd)>1e30)=NaN;
                    tmp.ystde(1:grid.nlon,1:grid.nlat) = tmp.dd;
                end
                netcdf.close(ncid);  
            else
                tmp.ystde(1:grid.nlon,1:grid.nlat)=0;
            end

%% LENS2 mean
            if tmp.fy <= cfg.max_iy
                cfg.lens2_fnm=[dirs.lens2dir, tmp.fs, ...
                    tmp.varname, '_y_ensmean_', tmp.fy_str, '.nc'];
                ncid=netcdf.open(cfg.lens2_fnm, 'NOWRITE');
                tmpvarid=netcdf.inqVarID(ncid,tmp.fvarname);
                if length(tmp.dimids)>3
                     tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 0], [grid.nlon grid.nlat cfg.vlayer_cnt 1]));
                     tmp.dd(abs(tmp.dd)>1e30)=NaN;
                     ddd=mean(tmp.dd,3,'omitnan'); % depth mean
                     tmp.ymean_lens2(1:grid.nlon,1:grid.nlat) = ddd; %depth averaged value
                else
                    tmp.dd=netcdf.getVar(ncid,tmpvarid);
                     tmp.dd(abs(tmp.dd)>1e30)=NaN;                        
                    tmp.ymean_lens2(1:grid.nlon,1:grid.nlat) =tmp.dd;
                end
                netcdf.close(ncid);
            else
                tmp.ymean_lens2(1:grid.nlon,1:grid.nlat) = NaN;                    
            end

%% LENS2 std
            if tmp.fy <= cfg.max_iy
                cfg.lens2_fnm=[dirs.lens2dir, tmp.fs, ...
                    tmp.varname, '_y_ensstd_', tmp.fy_str, '.nc'];
                if exist(cfg.lens2_fnm)~=0
                    ncid=netcdf.open(cfg.lens2_fnm, 'NOWRITE');
                    tmpvarid=netcdf.inqVarID(ncid,tmp.fvarname);
                    if length(tmp.dimids)>3
                         tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 0], [grid.nlon grid.nlat cfg.vlayer_cnt 1]));
                         tmp.dd(abs(tmp.dd)>1e30)=NaN;
                         ddd=mean(tmp.dd,3,'omitnan'); % depth mean
                         tmp.ystde_lens2(1:grid.nlon,1:grid.nlat) = ddd; %depth averaged value
                    else
                        tmp.dd=netcdf.getVar(ncid,tmpvarid);
                         tmp.dd(abs(tmp.dd)>1e30)=NaN;                        
                        tmp.ystde_lens2(1:grid.nlon,1:grid.nlat) =tmp.dd;
                    end
                    netcdf.close(ncid);
                else
                    tmp.ystde_lens2(1:grid.nlon,1:grid.nlat)=0;
                end
            else
                tmp.ystde_lens2(1:grid.nlon,1:grid.nlat) = NaN;                    
            end

%% OBS
            
            if tmp.fy <= cfg.max_iy
                for mi=1:12
                    tmp.m_str=num2str(mi,'%02i');
                    if strcmp(cfg.obs_name, 'ERA5')==1
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, cfg.obs_varname, '/monthly_reg_', cfg.obs_fname_module(2:4),tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
                    elseif (strcmp(cfg.obs_name, 'CMEMS')==1 && strcmp(tmp.varname, 'SOILWATER_10CM')==1)
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'SOIL_MOISTURE/COMBINED', '/monthly_reg_cam', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
                    elseif (strcmp(cfg.obs_name, 'CMEMS')==1 && strcmp(tmp.varname, 'photoC_TOT_zint_100m')==1)
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'monthly_reg_pop/PP/', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
                    elseif (strcmp(cfg.obs_name, 'CMEMS')==1 && strcmp(tmp.varname, 'photoC_TOT_zint')==1)
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'monthly_reg_pop/PP/', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
                    elseif (strcmp(cfg.obs_name, 'NOAA')==1 && strcmp(tmp.varname, 'TWS')==1)
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'TSW', '/monthly_reg_cam', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
                    elseif strcmp(cfg.obs_name, 'GPCC')==1 || strcmp(cfg.obs_name, 'ORNL_DAAC')==1
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'monthly_reg_cam' ,tmp.fs,cfg.obs_fname_mid,tmp.fy_str,tmp.m_str,'.nc'];
                    elseif (strcmp(cfg.obs_name, 'GLEAM')==1 && strcmp(tmp.varname, 'SOILWATER_10CM')==1)
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, cfg.obs_varname, tmp.fs, 'monthly_reg_cam', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
                    elseif (strcmp(cfg.obs_name, 'NOAA')==1 && strcmp(tmp.varname, 'TLAI')==1)
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, cfg.obs_varname, tmp.fs, 'monthly_reg_cam', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
                    elseif (strcmp(cfg.obs_name, 'MODIS')==1 && strcmp(tmp.varname, 'FAREA_BURNED')==1)
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, cfg.obs_varname, tmp.fs, 'monthly_reg_cam', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
                    elseif (strcmp(cfg.obs_name, 'AVHRR')==1 && strcmp(tmp.varname, 'FAREA_BURNED')==1)
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'AVHRR-LTDR', tmp.fs, 'monthly_reg_cam', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];  
                    elseif (strcmp(cfg.obs_name, 'GFED')==1 && strcmp(tmp.varname, 'COL_FIRE_CLOSS')==1)
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'FIRE_CLOSS', tmp.fs, 'monthly_reg_cam', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
                    elseif (strcmp(cfg.obs_name, 'VGPM')==1)
%                         cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'netcdf_regrid/s_vgpm', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'netcdf_regrid/ensmean', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];                                                
                    else
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'monthly_reg_', cfg.obs_fname_module(2:4),tmp.fs,cfg.obs_fname_mid,tmp.fy_str,tmp.m_str,'.nc'];
                    end
                    if exist(cfg.obs_fnm)~=0
                        ncid=netcdf.open(cfg.obs_fnm, 'NOWRITE');
                        tmpvarid=netcdf.inqVarID(ncid,cfg.obs_varname);
                        tmp.dd =  netcdf.getVar(ncid,tmpvarid);
                        tmp.dd=double(tmp.dd);
                        if (strcmp(cfg.obs_name, 'ERA5')==1 || strcmp(tmp.varname, 'TLAI')==1)
                            tmp.add_offset=netcdf.getAtt(ncid,tmpvarid,'add_offset');
                            tmp.scale_factor=netcdf.getAtt(ncid,tmpvarid,'scale_factor');
                            tmp.dd=tmp.dd.*tmp.scale_factor+tmp.add_offset;
                        end
                         tmp.dd(abs(tmp.dd)>1e30)=NaN;
                        if strcmp(cfg.obs_name, 'GPCC')==1
                            if strcmp(tmp.varname, 'PRECT')
                                tmp.dd=tmp.dd./1000.0./86400/eomday(tmp.fy,mi);
                            elseif strcmp(tmp.varname, 'RAIN')
                                tmp.dd=tmp.dd./86400/eomday(tmp.fy,mi);
                            end
                        elseif (strcmp(cfg.obs_name, 'CMEMS')==1 && strcmp(tmp.varname, 'SOILWATER_10CM')==1)
                            tmp.dd=tmp.dd.*1000.*(10./3);
                        elseif (strcmp(cfg.obs_name, 'CMEMS')==1 && strcmp(tmp.varname, 'SSH')==1)
                            tmp.dd=tmp.dd./100; % cm -> m
                        elseif (strcmp(cfg.obs_name, 'GLEAM')==1 && strcmp(tmp.varname, 'SOILWATER_10CM')==1)
                            tmp.dd=tmp.dd.*(1000)./10; % 1. m^3 -> kg, 2. 100cm(1m) -> 10cm,  m3/m3 -> 10cm(surface) soil kg/m2
                        elseif (strcmp(cfg.obs_name, 'NOAA')==1 && strcmp(tmp.varname, 'TWS')==1)
                            tmp.dd(tmp.dd<=0)=NaN;
                            tmp.dd(tmp.dd>740)=NaN;
                        elseif strcmp(tmp.varname, 'FAREA_BURNED')
%                                 tmp.dd=tmp.dd./86400./eomday(tmp.fy,mi)./grid.area;
                                tmp.dd=tmp.dd./86400./grid.area;
                        elseif strcmp(tmp.varname, 'COL_FIRE_CLOSS')
                            tmp.dd=tmp.dd./86400./eomday(tmp.fy,mi);
                        end
                        tmp.ydata_obs(1:grid.nlon,1:grid.nlat,mi) = tmp.dd;
                        netcdf.close(ncid);
                        
                    else
                        tmp.ydata_obs(1:grid.nlon,1:grid.nlat,mi) = NaN;                    
                    end
                end
                tmp.ymean_obs(1:grid.nlon,1:grid.nlat)=mean(tmp.ydata_obs,3);

            else
                tmp.ymean_obs(1:grid.nlon,1:grid.nlat) = NaN;          
            end
%             tmp.ydata_obs(:,:,mon)=ncread(cfg.obs_fnm, cfg.obs_varname);

%% ASSM
            if tmp.fy <= cfg.max_iy
                cfg.assm_fnm=[dirs.assmdir, tmp.fs, tmp.varname, '_ensmean_', tmp.fy_str, '.nc'];
                ncid=netcdf.open(cfg.assm_fnm, 'NOWRITE');
                tmpvarid=netcdf.inqVarID(ncid,tmp.fvarname);
                if length(tmp.dimids)>3
                    tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 0], [grid.nlon grid.nlat cfg.vlayer_cnt 1]));
%                          ddd=tmp.dd.*grid.dz_res; %% weight depth
%                          ddd=sum(ddd,3); % depth sum
%                          ddd=ddd./sum(grid.dz_res(:)); % divide total depth for avg
                     tmp.dd(abs(tmp.dd)>1e30)=NaN;
                     ddd=mean(tmp.dd,3,'omitnan'); % depth mean
                     tmp.ymean_assm(1:grid.nlon,1:grid.nlat) = ddd; %depth averaged value             
                else
                    tmp.ymean_assm(1:grid.nlon,1:grid.nlat) = netcdf.getVar(ncid,tmpvarid);
                end
                if(strcmp(tmp.varname(1:2),'mu')~=1 && strcmp(tmp.varname(1:2),'su')~=1)
                    data.units=netcdf.getAtt(ncid,tmpvarid,'units');
                else
                    data.units=' ';
                end
                netcdf.close(ncid);
            else
                tmp.ymean_assm(1:grid.nlon,1:grid.nlat) = NaN;                    
            end

%% ASSM std
            if tmp.fy <= cfg.max_iy
                cfg.assm_fnm=[dirs.assmdir, tmp.fs, tmp.varname, '_ensstd_', tmp.fy_str, '.nc'];
                if exist(cfg.assm_fnm)~=0
                    ncid=netcdf.open(cfg.assm_fnm, 'NOWRITE');
                    tmpvarid=netcdf.inqVarID(ncid,tmp.fvarname);
                    if length(tmp.dimids)>3
                        tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 0], [grid.nlon grid.nlat cfg.vlayer_cnt 1]));
    %                          ddd=tmp.dd.*grid.dz_res; %% weight depth
    %                          ddd=sum(ddd,3); % depth sum
    %                          ddd=ddd./sum(grid.dz_res(:)); % divide total depth for avg
                         tmp.dd(abs(tmp.dd)>1e30)=NaN;
                         ddd=mean(tmp.dd,3,'omitnan'); % depth mean
                         tmp.ystde_assm(1:grid.nlon,1:grid.nlat) = ddd; %depth averaged value             
                    else
                        tmp.ystde_assm(1:grid.nlon,1:grid.nlat) = netcdf.getVar(ncid,tmpvarid);
                    end
                    netcdf.close(ncid);
                else
                    tmp.ystde_assm(1:grid.nlon,1:grid.nlat)=0;
                end
            else
                tmp.ystde_assm(1:grid.nlon,1:grid.nlat) = NaN;                    
            end

%% observation mask (scarced data like marine npp)            
%             if (strcmp(cfg.var, 'sumChl') || strcmp(cfg.var, 'TS')) %for Chls, it uses all available data
%             if (strcmp(cfg.var, 'sumChl')) %for Chls, it uses all available data
            if (strcmp(cfg.var, 'sumChl')) %for Chls, it uses all available data
                tmp.obs_mask(:,:)=tmp.ymean_obs(:,:)./tmp.ymean_obs(:,:);
                tmp.ymean_mod_obs_masked(:,:) = tmp.ymean(:,:) .* tmp.obs_mask(:,:); % masking with available observation
                tmp.ymean_assm_obs_masked(:,:) = tmp.ymean_assm(:,:) .* tmp.obs_mask(:,:); % masking with available observation
                tmp.ymean_lens2_obs_masked(:,:) = tmp.ymean_lens2(:,:) .* tmp.obs_mask(:,:); % masking with available observation
                tmp.ystde_mod_obs_masked(:,:) = tmp.ystde(:,:) .* tmp.obs_mask(:,:); % masking with available observation
                tmp.ystde_assm_obs_masked(:,:) = tmp.ystde_assm(:,:) .* tmp.obs_mask(:,:); % masking with available observation
                tmp.ystde_lens2_obs_masked(:,:) = tmp.ystde_lens2(:,:) .* tmp.obs_mask(:,:); % masking with available observation
            end
    
    %             %% AR1 (integration)
    %             tmp.ydata_AR1(:,:,mon) = ...
    %                 Func_0030_AR1_prog(tmp.all_data_obs(:,:,(tmp.iind-1)*12+1), data.([tmp.varname, '_AC_lag1']), lyear*12+mon-1, 0);
            
    %         data.time=ncread(cfg.datafilename, 'time');
    %         tmp.ymean= mean(tmp.ydata,3);
    %         data.([tmp.varname, '_model', '_l', tmp.lyear_str])(:,:,iyear-min(cfg.iyears)+1+lyear)= tmp.ymean;
    %         tmp.ymean_obs= mean(tmp.ydata_obs,3);
    %         data.([tmp.varname, '_obs'])(:,:,iyear-min(cfg.iyears)+1)= tmp.ymean_obs;
    
    %         if ( strcmp(cfg.var, 'sumChl') || strcmp(cfg.var, 'TS')) %for Chls, it uses all available data
    %             tmp.ymean= mean(tmp.ydata,3, 'omitnan');
    %             tmp.ymean_obs= mean(tmp.ydata_obs,3, 'omitnan');
    %         else
           
            tmp.ymean_assm(tmp.ymean_assm>10e30)=NaN;
            tmp.ystde_assm(tmp.ystde_assm>10e30)=NaN;

            switch cfg.comp
                case 'ocn'
                    tmp.ymean = tmp.ymean .* grid.mask_ocn;
                    tmp.ystde = tmp.ystde .* grid.mask_ocn;
                    if isfield(tmp, 'ymean_mod_obs_masked')
                        tmp.ymean_mod_obs_masked = tmp.ymean_mod_obs_masked .* grid.mask_ocn;
                        tmp.ymean_assm_obs_masked = tmp.ymean_assm_obs_masked .* grid.mask_ocn;
                        tmp.ymean_lens2_obs_masked = tmp.ymean_lens2_obs_masked .* grid.mask_ocn;
                        tmp.ystde_mod_obs_masked = tmp.ystde_mod_obs_masked .* grid.mask_ocn;
                        tmp.ystde_assm_obs_masked = tmp.ystde_assm_obs_masked .* grid.mask_ocn;
                        tmp.ystde_lens2_obs_masked = tmp.ystde_lens2_obs_masked .* grid.mask_ocn;
                    end
                    tmp.ymean_obs = tmp.ymean_obs .* grid.mask_ocn;
                    tmp.ymean_assm = tmp.ymean_assm .* grid.mask_ocn;
                    tmp.ymean_lens2 = tmp.ymean_lens2 .* grid.mask_ocn;
%                     tmp.ystde_obs = tmp.ystde_obs .* grid.mask_ocn;
                    tmp.ystde_assm = tmp.ystde_assm .* grid.mask_ocn;
                    tmp.ystde_lens2 = tmp.ystde_lens2 .* grid.mask_ocn;    
            end
    %         end

%% save variables as structure
            data.([tmp.varname, '_model', '_l', tmp.lyear_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean;
            data.([tmp.varname, '_model_stde', '_l', tmp.lyear_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ystde;
            if isfield(tmp, 'ymean_mod_obs_masked')
                data.([tmp.varname, '_mod_obs_masked', '_l', tmp.lyear_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_mod_obs_masked;
                data.([tmp.varname, '_assm_obs_masked', '_l', tmp.lyear_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_assm_obs_masked;
                data.([tmp.varname, '_lens2_obs_masked', '_l', tmp.lyear_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_lens2_obs_masked;
                data.([tmp.varname, '_mod_stde_obs_masked', '_l', tmp.lyear_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ystde_mod_obs_masked;
                data.([tmp.varname, '_assm_stde_obs_masked', '_l', tmp.lyear_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ystde_assm_obs_masked;
                data.([tmp.varname, '_lens2_stde_obs_masked', '_l', tmp.lyear_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ystde_lens2_obs_masked;
            end
            data.([tmp.varname, '_obs'])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_obs;
            data.([tmp.varname, '_lens2', '_l', tmp.lyear_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_lens2;
            data.([tmp.varname, '_assm'])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_assm;
            data.([tmp.varname, '_lens2_stde', '_l', tmp.lyear_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ystde_lens2;
            data.([tmp.varname, '_assm_stde'])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ystde_assm;
            
            
    
        end  %% initialized year loop end
        if (strcmp(cfg.var, 'SST' ))
            data.([tmp.varname, '_model', '_l', tmp.lyear_str]) = data.([tmp.varname, '_model', '_l', tmp.lyear_str]) - 273.15;
            data.([tmp.varname, '_lens2', '_l', tmp.lyear_str]) = data.([tmp.varname, '_lens2', '_l', tmp.lyear_str]) - 273.15;
            data.([tmp.varname, '_assm']) = data.([tmp.varname, '_assm']) - 273.15;
            data.([tmp.varname, '_model', '_l', tmp.lyear_str])(data.([tmp.varname, '_model', '_l', tmp.lyear_str])<-30) = NaN;
            data.([tmp.varname, '_lens2', '_l', tmp.lyear_str])(data.([tmp.varname, '_lens2', '_l', tmp.lyear_str])<-30) = NaN;
            data.([tmp.varname, '_assm'])(data.([tmp.varname, '_assm'])<-30) = NaN;
            data.([tmp.varname, '_obs'])(data.([tmp.varname, '_obs'])==-999)=NaN;
        end
        disp('data read finished')
        fprintf('%7.1f sec\n', toc(lap_time) );
    
        for loni=1:size(data.([tmp.varname, '_obs']),1)
            for lati=1:size(data.([tmp.varname, '_obs']),2)
                if sum(isfinite(data.([tmp.varname, '_obs'])(loni,lati,:)))<floor(length(cfg.obs_iyears2).*0.8)
                    data.([tmp.varname, '_obs'])(loni,lati,:)=NaN;
                end
            end
        end


   %% get temporal std of ensmean
        data.([tmp.varname, '_assm_stdt'])= std(data.([tmp.varname, '_assm']),0,3,'omitnan');
        data.([tmp.varname, '_obs_stdt'])= std(data.([tmp.varname, '_obs']),0,3,'omitnan');
        data.([tmp.varname, '_model_stdt_l', tmp.lyear_str])= std(data.([tmp.varname, '_model_l', tmp.lyear_str]),0,3,'omitnan');
        data.([tmp.varname, '_lens2_stdt_l', tmp.lyear_str])= std(data.([tmp.varname, '_lens2_l', tmp.lyear_str]),0,3,'omitnan');

   %% get model error(bias) with ASSM 
        data.([tmp.varname, '_model_err', '_l', tmp.lyear_str]) = ...
            data.([tmp.varname, '_model', '_l', tmp.lyear_str]) - data.([tmp.varname, '_assm']);
        data.([tmp.varname, '_lens2_err', '_l', tmp.lyear_str]) = ...
            data.([tmp.varname, '_lens2', '_l', tmp.lyear_str]) - data.([tmp.varname, '_assm']);

   %% get model error(bias) with obs 
        data.([tmp.varname, '_model_err_obs', '_l', tmp.lyear_str]) = ...
            data.([tmp.varname, '_model', '_l', tmp.lyear_str]) - data.([tmp.varname, '_obs']);
        data.([tmp.varname, '_lens2_err_obs', '_l', tmp.lyear_str]) = ...
            data.([tmp.varname, '_lens2', '_l', tmp.lyear_str]) - data.([tmp.varname, '_obs']);
        data.([tmp.varname, '_assm_err_obs', '_l', tmp.lyear_str]) = ...
            data.([tmp.varname, '_assm']) - data.([tmp.varname, '_obs']);

   %% get RMSE with ASSM (drift-removed)
%         data.([tmp.varname, '_model_rmse', '_l', tmp.lyear_str]) = ...
%             sqrt(mean(data.([tmp.varname, '_model_err', '_l', tmp.lyear_str]).^2, 3, 'omitnan'));
%         data.([tmp.varname, '_lens2_rmse', '_l', tmp.lyear_str]) = ...
%             sqrt(mean(data.([tmp.varname, '_lens2_err', '_l', tmp.lyear_str]).^2, 3, 'omitnan'));
        data.([tmp.varname, '_model_rmse', '_l', tmp.lyear_str]) = ...
            sqrt( mean( ...
            (data.([tmp.varname, '_model_err', '_l', tmp.lyear_str]) ...
            - mean(data.([tmp.varname, '_model_err', '_l', tmp.lyear_str]), 3, 'omitnan')).^2, ...
            3, 'omitnan') );
        data.([tmp.varname, '_lens2_rmse', '_l', tmp.lyear_str]) = ...
            sqrt( mean( ...
            (data.([tmp.varname, '_lens2_err', '_l', tmp.lyear_str]) ...
            - mean(data.([tmp.varname, '_lens2_err', '_l', tmp.lyear_str]), 3, 'omitnan')).^2, ...
            3, 'omitnan') );
   %% get RMSE with obs (drift-removed)
        data.([tmp.varname, '_model_rmse_obs', '_l', tmp.lyear_str]) = ...
            sqrt( mean( ...
            (data.([tmp.varname, '_model_err_obs', '_l', tmp.lyear_str]) ...
            - mean(data.([tmp.varname, '_model_err_obs', '_l', tmp.lyear_str]), 3, 'omitnan')).^2, ...
            3, 'omitnan') );
        data.([tmp.varname, '_lens2_rmse_obs', '_l', tmp.lyear_str]) = ...
            sqrt( mean( ...
            (data.([tmp.varname, '_lens2_err_obs', '_l', tmp.lyear_str]) ...
            - mean(data.([tmp.varname, '_lens2_err_obs', '_l', tmp.lyear_str]), 3, 'omitnan')).^2, ...
            3, 'omitnan') );
        data.([tmp.varname, '_assm_rmse_obs', '_l', tmp.lyear_str]) = ...
            sqrt( mean( ...
            (data.([tmp.varname, '_assm_err_obs', '_l', tmp.lyear_str]) ...
            - mean(data.([tmp.varname, '_assm_err_obs', '_l', tmp.lyear_str]), 3, 'omitnan')).^2, ...
            3, 'omitnan') );

    %% get nRMSE with ASSM (drift-removed)
        data.([tmp.varname, '_model_nrmse', '_l', tmp.lyear_str])= ...
            data.([tmp.varname, '_model_rmse', '_l', tmp.lyear_str])./data.([tmp.varname, '_assm_stdt']);
        data.([tmp.varname, '_lens2_nrmse', '_l', tmp.lyear_str])= ...
            data.([tmp.varname, '_lens2_rmse', '_l', tmp.lyear_str])./data.([tmp.varname, '_assm_stdt']);

    %% get nRMSE with obs (drift-removed)
        data.([tmp.varname, '_model_nrmse_obs', '_l', tmp.lyear_str])= ...
            data.([tmp.varname, '_model_rmse_obs', '_l', tmp.lyear_str])./data.([tmp.varname, '_obs_stdt']);
        data.([tmp.varname, '_lens2_nrmse_obs', '_l', tmp.lyear_str])= ...
            data.([tmp.varname, '_lens2_rmse_obs', '_l', tmp.lyear_str])./data.([tmp.varname, '_obs_stdt']);
        data.([tmp.varname, '_assm_nrmse_obs', '_l', tmp.lyear_str])= ...
            data.([tmp.varname, '_assm_rmse_obs', '_l', tmp.lyear_str])./data.([tmp.varname, '_obs_stdt']);
    
    
    disp('RMSEs calculated')
    fprintf('%7.1f sec\n', toc(lap_time) );

    %% get hindcast climatology as the function of lead year
       switch cfg.obs_name
           case 'GPCC'
               cfg.clim_ys=1965;
               cfg.clim_ye=2019;
           otherwise
               cfg.clim_ys=1965;
               cfg.clim_ye=2020;
       end
       cfg.clim_tlen = (cfg.clim_ys-1959)-lyear:(cfg.clim_ye-2020)+cfg.len_t_y-lyear;
       cfg.clim_tlen2=length(cfg.clim_tlen);
       if lyear==0
            cfg.clim_tlen_AR1 = 1:cfg.len_t_y; % for AR1 integration, it should start from 1960 to integrate 1965, LY5 
       end

        data.([tmp.varname, '_model_clim', '_l', tmp.lyear_str]) = mean(data.([tmp.varname, '_model', '_l', tmp.lyear_str])(:,:,cfg.clim_tlen),3, 'omitnan'); %1964 ~ 2020
        data.([tmp.varname, '_assm_clim']) = mean(data.([tmp.varname, '_assm'])(:,:,cfg.clim_tlen),3, 'omitnan'); %1964 ~ 2020
        data.([tmp.varname, '_lens2_clim']) = mean(data.([tmp.varname, '_lens2', '_l', tmp.lyear_str])(:,:,cfg.clim_tlen),3, 'omitnan'); %1964 ~ 2020
        data.([tmp.varname, '_obs_clim']) = mean(data.([tmp.varname, '_obs'])(:,:,cfg.clim_tlen),3, 'omitnan'); %1964 ~ 2020
        
        if lyear==0
            tmp.assm_clim_AR1_period = mean(data.([tmp.varname, '_assm'])(:,:,cfg.clim_tlen_AR1),3, 'omitnan');
            tmp.obs_clim_AR1_period = mean(data.([tmp.varname, '_obs'])(:,:,cfg.clim_tlen_AR1),3, 'omitnan');
        end

   %% get anomaly values
        data2.([tmp.varname, '_model_ano', '_l', tmp.lyear_str]) = ...
            data.([tmp.varname, '_model', '_l', tmp.lyear_str])(:,:,cfg.clim_tlen) - data.([tmp.varname, '_model_clim', '_l', tmp.lyear_str]);
        data2.([tmp.varname, '_assm_ano', '_l', tmp.lyear_str]) = ...
            data.([tmp.varname, '_assm'])(:,:,cfg.clim_tlen) - data.([tmp.varname, '_assm_clim']);
        data2.([tmp.varname, '_lens2_ano', '_l', tmp.lyear_str]) = ...
            data.([tmp.varname, '_lens2', '_l', tmp.lyear_str])(:,:,cfg.clim_tlen) - data.([tmp.varname, '_lens2_clim']);
        data2.([tmp.varname, '_obs_ano', '_l', tmp.lyear_str]) = ...
            data.([tmp.varname, '_obs'])(:,:,cfg.clim_tlen) - data.([tmp.varname, '_obs_clim']);
        if lyear==0
            data2_l0_AR1.assm_ano_full = data.([tmp.varname, '_assm'])(:,:,cfg.clim_tlen_AR1) - tmp.assm_clim_AR1_period;
            data2_l0_AR1.obs_ano_full = data.([tmp.varname, '_obs'])(:,:,cfg.clim_tlen_AR1) - tmp.obs_clim_AR1_period;
        end

    %% get anomaly_detrended values
        data2.([tmp.varname, '_model_ano_det2', '_l', tmp.lyear_str]) = NaN(grid.nlon, grid.nlat, cfg.clim_tlen2);
        data2.([tmp.varname, '_assm_ano_det2', '_l', tmp.lyear_str]) = NaN(grid.nlon, grid.nlat, cfg.clim_tlen2);
        data2.([tmp.varname, '_lens2_ano_det2', '_l', tmp.lyear_str]) = NaN(grid.nlon, grid.nlat, cfg.clim_tlen2);
        data2.([tmp.varname, '_obs_ano_det2', '_l', tmp.lyear_str]) = NaN(grid.nlon, grid.nlat, cfg.clim_tlen2);

        data2.([tmp.varname, '_model_trend2', '_l', tmp.lyear_str]) = NaN(grid.nlon, grid.nlat);
        data2.([tmp.varname, '_assm_trend2', '_l', tmp.lyear_str]) = NaN(grid.nlon, grid.nlat);
        data2.([tmp.varname, '_lens2_trend2', '_l', tmp.lyear_str]) = NaN(grid.nlon, grid.nlat);
        data2.([tmp.varname, '_obs_trend2', '_l', tmp.lyear_str]) = NaN(grid.nlon, grid.nlat);

        data2.([tmp.varname, '_model_ano_det3', '_l', tmp.lyear_str]) = NaN(grid.nlon, grid.nlat, cfg.clim_tlen2);
        data2.([tmp.varname, '_assm_ano_det3', '_l', tmp.lyear_str]) = NaN(grid.nlon, grid.nlat, cfg.clim_tlen2);
        data2.([tmp.varname, '_lens2_ano_det3', '_l', tmp.lyear_str]) = NaN(grid.nlon, grid.nlat, cfg.clim_tlen2);
        data2.([tmp.varname, '_obs_ano_det3', '_l', tmp.lyear_str]) = NaN(grid.nlon, grid.nlat, cfg.clim_tlen2);

        data2.([tmp.varname, '_model_trend3', '_l', tmp.lyear_str]) = NaN(grid.nlon, grid.nlat);
        data2.([tmp.varname, '_assm_trend3', '_l', tmp.lyear_str]) = NaN(grid.nlon, grid.nlat);
        data2.([tmp.varname, '_lens2_trend3', '_l', tmp.lyear_str]) = NaN(grid.nlon, grid.nlat);
        data2.([tmp.varname, '_obs_trend3', '_l', tmp.lyear_str]) = NaN(grid.nlon, grid.nlat);

        for loni=1:grid.nlon
            for lati=1:grid.nlat
                 if (isnan(data.([tmp.varname, '_assm'])(loni,lati,1))~=1 && sum(data.([tmp.varname, '_model_l', tmp.lyear_str])(loni,lati,:),'omitnan')~=0)
                    tmp.data_obs = squeeze(data2.([tmp.varname, '_obs_ano', '_l', tmp.lyear_str])(loni,lati,:));
                    tmp.data_assm = squeeze(data2.([tmp.varname, '_assm_ano', '_l', tmp.lyear_str])(loni,lati,:));
                    tmp.data_lens2 = squeeze(data2.([tmp.varname, '_lens2_ano', '_l', tmp.lyear_str])(loni,lati,:));
                    tmp.data_model = squeeze(data2.([tmp.varname, '_model_ano', '_l', tmp.lyear_str])(loni,lati,:));
                    
                    [data2.([tmp.varname, '_obs_ano_det2', '_l', tmp.lyear_str])(loni,lati,:), ...
                        data2.([tmp.varname, '_obs_trend2', '_l', tmp.lyear_str])(loni,lati)] ...
                        = Func_0028_detrend_linear_1d(tmp.data_obs, 'omitnan');
                    [data2.([tmp.varname, '_assm_ano_det2', '_l', tmp.lyear_str])(loni,lati,:), ...
                        data2.([tmp.varname, '_assm_trend2', '_l', tmp.lyear_str])(loni,lati)] ...
                        = Func_0028_detrend_linear_1d(tmp.data_assm, 'omitnan');
                    [data2.([tmp.varname, '_lens2_ano_det2', '_l', tmp.lyear_str])(loni,lati,:), ...
                        data2.([tmp.varname, '_lens2_trend2', '_l', tmp.lyear_str])(loni,lati)] ...
                        = Func_0028_detrend_linear_1d(tmp.data_lens2, 'omitnan');
                    [data2.([tmp.varname, '_model_ano_det2', '_l', tmp.lyear_str])(loni,lati,:), ...
                        data2.([tmp.varname, '_model_trend2', '_l', tmp.lyear_str])(loni,lati)] ...
                        = Func_0028_detrend_linear_1d(tmp.data_model, 'omitnan');

                     tmp.data_assm2=tmp.data_assm;
                    tmp.data_assm2(isnan(tmp.data_obs))=NaN;
                    tmp.data_lens22=tmp.data_lens2;
                    tmp.data_lens22(isnan(tmp.data_obs))=NaN;
                    tmp.data_model2=tmp.data_model;
                    tmp.data_model2(isnan(tmp.data_obs))=NaN;

                    [data2.([tmp.varname, '_assm_ano_det3', '_l', tmp.lyear_str])(loni,lati,:), ...
                        data2.([tmp.varname, '_assm_trend3', '_l', tmp.lyear_str])(loni,lati)] ...
                        = Func_0028_detrend_linear_1d(tmp.data_assm2, 'omitnan');
                    [data2.([tmp.varname, '_lens2_ano_det3', '_l', tmp.lyear_str])(loni,lati,:), ...
                        data2.([tmp.varname, '_lens2_trend3', '_l', tmp.lyear_str])(loni,lati)] ...
                        = Func_0028_detrend_linear_1d(tmp.data_lens22, 'omitnan');
                    [data2.([tmp.varname, '_model_ano_det3', '_l', tmp.lyear_str])(loni,lati,:), ...
                        data2.([tmp.varname, '_model_trend3', '_l', tmp.lyear_str])(loni,lati)] ...
                        = Func_0028_detrend_linear_1d(tmp.data_model2, 'omitnan');
                 end
            end
        end
        data2.([tmp.varname, '_obs_ano_det2', '_l', tmp.lyear_str]) = ...
            data2.([tmp.varname, '_obs_ano_det2', '_l', tmp.lyear_str]) - ...
            mean(data2.([tmp.varname, '_obs_ano_det2', '_l', tmp.lyear_str]), 3, 'omitnan');
        data2.([tmp.varname, '_assm_ano_det2', '_l', tmp.lyear_str]) = ...
            data2.([tmp.varname, '_assm_ano_det2', '_l', tmp.lyear_str]) - ...
            mean(data2.([tmp.varname, '_assm_ano_det2', '_l', tmp.lyear_str]), 3, 'omitnan');
        data2.([tmp.varname, '_lens2_ano_det2', '_l', tmp.lyear_str]) = ...
            data2.([tmp.varname, '_lens2_ano_det2', '_l', tmp.lyear_str]) - ...
            mean(data2.([tmp.varname, '_lens2_ano_det2', '_l', tmp.lyear_str]), 3, 'omitnan');
        data2.([tmp.varname, '_model_ano_det2', '_l', tmp.lyear_str]) = ...
            data2.([tmp.varname, '_model_ano_det2', '_l', tmp.lyear_str]) - ...
            mean(data2.([tmp.varname, '_model_ano_det2', '_l', tmp.lyear_str]), 3, 'omitnan');

        data2.([tmp.varname, '_assm_ano_det3', '_l', tmp.lyear_str]) = ...
            data2.([tmp.varname, '_assm_ano_det3', '_l', tmp.lyear_str]) - ...
            mean(data2.([tmp.varname, '_assm_ano_det3', '_l', tmp.lyear_str]), 3, 'omitnan');
        data2.([tmp.varname, '_lens2_ano_det3', '_l', tmp.lyear_str]) = ...
            data2.([tmp.varname, '_lens2_ano_det3', '_l', tmp.lyear_str]) - ...
            mean(data2.([tmp.varname, '_lens2_ano_det3', '_l', tmp.lyear_str]), 3, 'omitnan');
        data2.([tmp.varname, '_model_ano_det3', '_l', tmp.lyear_str]) = ...
            data2.([tmp.varname, '_model_ano_det3', '_l', tmp.lyear_str]) - ...
            mean(data2.([tmp.varname, '_model_ano_det3', '_l', tmp.lyear_str]), 3, 'omitnan');

    %% get temporal std of ensmean
        data2.([tmp.varname, '_assm_ano_det2_stdt'])= std(data2.([tmp.varname, '_assm_ano_det2', '_l', tmp.lyear_str]),0,3,'omitnan');
        data2.([tmp.varname, '_obs_ano_det2_stdt'])= std(data2.([tmp.varname, '_obs_ano_det2', '_l', tmp.lyear_str]),0,3,'omitnan');
        data2.([tmp.varname, '_model_ano_det2_stdt_l', tmp.lyear_str])= std(data2.([tmp.varname, '_model_ano_det2', '_l', tmp.lyear_str]),0,3,'omitnan');
        data2.([tmp.varname, '_lens2_ano_det2_stdt_l', tmp.lyear_str])= std(data2.([tmp.varname, '_lens2_ano_det2', '_l', tmp.lyear_str]),0,3,'omitnan');

        data2.([tmp.varname, '_assm_ano_det3_stdt'])= std(data2.([tmp.varname, '_assm_ano_det3', '_l', tmp.lyear_str]),0,3,'omitnan');
        data2.([tmp.varname, '_model_ano_det3_stdt_l', tmp.lyear_str])= std(data2.([tmp.varname, '_model_ano_det3', '_l', tmp.lyear_str]),0,3,'omitnan');
        data2.([tmp.varname, '_lens2_ano_det3_stdt_l', tmp.lyear_str])= std(data2.([tmp.varname, '_lens2_ano_det3', '_l', tmp.lyear_str]),0,3,'omitnan');

   %% get model error(bias) with ASSM (between anomalies)
        data2.([tmp.varname, '_model_err2', '_l', tmp.lyear_str]) = ...
            data2.([tmp.varname, '_model_ano_det2', '_l', tmp.lyear_str]) - data2.([tmp.varname, '_assm_ano_det2', '_l', tmp.lyear_str]);
        data2.([tmp.varname, '_lens2_err2', '_l', tmp.lyear_str]) = ...
            data2.([tmp.varname, '_lens2_ano_det2', '_l', tmp.lyear_str]) - data2.([tmp.varname, '_assm_ano_det2', '_l', tmp.lyear_str]);
   %% get model error(bias) with obs (between anomalies)
        data2.([tmp.varname, '_model_err2_obs', '_l', tmp.lyear_str]) = ...
            data2.([tmp.varname, '_model_ano_det3', '_l', tmp.lyear_str]) - data2.([tmp.varname, '_obs_ano_det2', '_l', tmp.lyear_str]);
        data2.([tmp.varname, '_lens2_err2_obs', '_l', tmp.lyear_str]) = ...
            data2.([tmp.varname, '_lens2_ano_det3', '_l', tmp.lyear_str]) - data2.([tmp.varname, '_obs_ano_det2', '_l', tmp.lyear_str]);
        data2.([tmp.varname, '_assm_err2_obs', '_l', tmp.lyear_str]) = ...
            data2.([tmp.varname, '_assm_ano_det3', '_l', tmp.lyear_str]) - data2.([tmp.varname, '_obs_ano_det2', '_l', tmp.lyear_str]);

   %% get RMSE with ASSM based on anomaly (climatology removed, detrended, drift-removed)
        data2.([tmp.varname, '_model_rmse2', '_l', tmp.lyear_str]) = ...
            sqrt( mean(data2.([tmp.varname, '_model_err2', '_l', tmp.lyear_str]).^2, 3, 'omitnan') );
        data2.([tmp.varname, '_lens2_rmse2', '_l', tmp.lyear_str]) = ...
            sqrt( mean(data2.([tmp.varname, '_lens2_err2', '_l', tmp.lyear_str]).^2, 3, 'omitnan') );
   %% get RMSE with obs based on anomaly (climatology removed, detrended, drift-removed)
        data2.([tmp.varname, '_assm_rmse2_obs', '_l', tmp.lyear_str]) = ...
            sqrt( mean(data2.([tmp.varname, '_assm_err2_obs', '_l', tmp.lyear_str]).^2, 3, 'omitnan') );
        data2.([tmp.varname, '_model_rmse2_obs', '_l', tmp.lyear_str]) = ...
            sqrt( mean(data2.([tmp.varname, '_model_err2_obs', '_l', tmp.lyear_str]).^2, 3, 'omitnan') );
        data2.([tmp.varname, '_lens2_rmse2_obs', '_l', tmp.lyear_str]) = ...
            sqrt( mean(data2.([tmp.varname, '_lens2_err2_obs', '_l', tmp.lyear_str]).^2, 3, 'omitnan') );

    %% get nRMSE with ASSM based on anomaly (climatology removed, detrended, drift-removed)
        data2.([tmp.varname, '_model_nrmse2', '_l', tmp.lyear_str])= ...
            data2.([tmp.varname, '_model_rmse2', '_l', tmp.lyear_str])./data2.([tmp.varname, '_assm_ano_det2_stdt']);
        data2.([tmp.varname, '_lens2_nrmse2', '_l', tmp.lyear_str])= ...
            data2.([tmp.varname, '_lens2_rmse2', '_l', tmp.lyear_str])./data2.([tmp.varname, '_assm_ano_det2_stdt']);

    %% get nRMSE with obs based on anomaly (climatology removed, detrended, drift-removed)
        data2.([tmp.varname, '_model_nrmse2_obs', '_l', tmp.lyear_str])= ...
            data2.([tmp.varname, '_model_rmse2_obs', '_l', tmp.lyear_str])./data2.([tmp.varname, '_obs_ano_det2_stdt']);
        data2.([tmp.varname, '_lens2_nrmse2_obs', '_l', tmp.lyear_str])= ...
            data2.([tmp.varname, '_lens2_rmse2_obs', '_l', tmp.lyear_str])./data2.([tmp.varname, '_obs_ano_det2_stdt']);
        data2.([tmp.varname, '_assm_nrmse2_obs', '_l', tmp.lyear_str])= ...
            data2.([tmp.varname, '_assm_rmse2_obs', '_l', tmp.lyear_str])./data2.([tmp.varname, '_obs_ano_det2_stdt']);
    
    %% get MSSS (mean-square skill score, 1- (MSE_hcst/MSE_ref)
        data2.([tmp.varname, '_MSSS_assm_ano', '_l', tmp.lyear_str]) = ...
            1 -  data2.([tmp.varname, '_model_rmse2', '_l', tmp.lyear_str]).^2 ...
                 ./ data2.([tmp.varname, '_lens2_rmse2', '_l', tmp.lyear_str]).^2;
        data2.([tmp.varname, '_MSSS_obs_ano', '_l', tmp.lyear_str]) = ...
            1 -  data2.([tmp.varname, '_model_rmse2_obs', '_l', tmp.lyear_str]).^2 ...
                 ./ data2.([tmp.varname, '_lens2_rmse2_obs', '_l', tmp.lyear_str]).^2;
    
    disp('anomaly-related variables calculated')
    fprintf('%7.1f sec\n', toc(lap_time) );

    %% get svd modes (first 10 modes)
        tmp.svd_modes=5;
        % lens2 <-> hcst
        [data2.([tmp.varname, '_svd_lens2_model_lvmap_left', '_l', tmp.lyear_str]), ...
            data2.([tmp.varname, '_svd_lens2_model_pcs_left', '_l', tmp.lyear_str]), ...
            data2.([tmp.varname, '_svd_lens2_model_lvmap_right', '_l', tmp.lyear_str]), ...
            data2.([tmp.varname, '_svd_lens2_model_pcs_right', '_l', tmp.lyear_str]), ...
            data2.([tmp.varname, '_svd_lens2_model_lambda', '_l', tmp.lyear_str]), ...
            data2.([tmp.varname, '_svd_lens2_model_scf', '_l', tmp.lyear_str])] = ...
            mca( data2.([tmp.varname, '_lens2_ano', '_l', tmp.lyear_str]), ... 
                 data2.([tmp.varname, '_model_ano', '_l', tmp.lyear_str]), tmp.svd_modes);
        for orderi=1:tmp.svd_modes
            data2.([tmp.varname, '_svd_order_lens2_model', '_l', tmp.lyear_str]) = ...
                order(mean(abs(data2.([tmp.varname, '_svd_lens2_model_pcs_left', '_l', tmp.lyear_str])(orderi,:)),'all'));
            data2.([tmp.varname, '_svd_lens2_model_lvmap_left', '_l', tmp.lyear_str])(:,:,orderi) = ...
                data2.([tmp.varname, '_svd_lens2_model_lvmap_left', '_l', tmp.lyear_str])(:,:,orderi) ...
                .* 10.^data2.([tmp.varname, '_svd_order_lens2_model', '_l', tmp.lyear_str]);
            data2.([tmp.varname, '_svd_lens2_model_pcs_left', '_l', tmp.lyear_str])(orderi,:) = ...
                data2.([tmp.varname, '_svd_lens2_model_pcs_left', '_l', tmp.lyear_str])(orderi,:) ...
                .* 10.^(-data2.([tmp.varname, '_svd_order_lens2_model', '_l', tmp.lyear_str]));
            data2.([tmp.varname, '_svd_lens2_model_lvmap_right', '_l', tmp.lyear_str])(:,:,orderi) = ...
                data2.([tmp.varname, '_svd_lens2_model_lvmap_right', '_l', tmp.lyear_str])(:,:,orderi) ...
                .* 10.^data2.([tmp.varname, '_svd_order_lens2_model', '_l', tmp.lyear_str]);
            data2.([tmp.varname, '_svd_lens2_model_pcs_right', '_l', tmp.lyear_str])(orderi,:) = ...
                data2.([tmp.varname, '_svd_lens2_model_pcs_right', '_l', tmp.lyear_str])(orderi,:) ...
                .* 10.^(-data2.([tmp.varname, '_svd_order_lens2_model', '_l', tmp.lyear_str]));
        end
        tmp.pcs=data2.([tmp.varname, '_svd_lens2_model_pcs_right', '_l', tmp.lyear_str])(1,:);
        tmp.pcs=reshape(tmp.pcs,[1,1,length(tmp.pcs)]);
        data2.([tmp.varname, '_svd_lens2_model_right_1mode_recon', '_l', tmp.lyear_str]) = ...
            data2.([tmp.varname, '_svd_lens2_model_lvmap_right', '_l', tmp.lyear_str])(:,:,1).*tmp.pcs(1,1,:);

        % lens2 <-> assm
        [data2.([tmp.varname, '_svd_lens2_assm_lvmap_left', '_l', tmp.lyear_str]), ...
            data2.([tmp.varname, '_svd_lens2_assm_pcs_left', '_l', tmp.lyear_str]), ...
            data2.([tmp.varname, '_svd_lens2_assm_lvmap_right', '_l', tmp.lyear_str]), ...
            data2.([tmp.varname, '_svd_lens2_assm_pcs_right', '_l', tmp.lyear_str]), ...
            data2.([tmp.varname, '_svd_lens2_assm_lambda', '_l', tmp.lyear_str]), ...
            data2.([tmp.varname, '_svd_lens2_assm_scf', '_l', tmp.lyear_str])] = ...
            mca( data2.([tmp.varname, '_lens2_ano', '_l', tmp.lyear_str]), ... 
                 data2.([tmp.varname, '_assm_ano', '_l', tmp.lyear_str]), tmp.svd_modes);
        for orderi=1:tmp.svd_modes
            data2.([tmp.varname, '_svd_order_lens2_assm', '_l', tmp.lyear_str]) = ...
                order(mean(abs(data2.([tmp.varname, '_svd_lens2_assm_pcs_left', '_l', tmp.lyear_str])(orderi,:)),'all'));
            data2.([tmp.varname, '_svd_lens2_assm_lvmap_left', '_l', tmp.lyear_str])(:,:,orderi) = ...
                data2.([tmp.varname, '_svd_lens2_assm_lvmap_left', '_l', tmp.lyear_str])(:,:,orderi) ...
                .* 10.^data2.([tmp.varname, '_svd_order_lens2_assm', '_l', tmp.lyear_str]);
            data2.([tmp.varname, '_svd_lens2_assm_pcs_left', '_l', tmp.lyear_str])(orderi,:) = ...
                data2.([tmp.varname, '_svd_lens2_assm_pcs_left', '_l', tmp.lyear_str])(orderi,:) ...
                .* 10.^(-data2.([tmp.varname, '_svd_order_lens2_assm', '_l', tmp.lyear_str]));
            data2.([tmp.varname, '_svd_lens2_assm_lvmap_right', '_l', tmp.lyear_str])(:,:,orderi) = ...
                data2.([tmp.varname, '_svd_lens2_assm_lvmap_right', '_l', tmp.lyear_str])(:,:,orderi) ...
                .* 10.^data2.([tmp.varname, '_svd_order_lens2_assm', '_l', tmp.lyear_str]);
            data2.([tmp.varname, '_svd_lens2_assm_pcs_right', '_l', tmp.lyear_str])(orderi,:)= ...
                data2.([tmp.varname, '_svd_lens2_assm_pcs_right', '_l', tmp.lyear_str])(orderi,:) ...
                .* 10.^(-data2.([tmp.varname, '_svd_order_lens2_assm', '_l', tmp.lyear_str]));
        end
        tmp.pcs(1,1,:)=data2.([tmp.varname, '_svd_lens2_assm_pcs_right', '_l', tmp.lyear_str])(1,:);
        tmp.pcs=reshape(tmp.pcs,[1,1,length(tmp.pcs)]);
        data2.([tmp.varname, '_svd_lens2_assm_right_1mode_recon', '_l', tmp.lyear_str]) = ...
            data2.([tmp.varname, '_svd_lens2_assm_lvmap_right', '_l', tmp.lyear_str])(:,:,1).*tmp.pcs(1,1,:);

        % lens2 <-> obs
        if isnan(sum( data2.([tmp.varname, '_obs_ano', '_l', tmp.lyear_str])(:)))
            data2.([tmp.varname, '_svd_lens2_obs_lvmap_left', '_l', tmp.lyear_str]) = ...
                NaN(size(data2.([tmp.varname, '_svd_lens2_assm_lvmap_left', '_l', tmp.lyear_str])));
            data2.([tmp.varname, '_svd_lens2_obs_pcs_left', '_l', tmp.lyear_str]) = ...
                NaN(size(data2.([tmp.varname, '_svd_lens2_assm_pcs_left', '_l', tmp.lyear_str])));
            data2.([tmp.varname, '_svd_lens2_obs_lvmap_right', '_l', tmp.lyear_str]) = ...
                NaN(size(data2.([tmp.varname, '_svd_lens2_assm_lvmap_left', '_l', tmp.lyear_str])));
            data2.([tmp.varname, '_svd_lens2_obs_pcs_right', '_l', tmp.lyear_str]) = ...
                NaN(size(data2.([tmp.varname, '_svd_lens2_assm_pcs_left', '_l', tmp.lyear_str])));
            data2.([tmp.varname, '_svd_lens2_obs_lambda', '_l', tmp.lyear_str]) = ...
                NaN(size(data2.([tmp.varname, '_svd_lens2_model_lambda', '_l', tmp.lyear_str])));
            data2.([tmp.varname, '_svd_lens2_obs_scf', '_l', tmp.lyear_str]) = ...
                NaN(size(data2.([tmp.varname, '_svd_lens2_model_scf', '_l', tmp.lyear_str])));
        else
            [data2.([tmp.varname, '_svd_lens2_obs_lvmap_left', '_l', tmp.lyear_str]), ...
                data2.([tmp.varname, '_svd_lens2_obs_pcs_left', '_l', tmp.lyear_str]), ...
                data2.([tmp.varname, '_svd_lens2_obs_lvmap_right', '_l', tmp.lyear_str]), ...
                data2.([tmp.varname, '_svd_lens2_obs_pcs_right', '_l', tmp.lyear_str]), ...
                data2.([tmp.varname, '_svd_lens2_obs_lambda', '_l', tmp.lyear_str]), ...
                data2.([tmp.varname, '_svd_lens2_obs_scf', '_l', tmp.lyear_str])] = ...
                mca( data2.([tmp.varname, '_lens2_ano', '_l', tmp.lyear_str]), ... 
                     data2.([tmp.varname, '_obs_ano', '_l', tmp.lyear_str]), tmp.svd_modes);
            for orderi=1:tmp.svd_modes
                data2.([tmp.varname, '_svd_order_lens2_obs', '_l', tmp.lyear_str]) = ...
                    order(mean(abs(data2.([tmp.varname, '_svd_lens2_obs_pcs_left', '_l', tmp.lyear_str])(orderi,:)),'all'));
                data2.([tmp.varname, '_svd_lens2_obs_lvmap_left', '_l', tmp.lyear_str])(:,:,orderi) = ...
                    data2.([tmp.varname, '_svd_lens2_obs_lvmap_left', '_l', tmp.lyear_str])(:,:,orderi) ...
                    .* 10.^data2.([tmp.varname, '_svd_order_lens2_obs', '_l', tmp.lyear_str]);
                data2.([tmp.varname, '_svd_lens2_obs_pcs_left', '_l', tmp.lyear_str])(orderi,:) = ...
                    data2.([tmp.varname, '_svd_lens2_obs_pcs_left', '_l', tmp.lyear_str])(orderi,:) ...
                    .* 10.^(-data2.([tmp.varname, '_svd_order_lens2_obs', '_l', tmp.lyear_str]));
                data2.([tmp.varname, '_svd_lens2_obs_lvmap_right', '_l', tmp.lyear_str])(:,:,orderi) = ...
                    data2.([tmp.varname, '_svd_lens2_obs_lvmap_right', '_l', tmp.lyear_str])(:,:,orderi) ...
                    .* 10.^data2.([tmp.varname, '_svd_order_lens2_obs', '_l', tmp.lyear_str]);
                data2.([tmp.varname, '_svd_lens2_obs_pcs_right', '_l', tmp.lyear_str])(orderi,:) = ...
                    data2.([tmp.varname, '_svd_lens2_obs_pcs_right', '_l', tmp.lyear_str])(orderi,:)...
                    .* 10.^(-data2.([tmp.varname, '_svd_order_lens2_obs', '_l', tmp.lyear_str]));
            end
        end
      
        tmp.pcs(1,1,:)=data2.([tmp.varname, '_svd_lens2_obs_pcs_right', '_l', tmp.lyear_str])(1,:);
        tmp.pcs=reshape(tmp.pcs,[1,1,length(tmp.pcs)]);
        data2.([tmp.varname, '_svd_lens2_obs_right_1mode_recon', '_l', tmp.lyear_str]) = ...
            data2.([tmp.varname, '_svd_lens2_obs_lvmap_right', '_l', tmp.lyear_str])(:,:,1).*tmp.pcs(1,1,:);

    disp('svd finished')
    fprintf('%7.1f sec\n', toc(lap_time) );

   %% AR1 parameter estimation, getting prognostic value (assm based), should be estimated from anomaly
        data2.([tmp.varname, '_AR1', '_l', tmp.lyear_str])=NaN(grid.nlon, grid.nlat, length(cfg.clim_tlen));  %% AR1 initialization
        data2.([tmp.varname, '_obs_AR1', '_l', tmp.lyear_str])=NaN(grid.nlon, grid.nlat, length(cfg.clim_tlen));  %% AR1 initialization
        
%         if lyear==0
            for loni=1:grid.nlon
                for lati=1:grid.nlat
                    if isnan(data2.([tmp.varname, '_assm_ano', '_l', tmp.lyear_str])(loni,lati))~=1
%                         tmp.init_finite=squeeze(data2.([tmp.varname, '_assm_ano', '_l', tmp.lyear_str])(loni,lati,:));
                        tmp.init_finite=data2_l0_AR1.assm_ano_full(loni, lati, cfg.clim_ys-1960-lyear:cfg.clim_ye-1960-lyear);
                        tmp.init_finite=tmp.init_finite(isfinite(tmp.init_finite));
                        tmp.init_finite=tmp.init_finite-mean(tmp.init_finite);
                        
                        [tmp.coef, tmp.noise] = aryule(squeeze(tmp.init_finite), 1);
                        data2.([tmp.varname,'_assm_AR1_coef'])(loni,lati)=tmp.coef(2);
                        data2.([tmp.varname,'_assm_AR1_noise'])(loni,lati)=tmp.noise;
                    else
                        data2.([tmp.varname,'_assm_AR1_coef'])(loni,lati)=NaN;
                        data2.([tmp.varname,'_assm_AR1_noise'])(loni,lati)=NaN;
                    end
                end
            end
%             data2_l0=data2;
%         end
        for loni=1:grid.nlon
            for lati=1:grid.nlat
                if isnan(data2.([tmp.varname, '_assm_ano', '_l', tmp.lyear_str])(loni,lati))~=1
%                     prog=data2_l0.([tmp.varname,'_assm_ano', '_l', tmp.lyear_str])(loni,lati,:); %initialize
%                     prog=tmp.assm_ano(loni,lati,1:end-1); %initialize  
                      prog=data2_l0_AR1.assm_ano_full(loni, lati, cfg.clim_ys-1960-lyear:cfg.clim_ye-1960-lyear); %initialize
                    for i=1:lyear+1 % integration of autoregressive model
                        prog = -prog .* (data2.([tmp.varname,'_assm_AR1_coef'])(loni,lati)) ...
                            + (data2.([tmp.varname,'_assm_AR1_noise'])(loni,lati));
                    end
                    data2.([tmp.varname,'_AR1_l',tmp.lyear_str])(loni,lati,:)=prog;
                end
            end
        end

   %% AR1 parameter estimation, getting prognostic value (obs based)
%         if lyear==0 
            for loni=1:grid.nlon
                for lati=1:grid.nlat
                    if isnan(data2.([tmp.varname, '_obs_ano', '_l', tmp.lyear_str])(loni,lati))~=1
%                         tmp.init_finite=squeeze(data2.([tmp.varname, '_obs_ano', '_l', tmp.lyear_str])(loni,lati,:));
%                         tmp.init_finite=tmp.obs_ano(loni,lati,1:end-1);
                        tmp.init_finite=data2_l0_AR1.obs_ano_full(loni, lati, cfg.clim_ys-1960-lyear:cfg.clim_ye-1960-lyear);
                        tmp.init_finite=tmp.init_finite(isfinite(tmp.init_finite));
                        tmp.init_finite= tmp.init_finite-mean(tmp.init_finite);
                        [tmp.coef, tmp.noise] = aryule(squeeze(tmp.init_finite), 1);
                        data2.([tmp.varname,'_obs_AR1_coef'])(loni,lati)=tmp.coef(2);
                        data2.([tmp.varname,'_obs_AR1_noise'])(loni,lati)=tmp.noise;
                    else
                        data2.([tmp.varname,'_obs_AR1_coef'])(loni,lati)=NaN;
                        data2.([tmp.varname,'_obs_AR1_noise'])(loni,lati)=NaN;
                    end
                end
            end
%             data2_l0=data2;
%         end
        for loni=1:grid.nlon
            for lati=1:grid.nlat
                if isnan(data2.([tmp.varname, '_obs_ano', '_l', tmp.lyear_str])(loni,lati))~=1
%                     prog=data2_l0.([tmp.varname,'_obs_ano', '_l', tmp.lyear_str])(loni,lati,:);
%                     prog=tmp.obs_ano(loni,lati,1:end-1);        
                    prog=data2_l0_AR1.obs_ano_full(loni, lati, cfg.clim_ys-1960-lyear:cfg.clim_ye-1960-lyear); 
                    for i=1:lyear+1
                        prog = -prog .* (data2.([tmp.varname,'_obs_AR1_coef'])(loni,lati)) ...
                            + (data2.([tmp.varname,'_obs_AR1_noise'])(loni,lati));
                    end
                    data2.([tmp.varname,'_obs_AR1_l',tmp.lyear_str])(loni,lati,:)=prog;
                end
            end
        end

    disp('AR1 estimation, integration finished')
    fprintf('%7.1f sec\n', toc(lap_time) );
        %% get correlation coefficient
        for loni=1:grid.nlon
            for lati=1:grid.nlat
                 if (isnan(data.([tmp.varname, '_assm'])(loni,lati,1))~=1 && sum(data.([tmp.varname, '_model_l', tmp.lyear_str])(loni,lati,:), 'omitnan')~=0)
    
                 %% corr assm
                     tmp.data_assm = data.([tmp.varname, '_assm'])(loni,lati,:);
                     tmp.data = squeeze(data.([tmp.varname, '_model', '_l', tmp.lyear_str])(loni,lati,:));
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_assm)), tmp.data_assm(isfinite(tmp.data_assm)));
                     [tmp.data_det, tmp.data_trend] = Func_0028_detrend_linear_1d(squeeze(data.([tmp.varname, '_model', '_l', tmp.lyear_str])(loni,lati,:)), 'omitnan');                         
                     [tmp.data_assm_det, tmp.data_assm_trend] = Func_0028_detrend_linear_1d(squeeze(data.([tmp.varname, '_assm'])(loni,lati,:)), 'omitnan');
                     [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data_det(isfinite(tmp.data_assm_det)), tmp.data_assm_det(isfinite(tmp.data_assm_det)));

                     data.([tmp.varname, '_corr_assm', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_assm_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);
                     data.([tmp.varname, '_corr_assm_det', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det(1,2);
                     data.([tmp.varname, '_corr_assm_det_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det_p(1,2);
                     
                     data.([tmp.varname, '_assm_trend', '_l', tmp.lyear_str])(loni,lati)=tmp.data_assm_trend;
                     data.([tmp.varname, '_model_trend', '_l', tmp.lyear_str])(loni,lati)=tmp.data_trend;

                 %% hcst-lens2 ~ assm
                     tmp.data_hcst_lens2 = squeeze(data.([tmp.varname, '_model', '_l', tmp.lyear_str])(loni,lati,:)) - ...
                                    squeeze(data.([tmp.varname, '_lens2', '_l', tmp.lyear_str])(loni,lati,:));
                     [tmp.corr_int, tmp.corr_int_p]=corrcoef(tmp.data_hcst_lens2(isfinite(tmp.data_assm)), tmp.data_assm(isfinite(tmp.data_assm)));
    
                     
                     data.([tmp.varname, '_corr_assm_int', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_int(1,2);
                     data.([tmp.varname, '_corr_assm_int_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_int_p(1,2);

                 %% hcst-lens2 ~ assm-lens2
                     tmp.data_hcst_lens2 = squeeze(data.([tmp.varname, '_model', '_l', tmp.lyear_str])(loni,lati,:)) - ...
                                    squeeze(data.([tmp.varname, '_lens2', '_l', tmp.lyear_str])(loni,lati,:));
                     tmp.data_assm_lens2 = squeeze(data.([tmp.varname, '_assm'])(loni,lati,:)) - ...
                                    squeeze(data.([tmp.varname, '_lens2', '_l', tmp.lyear_str])(loni,lati,:));
                     [tmp.corr_int, tmp.corr_int_p]=corrcoef(tmp.data_hcst_lens2(isfinite(tmp.data_assm_lens2)), tmp.data_assm_lens2(isfinite(tmp.data_assm_lens2)));
    
                     
                     data.([tmp.varname, '_corr_assm_int2', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_int(1,2);
                     data.([tmp.varname, '_corr_assm_int2_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_int_p(1,2);
    
                 %% corr lens2 ~ assm
                     tmp.lens2 = squeeze(data.([tmp.varname, '_lens2', '_l', tmp.lyear_str])(loni,lati,:));
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.lens2(isfinite(tmp.data_assm)), tmp.data_assm(isfinite(tmp.data_assm)));
                     [tmp.data_lens2_det, tmp.data_lens2_trend] = Func_0028_detrend_linear_1d(squeeze(data.([tmp.varname, '_lens2', '_l', tmp.lyear_str])(loni,lati,:)), 'omitnan');                         
                     [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data_lens2_det(isfinite(tmp.data_assm_det)), tmp.data_assm_det(isfinite(tmp.data_assm_det)));

                     data.([tmp.varname, '_corr_assm_lens2', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_assm_lens2_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);
                     data.([tmp.varname, '_corr_assm_lens2_det', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det(1,2);
                     data.([tmp.varname, '_corr_assm_lens2_det_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det_p(1,2);

                     data.([tmp.varname, '_lens2_trend', '_l', tmp.lyear_str])(loni,lati)=tmp.data_lens2_trend;

%                  %% corr AR1(assm) ~ assm
%                     tmp.AR1 = squeeze(data2.([tmp.varname, '_AR1', '_l', tmp.lyear_str])(loni,lati,:));
%                     tmp.AR1_mask=squeeze(tmp.AR1./tmp.AR1);
%                     tmp.data_assm_mask=squeeze(tmp.data_assm./tmp.data_assm);
%                      [tmp.corr, tmp.corr_p]=corrcoef(tmp.AR1.*tmp.data_assm_mask.*tmp.AR1_mask, ...
%                          squeeze(tmp.data_assm).*tmp.data_assm_mask.*tmp.AR1_mask, 'Rows', 'complete');
%                      data.([tmp.varname, '_corr_assm_AR1', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
%                      data.([tmp.varname, '_corr_assm_AR1_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);
    
                
                 %% corr obs ~ hcst
                     if isfield(tmp, 'ydata_mod_obs_masked')
                        tmp.data = squeeze(data.([tmp.varname, '_mod_obs_masked', '_l', tmp.lyear_str])(loni,lati,:));
                     else
                        tmp.data = squeeze(data.([tmp.varname, '_model', '_l', tmp.lyear_str])(loni,lati,:));
                     end
                     tmp.data_obs = squeeze(data.([tmp.varname, '_obs'])(loni,lati,:));
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
                     if isfield(tmp, 'ydata_mod_obs_masked')
                        tmp.data_det = Func_0028_detrend_linear_1d(squeeze(data.([tmp.varname, '_mod_obs_masked', '_l', tmp.lyear_str])(loni,lati,:)), 'omitnan');                                              
                     else
                        tmp.data_det = Func_0028_detrend_linear_1d(squeeze(data.([tmp.varname, '_model', '_l', tmp.lyear_str])(loni,lati,:)), 'omitnan');                         
                     end
                     [tmp.data_obs_det, tmp.data_obs_trend] = Func_0028_detrend_linear_1d(data.([tmp.varname, '_obs'])(loni,lati,:), 'omitnan');
                     [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data_det(isfinite(tmp.data_obs_det)), tmp.data_obs_det(isfinite(tmp.data_obs_det)));
    %                  [tmp.corr, tmp.corr_p]=corrcoef(squeeze(data.([tmp.varname, '_model', '_l', tmp.lyear_str])(loni,lati,1+lyear:end-cfg.proj_year+1)), ...
    %                      squeeze(data.([tmp.varname, '_obs'])(loni,lati,1+lyear:end)));

                     if sum(squeeze(isfinite(tmp.data_obs)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                        tmp.corr_det=NaN(2,2);
                        tmp.corr_det_p=NaN(2,2);
                        tmp.data_obs_trend=NaN;
                     end

                     data.([tmp.varname, '_corr_obs', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_obs_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);
                     data.([tmp.varname, '_corr_obs_det', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det(1,2);
                     data.([tmp.varname, '_corr_obs_det_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det_p(1,2);

                     data.([tmp.varname, '_obs_trend', '_l', tmp.lyear_str])(loni,lati)=tmp.data_obs_trend;

                 %% corr obs ~ assm
                    if isfield(tmp, 'ydata_mod_obs_masked')
                        tmp.data_assm = squeeze(data.([tmp.varname, '_assm_obs_masked', '_l', tmp.lyear_str])(loni,lati,:));
                     else
                        tmp.data_assm = squeeze(data.([tmp.varname, '_assm'])(loni,lati,:));
                     end
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data_assm(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
                     
                     if sum(squeeze(isfinite(tmp.data_obs)))<10
                         tmp.data_assm_det = NaN(size(tmp.data));
                         tmp.data_assm_trend = NaN(size(tmp.data));
                     else
                         [tmp.data_assm_det, tmp.data_assm_trend] = Func_0028_detrend_linear_1d(tmp.data_assm, 'omitnan');
                     end
                     [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data_assm_det(isfinite(tmp.data_obs_det)), tmp.data_obs_det(isfinite(tmp.data_obs_det)));
    %                  [tmp.corr, tmp.corr_p]=corrcoef(squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,1+lyear:end-cfg.proj_year+1)), ...
    %                      squeeze(data.([tmp.varname, '_obs'])(loni,lati,1+lyear:end)));

                     if sum(squeeze(isfinite(tmp.data_obs)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                        tmp.corr_det=NaN(2,2);
                        tmp.corr_det_p=NaN(2,2);
                     end

                     data.([tmp.varname, '_corr_obs_assm', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_obs_assm_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);
                     data.([tmp.varname, '_corr_obs_assm_det', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det(1,2);
                     data.([tmp.varname, '_corr_obs_assm_det_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det_p(1,2);

                 %% corr obs ~ lens2
                    if isfield(tmp, 'ydata_mod_obs_masked')
                        tmp.data_lens2 = squeeze(data.([tmp.varname, '_lens2_obs_masked', '_l', tmp.lyear_str])(loni,lati,:));
                     else
                        tmp.data_lens2 = squeeze(data.([tmp.varname, '_lens2', '_l', tmp.lyear_str])(loni,lati,:));
                     end
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data_lens2(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
                     
                     if sum(squeeze(isfinite(tmp.data_obs)))<10
                         tmp.data_lens2_det = NaN(size(tmp.data));
                         tmp.data_lens2_trend = NaN(size(tmp.data));
                     else
                         [tmp.data_lens2_det, tmp.data_lens2_trend] = Func_0028_detrend_linear_1d(tmp.data_lens2, 'omitnan');
                     end
                     [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data_lens2_det(isfinite(tmp.data_obs_det)), tmp.data_obs_det(isfinite(tmp.data_obs_det)));

                     if sum(squeeze(isfinite(tmp.data_obs)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                        tmp.corr_det=NaN(2,2);
                        tmp.corr_det_p=NaN(2,2);
                        tmp.data_lens2_trend=NaN;
                     end

                     data.([tmp.varname, '_corr_obs_lens2', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_obs_lens2_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);
                     data.([tmp.varname, '_corr_obs_lens2_det', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det(1,2);
                     data.([tmp.varname, '_corr_obs_lens2_det_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det_p(1,2);


                 %% corr obs ~ hcst-lens2
                    tmp.data_hcst_lens2 = squeeze(tmp.data) - squeeze(tmp.data_lens2);
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data_hcst_lens2(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
                     if sum(squeeze(isfinite(tmp.data_obs)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                     end
                     data.([tmp.varname, '_corr_obs_int', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_obs_int_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);
                 
                %% corr obs-lens2 ~ hcst-lens2
                    tmp.data_hcst_lens2 = squeeze(tmp.data) - squeeze(tmp.data_lens2);
                     tmp.data_obs_lens2 = squeeze(tmp.data_obs) - squeeze(tmp.data_lens2);
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data_hcst_lens2(isfinite(tmp.data_obs_lens2)), tmp.data_obs_lens2(isfinite(tmp.data_obs_lens2)));
                    
                     if sum(squeeze(isfinite(tmp.data_obs)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                     end
                     data.([tmp.varname, '_corr_obs_int2', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_obs_int2_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);

%                 %% corr AR1(obs) ~ obs
% %                     tmp.obs_AR1 = squeeze(data.([tmp.varname, '_obs_AR1', '_l', tmp.lyear_str])(loni,lati,:));
% %                      [tmp.corr, tmp.corr_p]=corrcoef(tmp.obs_AR1(isfinite(tmp.data_assm).*isfinite(tmp.obs_AR1)), ...
% %                          tmp.data_obs(isfinite(tmp.data_assm).*isfinite(tmp.obs_AR1)));
%                     tmp.obs_AR1 = squeeze(data2.([tmp.varname, '_obs_AR1', '_l', tmp.lyear_str])(loni,lati,:));
%                     tmp.obs_AR1_mask=squeeze(tmp.obs_AR1./tmp.obs_AR1);
%                     tmp.data_obs_mask=squeeze(tmp.data_obs./tmp.data_obs);
%                      [tmp.corr, tmp.corr_p]=corrcoef(tmp.obs_AR1.*tmp.data_obs_mask.*tmp.obs_AR1_mask, ...
%                          squeeze(tmp.data_obs).*tmp.data_obs_mask.*tmp.obs_AR1_mask, 'Rows', 'complete');
% 
%                      data.([tmp.varname, '_corr_obs_AR1', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
%                      data.([tmp.varname, '_corr_obs_AR1_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);
    
                 else
                    %% obs (NaN)
                     data.([tmp.varname, '_corr_obs', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_det', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_det_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_obs_trend', '_l', tmp.lyear_str])(loni,lati)=NaN;
    
                    %% AR1 (NaN)
                     data.([tmp.varname, '_corr_assm_AR1', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_AR1_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_AR1', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_AR1_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
    
                    %% ASSM (NaN)
                     data.([tmp.varname, '_corr_assm', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_det', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_det_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_int', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_int_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_int2', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_int2_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_assm_trend', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_model_trend', '_l', tmp.lyear_str])(loni,lati)=NaN;
    
                    %% LENS2 (NaN)
                     data.([tmp.varname, '_corr_assm_lens2', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_lens2_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_lens2_det', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_lens2_det_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_lens2_trend', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     %% obs ~ ASSM (NaN)
                     data.([tmp.varname, '_corr_obs_assm', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_assm_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_assm_det', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_assm_det_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     
                     %% obs ~ LENS2 (NaN)
                     data.([tmp.varname, '_corr_obs_lens2', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_lens2_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_lens2_det', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_lens2_det_p', '_l', tmp.lyear_str])(loni,lati)=NaN;

                     %% obs ~ HCST-LENS2 (NaN)
                     data.([tmp.varname, '_corr_obs_int', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_int_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_int2', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_int2_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                 end
            end
        end
        disp('correlation1 finished')
        fprintf('%7.1f sec\n', toc(lap_time) );

        %% get correlation coefficient2 (anomaly based)
        for loni=1:grid.nlon
            for lati=1:grid.nlat
%                  if (~isnan(data2.([tmp.varname, '_assm_ano', '_l', tmp.lyear_str])(loni,lati,1)) ...
%                          && sum(data2.([tmp.varname, '_model_ano_l', tmp.lyear_str])(loni,lati,:),'omitnan')~=0)
                if (~isnan(data2.([tmp.varname, '_assm_ano', '_l', tmp.lyear_str])(loni,lati,1)))

                 %% corr assm (anomaly based)
                     tmp.data2_assm = data2.([tmp.varname, '_assm_ano','_l',tmp.lyear_str])(loni,lati,:);
                     tmp.data2 = squeeze(data2.([tmp.varname, '_model_ano', '_l', tmp.lyear_str])(loni,lati,:));
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data2(isfinite(tmp.data2_assm)), tmp.data2_assm(isfinite(tmp.data2_assm)));
                     tmp.data2_det = data2.([tmp.varname, '_model_ano_det2', '_l', tmp.lyear_str])(loni,lati,:);
                     tmp.data2_assm_det = data2.([tmp.varname, '_assm_ano_det2', '_l', tmp.lyear_str])(loni,lati,:);
                     [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data2_det(isfinite(tmp.data2_assm_det)), tmp.data2_assm_det(isfinite(tmp.data2_assm_det)));
                                    
                     data2.([tmp.varname, '_corr_assm_ano', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
                     data2.([tmp.varname, '_corr_assm_ano_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);
                     data2.([tmp.varname, '_corr_assm_det_ano', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det(1,2);
                     data2.([tmp.varname, '_corr_assm_det_ano_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det_p(1,2);


                 %% hcst-lens2 ~ assm (anomaly based)
                     tmp.data2_hcst_lens2 = squeeze(data2.([tmp.varname, '_model_ano', '_l', tmp.lyear_str])(loni,lati,:)) - ...
                                    squeeze(data2.([tmp.varname, '_lens2_ano', '_l', tmp.lyear_str])(loni,lati,:));
                     [tmp.corr_int, tmp.corr_int_p]=corrcoef(tmp.data2_hcst_lens2(isfinite(tmp.data2_assm)), tmp.data2_assm(isfinite(tmp.data2_assm)));
    
                     
                     data2.([tmp.varname, '_corr_assm_int_ano', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_int(1,2);
                     data2.([tmp.varname, '_corr_assm_int_ano_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_int_p(1,2);
                 
                 %% hcst-lens2 ~ assm-lens2 (anomaly based)
                     tmp.data2_hcst_lens2 = squeeze(data2.([tmp.varname, '_model_ano', '_l', tmp.lyear_str])(loni,lati,:)) - ...
                                    squeeze(data2.([tmp.varname, '_lens2_ano', '_l', tmp.lyear_str])(loni,lati,:));
                     tmp.data2_assm_lens2 = squeeze(data2.([tmp.varname, '_assm_ano', '_l', tmp.lyear_str])(loni,lati,:)) - ...
                                    squeeze(data2.([tmp.varname, '_lens2_ano', '_l', tmp.lyear_str])(loni,lati,:));

                     [tmp.corr_int, tmp.corr_int_p]=corrcoef(tmp.data2_hcst_lens2(isfinite(tmp.data2_assm_lens2)), tmp.data2_assm_lens2(isfinite(tmp.data2_assm_lens2)));
    
                     
                     data2.([tmp.varname, '_corr_assm_int2_ano', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_int(1,2);
                     data2.([tmp.varname, '_corr_assm_int2_ano_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_int_p(1,2);
                
                 %% hcst-svd 1st mode ~ assm (anomaly based)
                    tmp.data2_hcst_svd = squeeze(data2.([tmp.varname, '_model_ano', '_l', tmp.lyear_str])(loni,lati,:)) - ...
                                    squeeze(data2.([tmp.varname, '_svd_lens2_model_right_1mode_recon', '_l', tmp.lyear_str])(loni,lati,:));
                    [tmp.corr_int2, tmp.corr_int2_p]=corrcoef(tmp.data2_hcst_svd(isfinite(tmp.data2_assm)), tmp.data2_assm(isfinite(tmp.data2_assm)));
                    data2.([tmp.varname, '_corr_assm_int_svd', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_int2(1,2);
                    data2.([tmp.varname, '_corr_assm_int_svd_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_int2_p(1,2);

                %% hcst-svd 1st mode ~ assm-svd 1st mode (anomaly based)
                    tmp.data2_hcst_svd = squeeze(data2.([tmp.varname, '_model_ano', '_l', tmp.lyear_str])(loni,lati,:)) - ...
                                    squeeze(data2.([tmp.varname, '_svd_lens2_model_right_1mode_recon', '_l', tmp.lyear_str])(loni,lati,:));
                    tmp.data2_assm_svd = squeeze(data2.([tmp.varname, '_assm_ano', '_l', tmp.lyear_str])(loni,lati,:)) - ...
                                    squeeze(data2.([tmp.varname, '_svd_lens2_assm_right_1mode_recon', '_l', tmp.lyear_str])(loni,lati,:));
                    [tmp.corr_int2, tmp.corr_int2_p]=corrcoef(tmp.data2_hcst_svd(isfinite(tmp.data2_assm_svd)), tmp.data2_assm_svd(isfinite(tmp.data2_assm_svd)));
                    data2.([tmp.varname, '_corr_assm_int2_svd', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_int2(1,2);
                    data2.([tmp.varname, '_corr_assm_int2_svd_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_int2_p(1,2);
                

                 %% corr lens2 ~ assm (anomaly based)
                     tmp.lens2 = squeeze(data2.([tmp.varname, '_lens2_ano', '_l', tmp.lyear_str])(loni,lati,:));
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.lens2(isfinite(tmp.data2_assm)), tmp.data2_assm(isfinite(tmp.data2_assm)));
                     tmp.data2_lens2_det = data2.([tmp.varname, '_lens2_ano_det2', '_l', tmp.lyear_str])(loni,lati,:);
                     [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data2_lens2_det(isfinite(tmp.data2_assm_det)), tmp.data2_assm_det(isfinite(tmp.data2_assm_det)));

                     data2.([tmp.varname, '_corr_assm_lens2_ano', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
                     data2.([tmp.varname, '_corr_assm_lens2_ano_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);
                     data2.([tmp.varname, '_corr_assm_lens2_det_ano', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det(1,2);
                     data2.([tmp.varname, '_corr_assm_lens2_det_ano_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det_p(1,2);

                 %% corr AR1(assm) ~ assm (anomaly based)
                
                    tmp.AR1 = squeeze(data2.([tmp.varname, '_AR1', '_l', tmp.lyear_str])(loni,lati,:));
                    tmp.AR1_mask=squeeze(tmp.AR1./tmp.AR1);
                    tmp.data2_assm_mask=squeeze(tmp.data2_assm./tmp.data2_assm);
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.AR1.*tmp.data2_assm_mask.*tmp.AR1_mask, ...
                         squeeze(tmp.data2_assm).*tmp.data2_assm_mask.*tmp.AR1_mask, 'Rows', 'complete');
                     data2.([tmp.varname, '_corr_assm_AR1_ano', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
                     data2.([tmp.varname, '_corr_assm_AR1_ano_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);
    
                
                 %% corr obs ~ hcst (anomaly based)
                     if isfield(tmp, 'ydata_mod_obs_masked')
                        tmp.data2 = squeeze(data.([tmp.varname, '_mod_obs_masked', '_l', tmp.lyear_str])(loni,lati,:)); %% not anomaly, caution!
                     else
                        tmp.data2 = squeeze(data2.([tmp.varname, '_model_ano', '_l', tmp.lyear_str])(loni,lati,:));
                     end
                     tmp.data2_obs = squeeze(data2.([tmp.varname, '_obs_ano', '_l', tmp.lyear_str])(loni,lati,:));
                     if sum(isfinite(tmp.data2_obs)) ==1
                         tmp.data2_obs=NaN(size(tmp.data2_obs));
                     end
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data2(isfinite(tmp.data2_obs)), tmp.data2_obs(isfinite(tmp.data2_obs)));
                     tmp.data2_obs_det = data2.([tmp.varname, '_obs_ano_det2', '_l', tmp.lyear_str])(loni,lati,:);
                     tmp.data2_det3 = data2.([tmp.varname, '_model_ano_det3', '_l', tmp.lyear_str])(loni,lati,:);
                     tmp.data2_assm_det3 = data2.([tmp.varname, '_assm_ano_det3', '_l', tmp.lyear_str])(loni,lati,:);
                     tmp.data2_lens2_det3 = data2.([tmp.varname, '_lens2_ano_det3', '_l', tmp.lyear_str])(loni,lati,:);
                     [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data2_det3(isfinite(tmp.data2_obs_det)), tmp.data2_obs_det(isfinite(tmp.data2_obs_det)));

                     if sum(squeeze(isfinite(tmp.data2_obs)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                        tmp.corr_det=NaN(2,2);
                        tmp.corr_det_p=NaN(2,2);
                        tmp.data2_obs_trend=NaN;
                     end

                     data2.([tmp.varname, '_corr_obs_ano', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
                     data2.([tmp.varname, '_corr_obs_ano_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);
                     data2.([tmp.varname, '_corr_obs_det_ano', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det(1,2);
                     data2.([tmp.varname, '_corr_obs_det_ano_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det_p(1,2);    

                 %% corr obs ~ assm (anomaly based)
                    if isfield(tmp, 'ydata_mod_obs_masked')
                        tmp.data2_assm = squeeze(data.([tmp.varname, '_assm_obs_masked', '_l', tmp.lyear_str])(loni,lati,:));
                     else
                        tmp.data2_assm = squeeze(data2.([tmp.varname, '_assm_ano', '_l', tmp.lyear_str])(loni,lati,:));
                     end
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data2_assm(isfinite(tmp.data2_obs)), tmp.data2_obs(isfinite(tmp.data2_obs)));
                     
                     [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data2_assm_det3(isfinite(tmp.data2_obs_det)), tmp.data2_obs_det(isfinite(tmp.data2_obs_det)));

                     if sum(squeeze(isfinite(tmp.data2_obs)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                        tmp.corr_det=NaN(2,2);
                        tmp.corr_det_p=NaN(2,2);
                     end

                     data2.([tmp.varname, '_corr_obs_assm_ano', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
                     data2.([tmp.varname, '_corr_obs_assm_ano_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);
                     data2.([tmp.varname, '_corr_obs_assm_det_ano', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det(1,2);
                     data2.([tmp.varname, '_corr_obs_assm_det_ano_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det_p(1,2);

                 %% corr obs ~ lens2 (anomaly based)
                    if isfield(tmp, 'ydata_mod_obs_masked')
                        tmp.data2_lens2 = squeeze(data.([tmp.varname, '_lens2_obs_masked', '_l', tmp.lyear_str])(loni,lati,:));
                     else
                        tmp.data2_lens2 = squeeze(data2.([tmp.varname, '_lens2_ano', '_l', tmp.lyear_str])(loni,lati,:));
                     end
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data2_lens2(isfinite(tmp.data2_obs)), tmp.data2_obs(isfinite(tmp.data2_obs)));
                     
                     [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data2_lens2_det3(isfinite(tmp.data2_obs_det)), tmp.data2_obs_det(isfinite(tmp.data2_obs_det)));

                     if sum(squeeze(isfinite(tmp.data2_obs)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                        tmp.corr_det=NaN(2,2);
                        tmp.corr_det_p=NaN(2,2);
                        tmp.data2_lens2_trend=NaN;
                     end

                     data2.([tmp.varname, '_corr_obs_lens2_ano', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
                     data2.([tmp.varname, '_corr_obs_lens2_ano_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);
                     data2.([tmp.varname, '_corr_obs_lens2_det_ano', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det(1,2);
                     data2.([tmp.varname, '_corr_obs_lens2_det_ano_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det_p(1,2);


                 %% corr obs ~ hcst-lens2 (anomaly based)
                    tmp.data2_hcst_lens2 = squeeze(tmp.data2) - squeeze(tmp.data2_lens2);
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data2_hcst_lens2(isfinite(tmp.data2_obs)), tmp.data2_obs(isfinite(tmp.data2_obs)));
                     if sum(squeeze(isfinite(tmp.data2_obs)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                     end
                     data2.([tmp.varname, '_corr_obs_int_ano', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
                     data2.([tmp.varname, '_corr_obs_int_ano_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);
                 
                 %% corr obs-lens2 ~ hcst-lens2 (anomaly based)
                    tmp.data2_hcst_lens2 = squeeze(tmp.data2) - squeeze(tmp.data2_lens2);
                    tmp.data2_obs_lens2 = squeeze(tmp.data2_obs) - squeeze(tmp.data2_lens2);
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data2_hcst_lens2(isfinite(tmp.data2_obs_lens2)), tmp.data2_obs_lens2(isfinite(tmp.data2_obs_lens2)));
                     if sum(squeeze(isfinite(tmp.data2_obs_lens2)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                     end
                     data2.([tmp.varname, '_corr_obs_int2_ano', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
                     data2.([tmp.varname, '_corr_obs_int2_ano_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);

                 %% hcst-svd 1st mode ~ obs (anomaly based)
                    tmp.data2_hcst_svd = squeeze(data2.([tmp.varname, '_model_ano', '_l', tmp.lyear_str])(loni,lati,:)) - ...
                                    squeeze(data2.([tmp.varname, '_svd_lens2_model_right_1mode_recon', '_l', tmp.lyear_str])(loni,lati,:));
                    [tmp.corr_int2, tmp.corr_int2_p]=corrcoef(tmp.data2_hcst_svd(isfinite(tmp.data2_obs)), tmp.data2_obs(isfinite(tmp.data2_obs)));
                    data2.([tmp.varname, '_corr_obs_int_svd', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_int2(1,2);
                    data2.([tmp.varname, '_corr_obs_int_svd_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_int2_p(1,2);

                %% hcst-svd 1st mode ~ obs-svd 1st mode (anomaly based)
                    tmp.data2_hcst_svd = squeeze(data2.([tmp.varname, '_model_ano', '_l', tmp.lyear_str])(loni,lati,:)) - ...
                                    squeeze(data2.([tmp.varname, '_svd_lens2_model_right_1mode_recon', '_l', tmp.lyear_str])(loni,lati,:));
                    tmp.data2_obs_svd = squeeze(data2.([tmp.varname, '_obs_ano', '_l', tmp.lyear_str])(loni,lati,:)) - ...
                                    squeeze(data2.([tmp.varname, '_svd_lens2_obs_right_1mode_recon', '_l', tmp.lyear_str])(loni,lati,:));
                    [tmp.corr_int2, tmp.corr_int2_p]=corrcoef(tmp.data2_hcst_svd(isfinite(tmp.data2_obs_svd)), tmp.data2_obs_svd(isfinite(tmp.data2_obs_svd)));
                    data2.([tmp.varname, '_corr_obs_int2_svd', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_int2(1,2);
                    data2.([tmp.varname, '_corr_obs_int2_svd_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_int2_p(1,2);

                %% corr AR1(obs) ~ obs
                    tmp.obs_AR1 = squeeze(data2.([tmp.varname, '_obs_AR1', '_l', tmp.lyear_str])(loni,lati,:));
                    tmp.obs_AR1_mask=squeeze(tmp.obs_AR1./tmp.obs_AR1);
                    tmp.data2_obs_mask=squeeze(tmp.data2_obs./tmp.data2_obs);
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.obs_AR1.*tmp.data2_obs_mask.*tmp.obs_AR1_mask, ...
                         squeeze(tmp.data2_obs).*tmp.data2_obs_mask.*tmp.obs_AR1_mask, 'Rows', 'complete');

                     data2.([tmp.varname, '_corr_obs_AR1', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
                     data2.([tmp.varname, '_corr_obs_AR1_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);
    
                 else
                    %% obs (NaN)
                     data2.([tmp.varname, '_corr_obs_ano', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_ano_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_det_ano', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_det_ano_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
    
                    %% AR1 (NaN)
                     data2.([tmp.varname, '_corr_assm_AR1', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_AR1_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_AR1', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_AR1_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_AR1_ano', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_AR1_ano_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
    
    
                    %% ASSM (NaN)
                     data2.([tmp.varname, '_corr_assm_ano', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_ano_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_det_ano', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_det_ano_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_int_ano', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_int_ano_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_int2_ano', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_int2_ano_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_int_svd', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_int_svd_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_int2_svd', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_int2_svd_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
    
                    %% LENS2 (NaN)
                     data2.([tmp.varname, '_corr_assm_lens2_ano', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_lens2_ano_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_lens2_det_ano', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_lens2_det_ano_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     %% obs ~ ASSM (NaN)
                     data2.([tmp.varname, '_corr_obs_assm_ano', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_assm_ano_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_assm_det_ano', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_assm_det_ano_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     
                     %% obs ~ LENS2 (NaN)
                     data2.([tmp.varname, '_corr_obs_lens2_ano', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_lens2_ano_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_lens2_det_ano', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_lens2_det_ano_p', '_l', tmp.lyear_str])(loni,lati)=NaN;

                     %% obs ~ HCST-LENS2 (NaN)
                     data2.([tmp.varname, '_corr_obs_int_ano', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_int_ano_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_int2_ano', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_int2_ano_p', '_l', tmp.lyear_str])(loni,lati)=NaN;

                     %% obs ~ HCST-svd (NaN)
                     data2.([tmp.varname, '_corr_obs_int_svd', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_int_svd_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_int2_svd', '_l', tmp.lyear_str])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_int2_svd_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                 end
            end
        end
        disp('correlation2(anomaly based) finished')
        fprintf('%7.1f sec\n', toc(lap_time) );

        if isfield(tmp,'ymean_mod_obs_masked')
            tmp=rmfield(tmp, 'ymean_mod_obs_masked');
        end
        
        if ~exist(dirs.matroot,'dir'), mkdir(dirs.matroot); end
        fig_cfg.mat_name=[dirs.matroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l', tmp.lyear_str, 'y.mat'];
        save(fig_cfg.mat_name, 'data', 'data2')
        clear data data2
        fprintf('%7.1f sec\n', toc(lap_time) );

    end

    


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
            obsname_simple='ersst_reg_cesm2.v5.';
        case 'PRECT'
            obsname_simple='GPCC_reg_cesm2.v5.';
        case 'RAIN'
            obsname_simple='GPCC_reg_cesm2.v5.';
        case 'PSL'
            obsname_simple='ERA5_msl_reg_cesm2.';
        case 'SOILWATER_10CM'
%             obsname_simple='SM_reg_cesm2.';
            obsname_simple='GLEAM_reg_cesm2.v5.';
        case 'TWS'
            obsname_simple='TSW_reg_cesm2.';
        case 'SSH'
            obsname_simple='CMEMS_reg_cesm2.';
        case 'TS'
%             obsname_simple='HadCRUT5_reg_cesm2.';
            obsname_simple='ERA5_t2m_reg_cesm2.';
        case 'sumChl'
            obsname_simple='OC_CCI_reg_cesm2.';
        case 'TLAI'
            obsname_simple='LAI_reg_cesm2.'; % NOAA LAI
        case 'FAREA_BURNED'
            obsname_simple='burned_area_reg_cesm2.'; % MODIS Fire_cci v5.1
            obsname_simple='AVHRR-LTDR_reg_cesm2.'; % AVHRR-LTDR
        case 'COL_FIRE_CLOSS'
            obsname_simple='FIRE_CLOSS_reg_cesm2.'; % GFED
        case 'photoC_TOT_zint'
%             obsname_simple='s_vgpm_reg_cesm2.'; % VGPM
%             obsname_simple='ensmean_reg_cesm2.'; % VGPM   
            obsname_simple='CMEMS_reg_cesm2.'; %Globcolour;
        case 'photoC_TOT_zint_100m'
%             obsname_simple='s_vgpm_reg_cesm2.'; % VGPM
%             obsname_simple='ensmean_reg_cesm2.'; % VGPM
            obsname_simple='CMEMS_reg_cesm2.'; %Globcolour;
        case 'GPP'
%             obsname_simple='ORNL_DAAC_reg_cesm2.'; % ORNL_DAAC
            obsname_simple='VODCA2GPP_reg_cesm2.'; % VODCA2GPP
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

