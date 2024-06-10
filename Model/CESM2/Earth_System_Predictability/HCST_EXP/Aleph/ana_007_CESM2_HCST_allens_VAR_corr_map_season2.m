% %  Created 22-Sep-2023 by Yong-Yub Kim
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

%% model configuration

cfg.vars = {'TS', 'PSL', 'PRECT', 'SST'};
% cfg.vars = {'PRECT', 'PSL'};
% cfg.vars = {'SST'};
% cfg.vars = {'NPP', 'GPP', 'RAIN', 'TOTVEGC', 'FIRE'};
% cfg.vars = {'TWS'};
% cfg.vars = {'SOILWATER_10CM'};
cfg.vars = {'TS', 'PRECT', 'SST', 'PSL'};
cfg.vars = {'TS', 'PRECT', 'SST', 'PSL'};
cfg.vars={'TWS'};
cfg.vars = {'TWS', 'SOILWATER_10CM', 'TLAI', 'FAREA_BURNED', 'COL_FIRE_CLOSS'};
cfg.vars = {'TS'};
% cfg.vars= {'photoC_TOT_zint', 'photoC_TOT_zint_100m', 'SSH'};
% tmp.dimids= [1, 2, 4];
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
    cfg.obs_iyears=1970:2020;
%     cfg.obs_iyears=1960:2020;

    cfg.obs_iyears2=f_obs_iyears(cfg.var);
    
    % dirs.hcstroot=['/mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive/', cfg.comp, '/', cfg.var];
    % dirs.obsroot=['/Volumes/kyy_raid/kimyy/Observation/', cfg.obs_name, '/monthly_reg_', cfg.obs_fname_module(2:4)];
    % dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', cfg.var];
    % dirs.lens2root=['/mnt/lustre/proj/earth.system.predictability/LENS2/archive_analysis/', cfg.comp, '/', cfg.var];
    % dirs.assmroot=['/mnt/lustre/proj/earth.system.predictability/ASSM_EXP/archive_analysis/', cfg.comp, '/', cfg.var];
    
    dirs.hcstroot=['/mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive_seasonal_transfer/', cfg.comp, '/', cfg.var];
    dirs.obsroot=['/mnt/lustre/proj/kimyy/Observation/', cfg.obs_name];
    dirs.matroot=['/mnt/lustre/proj/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp,'/', cfg.var];
    dirs.lens2root=['/mnt/lustre/proj/earth.system.predictability/LENS2/archive_seasonal_transfer/', cfg.comp, '/', cfg.var];
    dirs.assmroot=['/mnt/lustre/proj/earth.system.predictability/ASSM_EXP/archive_seasonal_transfer/', cfg.comp, '/', cfg.var];    
    
    cfg.iyears=cfg.obs_iyears;
    cfg.gnm='f09_g17';
%     cfg.proj_year=1;
    cfg.proj_year=5;
%     cfg.season = {'JFM1', 'AMJ1', 'JAS1', 'OND1'};
    cfg.season = {'MAM1', 'JJA1', 'SON1', 'DJF1', 'MAM2', 'JJA2', 'SON2'};

    cfg.max_iy=max(cfg.iyears);
    cfg.min_iy=min(cfg.iyears);
    cfg.len_t_y = length(cfg.iyears);
    
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
    %% lead season loop
    for lss=1:length(cfg.season)
        tmp.season=cfg.season{lss};
%         tmp.mons = f_season_mons(tmp.season);
        cfg.casename_m=['all'];

    
        dirs.datadir= [dirs.hcstroot, filesep, 'ens_all', filesep];
        dirs.assmdir= [dirs.assmroot, filesep, 'ens_all', filesep];
        dirs.lens2dir= [dirs.lens2root, filesep, 'ens_smbb', filesep];
        fprintf('%02d_%s_%s  ',tmp.season, ',', cfg.casename_m,'_', tmp.varname); lap_time = tic;
        
%% variables initialization
        data.([tmp.varname, '_model', '_', tmp.season])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data.([tmp.varname, '_lens2', '_', tmp.season])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data.([tmp.varname, '_bias', '_', tmp.season])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data.([tmp.varname, '_obs'])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data.([tmp.varname, '_assm'])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        

    %     data.([tmp.varname, '_assm'])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        
%% read variables
        for iyear=min(cfg.iyears):max(cfg.iyears)
            tmp.iyear_str=num2str(iyear, '%04i');
            cfg.casename=['ensmean_', cfg.casename_m, '_i', tmp.iyear_str];
           

%% yearly filename    
%% HCST mean
            cfg.mod_fnm=[dirs.datadir, tmp.fs, 'ens_all_i',tmp.iyear_str, tmp.fs, tmp.season, tmp.fs, ...
            tmp.varname, '_', cfg.gnm, '.hcst.', cfg.casename, cfg.obs_fname_module, tmp.iyear_str, '_', tmp.season, '.nc'];
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
            cfg.mod_fnm=[dirs.datadir, tmp.fs, 'ens_all_i',tmp.iyear_str, tmp.fs, tmp.season, tmp.fs, ...
                tmp.varname, '_', cfg.gnm, '.hcst.', cfg.casename_std, cfg.obs_fname_module,tmp.iyear_str,'_',tmp.season '.nc'];
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

%% LENS2 mean
%             if iyear <= cfg.max_iy
                cfg.lens2_fnm=[dirs.lens2dir, tmp.fs, tmp.season, tmp.fs, ...
                    tmp.varname, '_', tmp.season, '_ensmean_', tmp.iyear_str, '.nc'];

                if exist(cfg.lens2_fnm)~=0
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
%             else
%                 tmp.ymean_lens2(1:grid.nlon,1:grid.nlat) = NaN;                    
%             end

%% LENS2 std
%             if iyear <= cfg.max_iy
                cfg.lens2_fnm=[dirs.lens2dir, tmp.fs, tmp.season, tmp.fs, ...
                    tmp.varname, '_', tmp.season, '_ensstd_', tmp.iyear_str, '.nc'];
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
                    tmp.ystde_lens2(1:grid.nlon,1:grid.nlat) = NaN;  
                end
%             else
%                 tmp.ystde_lens2(1:grid.nlon,1:grid.nlat) = NaN;                    
%             end

%% OBS
               tmp.mons = f_season_mons(tmp.season);
                for mi=1:length(tmp.mons)
                    tmp.mon=tmp.mons(mi);
                   if tmp.mon>=12 && mod(tmp.mon,12)~=0
                        tmp.mon_str=num2str(tmp.mon-12*floor(tmp.mon/12), '%02i');
                        tmp.fy = iyear+floor(tmp.mon/12);
                        tmp.fy_str=num2str(iyear+floor(tmp.mon/12), '%04i');
                    elseif tmp.mon>=12 && mod(tmp.mon,12) ==0
                        tmp.mon_str='12';
                        tmp.fy=iyear+floor(tmp.mon/12)-1;
                        tmp.fy_str=num2str(iyear+floor(tmp.mon/12)-1, '%04i');
                    else
                        tmp.mon_str=num2str(tmp.mon, '%02i');
                        tmp.fy=iyear;
                        tmp.fy_str=num2str(iyear, '%04i');
                   end
                    tmp.m_str=tmp.mon_str;
                    if strcmp(cfg.obs_name, 'ERA5')==1
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, cfg.obs_varname, '/monthly_reg_', cfg.obs_fname_module(2:4),tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
                    elseif (strcmp(cfg.obs_name, 'CMEMS')==1 && strcmp(tmp.varname, 'SOILWATER_10CM')==1)
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'SOIL_MOISTURE/COMBINED', '/monthly_reg_cam', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
                    elseif (strcmp(cfg.obs_name, 'NOAA')==1 && strcmp(tmp.varname, 'TWS')==1)
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'TSW', '/monthly_reg_cam', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
                    elseif strcmp(cfg.obs_name, 'GPCC')==1
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'monthly_reg_cam' ,tmp.fs,cfg.obs_fname_mid,tmp.fy_str,tmp.m_str,'.nc'];
                    elseif (strcmp(cfg.obs_name, 'GLEAM')==1 && strcmp(tmp.varname, 'SOILWATER_10CM')==1)
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, cfg.obs_varname, tmp.fs, 'monthly_reg_cam', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
                    elseif (strcmp(cfg.obs_name, 'NOAA')==1 && strcmp(tmp.varname, 'TLAI')==1)
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, cfg.obs_varname, tmp.fs, 'monthly_reg_cam', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
                    elseif (strcmp(cfg.obs_name, 'MODIS')==1 && strcmp(tmp.varname, 'FAREA_BURNED')==1)
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, cfg.obs_varname, tmp.fs, 'monthly_reg_cam', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
                    elseif (strcmp(cfg.obs_name, 'GFED')==1 && strcmp(tmp.varname, 'COL_FIRE_CLOSS')==1)
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'FIRE_CLOSS', tmp.fs, 'monthly_reg_cam', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
                    elseif (strcmp(cfg.obs_name, 'VGPM')==1)
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'netcdf_regrid/s_vgpm', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
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
                                tmp.dd=tmp.dd./1000.0./86400/eomday(iyear,mi);
                            elseif strcmp(tmp.varname, 'RAIN')
                                tmp.dd=tmp.dd./86400/eomday(iyear,mi);
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
                                tmp.dd=tmp.dd./86400./eomday(tmp.fy,mi);
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
  
%             tmp.ydata_obs(:,:,mon)=ncread(cfg.obs_fnm, cfg.obs_varname);

%% ASSM
%             if iyear <= cfg.max_iy
                cfg.assm_fnm=[dirs.assmdir, tmp.fs, tmp.season, tmp.fs, tmp.varname, '_ensmean_', tmp.iyear_str, '_', tmp.season, '.nc'];
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
                data.units=netcdf.getAtt(ncid,tmpvarid,'units');
                netcdf.close(ncid);
%             else
%                 tmp.ymean_assm(1:grid.nlon,1:grid.nlat) = NaN;                    
%             end

%% ASSM std
%             if iyear <= cfg.max_iy
                cfg.assm_fnm=[dirs.assmdir, tmp.fs, tmp.season, tmp.fs, tmp.varname, '_ensstd_', tmp.iyear_str, '_', tmp.season, '.nc'];
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
%             else
%                 tmp.ystde_assm(1:grid.nlon,1:grid.nlat) = NaN;                    
%             end

%% observation mask (scarced data like marine npp)            
%             if (strcmp(cfg.var, 'sumChl') || strcmp(cfg.var, 'TS')) %for Chls, it uses all available data
            if (strcmp(cfg.var, 'sumChl')) %for Chls, it uses all available data
%             if (strcmp(cfg.var, 'sumChl') || strcmp(cfg.var, 'SOILWATER_10CM')) %for Chls, it uses all available data
                tmp.obs_mask(:,:)=tmp.ymean_obs(:,:)./tmp.ymean_obs(:,:);
                tmp.ymean_mod_obs_masked(:,:) = tmp.ymean(:,:) .* tmp.obs_mask(:,:); % masking with available observation
                tmp.ymean_assm_obs_masked(:,:) = tmp.ymean_assm(:,:) .* tmp.obs_mask(:,:); % masking with available observation
                tmp.ymean_lens2_obs_masked(:,:) = tmp.ymean_lens2(:,:) .* tmp.obs_mask(:,:); % masking with available observation
                tmp.ystde_mod_obs_masked(:,:) = tmp.ystde(:,:) .* tmp.obs_mask(:,:); % masking with available observation
                tmp.ystde_assm_obs_masked(:,:) = tmp.ystde_assm(:,:) .* tmp.obs_mask(:,:); % masking with available observation
                tmp.ystde_lens2_obs_masked(:,:) = tmp.ystde_lens2(:,:) .* tmp.obs_mask(:,:); % masking with available observation
            end
    
           
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
                    tmp.ystde_obs = tmp.ystde_obs .* grid.mask_ocn;
                    tmp.ystde_assm = tmp.ystde_assm .* grid.mask_ocn;
                    tmp.ystde_lens2 = tmp.ystde_lens2 .* grid.mask_ocn;    
            end
    %         end

%% save variables as structure
            data.([tmp.varname, '_model', '_', tmp.season])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean;
            data.([tmp.varname, '_model_stde', '_', tmp.season])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ystde;
            if isfield(tmp, 'ymean_mod_obs_masked')
                data.([tmp.varname, '_mod_obs_masked', '_', tmp.season])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_mod_obs_masked;
                data.([tmp.varname, '_assm_obs_masked', '_', tmp.season])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_assm_obs_masked;
                data.([tmp.varname, '_lens2_obs_masked', '_', tmp.season])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_lens2_obs_masked;
                data.([tmp.varname, '_mod_stde_obs_masked', '_', tmp.season])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ystde_mod_obs_masked;
                data.([tmp.varname, '_assm_stde_obs_masked', '_', tmp.season])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ystde_assm_obs_masked;
                data.([tmp.varname, '_lens2_stde_obs_masked', '_', tmp.season])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ystde_lens2_obs_masked;
            end
            data.([tmp.varname, '_obs'])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_obs;
            data.([tmp.varname, '_lens2', '_', tmp.season])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_lens2;
            data.([tmp.varname, '_assm'])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_assm;
            data.([tmp.varname, '_lens2_stde', '_', tmp.season])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ystde_lens2;
            data.([tmp.varname, '_assm_stde'])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ystde_assm;
            
            
    
        end  %% initialized year loop end
        if (strcmp(cfg.var, 'SST' ))
            data.([tmp.varname, '_model', '_', tmp.season]) = data.([tmp.varname, '_model', '_', tmp.season]) - 273.15;
            data.([tmp.varname, '_lens2', '_', tmp.season]) = data.([tmp.varname, '_lens2', '_', tmp.season]) - 273.15;
            data.([tmp.varname, '_assm']) = data.([tmp.varname, '_assm']) - 273.15;
            data.([tmp.varname, '_model', '_', tmp.season])(data.([tmp.varname, '_model', '_', tmp.season])<-30) = NaN;
            data.([tmp.varname, '_lens2', '_', tmp.season])(data.([tmp.varname, '_lens2', '_', tmp.season])<-30) = NaN;
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
        data.([tmp.varname, '_model_stdt_', tmp.season])= std(data.([tmp.varname, '_model_', tmp.season]),0,3,'omitnan');
        data.([tmp.varname, '_lens2_stdt_', tmp.season])= std(data.([tmp.varname, '_lens2_', tmp.season]),0,3,'omitnan');

   %% get model error(bias) with ASSM 
        data.([tmp.varname, '_model_err', '_', tmp.season]) = ...
            data.([tmp.varname, '_model', '_', tmp.season]) - data.([tmp.varname, '_assm']);
        data.([tmp.varname, '_lens2_err', '_', tmp.season]) = ...
            data.([tmp.varname, '_lens2', '_', tmp.season]) - data.([tmp.varname, '_assm']);

   %% get model error(bias) with obs 
        data.([tmp.varname, '_model_err_obs', '_', tmp.season]) = ...
            data.([tmp.varname, '_model', '_', tmp.season]) - data.([tmp.varname, '_obs']);
        data.([tmp.varname, '_lens2_err_obs', '_', tmp.season]) = ...
            data.([tmp.varname, '_lens2', '_', tmp.season]) - data.([tmp.varname, '_obs']);
        data.([tmp.varname, '_assm_err_obs', '_', tmp.season]) = ...
            data.([tmp.varname, '_assm']) - data.([tmp.varname, '_obs']);

   %% get RMSE with ASSM (drift-removed)
%         data.([tmp.varname, '_model_rmse', '_', tmp.season]) = ...
%             sqrt(mean(data.([tmp.varname, '_model_err', '_', tmp.season]).^2, 3, 'omitnan'));
%         data.([tmp.varname, '_lens2_rmse', '_', tmp.season]) = ...
%             sqrt(mean(data.([tmp.varname, '_lens2_err', '_', tmp.season]).^2, 3, 'omitnan'));
        data.([tmp.varname, '_model_rmse', '_', tmp.season]) = ...
            sqrt( mean( ...
            (data.([tmp.varname, '_model_err', '_', tmp.season]) ...
            - mean(data.([tmp.varname, '_model_err', '_', tmp.season]), 3, 'omitnan')).^2, ...
            3, 'omitnan') );
        data.([tmp.varname, '_lens2_rmse', '_', tmp.season]) = ...
            sqrt( mean( ...
            (data.([tmp.varname, '_lens2_err', '_', tmp.season]) ...
            - mean(data.([tmp.varname, '_lens2_err', '_', tmp.season]), 3, 'omitnan')).^2, ...
            3, 'omitnan') );
   %% get RMSE with obs (drift-removed)
        data.([tmp.varname, '_model_rmse_obs', '_', tmp.season]) = ...
            sqrt( mean( ...
            (data.([tmp.varname, '_model_err_obs', '_', tmp.season]) ...
            - mean(data.([tmp.varname, '_model_err_obs', '_', tmp.season]), 3, 'omitnan')).^2, ...
            3, 'omitnan') );
        data.([tmp.varname, '_lens2_rmse_obs', '_', tmp.season]) = ...
            sqrt( mean( ...
            (data.([tmp.varname, '_lens2_err_obs', '_', tmp.season]) ...
            - mean(data.([tmp.varname, '_lens2_err_obs', '_', tmp.season]), 3, 'omitnan')).^2, ...
            3, 'omitnan') );
        data.([tmp.varname, '_assm_rmse_obs', '_', tmp.season]) = ...
            sqrt( mean( ...
            (data.([tmp.varname, '_assm_err_obs', '_', tmp.season]) ...
            - mean(data.([tmp.varname, '_assm_err_obs', '_', tmp.season]), 3, 'omitnan')).^2, ...
            3, 'omitnan') );

    %% get nRMSE with ASSM (drift-removed)
        data.([tmp.varname, '_model_nrmse', '_', tmp.season])= ...
            data.([tmp.varname, '_model_rmse', '_', tmp.season])./data.([tmp.varname, '_assm_stdt']);
        data.([tmp.varname, '_lens2_nrmse', '_', tmp.season])= ...
            data.([tmp.varname, '_lens2_rmse', '_', tmp.season])./data.([tmp.varname, '_assm_stdt']);

    %% get nRMSE with obs (drift-removed)
        data.([tmp.varname, '_model_nrmse_obs', '_', tmp.season])= ...
            data.([tmp.varname, '_model_rmse_obs', '_', tmp.season])./data.([tmp.varname, '_obs_stdt']);
        data.([tmp.varname, '_lens2_nrmse_obs', '_', tmp.season])= ...
            data.([tmp.varname, '_lens2_rmse_obs', '_', tmp.season])./data.([tmp.varname, '_obs_stdt']);
        data.([tmp.varname, '_assm_nrmse_obs', '_', tmp.season])= ...
            data.([tmp.varname, '_assm_rmse_obs', '_', tmp.season])./data.([tmp.varname, '_obs_stdt']);
    
    
    disp('RMSEs calculated')
    fprintf('%7.1f sec\n', toc(lap_time) );

    %% get hindcast climatology as the function of lead year
       switch cfg.obs_name
           case 'GPCC'
               cfg.clim_ys=1961;
               cfg.clim_ye=2019;
           otherwise
               cfg.clim_ys=1961;
               cfg.clim_ye=2020;
       end
       cfg.clim_tlen = (cfg.clim_ys-1959):(cfg.clim_ye-2020)+cfg.len_t_y;
       cfg.clim_tlen2=length(cfg.clim_tlen);
       
       cfg.clim_tlen_AR1 = 1:cfg.len_t_y; % for AR1 integration, it should start from 1960 to integrate 1965, LY5 
       

        data.([tmp.varname, '_model_clim', '_', tmp.season]) = mean(data.([tmp.varname, '_model', '_', tmp.season])(:,:,cfg.clim_tlen),3, 'omitnan'); %1964 ~ 2020
        data.([tmp.varname, '_assm_clim']) = mean(data.([tmp.varname, '_assm'])(:,:,cfg.clim_tlen),3, 'omitnan'); %1964 ~ 2020
        data.([tmp.varname, '_lens2_clim']) = mean(data.([tmp.varname, '_lens2', '_', tmp.season])(:,:,cfg.clim_tlen),3, 'omitnan'); %1964 ~ 2020
        data.([tmp.varname, '_obs_clim']) = mean(data.([tmp.varname, '_obs'])(:,:,cfg.clim_tlen),3, 'omitnan'); %1964 ~ 2020
        
        tmp.assm_clim_AR1_period = mean(data.([tmp.varname, '_assm'])(:,:,cfg.clim_tlen_AR1),3, 'omitnan');
        tmp.obs_clim_AR1_period = mean(data.([tmp.varname, '_obs'])(:,:,cfg.clim_tlen_AR1),3, 'omitnan');

   %% get anomaly values
        data2.([tmp.varname, '_model_ano', '_', tmp.season]) = ...
            data.([tmp.varname, '_model', '_', tmp.season])(:,:,cfg.clim_tlen) - data.([tmp.varname, '_model_clim', '_', tmp.season]);
        data2.([tmp.varname, '_assm_ano', '_', tmp.season]) = ...
            data.([tmp.varname, '_assm'])(:,:,cfg.clim_tlen) - data.([tmp.varname, '_assm_clim']);
        data2.([tmp.varname, '_lens2_ano', '_', tmp.season]) = ...
            data.([tmp.varname, '_lens2', '_', tmp.season])(:,:,cfg.clim_tlen) - data.([tmp.varname, '_lens2_clim']);
        data2.([tmp.varname, '_obs_ano', '_', tmp.season]) = ...
            data.([tmp.varname, '_obs'])(:,:,cfg.clim_tlen) - data.([tmp.varname, '_obs_clim']);
        
        data2_l0_AR1.assm_ano_full = data.([tmp.varname, '_assm'])(:,:,cfg.clim_tlen_AR1) - tmp.assm_clim_AR1_period;
        data2_l0_AR1.obs_ano_full = data.([tmp.varname, '_obs'])(:,:,cfg.clim_tlen_AR1) - tmp.obs_clim_AR1_period;

    %% get anomaly_detrended values
        data2.([tmp.varname, '_model_ano_det2', '_', tmp.season]) = NaN(grid.nlon, grid.nlat, cfg.clim_tlen2);
        data2.([tmp.varname, '_assm_ano_det2', '_', tmp.season]) = NaN(grid.nlon, grid.nlat, cfg.clim_tlen2);
        data2.([tmp.varname, '_lens2_ano_det2', '_', tmp.season]) = NaN(grid.nlon, grid.nlat, cfg.clim_tlen2);
        data2.([tmp.varname, '_obs_ano_det2', '_', tmp.season]) = NaN(grid.nlon, grid.nlat, cfg.clim_tlen2);

        data2.([tmp.varname, '_model_trend2', '_', tmp.season]) = NaN(grid.nlon, grid.nlat);
        data2.([tmp.varname, '_assm_trend2', '_', tmp.season]) = NaN(grid.nlon, grid.nlat);
        data2.([tmp.varname, '_lens2_trend2', '_', tmp.season]) = NaN(grid.nlon, grid.nlat);
        data2.([tmp.varname, '_obs_trend2', '_', tmp.season]) = NaN(grid.nlon, grid.nlat);

        for loni=1:grid.nlon
            for lati=1:grid.nlat
                 if (isnan(data.([tmp.varname, '_assm'])(loni,lati,1))~=1 & nansum(data.([tmp.varname, '_model_', tmp.season])(loni,lati,:))~=0)
                    tmp.data_obs = squeeze(data2.([tmp.varname, '_obs_ano', '_', tmp.season])(loni,lati,:));
                    tmp.data_assm = squeeze(data2.([tmp.varname, '_assm_ano', '_', tmp.season])(loni,lati,:));
                    tmp.data_lens2 = squeeze(data2.([tmp.varname, '_lens2_ano', '_', tmp.season])(loni,lati,:));
                    tmp.data_model = squeeze(data2.([tmp.varname, '_model_ano', '_', tmp.season])(loni,lati,:));
                    
                    [data2.([tmp.varname, '_obs_ano_det2', '_', tmp.season])(loni,lati,:), ...
                        data2.([tmp.varname, '_obs_trend2', '_', tmp.season])(loni,lati)] ...
                        = Func_0028_detrend_linear_1d(tmp.data_obs, 'omitnan');
                    [data2.([tmp.varname, '_assm_ano_det2', '_', tmp.season])(loni,lati,:), ...
                        data2.([tmp.varname, '_assm_trend2', '_', tmp.season])(loni,lati)] ...
                        = Func_0028_detrend_linear_1d(tmp.data_assm, 'omitnan');
                    [data2.([tmp.varname, '_lens2_ano_det2', '_', tmp.season])(loni,lati,:), ...
                        data2.([tmp.varname, '_lens2_trend2', '_', tmp.season])(loni,lati)] ...
                        = Func_0028_detrend_linear_1d(tmp.data_lens2, 'omitnan');
                    [data2.([tmp.varname, '_model_ano_det2', '_', tmp.season])(loni,lati,:), ...
                        data2.([tmp.varname, '_model_trend2', '_', tmp.season])(loni,lati)] ...
                        = Func_0028_detrend_linear_1d(tmp.data_model, 'omitnan');
                 end
            end
        end
        data2.([tmp.varname, '_obs_ano_det2', '_', tmp.season]) = ...
            data2.([tmp.varname, '_obs_ano_det2', '_', tmp.season]) - ...
            mean(data2.([tmp.varname, '_obs_ano_det2', '_', tmp.season]), 3, 'omitnan');
        data2.([tmp.varname, '_assm_ano_det2', '_', tmp.season]) = ...
            data2.([tmp.varname, '_assm_ano_det2', '_', tmp.season]) - ...
            mean(data2.([tmp.varname, '_assm_ano_det2', '_', tmp.season]), 3, 'omitnan');
        data2.([tmp.varname, '_lens2_ano_det2', '_', tmp.season]) = ...
            data2.([tmp.varname, '_lens2_ano_det2', '_', tmp.season]) - ...
            mean(data2.([tmp.varname, '_lens2_ano_det2', '_', tmp.season]), 3, 'omitnan');
        data2.([tmp.varname, '_model_ano_det2', '_', tmp.season]) = ...
            data2.([tmp.varname, '_model_ano_det2', '_', tmp.season]) - ...
            mean(data2.([tmp.varname, '_model_ano_det2', '_', tmp.season]), 3, 'omitnan');

    %% get temporal std of ensmean
        data2.([tmp.varname, '_assm_ano_det2_stdt'])= std(data2.([tmp.varname, '_assm_ano_det2', '_', tmp.season]),0,3,'omitnan');
        data2.([tmp.varname, '_obs_ano_det2_stdt'])= std(data2.([tmp.varname, '_obs_ano_det2', '_', tmp.season]),0,3,'omitnan');
        data2.([tmp.varname, '_model_ano_det2_stdt_', tmp.season])= std(data2.([tmp.varname, '_model_ano_det2', '_', tmp.season]),0,3,'omitnan');
        data2.([tmp.varname, '_lens2_ano_det2_stdt_', tmp.season])= std(data2.([tmp.varname, '_lens2_ano_det2', '_', tmp.season]),0,3,'omitnan');

   %% get model error(bias) with ASSM (between anomalies)
        data2.([tmp.varname, '_model_err2', '_', tmp.season]) = ...
            data2.([tmp.varname, '_model_ano_det2', '_', tmp.season]) - data2.([tmp.varname, '_assm_ano_det2', '_', tmp.season]);
        data2.([tmp.varname, '_lens2_err2', '_', tmp.season]) = ...
            data2.([tmp.varname, '_lens2_ano_det2', '_', tmp.season]) - data2.([tmp.varname, '_assm_ano_det2', '_', tmp.season]);
   %% get model error(bias) with obs (between anomalies)
        data2.([tmp.varname, '_model_err2_obs', '_', tmp.season]) = ...
            data2.([tmp.varname, '_model_ano_det2', '_', tmp.season]) - data2.([tmp.varname, '_obs_ano_det2', '_', tmp.season]);
        data2.([tmp.varname, '_lens2_err2_obs', '_', tmp.season]) = ...
            data2.([tmp.varname, '_lens2_ano_det2', '_', tmp.season]) - data2.([tmp.varname, '_obs_ano_det2', '_', tmp.season]);
        data2.([tmp.varname, '_assm_err2_obs', '_', tmp.season]) = ...
            data2.([tmp.varname, '_assm_ano_det2', '_', tmp.season]) - data2.([tmp.varname, '_obs_ano_det2', '_', tmp.season]);

   %% get RMSE with ASSM based on anomaly (climatology removed, detrended, drift-removed)
        data2.([tmp.varname, '_model_rmse2', '_', tmp.season]) = ...
            sqrt( mean(data2.([tmp.varname, '_model_err2', '_', tmp.season]).^2, 3, 'omitnan') );
        data2.([tmp.varname, '_lens2_rmse2', '_', tmp.season]) = ...
            sqrt( mean(data2.([tmp.varname, '_lens2_err2', '_', tmp.season]).^2, 3, 'omitnan') );
   %% get RMSE with obs based on anomaly (climatology removed, detrended, drift-removed)
        data2.([tmp.varname, '_assm_rmse2_obs', '_', tmp.season]) = ...
            sqrt( mean(data2.([tmp.varname, '_assm_err2_obs', '_', tmp.season]).^2, 3, 'omitnan') );
        data2.([tmp.varname, '_model_rmse2_obs', '_', tmp.season]) = ...
            sqrt( mean(data2.([tmp.varname, '_model_err2_obs', '_', tmp.season]).^2, 3, 'omitnan') );
        data2.([tmp.varname, '_lens2_rmse2_obs', '_', tmp.season]) = ...
            sqrt( mean(data2.([tmp.varname, '_lens2_err2_obs', '_', tmp.season]).^2, 3, 'omitnan') );

    %% get nRMSE with ASSM based on anomaly (climatology removed, detrended, drift-removed)
        data2.([tmp.varname, '_model_nrmse2', '_', tmp.season])= ...
            data2.([tmp.varname, '_model_rmse2', '_', tmp.season])./data2.([tmp.varname, '_assm_ano_det2_stdt']);
        data2.([tmp.varname, '_lens2_nrmse2', '_', tmp.season])= ...
            data2.([tmp.varname, '_lens2_rmse2', '_', tmp.season])./data2.([tmp.varname, '_assm_ano_det2_stdt']);

    %% get nRMSE with obs based on anomaly (climatology removed, detrended, drift-removed)
        data2.([tmp.varname, '_model_nrmse2_obs', '_', tmp.season])= ...
            data2.([tmp.varname, '_model_rmse2_obs', '_', tmp.season])./data2.([tmp.varname, '_obs_ano_det2_stdt']);
        data2.([tmp.varname, '_lens2_nrmse2_obs', '_', tmp.season])= ...
            data2.([tmp.varname, '_lens2_rmse2_obs', '_', tmp.season])./data2.([tmp.varname, '_obs_ano_det2_stdt']);
        data2.([tmp.varname, '_assm_nrmse2_obs', '_', tmp.season])= ...
            data2.([tmp.varname, '_assm_rmse2_obs', '_', tmp.season])./data2.([tmp.varname, '_obs_ano_det2_stdt']);
    
    %% get MSSS (mean-square skill score, 1- (MSE_hcst/MSE_ref)
        data2.([tmp.varname, '_MSSS_assm_ano', '_', tmp.season]) = ...
            1 -  data2.([tmp.varname, '_model_rmse2', '_', tmp.season]).^2 ...
                 ./ data2.([tmp.varname, '_lens2_rmse2', '_', tmp.season]).^2;
        data2.([tmp.varname, '_MSSS_obs_ano', '_', tmp.season]) = ...
            1 -  data2.([tmp.varname, '_model_rmse2_obs', '_', tmp.season]).^2 ...
                 ./ data2.([tmp.varname, '_lens2_rmse2_obs', '_', tmp.season]).^2;
    
    disp('anomaly-related variables calculated')
    fprintf('%7.1f sec\n', toc(lap_time) );

    %% get svd modes (first 10 modes)
        tmp.svd_modes=5;
        % lens2 <-> hcst
        [data2.([tmp.varname, '_svd_lens2_model_lvmap_left', '_', tmp.season]), ...
            data2.([tmp.varname, '_svd_lens2_model_pcs_left', '_', tmp.season]), ...
            data2.([tmp.varname, '_svd_lens2_model_lvmap_right', '_', tmp.season]), ...
            data2.([tmp.varname, '_svd_lens2_model_pcs_right', '_', tmp.season]), ...
            data2.([tmp.varname, '_svd_lens2_model_lambda', '_', tmp.season]), ...
            data2.([tmp.varname, '_svd_lens2_model_scf', '_', tmp.season])] = ...
            mca( data2.([tmp.varname, '_lens2_ano', '_', tmp.season]), ... 
                 data2.([tmp.varname, '_model_ano', '_', tmp.season]), tmp.svd_modes);
        for orderi=1:tmp.svd_modes
            data2.([tmp.varname, '_svd_order_lens2_model', '_', tmp.season]) = ...
                order(mean(abs(data2.([tmp.varname, '_svd_lens2_model_pcs_left', '_', tmp.season])(orderi,:)),'all'));
            data2.([tmp.varname, '_svd_lens2_model_lvmap_left', '_', tmp.season])(:,:,orderi) = ...
                data2.([tmp.varname, '_svd_lens2_model_lvmap_left', '_', tmp.season])(:,:,orderi) ...
                .* 10.^data2.([tmp.varname, '_svd_order_lens2_model', '_', tmp.season]);
            data2.([tmp.varname, '_svd_lens2_model_pcs_left', '_', tmp.season])(orderi,:) = ...
                data2.([tmp.varname, '_svd_lens2_model_pcs_left', '_', tmp.season])(orderi,:) ...
                .* 10.^(-data2.([tmp.varname, '_svd_order_lens2_model', '_', tmp.season]));
            data2.([tmp.varname, '_svd_lens2_model_lvmap_right', '_', tmp.season])(:,:,orderi) = ...
                data2.([tmp.varname, '_svd_lens2_model_lvmap_right', '_', tmp.season])(:,:,orderi) ...
                .* 10.^data2.([tmp.varname, '_svd_order_lens2_model', '_', tmp.season]);
            data2.([tmp.varname, '_svd_lens2_model_pcs_right', '_', tmp.season])(orderi,:) = ...
                data2.([tmp.varname, '_svd_lens2_model_pcs_right', '_', tmp.season])(orderi,:) ...
                .* 10.^(-data2.([tmp.varname, '_svd_order_lens2_model', '_', tmp.season]));
        end
        tmp.pcs=data2.([tmp.varname, '_svd_lens2_model_pcs_right', '_', tmp.season])(1,:);
        tmp.pcs=reshape(tmp.pcs,[1,1,length(tmp.pcs)]);
        data2.([tmp.varname, '_svd_lens2_model_right_1mode_recon', '_', tmp.season]) = ...
            data2.([tmp.varname, '_svd_lens2_model_lvmap_right', '_', tmp.season])(:,:,1).*tmp.pcs(1,1,:);

        % lens2 <-> assm
        [data2.([tmp.varname, '_svd_lens2_assm_lvmap_left', '_', tmp.season]), ...
            data2.([tmp.varname, '_svd_lens2_assm_pcs_left', '_', tmp.season]), ...
            data2.([tmp.varname, '_svd_lens2_assm_lvmap_right', '_', tmp.season]), ...
            data2.([tmp.varname, '_svd_lens2_assm_pcs_right', '_', tmp.season]), ...
            data2.([tmp.varname, '_svd_lens2_assm_lambda', '_', tmp.season]), ...
            data2.([tmp.varname, '_svd_lens2_assm_scf', '_', tmp.season])] = ...
            mca( data2.([tmp.varname, '_lens2_ano', '_', tmp.season]), ... 
                 data2.([tmp.varname, '_assm_ano', '_', tmp.season]), tmp.svd_modes);
        for orderi=1:tmp.svd_modes
            data2.([tmp.varname, '_svd_order_lens2_assm', '_', tmp.season]) = ...
                order(mean(abs(data2.([tmp.varname, '_svd_lens2_assm_pcs_left', '_', tmp.season])(orderi,:)),'all'));
            data2.([tmp.varname, '_svd_lens2_assm_lvmap_left', '_', tmp.season])(:,:,orderi) = ...
                data2.([tmp.varname, '_svd_lens2_assm_lvmap_left', '_', tmp.season])(:,:,orderi) ...
                .* 10.^data2.([tmp.varname, '_svd_order_lens2_assm', '_', tmp.season]);
            data2.([tmp.varname, '_svd_lens2_assm_pcs_left', '_', tmp.season])(orderi,:) = ...
                data2.([tmp.varname, '_svd_lens2_assm_pcs_left', '_', tmp.season])(orderi,:) ...
                .* 10.^(-data2.([tmp.varname, '_svd_order_lens2_assm', '_', tmp.season]));
            data2.([tmp.varname, '_svd_lens2_assm_lvmap_right', '_', tmp.season])(:,:,orderi) = ...
                data2.([tmp.varname, '_svd_lens2_assm_lvmap_right', '_', tmp.season])(:,:,orderi) ...
                .* 10.^data2.([tmp.varname, '_svd_order_lens2_assm', '_', tmp.season]);
            data2.([tmp.varname, '_svd_lens2_assm_pcs_right', '_', tmp.season])(orderi,:)= ...
                data2.([tmp.varname, '_svd_lens2_assm_pcs_right', '_', tmp.season])(orderi,:) ...
                .* 10.^(-data2.([tmp.varname, '_svd_order_lens2_assm', '_', tmp.season]));
        end
        tmp.pcs(1,1,:)=data2.([tmp.varname, '_svd_lens2_assm_pcs_right', '_', tmp.season])(1,:);
        tmp.pcs=reshape(tmp.pcs,[1,1,length(tmp.pcs)]);
        data2.([tmp.varname, '_svd_lens2_assm_right_1mode_recon', '_', tmp.season]) = ...
            data2.([tmp.varname, '_svd_lens2_assm_lvmap_right', '_', tmp.season])(:,:,1).*tmp.pcs(1,1,:);

        % lens2 <-> obs
         if isnan(sum( data2.([tmp.varname, '_obs_ano', '_', tmp.season])(:)))
            data2.([tmp.varname, '_svd_lens2_obs_lvmap_left', '_', tmp.season]) = ...
                NaN(size(data2.([tmp.varname, '_svd_lens2_assm_lvmap_left', '_', tmp.season])));
            data2.([tmp.varname, '_svd_lens2_obs_pcs_left', '_', tmp.season]) = ...
                NaN(size(data2.([tmp.varname, '_svd_lens2_assm_pcs_left', '_', tmp.season])));
            data2.([tmp.varname, '_svd_lens2_obs_lvmap_right', '_', tmp.season]) = ...
                NaN(size(data2.([tmp.varname, '_svd_lens2_assm_lvmap_left', '_', tmp.season])));
            data2.([tmp.varname, '_svd_lens2_obs_pcs_right', '_', tmp.season]) = ...
                NaN(size(data2.([tmp.varname, '_svd_lens2_assm_pcs_left', '_', tmp.season])));
            data2.([tmp.varname, '_svd_lens2_obs_lambda', '_', tmp.season]) = ...
                NaN(size(data2.([tmp.varname, '_svd_lens2_model_lambda', '_', tmp.season])));
            data2.([tmp.varname, '_svd_lens2_obs_scf', '_', tmp.season]) = ...
                NaN(size(data2.([tmp.varname, '_svd_lens2_model_scf', '_', tmp.season])));
        else
            [data2.([tmp.varname, '_svd_lens2_obs_lvmap_left', '_', tmp.season]), ...
                data2.([tmp.varname, '_svd_lens2_obs_pcs_left', '_', tmp.season]), ...
                data2.([tmp.varname, '_svd_lens2_obs_lvmap_right', '_', tmp.season]), ...
                data2.([tmp.varname, '_svd_lens2_obs_pcs_right', '_', tmp.season]), ...
                data2.([tmp.varname, '_svd_lens2_obs_lambda', '_', tmp.season]), ...
                data2.([tmp.varname, '_svd_lens2_obs_scf', '_', tmp.season])] = ...
                mca( data2.([tmp.varname, '_lens2_ano', '_', tmp.season]), ... 
                     data2.([tmp.varname, '_obs_ano', '_', tmp.season]), tmp.svd_modes);
            for orderi=1:tmp.svd_modes
                data2.([tmp.varname, '_svd_order_lens2_obs', '_', tmp.season]) = ...
                    order(mean(abs(data2.([tmp.varname, '_svd_lens2_obs_pcs_left', '_', tmp.season])(orderi,:)),'all'));
                data2.([tmp.varname, '_svd_lens2_obs_lvmap_left', '_', tmp.season])(:,:,orderi) = ...
                    data2.([tmp.varname, '_svd_lens2_obs_lvmap_left', '_', tmp.season])(:,:,orderi) ...
                    .* 10.^data2.([tmp.varname, '_svd_order_lens2_obs', '_', tmp.season]);
                data2.([tmp.varname, '_svd_lens2_obs_pcs_left', '_', tmp.season])(orderi,:) = ...
                    data2.([tmp.varname, '_svd_lens2_obs_pcs_left', '_', tmp.season])(orderi,:) ...
                    .* 10.^(-data2.([tmp.varname, '_svd_order_lens2_obs', '_', tmp.season]));
                data2.([tmp.varname, '_svd_lens2_obs_lvmap_right', '_', tmp.season])(:,:,orderi) = ...
                    data2.([tmp.varname, '_svd_lens2_obs_lvmap_right', '_', tmp.season])(:,:,orderi) ...
                    .* 10.^data2.([tmp.varname, '_svd_order_lens2_obs', '_', tmp.season]);
                data2.([tmp.varname, '_svd_lens2_obs_pcs_right', '_', tmp.season])(orderi,:) = ...
                    data2.([tmp.varname, '_svd_lens2_obs_pcs_right', '_', tmp.season])(orderi,:)...
                    .* 10.^(-data2.([tmp.varname, '_svd_order_lens2_obs', '_', tmp.season]));
            end
        end


%         [data2.([tmp.varname, '_svd_lens2_obs_lvmap_left', '_', tmp.season]), ...
%             data2.([tmp.varname, '_svd_lens2_obs_pcs_left', '_', tmp.season]), ...
%             data2.([tmp.varname, '_svd_lens2_obs_lvmap_right', '_', tmp.season]), ...
%             data2.([tmp.varname, '_svd_lens2_obs_pcs_right', '_', tmp.season]), ...
%             data2.([tmp.varname, '_svd_lens2_obs_lambda', '_', tmp.season]), ...
%             data2.([tmp.varname, '_svd_lens2_obs_scf', '_', tmp.season])] = ...
%             mca( data2.([tmp.varname, '_lens2_ano', '_', tmp.season]), ... 
%                  data2.([tmp.varname, '_obs_ano', '_', tmp.season]), tmp.svd_modes);
%         for orderi=1:tmp.svd_modes
%             data2.([tmp.varname, '_svd_order_lens2_obs', '_', tmp.season]) = ...
%                 order(mean(abs(data2.([tmp.varname, '_svd_lens2_obs_pcs_left', '_', tmp.season])(orderi,:)),'all'));
%             data2.([tmp.varname, '_svd_lens2_obs_lvmap_left', '_', tmp.season])(:,:,orderi) = ...
%                 data2.([tmp.varname, '_svd_lens2_obs_lvmap_left', '_', tmp.season])(:,:,orderi) ...
%                 .* 10.^data2.([tmp.varname, '_svd_order_lens2_obs', '_', tmp.season]);
%             data2.([tmp.varname, '_svd_lens2_obs_pcs_left', '_', tmp.season])(orderi,:) = ...
%                 data2.([tmp.varname, '_svd_lens2_obs_pcs_left', '_', tmp.season])(orderi,:) ...
%                 .* 10.^(-data2.([tmp.varname, '_svd_order_lens2_obs', '_', tmp.season]));
%             data2.([tmp.varname, '_svd_lens2_obs_lvmap_right', '_', tmp.season])(:,:,orderi) = ...
%                 data2.([tmp.varname, '_svd_lens2_obs_lvmap_right', '_', tmp.season])(:,:,orderi) ...
%                 .* 10.^data2.([tmp.varname, '_svd_order_lens2_obs', '_', tmp.season]);
%             data2.([tmp.varname, '_svd_lens2_obs_pcs_right', '_', tmp.season])(orderi,:) = ...
%                 data2.([tmp.varname, '_svd_lens2_obs_pcs_right', '_', tmp.season])(orderi,:)...
%                 .* 10.^(-data2.([tmp.varname, '_svd_order_lens2_obs', '_', tmp.season]));
%         end





        tmp.pcs(1,1,:)=data2.([tmp.varname, '_svd_lens2_obs_pcs_right', '_', tmp.season])(1,:);
        tmp.pcs=reshape(tmp.pcs,[1,1,length(tmp.pcs)]);
        data2.([tmp.varname, '_svd_lens2_obs_right_1mode_recon', '_', tmp.season]) = ...
            data2.([tmp.varname, '_svd_lens2_obs_lvmap_right', '_', tmp.season])(:,:,1).*tmp.pcs(1,1,:);

    disp('svd finished')
    fprintf('%7.1f sec\n', toc(lap_time) );

   %% AR1 parameter estimation, getting prognostic value (assm based), should be estimated from anomaly
        data2.([tmp.varname, '_AR1', '_', tmp.season])=NaN(grid.nlon, grid.nlat, length(cfg.clim_tlen));  %% AR1 initialization
        data2.([tmp.varname, '_obs_AR1', '_', tmp.season])=NaN(grid.nlon, grid.nlat, length(cfg.clim_tlen));  %% AR1 initialization
        
%         if lyear==0
            for loni=1:grid.nlon
                for lati=1:grid.nlat
                    if isnan(data2.([tmp.varname, '_assm_ano', '_', tmp.season])(loni,lati))~=1
%                         tmp.init_finite=squeeze(data2.([tmp.varname, '_assm_ano', '_', tmp.season])(loni,lati,:));
                        tmp.init_finite=data2_l0_AR1.assm_ano_full(loni, lati, cfg.clim_ys-1960:cfg.clim_ye-1960);
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
                if isnan(data2.([tmp.varname, '_assm_ano', '_', tmp.season])(loni,lati))~=1
%                     prog=data2_l0.([tmp.varname,'_assm_ano', '_', tmp.season])(loni,lati,:); %initialize
%                     prog=tmp.assm_ano(loni,lati,1:end-1); %initialize  
                      prog=data2_l0_AR1.assm_ano_full(loni, lati, cfg.clim_ys-1960:cfg.clim_ye-1960); %initialize
                    for i=1:1 % integration of autoregressive model % seasonal -> ly=1
                        prog = -prog .* (data2.([tmp.varname,'_assm_AR1_coef'])(loni,lati)) ...
                            + (data2.([tmp.varname,'_assm_AR1_noise'])(loni,lati));
                    end
                    data2.([tmp.varname,'_AR1_', tmp.season])(loni,lati,:)=prog;
                end
            end
        end

   %% AR1 parameter estimation, getting prognostic value (obs based)
%         if lyear==0 
            for loni=1:grid.nlon
                for lati=1:grid.nlat
                    if isnan(data2.([tmp.varname, '_obs_ano', '_', tmp.season])(loni,lati))~=1
%                         tmp.init_finite=squeeze(data2.([tmp.varname, '_obs_ano', '_', tmp.season])(loni,lati,:));
%                         tmp.init_finite=tmp.obs_ano(loni,lati,1:end-1);
                        tmp.init_finite=data2_l0_AR1.obs_ano_full(loni, lati, cfg.clim_ys-1960:cfg.clim_ye-1960);
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
                if isnan(data2.([tmp.varname, '_obs_ano', '_', tmp.season])(loni,lati))~=1
%                     prog=data2_l0.([tmp.varname,'_obs_ano', '_', tmp.season])(loni,lati,:);
%                     prog=tmp.obs_ano(loni,lati,1:end-1);        
                    prog=data2_l0_AR1.obs_ano_full(loni, lati, cfg.clim_ys-1960:cfg.clim_ye-1960); 
                    for i=1:1 % integration of autoregressive model % seasonal -> ly=1
                        prog = -prog .* (data2.([tmp.varname,'_obs_AR1_coef'])(loni,lati)) ...
                            + (data2.([tmp.varname,'_obs_AR1_noise'])(loni,lati));
                    end
                    data2.([tmp.varname,'_obs_AR1_', tmp.season])(loni,lati,:)=prog;
                end
            end
        end

    disp('AR1 estimation, integration finished')
    fprintf('%7.1f sec\n', toc(lap_time) );
        %% get correlation coefficient
        for loni=1:grid.nlon
            for lati=1:grid.nlat
                 if (isnan(data.([tmp.varname, '_assm'])(loni,lati,1))~=1 & nansum(data.([tmp.varname, '_model_', tmp.season])(loni,lati,:))~=0)
    
                 %% corr assm
                     tmp.data_assm = data.([tmp.varname, '_assm'])(loni,lati,:);
                     tmp.data = squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,:));
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_assm)), tmp.data_assm(isfinite(tmp.data_assm)));
                     [tmp.data_det, tmp.data_trend] = Func_0028_detrend_linear_1d(squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,:)), 'omitnan');                         
                     [tmp.data_assm_det, tmp.data_assm_trend] = Func_0028_detrend_linear_1d(squeeze(data.([tmp.varname, '_assm'])(loni,lati,:)), 'omitnan');
                     [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data_det(isfinite(tmp.data_assm_det)), tmp.data_assm_det(isfinite(tmp.data_assm_det)));

                     data.([tmp.varname, '_corr_assm', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_assm_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
                     data.([tmp.varname, '_corr_assm_det', '_', tmp.season])(loni,lati)=tmp.corr_det(1,2);
                     data.([tmp.varname, '_corr_assm_det_p', '_', tmp.season])(loni,lati)=tmp.corr_det_p(1,2);
                     
                     data.([tmp.varname, '_assm_trend', '_', tmp.season])(loni,lati)=tmp.data_assm_trend;
                     data.([tmp.varname, '_model_trend', '_', tmp.season])(loni,lati)=tmp.data_trend;

                 %% hcst-lens2 ~ assm
                     tmp.data_hcst_lens2 = squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,:)) - ...
                                    squeeze(data.([tmp.varname, '_lens2', '_', tmp.season])(loni,lati,:));
                     [tmp.corr_int, tmp.corr_int_p]=corrcoef(tmp.data_hcst_lens2(isfinite(tmp.data_assm)), tmp.data_assm(isfinite(tmp.data_assm)));
    
                     
                     data.([tmp.varname, '_corr_assm_int', '_', tmp.season])(loni,lati)=tmp.corr_int(1,2);
                     data.([tmp.varname, '_corr_assm_int_p', '_', tmp.season])(loni,lati)=tmp.corr_int_p(1,2);

                 %% hcst-lens2 ~ assm-lens2
                     tmp.data_hcst_lens2 = squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,:)) - ...
                                    squeeze(data.([tmp.varname, '_lens2', '_', tmp.season])(loni,lati,:));
                     tmp.data_assm_lens2 = squeeze(data.([tmp.varname, '_assm'])(loni,lati,:)) - ...
                                    squeeze(data.([tmp.varname, '_lens2', '_', tmp.season])(loni,lati,:));
                     [tmp.corr_int, tmp.corr_int_p]=corrcoef(tmp.data_hcst_lens2(isfinite(tmp.data_assm_lens2)), tmp.data_assm_lens2(isfinite(tmp.data_assm_lens2)));
    
                     
                     data.([tmp.varname, '_corr_assm_int2', '_', tmp.season])(loni,lati)=tmp.corr_int(1,2);
                     data.([tmp.varname, '_corr_assm_int2_p', '_', tmp.season])(loni,lati)=tmp.corr_int_p(1,2);
    
                 %% corr lens2 ~ assm
                     tmp.lens2 = squeeze(data.([tmp.varname, '_lens2', '_', tmp.season])(loni,lati,:));
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.lens2(isfinite(tmp.data_assm)), tmp.data_assm(isfinite(tmp.data_assm)));
                     [tmp.data_lens2_det, tmp.data_lens2_trend] = Func_0028_detrend_linear_1d(squeeze(data.([tmp.varname, '_lens2', '_', tmp.season])(loni,lati,:)), 'omitnan');                         
                     [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data_lens2_det(isfinite(tmp.data_assm_det)), tmp.data_assm_det(isfinite(tmp.data_assm_det)));

                     data.([tmp.varname, '_corr_assm_lens2', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_assm_lens2_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
                     data.([tmp.varname, '_corr_assm_lens2_det', '_', tmp.season])(loni,lati)=tmp.corr_det(1,2);
                     data.([tmp.varname, '_corr_assm_lens2_det_p', '_', tmp.season])(loni,lati)=tmp.corr_det_p(1,2);

                     data.([tmp.varname, '_lens2_trend', '_', tmp.season])(loni,lati)=tmp.data_lens2_trend;

%                  %% corr AR1(assm) ~ assm
%                     tmp.AR1 = squeeze(data2.([tmp.varname, '_AR1', '_', tmp.season])(loni,lati,:));
%                     tmp.AR1_mask=squeeze(tmp.AR1./tmp.AR1);
%                     tmp.data_assm_mask=squeeze(tmp.data_assm./tmp.data_assm);
%                      [tmp.corr, tmp.corr_p]=corrcoef(tmp.AR1.*tmp.data_assm_mask.*tmp.AR1_mask, ...
%                          squeeze(tmp.data_assm).*tmp.data_assm_mask.*tmp.AR1_mask, 'Rows', 'complete');
%                      data.([tmp.varname, '_corr_assm_AR1', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
%                      data.([tmp.varname, '_corr_assm_AR1_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
    
                
                 %% corr obs ~ hcst
                     if isfield(tmp, 'ydata_mod_obs_masked')
                        tmp.data = squeeze(data.([tmp.varname, '_mod_obs_masked', '_', tmp.season])(loni,lati,:));
                     else
                        tmp.data = squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,:));
                     end
                     tmp.data_obs = squeeze(data.([tmp.varname, '_obs'])(loni,lati,:));
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
                     if isfield(tmp, 'ydata_mod_obs_masked')
                        tmp.data_det = Func_0028_detrend_linear_1d(squeeze(data.([tmp.varname, '_mod_obs_masked', '_', tmp.season])(loni,lati,:)), 'omitnan');                                              
                     else
                        tmp.data_det = Func_0028_detrend_linear_1d(squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,:)), 'omitnan');                         
                     end
                     [tmp.data_obs_det, tmp.data_obs_trend] = Func_0028_detrend_linear_1d(data.([tmp.varname, '_obs'])(loni,lati,:), 'omitnan');
                     [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data_det(isfinite(tmp.data_obs_det)), tmp.data_obs_det(isfinite(tmp.data_obs_det)));
    %                  [tmp.corr, tmp.corr_p]=corrcoef(squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,1+lyear:end-cfg.proj_year+1)), ...
    %                      squeeze(data.([tmp.varname, '_obs'])(loni,lati,1+lyear:end)));

                     if sum(squeeze(isfinite(tmp.data_obs)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                        tmp.corr_det=NaN(2,2);
                        tmp.corr_det_p=NaN(2,2);
                        tmp.data_obs_trend=NaN;
                     end

                     data.([tmp.varname, '_corr_obs', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_obs_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
                     data.([tmp.varname, '_corr_obs_det', '_', tmp.season])(loni,lati)=tmp.corr_det(1,2);
                     data.([tmp.varname, '_corr_obs_det_p', '_', tmp.season])(loni,lati)=tmp.corr_det_p(1,2);

                     data.([tmp.varname, '_obs_trend', '_', tmp.season])(loni,lati)=tmp.data_obs_trend;

                 %% corr obs ~ assm
                    if isfield(tmp, 'ydata_mod_obs_masked')
                        tmp.data_assm = squeeze(data.([tmp.varname, '_assm_obs_masked', '_', tmp.season])(loni,lati,:));
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

                     data.([tmp.varname, '_corr_obs_assm', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_obs_assm_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
                     data.([tmp.varname, '_corr_obs_assm_det', '_', tmp.season])(loni,lati)=tmp.corr_det(1,2);
                     data.([tmp.varname, '_corr_obs_assm_det_p', '_', tmp.season])(loni,lati)=tmp.corr_det_p(1,2);

                 %% corr obs ~ lens2
                    if isfield(tmp, 'ydata_mod_obs_masked')
                        tmp.data_lens2 = squeeze(data.([tmp.varname, '_lens2_obs_masked', '_', tmp.season])(loni,lati,:));
                     else
                        tmp.data_lens2 = squeeze(data.([tmp.varname, '_lens2', '_', tmp.season])(loni,lati,:));
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

                     data.([tmp.varname, '_corr_obs_lens2', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_obs_lens2_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
                     data.([tmp.varname, '_corr_obs_lens2_det', '_', tmp.season])(loni,lati)=tmp.corr_det(1,2);
                     data.([tmp.varname, '_corr_obs_lens2_det_p', '_', tmp.season])(loni,lati)=tmp.corr_det_p(1,2);


                 %% corr obs ~ hcst-lens2
                    tmp.data_hcst_lens2 = squeeze(tmp.data) - squeeze(tmp.data_lens2);
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data_hcst_lens2(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
                     if sum(squeeze(isfinite(tmp.data_obs)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                     end
                     data.([tmp.varname, '_corr_obs_int', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_obs_int_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
                 
                %% corr obs-lens2 ~ hcst-lens2
                    tmp.data_hcst_lens2 = squeeze(tmp.data) - squeeze(tmp.data_lens2);
                     tmp.data_obs_lens2 = squeeze(tmp.data_obs) - squeeze(tmp.data_lens2);
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data_hcst_lens2(isfinite(tmp.data_obs_lens2)), tmp.data_obs_lens2(isfinite(tmp.data_obs_lens2)));
                    
                     if sum(squeeze(isfinite(tmp.data_obs)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                     end
                     data.([tmp.varname, '_corr_obs_int2', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_obs_int2_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);

%                 %% corr AR1(obs) ~ obs
% %                     tmp.obs_AR1 = squeeze(data.([tmp.varname, '_obs_AR1', '_', tmp.season])(loni,lati,:));
% %                      [tmp.corr, tmp.corr_p]=corrcoef(tmp.obs_AR1(isfinite(tmp.data_assm).*isfinite(tmp.obs_AR1)), ...
% %                          tmp.data_obs(isfinite(tmp.data_assm).*isfinite(tmp.obs_AR1)));
%                     tmp.obs_AR1 = squeeze(data2.([tmp.varname, '_obs_AR1', '_', tmp.season])(loni,lati,:));
%                     tmp.obs_AR1_mask=squeeze(tmp.obs_AR1./tmp.obs_AR1);
%                     tmp.data_obs_mask=squeeze(tmp.data_obs./tmp.data_obs);
%                      [tmp.corr, tmp.corr_p]=corrcoef(tmp.obs_AR1.*tmp.data_obs_mask.*tmp.obs_AR1_mask, ...
%                          squeeze(tmp.data_obs).*tmp.data_obs_mask.*tmp.obs_AR1_mask, 'Rows', 'complete');
% 
%                      data.([tmp.varname, '_corr_obs_AR1', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
%                      data.([tmp.varname, '_corr_obs_AR1_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
    
                 else
                    %% obs (NaN)
                     data.([tmp.varname, '_corr_obs', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_p', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_det', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_det_p', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_obs_trend', '_', tmp.season])(loni,lati)=NaN;
    
                    %% AR1 (NaN)
                     data.([tmp.varname, '_corr_assm_AR1', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_AR1_p', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_AR1', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_AR1_p', '_', tmp.season])(loni,lati)=NaN;
    
                    %% ASSM (NaN)
                     data.([tmp.varname, '_corr_assm', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_p', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_det', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_det_p', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_int', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_int_p', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_int2', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_int2_p', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_assm_trend', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_model_trend', '_', tmp.season])(loni,lati)=NaN;
    
                    %% LENS2 (NaN)
                     data.([tmp.varname, '_corr_assm_lens2', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_lens2_p', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_lens2_det', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_lens2_det_p', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_lens2_trend', '_', tmp.season])(loni,lati)=NaN;
                     %% obs ~ ASSM (NaN)
                     data.([tmp.varname, '_corr_obs_assm', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_assm_p', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_assm_det', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_assm_det_p', '_', tmp.season])(loni,lati)=NaN;
                     
                     %% obs ~ LENS2 (NaN)
                     data.([tmp.varname, '_corr_obs_lens2', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_lens2_p', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_lens2_det', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_lens2_det_p', '_', tmp.season])(loni,lati)=NaN;

                     %% obs ~ HCST-LENS2 (NaN)
                     data.([tmp.varname, '_corr_obs_int', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_int_p', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_int2', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_int2_p', '_', tmp.season])(loni,lati)=NaN;
                 end
            end
        end
        disp('correlation1 finished')
        fprintf('%7.1f sec\n', toc(lap_time) );

        %% get correlation coefficient2 (anomaly based)
        for loni=1:grid.nlon
            for lati=1:grid.nlat
%                  if (~isnan(data2.([tmp.varname, '_assm_ano', '_', tmp.season])(loni,lati,1)) ...
%                          && sum(data2.([tmp.varname, '_model_ano_', tmp.season])(loni,lati,:),'omitnan')~=0)
                if (~isnan(data2.([tmp.varname, '_assm_ano', '_', tmp.season])(loni,lati,1)))

                 %% corr assm (anomaly based)
                     tmp.data2_assm = data2.([tmp.varname, '_assm_ano','_', tmp.season])(loni,lati,:);
                     tmp.data2 = squeeze(data2.([tmp.varname, '_model_ano', '_', tmp.season])(loni,lati,:));
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data2(isfinite(tmp.data2_assm)), tmp.data2_assm(isfinite(tmp.data2_assm)));
                     tmp.data2_det = data2.([tmp.varname, '_model_ano_det2', '_', tmp.season])(loni,lati,:);
                     tmp.data2_assm_det = data2.([tmp.varname, '_assm_ano_det2', '_', tmp.season])(loni,lati,:);
                     [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data2_det(isfinite(tmp.data2_assm_det)), tmp.data2_assm_det(isfinite(tmp.data2_assm_det)));
                                    
                     data2.([tmp.varname, '_corr_assm_ano', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data2.([tmp.varname, '_corr_assm_ano_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
                     data2.([tmp.varname, '_corr_assm_det_ano', '_', tmp.season])(loni,lati)=tmp.corr_det(1,2);
                     data2.([tmp.varname, '_corr_assm_det_ano_p', '_', tmp.season])(loni,lati)=tmp.corr_det_p(1,2);


                 %% hcst-lens2 ~ assm (anomaly based)
                     tmp.data2_hcst_lens2 = squeeze(data2.([tmp.varname, '_model_ano', '_', tmp.season])(loni,lati,:)) - ...
                                    squeeze(data2.([tmp.varname, '_lens2_ano', '_', tmp.season])(loni,lati,:));
                     [tmp.corr_int, tmp.corr_int_p]=corrcoef(tmp.data2_hcst_lens2(isfinite(tmp.data2_assm)), tmp.data2_assm(isfinite(tmp.data2_assm)));
    
                     
                     data2.([tmp.varname, '_corr_assm_int_ano', '_', tmp.season])(loni,lati)=tmp.corr_int(1,2);
                     data2.([tmp.varname, '_corr_assm_int_ano_p', '_', tmp.season])(loni,lati)=tmp.corr_int_p(1,2);
                 
                 %% hcst-lens2 ~ assm-lens2 (anomaly based)
                     tmp.data2_hcst_lens2 = squeeze(data2.([tmp.varname, '_model_ano', '_', tmp.season])(loni,lati,:)) - ...
                                    squeeze(data2.([tmp.varname, '_lens2_ano', '_', tmp.season])(loni,lati,:));
                     tmp.data2_assm_lens2 = squeeze(data2.([tmp.varname, '_assm_ano', '_', tmp.season])(loni,lati,:)) - ...
                                    squeeze(data2.([tmp.varname, '_lens2_ano', '_', tmp.season])(loni,lati,:));

                     [tmp.corr_int, tmp.corr_int_p]=corrcoef(tmp.data2_hcst_lens2(isfinite(tmp.data2_assm_lens2)), tmp.data2_assm_lens2(isfinite(tmp.data2_assm_lens2)));
    
                     
                     data2.([tmp.varname, '_corr_assm_int2_ano', '_', tmp.season])(loni,lati)=tmp.corr_int(1,2);
                     data2.([tmp.varname, '_corr_assm_int2_ano_p', '_', tmp.season])(loni,lati)=tmp.corr_int_p(1,2);
                
                 %% hcst-svd 1st mode ~ assm (anomaly based)
                    tmp.data2_hcst_svd = squeeze(data2.([tmp.varname, '_model_ano', '_', tmp.season])(loni,lati,:)) - ...
                                    squeeze(data2.([tmp.varname, '_svd_lens2_model_right_1mode_recon', '_', tmp.season])(loni,lati,:));
                    [tmp.corr_int2, tmp.corr_int2_p]=corrcoef(tmp.data2_hcst_svd(isfinite(tmp.data2_assm)), tmp.data2_assm(isfinite(tmp.data2_assm)));
                    data2.([tmp.varname, '_corr_assm_int_svd', '_', tmp.season])(loni,lati)=tmp.corr_int2(1,2);
                    data2.([tmp.varname, '_corr_assm_int_svd_p', '_', tmp.season])(loni,lati)=tmp.corr_int2_p(1,2);

                %% hcst-svd 1st mode ~ assm-svd 1st mode (anomaly based)
                    tmp.data2_hcst_svd = squeeze(data2.([tmp.varname, '_model_ano', '_', tmp.season])(loni,lati,:)) - ...
                                    squeeze(data2.([tmp.varname, '_svd_lens2_model_right_1mode_recon', '_', tmp.season])(loni,lati,:));
                    tmp.data2_assm_svd = squeeze(data2.([tmp.varname, '_assm_ano', '_', tmp.season])(loni,lati,:)) - ...
                                    squeeze(data2.([tmp.varname, '_svd_lens2_assm_right_1mode_recon', '_', tmp.season])(loni,lati,:));
                    [tmp.corr_int2, tmp.corr_int2_p]=corrcoef(tmp.data2_hcst_svd(isfinite(tmp.data2_assm_svd)), tmp.data2_assm_svd(isfinite(tmp.data2_assm_svd)));
                    data2.([tmp.varname, '_corr_assm_int2_svd', '_', tmp.season])(loni,lati)=tmp.corr_int2(1,2);
                    data2.([tmp.varname, '_corr_assm_int2_svd_p', '_', tmp.season])(loni,lati)=tmp.corr_int2_p(1,2);
                

                 %% corr lens2 ~ assm (anomaly based)
                     tmp.lens2 = squeeze(data2.([tmp.varname, '_lens2_ano', '_', tmp.season])(loni,lati,:));
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.lens2(isfinite(tmp.data2_assm)), tmp.data2_assm(isfinite(tmp.data2_assm)));
                     tmp.data2_lens2_det = data2.([tmp.varname, '_lens2_ano_det2', '_', tmp.season])(loni,lati,:);
                     [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data2_lens2_det(isfinite(tmp.data2_assm_det)), tmp.data2_assm_det(isfinite(tmp.data2_assm_det)));

                     data2.([tmp.varname, '_corr_assm_lens2_ano', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data2.([tmp.varname, '_corr_assm_lens2_ano_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
                     data2.([tmp.varname, '_corr_assm_lens2_det_ano', '_', tmp.season])(loni,lati)=tmp.corr_det(1,2);
                     data2.([tmp.varname, '_corr_assm_lens2_det_ano_p', '_', tmp.season])(loni,lati)=tmp.corr_det_p(1,2);

                 %% corr AR1(assm) ~ assm (anomaly based)
                
                    tmp.AR1 = squeeze(data2.([tmp.varname, '_AR1', '_', tmp.season])(loni,lati,:));
                    tmp.AR1_mask=squeeze(tmp.AR1./tmp.AR1);
                    tmp.data2_assm_mask=squeeze(tmp.data2_assm./tmp.data2_assm);
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.AR1.*tmp.data2_assm_mask.*tmp.AR1_mask, ...
                         squeeze(tmp.data2_assm).*tmp.data2_assm_mask.*tmp.AR1_mask, 'Rows', 'complete');
                     data2.([tmp.varname, '_corr_assm_AR1_ano', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data2.([tmp.varname, '_corr_assm_AR1_ano_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
    
                
                 %% corr obs ~ hcst (anomaly based)
                     if isfield(tmp, 'ydata_mod_obs_masked')
                        tmp.data2 = squeeze(data.([tmp.varname, '_mod_obs_masked', '_', tmp.season])(loni,lati,:)); %% not anomaly, caution!
                     else
                        tmp.data2 = squeeze(data2.([tmp.varname, '_model_ano', '_', tmp.season])(loni,lati,:));
                     end
                     tmp.data2_obs = squeeze(data2.([tmp.varname, '_obs_ano', '_', tmp.season])(loni,lati,:));
                     if sum(isfinite(tmp.data2_obs)) ==1
                         tmp.data2_obs=NaN(size(tmp.data2_obs));
                     end
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data2(isfinite(tmp.data2_obs)), tmp.data2_obs(isfinite(tmp.data2_obs)));
                     tmp.data2_obs_det = data2.([tmp.varname, '_obs_ano_det2', '_', tmp.season])(loni,lati,:);
                     [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data2_det(isfinite(tmp.data2_obs_det)), tmp.data2_obs_det(isfinite(tmp.data2_obs_det)));

                     if sum(squeeze(isfinite(tmp.data2_obs)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                        tmp.corr_det=NaN(2,2);
                        tmp.corr_det_p=NaN(2,2);
                        tmp.data2_obs_trend=NaN;
                     end

                     data2.([tmp.varname, '_corr_obs_ano', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data2.([tmp.varname, '_corr_obs_ano_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
                     data2.([tmp.varname, '_corr_obs_det_ano', '_', tmp.season])(loni,lati)=tmp.corr_det(1,2);
                     data2.([tmp.varname, '_corr_obs_det_ano_p', '_', tmp.season])(loni,lati)=tmp.corr_det_p(1,2);    

                 %% corr obs ~ assm (anomaly based)
                    if isfield(tmp, 'ydata_mod_obs_masked')
                        tmp.data2_assm = squeeze(data.([tmp.varname, '_assm_obs_masked', '_', tmp.season])(loni,lati,:));
                     else
                        tmp.data2_assm = squeeze(data2.([tmp.varname, '_assm_ano', '_', tmp.season])(loni,lati,:));
                     end
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data2_assm(isfinite(tmp.data2_obs)), tmp.data2_obs(isfinite(tmp.data2_obs)));
                     
                     [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data2_assm_det(isfinite(tmp.data2_obs_det)), tmp.data2_obs_det(isfinite(tmp.data2_obs_det)));

                     if sum(squeeze(isfinite(tmp.data2_obs)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                        tmp.corr_det=NaN(2,2);
                        tmp.corr_det_p=NaN(2,2);
                     end

                     data2.([tmp.varname, '_corr_obs_assm_ano', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data2.([tmp.varname, '_corr_obs_assm_ano_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
                     data2.([tmp.varname, '_corr_obs_assm_det_ano', '_', tmp.season])(loni,lati)=tmp.corr_det(1,2);
                     data2.([tmp.varname, '_corr_obs_assm_det_ano_p', '_', tmp.season])(loni,lati)=tmp.corr_det_p(1,2);

                 %% corr obs ~ lens2 (anomaly based)
                    if isfield(tmp, 'ydata_mod_obs_masked')
                        tmp.data2_lens2 = squeeze(data.([tmp.varname, '_lens2_obs_masked', '_', tmp.season])(loni,lati,:));
                     else
                        tmp.data2_lens2 = squeeze(data2.([tmp.varname, '_lens2_ano', '_', tmp.season])(loni,lati,:));
                     end
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data2_lens2(isfinite(tmp.data2_obs)), tmp.data2_obs(isfinite(tmp.data2_obs)));
                     
                     [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data2_lens2_det(isfinite(tmp.data2_obs_det)), tmp.data2_obs_det(isfinite(tmp.data2_obs_det)));

                     if sum(squeeze(isfinite(tmp.data2_obs)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                        tmp.corr_det=NaN(2,2);
                        tmp.corr_det_p=NaN(2,2);
                        tmp.data2_lens2_trend=NaN;
                     end

                     data2.([tmp.varname, '_corr_obs_lens2_ano', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data2.([tmp.varname, '_corr_obs_lens2_ano_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
                     data2.([tmp.varname, '_corr_obs_lens2_det_ano', '_', tmp.season])(loni,lati)=tmp.corr_det(1,2);
                     data2.([tmp.varname, '_corr_obs_lens2_det_ano_p', '_', tmp.season])(loni,lati)=tmp.corr_det_p(1,2);


                 %% corr obs ~ hcst-lens2 (anomaly based)
                    tmp.data2_hcst_lens2 = squeeze(tmp.data2) - squeeze(tmp.data2_lens2);
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data2_hcst_lens2(isfinite(tmp.data2_obs)), tmp.data2_obs(isfinite(tmp.data2_obs)));
                     if sum(squeeze(isfinite(tmp.data2_obs)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                     end
                     data2.([tmp.varname, '_corr_obs_int_ano', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data2.([tmp.varname, '_corr_obs_int_ano_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
                 
                 %% corr obs-lens2 ~ hcst-lens2 (anomaly based)
                    tmp.data2_hcst_lens2 = squeeze(tmp.data2) - squeeze(tmp.data2_lens2);
                    tmp.data2_obs_lens2 = squeeze(tmp.data2_obs) - squeeze(tmp.data2_lens2);
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data2_hcst_lens2(isfinite(tmp.data2_obs_lens2)), tmp.data2_obs_lens2(isfinite(tmp.data2_obs_lens2)));
                     if sum(squeeze(isfinite(tmp.data2_obs_lens2)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                     end
                     data2.([tmp.varname, '_corr_obs_int2_ano', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data2.([tmp.varname, '_corr_obs_int2_ano_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);

                 %% hcst-svd 1st mode ~ obs (anomaly based)
                    tmp.data2_hcst_svd = squeeze(data2.([tmp.varname, '_model_ano', '_', tmp.season])(loni,lati,:)) - ...
                                    squeeze(data2.([tmp.varname, '_svd_lens2_model_right_1mode_recon', '_', tmp.season])(loni,lati,:));
                    [tmp.corr_int2, tmp.corr_int2_p]=corrcoef(tmp.data2_hcst_svd(isfinite(tmp.data2_obs)), tmp.data2_obs(isfinite(tmp.data2_obs)));
                    data2.([tmp.varname, '_corr_obs_int_svd', '_', tmp.season])(loni,lati)=tmp.corr_int2(1,2);
                    data2.([tmp.varname, '_corr_obs_int_svd_p', '_', tmp.season])(loni,lati)=tmp.corr_int2_p(1,2);

                %% hcst-svd 1st mode ~ obs-svd 1st mode (anomaly based)
                    tmp.data2_hcst_svd = squeeze(data2.([tmp.varname, '_model_ano', '_', tmp.season])(loni,lati,:)) - ...
                                    squeeze(data2.([tmp.varname, '_svd_lens2_model_right_1mode_recon', '_', tmp.season])(loni,lati,:));
                    tmp.data2_obs_svd = squeeze(data2.([tmp.varname, '_obs_ano', '_', tmp.season])(loni,lati,:)) - ...
                                    squeeze(data2.([tmp.varname, '_svd_lens2_obs_right_1mode_recon', '_', tmp.season])(loni,lati,:));
                    [tmp.corr_int2, tmp.corr_int2_p]=corrcoef(tmp.data2_hcst_svd(isfinite(tmp.data2_obs_svd)), tmp.data2_obs_svd(isfinite(tmp.data2_obs_svd)));
                    data2.([tmp.varname, '_corr_obs_int2_svd', '_', tmp.season])(loni,lati)=tmp.corr_int2(1,2);
                    data2.([tmp.varname, '_corr_obs_int2_svd_p', '_', tmp.season])(loni,lati)=tmp.corr_int2_p(1,2);

                %% corr AR1(obs) ~ obs
                    tmp.obs_AR1 = squeeze(data2.([tmp.varname, '_obs_AR1', '_', tmp.season])(loni,lati,:));
                    tmp.obs_AR1_mask=squeeze(tmp.obs_AR1./tmp.obs_AR1);
                    tmp.data2_obs_mask=squeeze(tmp.data2_obs./tmp.data2_obs);
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.obs_AR1.*tmp.data2_obs_mask.*tmp.obs_AR1_mask, ...
                         squeeze(tmp.data2_obs).*tmp.data2_obs_mask.*tmp.obs_AR1_mask, 'Rows', 'complete');

                     data2.([tmp.varname, '_corr_obs_AR1', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data2.([tmp.varname, '_corr_obs_AR1_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
    
                 else
                    %% obs (NaN)
                     data2.([tmp.varname, '_corr_obs_ano', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_ano_p', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_det_ano', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_det_ano_p', '_', tmp.season])(loni,lati)=NaN;
    
                    %% AR1 (NaN)
                     data2.([tmp.varname, '_corr_assm_AR1', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_AR1_p', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_AR1', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_AR1_p', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_AR1_ano', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_AR1_ano_p', '_', tmp.season])(loni,lati)=NaN;
    
                    %% ASSM (NaN)
                     data2.([tmp.varname, '_corr_assm_ano', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_ano_p', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_det_ano', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_det_ano_p', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_int_ano', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_int_ano_p', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_int2_ano', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_int2_ano_p', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_int_svd', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_int_svd_p', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_int2_svd', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_int2_svd_p', '_', tmp.season])(loni,lati)=NaN;
    
                    %% LENS2 (NaN)
                     data2.([tmp.varname, '_corr_assm_lens2_ano', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_lens2_ano_p', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_lens2_det_ano', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_assm_lens2_det_ano_p', '_', tmp.season])(loni,lati)=NaN;
                     %% obs ~ ASSM (NaN)
                     data2.([tmp.varname, '_corr_obs_assm_ano', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_assm_ano_p', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_assm_det_ano', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_assm_det_ano_p', '_', tmp.season])(loni,lati)=NaN;
                     
                     %% obs ~ LENS2 (NaN)
                     data2.([tmp.varname, '_corr_obs_lens2_ano', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_lens2_ano_p', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_lens2_det_ano', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_lens2_det_ano_p', '_', tmp.season])(loni,lati)=NaN;

                     %% obs ~ HCST-LENS2 (NaN)
                     data2.([tmp.varname, '_corr_obs_int_ano', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_int_ano_p', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_int2_ano', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_int2_ano_p', '_', tmp.season])(loni,lati)=NaN;

                     %% obs ~ HCST-svd (NaN)
                     data2.([tmp.varname, '_corr_obs_int_svd', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_int_svd_p', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_int2_svd', '_', tmp.season])(loni,lati)=NaN;
                     data2.([tmp.varname, '_corr_obs_int2_svd_p', '_', tmp.season])(loni,lati)=NaN;
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
            '_y', num2str(cfg.min_iy, '%04i'), '_y', num2str(cfg.max_iy, '%04i'), ...
            '_', tmp.season, '.mat'];
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
        case 'COL_FIRE_CLOSS'
            obsname_simple='GFED'; % GFED
        case 'photoC_TOT_zint'
            obsname_simple='VGPM'; % VGPM
        case 'photoC_TOT_zint_100m'
            obsname_simple='VGPM'; % VGPM
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
        case 'COL_FIRE_CLOSS'
            obsname_simple='FIRE_CLOSS_reg_cesm2.'; % GFED
        case 'photoC_TOT_zint'
            obsname_simple='s_vgpm_reg_cesm2.'; % VGPM
        case 'photoC_TOT_zint_100m'
            obsname_simple='s_vgpm_reg_cesm2.'; % VGPM
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
            obsname_simple='burned_area'; % MODIS Fire_cci v5.1
        case 'COL_FIRE_CLOSS'
            obsname_simple='C'; % GFED
        case 'photoC_TOT_zint'
            obsname_simple='npp'; % VGPM
        case 'photoC_TOT_zint_100m'
            obsname_simple='npp'; % VGPM
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
            obsname_simple=2001:2020; % MODIS Fire_cci v5.1
        case 'COL_FIRE_CLOSS'
            obsname_simple=1997:2020; % GFED
        case 'photoC_TOT_zint'
            obsname_simple=2003:2020; % VGPM
        case 'photoC_TOT_zint_100m'
            obsname_simple=2003:2020; % VGPM
        otherwise
            obsname_simple=1960:2020;
    end
end

function mons = f_season_mons(season)
    switch season
        case 'INI'
            mons = [1];
        case 'FMA'
            mons = [2,3,4];
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
        case 'JFM1'
            mons = [1,2,3];
        case 'MAM1'
            mons = [3,4,5];
        case 'AMJ1'
            mons = [4,5,6]; 
        case 'JJA1'
            mons = [6,7,8];
        case 'JAS1'
            mons = [7,8,9];
        case 'SON1'
            mons = [9,10,11];
        case 'OND1'
            mons = [10,11,12];
        case 'DJF1'
            mons = [12,13,14];
        case 'MAM2'
            mons = [15,16,17];
        case 'AMJ2'
            mons = [16,17,18];
        case 'JJA2'
            mons = [18,19,20];
        case 'JAS2'
            mons = [19,20,21];
        case 'SON2'
            mons = [21,22,23];
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