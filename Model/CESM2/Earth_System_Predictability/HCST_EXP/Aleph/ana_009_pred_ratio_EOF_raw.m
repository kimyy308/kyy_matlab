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
cfg.vars={'SST', 'PRECT'};
cfg.vars={'photoC_TOT_zint_100m'};
% tmp.dimids= [1, 2, 4];
cfg.vlayer=1; % surf, vertical slice 
% cfg.vlayer=1:10; % 10layer. don't put more than 15
cfg.vlayer=1; % 100m, vertical slice 
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
    cfg.obs_iyears=1960:2016;
    cfg.obs_iyears2=f_obs_iyears(cfg.var);
  
    % dirs.hcstroot=['/mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive/', cfg.comp, '/', cfg.var];
    % dirs.obsroot=['/Volumes/kyy_raid/kimyy/Observation/', cfg.obs_name, '/monthly_reg_', cfg.obs_fname_module(2:4)];
    % dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', cfg.var];
    % dirs.lens2root=['/mnt/lustre/proj/earth.system.predictability/LENS2/archive_analysis/', cfg.comp, '/', cfg.var];
    % dirs.assmroot=['/mnt/lustre/proj/earth.system.predictability/ASSM_EXP/archive_analysis/', cfg.comp, '/', cfg.var];
    
    dirs.hcstroot=['/mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive_transfer/', cfg.comp, '/', cfg.var];
    dirs.obsroot=['/mnt/lustre/proj/kimyy/Observation/', cfg.obs_name];
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
    for lmonth=1:cfg.proj_year*12
        tmp.lmonth_str=num2str(lmonth, '%02i');
        cfg.casename_m=['all'];
    
        dirs.datadir= [dirs.hcstroot, filesep, 'ens_all', filesep];
        dirs.assmdir= [dirs.assmroot, filesep, 'ens_all', filesep];
        dirs.lens2dir= [dirs.lens2root, filesep, 'ens_smbb', filesep];
        fprintf('%02d_%s_%s  ',lmonth, ',', cfg.casename_m,'_', tmp.varname); lap_time = tic;
    
        
%% variables initialization
        data.([tmp.varname, '_model', '_l', tmp.lmonth_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data.([tmp.varname, '_model_inc', '_l', tmp.lmonth_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data.([tmp.varname, '_model_ratio', '_l', tmp.lmonth_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
% % % % %         data.([tmp.varname, '_model_rmse', '_l', tmp.lmonth_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
% % % % %         data.([tmp.varname, '_model_rmse_ratio', '_l', tmp.lmonth_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);

% % % % %         if lmonth<13
% % % % %             data.([tmp.varname, '_lens2', '_l', tmp.lmonth_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
% % % % %             data.([tmp.varname, '_assm', '_l', tmp.lmonth_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
% % % % %             data.([tmp.varname, '_lens2_inc', '_l', tmp.lmonth_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
% % % % %             data.([tmp.varname, '_assm_inc', '_l', tmp.lmonth_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
% % % % %         end

    %     data.([tmp.varname, '_assm'])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        
%% read variables
        for iyear=min(cfg.iyears):max(cfg.iyears)
            tmp.iyear_str=num2str(iyear, '%04i');
            cfg.casename=['ensmean_', cfg.casename_m, '_i', tmp.iyear_str];
            tmp.fy=iyear+floor((lmonth-1)/12);
            tmp.fy_str=num2str(tmp.fy, '%04i');
            tmp.fm=lmonth - floor((lmonth-1)/12)*12;
            tmp.fm_str = num2str(tmp.fm, '%02i');

%% HCST mean
            cfg.mod_fnm=[dirs.datadir, tmp.fs, 'ens_all_i',tmp.iyear_str, tmp.fs, ...
            tmp.varname, '_', cfg.gnm, '.hcst.', cfg.casename, cfg.obs_fname_module, tmp.fy_str, '-', tmp.fm_str, '.nc'];
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
                tmp.varname, '_', cfg.gnm, '.hcst.', cfg.casename_std, cfg.obs_fname_module, tmp.fy_str, '-', tmp.fm_str, '.nc'];
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

% % % % % %% LENS2 mean
% % % % %             if tmp.fy <= cfg.max_iy
% % % % %                 cfg.lens2_fnm=[dirs.lens2dir, tmp.fs, ...
% % % % %                     tmp.varname, '_y_ensmean_', tmp.fy_str, '.nc'];
% % % % %                 ncid=netcdf.open(cfg.lens2_fnm, 'NOWRITE');
% % % % %                 tmpvarid=netcdf.inqVarID(ncid,tmp.fvarname);
% % % % %                 if length(tmp.dimids)>3
% % % % %                      tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 0], [grid.nlon grid.nlat cfg.vlayer_cnt 1]));
% % % % %                      tmp.dd(abs(tmp.dd)>1e30)=NaN;
% % % % %                      ddd=mean(tmp.dd,3,'omitnan'); % depth mean
% % % % %                      tmp.ymean_lens2(1:grid.nlon,1:grid.nlat) = ddd; %depth averaged value
% % % % %                 else
% % % % %                     tmp.dd=netcdf.getVar(ncid,tmpvarid);
% % % % %                      tmp.dd(abs(tmp.dd)>1e30)=NaN;                        
% % % % %                     tmp.ymean_lens2(1:grid.nlon,1:grid.nlat) =tmp.dd;
% % % % %                 end
% % % % %                 netcdf.close(ncid);
% % % % %             else
% % % % %                 tmp.ymean_lens2(1:grid.nlon,1:grid.nlat) = NaN;                    
% % % % %             end
% % % % % 
% % % % % %% LENS2 std
% % % % %             if tmp.fy <= cfg.max_iy
% % % % %                 cfg.lens2_fnm=[dirs.lens2dir, tmp.fs, ...
% % % % %                     tmp.varname, '_y_ensstd_', tmp.fy_str, '.nc'];
% % % % %                 if exist(cfg.lens2_fnm)~=0
% % % % %                     ncid=netcdf.open(cfg.lens2_fnm, 'NOWRITE');
% % % % %                     tmpvarid=netcdf.inqVarID(ncid,tmp.fvarname);
% % % % %                     if length(tmp.dimids)>3
% % % % %                          tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 0], [grid.nlon grid.nlat cfg.vlayer_cnt 1]));
% % % % %                          tmp.dd(abs(tmp.dd)>1e30)=NaN;
% % % % %                          ddd=mean(tmp.dd,3,'omitnan'); % depth mean
% % % % %                          tmp.ystde_lens2(1:grid.nlon,1:grid.nlat) = ddd; %depth averaged value
% % % % %                     else
% % % % %                         tmp.dd=netcdf.getVar(ncid,tmpvarid);
% % % % %                          tmp.dd(abs(tmp.dd)>1e30)=NaN;                        
% % % % %                         tmp.ystde_lens2(1:grid.nlon,1:grid.nlat) =tmp.dd;
% % % % %                     end
% % % % %                     netcdf.close(ncid);
% % % % %                 else
% % % % %                     tmp.ystde_lens2(1:grid.nlon,1:grid.nlat)=0;
% % % % %                 end
% % % % %             else
% % % % %                 tmp.ystde_lens2(1:grid.nlon,1:grid.nlat) = NaN;                    
% % % % %             end
% % % % % 
% % % % % %% ASSM
% % % % %             if tmp.fy <= cfg.max_iy
% % % % %                 cfg.assm_fnm=[dirs.assmdir, tmp.fs, tmp.varname, '_ensmean_', tmp.fy_str, '.nc'];
% % % % %                 ncid=netcdf.open(cfg.assm_fnm, 'NOWRITE');
% % % % %                 tmpvarid=netcdf.inqVarID(ncid,tmp.fvarname);
% % % % %                 if length(tmp.dimids)>3
% % % % %                     tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 0], [grid.nlon grid.nlat cfg.vlayer_cnt 1]));
% % % % % %                          ddd=tmp.dd.*grid.dz_res; %% weight depth
% % % % % %                          ddd=sum(ddd,3); % depth sum
% % % % % %                          ddd=ddd./sum(grid.dz_res(:)); % divide total depth for avg
% % % % %                      tmp.dd(abs(tmp.dd)>1e30)=NaN;
% % % % %                      ddd=mean(tmp.dd,3,'omitnan'); % depth mean
% % % % %                      tmp.ymean_assm(1:grid.nlon,1:grid.nlat) = ddd; %depth averaged value             
% % % % %                 else
% % % % %                     tmp.ymean_assm(1:grid.nlon,1:grid.nlat) = netcdf.getVar(ncid,tmpvarid);
% % % % %                 end
% % % % %                 if(strcmp(tmp.varname(1:2),'mu')~=1 && strcmp(tmp.varname(1:2),'su')~=1)
% % % % %                     data.units=netcdf.getAtt(ncid,tmpvarid,'units');
% % % % %                 else
% % % % %                     data.units=' ';
% % % % %                 end
% % % % %                 netcdf.close(ncid);
% % % % %             else
% % % % %                 tmp.ymean_assm(1:grid.nlon,1:grid.nlat) = NaN;                    
% % % % %             end
% % % % % 
% % % % % %% ASSM std
% % % % %             if tmp.fy <= cfg.max_iy
% % % % %                 cfg.assm_fnm=[dirs.assmdir, tmp.fs, tmp.varname, '_ensstd_', tmp.fy_str, '.nc'];
% % % % %                 if exist(cfg.assm_fnm)~=0
% % % % %                     ncid=netcdf.open(cfg.assm_fnm, 'NOWRITE');
% % % % %                     tmpvarid=netcdf.inqVarID(ncid,tmp.fvarname);
% % % % %                     if length(tmp.dimids)>3
% % % % %                         tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 0], [grid.nlon grid.nlat cfg.vlayer_cnt 1]));
% % % % %     %                          ddd=tmp.dd.*grid.dz_res; %% weight depth
% % % % %     %                          ddd=sum(ddd,3); % depth sum
% % % % %     %                          ddd=ddd./sum(grid.dz_res(:)); % divide total depth for avg
% % % % %                          tmp.dd(abs(tmp.dd)>1e30)=NaN;
% % % % %                          ddd=mean(tmp.dd,3,'omitnan'); % depth mean
% % % % %                          tmp.ystde_assm(1:grid.nlon,1:grid.nlat) = ddd; %depth averaged value             
% % % % %                     else
% % % % %                         tmp.ystde_assm(1:grid.nlon,1:grid.nlat) = netcdf.getVar(ncid,tmpvarid);
% % % % %                     end
% % % % %                     netcdf.close(ncid);
% % % % %                 else
% % % % %                     tmp.ystde_assm(1:grid.nlon,1:grid.nlat)=0;
% % % % %                 end
% % % % %             else
% % % % %                 tmp.ystde_assm(1:grid.nlon,1:grid.nlat) = NaN;                    
% % % % %             end
           
% % % % %             tmp.ymean_assm(tmp.ymean_assm>10e30)=NaN;
% % % % %             tmp.ystde_assm(tmp.ystde_assm>10e30)=NaN;

            switch cfg.comp
                case 'ocn'
                    tmp.ymean = tmp.ymean .* grid.mask_ocn;
                    tmp.ystde = tmp.ystde .* grid.mask_ocn;
                    if isfield(tmp, 'ymean_mod_obs_masked')
                        tmp.ymean_mod_obs_masked = tmp.ymean_mod_obs_masked .* grid.mask_ocn;
                        tmp.ystde_mod_obs_masked = tmp.ystde_mod_obs_masked .* grid.mask_ocn;
% % % % %                         tmp.ymean_assm_obs_masked = tmp.ymean_assm_obs_masked .* grid.mask_ocn;
% % % % %                         tmp.ymean_lens2_obs_masked = tmp.ymean_lens2_obs_masked .* grid.mask_ocn;
% % % % %                         tmp.ystde_assm_obs_masked = tmp.ystde_assm_obs_masked .* grid.mask_ocn;
% % % % %                         tmp.ystde_lens2_obs_masked = tmp.ystde_lens2_obs_masked .* grid.mask_ocn;
                    end
% % % % %                     tmp.ymean_assm = tmp.ymean_assm .* grid.mask_ocn;
% % % % %                     tmp.ymean_lens2 = tmp.ymean_lens2 .* grid.mask_ocn;
% % % % %                     tmp.ystde_assm = tmp.ystde_assm .* grid.mask_ocn;
% % % % %                     tmp.ystde_lens2 = tmp.ystde_lens2 .* grid.mask_ocn;    
            end
    %         end

%% save variables as structure
            data.([tmp.varname, '_model', '_l', tmp.lmonth_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean;
            data.([tmp.varname, '_model_stde', '_l', tmp.lmonth_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ystde;
            
% % % % %             data.([tmp.varname, '_lens2', '_l', tmp.lmonth_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_lens2;
% % % % %             data.([tmp.varname, '_assm'])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_assm;
% % % % %             data.([tmp.varname, '_lens2_stde', '_l', tmp.lmonth_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ystde_lens2;
% % % % %             data.([tmp.varname, '_assm_stde'])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ystde_assm;
            
            ind_t=iyear-min(cfg.iyears)+1;
            data.([tmp.varname, '_model_inc', '_l', tmp.lmonth_str])(1:grid.nlon,1:grid.nlat,ind_t)= ...
                abs(tmp.ymean - data.([tmp.varname, '_model', '_l01'])(1:grid.nlon,1:grid.nlat,ind_t));
            data.([tmp.varname, '_model_ratio', '_l', tmp.lmonth_str])(1:grid.nlon,1:grid.nlat,ind_t)= ...
                data.([tmp.varname, '_model_inc', '_l', tmp.lmonth_str])(1:grid.nlon,1:grid.nlat,ind_t) ...
                ./ data.([tmp.varname, '_model_stde', '_l', tmp.lmonth_str])(1:grid.nlon,1:grid.nlat,ind_t);
% % % % %             % rmse
% % % % %             for lm2=1:lmonth
% % % % %                 tmp.lm2_str=num2str(lm2, '%02i');
% % % % %                 data.([tmp.varname, '_model_rmse', '_l', tmp.lmonth_str])(1:grid.nlon,1:grid.nlat,ind_t) = ...
% % % % %                     data.([tmp.varname, '_model_rmse', '_l', tmp.lmonth_str])(1:grid.nlon,1:grid.nlat,ind_t) + ...
% % % % %                 data.([tmp.varname, '_model', '_l', tmp.lmonth_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1) - ...
% % % % %                 data.([tmp.varname, '_model', '_l', tmp.lmonth_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)
% % % % %             
% % % % %             end
    
        end  %% initialized year loop end
        data.([tmp.varname, '_model_inc', '_tm_l', tmp.lmonth_str])= ...
            mean( data.([tmp.varname, '_model_inc', '_l', tmp.lmonth_str]), 3);
        data.([tmp.varname, '_model_ratio', '_tm_l', tmp.lmonth_str])= ...
            mean( data.([tmp.varname, '_model_ratio', '_l', tmp.lmonth_str]), 3);
        data.([tmp.varname, '_model_stde', '_tm_l', tmp.lmonth_str])= ...
            mean( data.([tmp.varname, '_model_stde', '_l', tmp.lmonth_str]), 3);

        disp(['lmonth: ', num2str(lmonth), ', data read finished'])
        fprintf('%7.1f sec\n', toc(lap_time) );
    

    end
    
    for lm=1:cfg.proj_year*12       
        tmp.lmonth_str=num2str(lm,'%02i');
        data_lm.ratio_for_EOF(:,:,lm)=data.([tmp.varname, '_model_ratio', '_tm_l', tmp.lmonth_str]);
        data_lm.inc_for_EOF(:,:,lm)=data.([tmp.varname, '_model_inc', '_tm_l', tmp.lmonth_str]);
        data_lm.stde_for_EOF(:,:,lm)=data.([tmp.varname, '_model_stde', '_tm_l', tmp.lmonth_str]);
    end
% % % % %     tmp.inc_stdt=std(data_lm.inc_for_EOF,0,3);

% % % % %     %% EOF
% % % % %     tmp.mode_want=3;
% % % % %     [EOF.lv, EOF.pc, EOF.var_exp] = ...
% % % % %         Func_0024_EOF_3d(tmp.ratio_for_EOF, tmp.mode_want);
% % % % % 
% % % % %     plot(EOF.pc(:,1));
% % % % %     pcolor(EOF.lv(:,:,1)'); shading flat; colorbar;
% % % % % 
% % % % %      plot(EOF.pc(:,2));
% % % % %     pcolor(EOF.lv(:,:,2)'); shading flat; colorbar;
% % % % % 
% % % % %     plot(EOF.pc(:,3));
% % % % %     pcolor(EOF.lv(:,:,3)'); shading flat; colorbar;
% % % % %     
% % % % %     
% % % % % 
% % % % %     %% valid grids only (value >=1 at 2nd time) % not effective
% % % % %     grid.ratio_valid_mask=tmp.ratio_for_EOF(:,:,2);
% % % % %     grid.ratio_valid_mask(grid.ratio_valid_mask<1)=NaN;
% % % % %     grid.ratio_valid_mask(grid.ratio_valid_mask>=1)=1;
% % % % %     pcolor(grid.ratio_valid_mask'); shading flat; colorbar;
% % % % %     tmp.mode_want=3;
% % % % %     [EOF.lv, EOF.pc, EOF.var_exp] = ...
% % % % %         Func_0024_EOF_3d(tmp.ratio_for_EOF.*grid.ratio_valid_mask, tmp.mode_want);
% % % % % 
% % % % %     plot(EOF.pc(:,1));
% % % % %     pcolor(EOF.lv(:,:,1)'); shading flat; colorbar;
% % % % % 
% % % % %      plot(EOF.pc(:,2));
% % % % %     pcolor(EOF.lv(:,:,2)'); shading flat; colorbar;
% % % % % 
% % % % %     plot(EOF.pc(:,3));
% % % % %     pcolor(EOF.lv(:,:,3)'); shading flat; colorbar;
% % % % % 
% % % % %     %% std(inc)_t < stde at first step = NaN
% % % % %     grid.ratio_valid_mask=zeros(size(tmp.inc_for_EOF(:,:,1)))+1;
% % % % %     grid.ratio_valid_mask(tmp.inc_stdt<data.([tmp.varname, '_model_stde', '_tm_l01']))=NaN;
% % % % %     pcolor(grid.ratio_valid_mask'); shading flat; colorbar;
% % % % %     tmp.mode_want=3;
% % % % %     [EOF.lv, EOF.pc, EOF.var_exp] = ...
% % % % %         Func_0024_EOF_3d(tmp.ratio_for_EOF.*grid.ratio_valid_mask, tmp.mode_want);
% % % % % 
% % % % %     plot(EOF.pc(:,1));
% % % % %     pcolor(EOF.lv(:,:,1)'); shading flat; colorbar;
% % % % % 
% % % % % 
% % % % %     %% movmean 12
% % % % %     lyt=0.5:1/12:4.5;
% % % % %     tmp.mode_want=3;
% % % % %     [EOF.lv, EOF.pc, EOF.var_exp] = ...
% % % % %         Func_0024_EOF_3d(movmean(tmp.ratio_for_EOF(:,85:105,:),12,3,'Endpoints', 'discard'), tmp.mode_want);
% % % % % 
% % % % %     figure(1); plot(lyt, EOF.pc(:,1));
% % % % %     figure(2); pcolor(EOF.lv(:,:,1)'); shading flat; colorbar;
% % % % % 
% % % % %     figure(1);plot(lyt, EOF.pc(:,2));
% % % % %     figure(2); pcolor(EOF.lv(:,:,2)'); shading flat; colorbar;
% % % % % 
% % % % %     figure(1); plot(lyt, EOF.pc(:,3));
% % % % %     figure(2); pcolor(EOF.lv(:,:,3)'); shading flat; colorbar;
% % % % % 
% % % % %     %% for stde
% % % % %      tmp.mode_want=3;
% % % % %     [EOF.lv, EOF.pc, EOF.var_exp] = ...
% % % % %         Func_0024_EOF_3d(movmean(tmp.stde_for_EOF(:,:,:),12,3,'Endpoints', 'discard'), tmp.mode_want);
% % % % % 
% % % % %     figure(1); plot(lyt, EOF.pc(:,1));
% % % % %     figure(2); pcolor(EOF.lv(:,:,1)'); shading flat; colorbar;
    
    %% for inc
     tmp.mode_want=3;
    [EOF.lv, EOF.pc, EOF.var_exp] = ...
        Func_0024_EOF_3d(movmean(data_lm.inc_for_EOF(:,:,:),12,3,'Endpoints', 'discard'), tmp.mode_want);

    figure(1); plot(lyt, EOF.pc(:,1));
    figure(2); pcolor(EOF.lv(:,:,1)'); shading flat; colorbar; caxis([-0.01 0])

    
%         if isfield(tmp,'ymean_mod_obs_masked')
%             tmp=rmfield(tmp, 'ymean_mod_obs_masked');
%         end
        if ~exist(dirs.matroot,'dir'), mkdir(dirs.matroot); end
        fig_cfg.mat_name=[dirs.matroot, filesep, 'pred_ratio_hcst_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l', tmp.lmonth_str, 'm.mat'];
        save(fig_cfg.mat_name, 'data_lm')
        clear data data_lm
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

