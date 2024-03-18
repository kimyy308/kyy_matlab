% %  Created 11-May-2023 by Yong-Yub Kim
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

%% model configuration
% cfg.var='TS'; %SST PRECT PSL TS SSH sumChl
% cfg.vars = {'SSH', 'sumChl'};
cfg.vars = {'SST', 'PRECT', 'PSL', 'TS', 'SSH', 'BSF'};
% cfg.vars = {'SSH'};
% cfg.vars = {'SST', 'PRECT', 'PSL', 'TS', 'sumChl', 'photoC_TOT_zint'};
% cfg.vars = { 'sumChl', 'photoC_TOT_zint',  'SSH'};
% cfg.vars = {  'SST', 'PRECT', 'PSL'};
% % 
% cfg.vars = { 'PAR_avg', 'SiO3', 'Fe', 'PO4', 'NO3'};
cfg.vars = {'SiO3', 'Fe', 'PO4', 'NO3', 'diatC', 'diazC', 'spC', 'SSH', 'BSF'};

% cfg.vars = {'diatChl', 'diazChl', 'spChl'};

% cfg.vars = {'diatChl'};
% cfg.vars = {'spChl'};
% cfg.vars={'WVEL2};
cfg.vars = {'SST'};
cfg.vars = {'diatChl', 'diazChl', 'spChl'};

cfg.vars={'diat_agg_zint_100m', 'diat_Fe_lim_Cweight_avg_100m', ...
    'diat_Fe_lim_surf', 'diat_light_lim_Cweight_avg_100m', 'diat_light_lim_surf', ...
    'diat_loss_zint_100m', 'diat_N_lim_Cweight_avg_100m', 'diat_N_lim_surf', ...
    'diat_P_lim_Cweight_avg_100m', 'diat_P_lim_surf', 'diat_SiO3_lim_Cweight_avg_100m', ...
    'diat_SiO3_lim_surf', 'diaz_agg_zint_100m', 'diaz_Fe_lim_Cweight_avg_100m', ...
    'diaz_Fe_lim_surf', 'diaz_light_lim_Cweight_avg_100m', 'diaz_light_lim_surf', ...
    'diaz_loss_zint_100m', 'diaz_P_lim_Cweight_avg_100m', 'diaz_P_lim_surf', 'dustToSed', ...
    'graze_diat_zint_100m', 'graze_diat_zoo_zint_100m', 'graze_diaz_zint_100m', ...
    'graze_diaz_zoo_zint_100m', 'graze_sp_zint_100m', 'graze_sp_zoo_zint_100m', ...
    'LWUP_F', 'O2_ZMIN_DEPTH', 'photoC_diat_zint_100m', 'photoC_diaz_zint_100m', ...
    'photoC_NO3_TOT_zint_100m', 'photoC_sp_zint_100m', ...
    'sp_agg_zint_100m', 'sp_Fe_lim_Cweight_avg_100m', 'sp_Fe_lim_surf', ...
    'sp_light_lim_Cweight_avg_100m', 'sp_light_lim_surf', 'sp_loss_zint_100m', ...
    'sp_N_lim_Cweight_avg_100m', 'sp_N_lim_surf', 'sp_P_lim_Cweight_avg_100m', ...
    'sp_P_lim_surf', 'zoo_loss_zint_100m, photoC_TOT_zint_100m'};

cfg.vars = {'sutmChl'};
cfg.vars = {'TS'};

% cfg.vlayer=1; % surf 
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
%     cfg.proj_year=5;
    cfg.season = {'AMJ', 'JAS', 'OND', 'JFM'};
    
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
        case 'ocn'
            grid.tlong=ncread(tmp.gridname, 'TLONG');
            grid.tlat=ncread(tmp.gridname, 'TLAT');
            grid.mask_ocn=ncread(tmp.maskname, 'open_ocean');
            grid.mask_ocn(grid.mask_ocn<-10e10)=NaN;
            grid.mask_ocn=grid.mask_ocn./grid.mask_ocn;
            grid.z_t=ncread(tmp.gridname, 'z_t')./100; % meter
            grid.dz=ncread(tmp.gridname, 'dz')./100; % meter
            grid.dz_res=reshape(grid.dz(cfg.vlayer), [1 1 length(grid.dz(cfg.vlayer))]);
        case 'atm'
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

    clear tmp.ydata tmp.ydata_lens2 tmp.ydata_obs tmp.ydata_assm

    for lss=1:length(cfg.season)
        tmp.season=cfg.season{lss};
        tmp.mons = f_season_mons(tmp.season);
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
        data.([tmp.varname, '_AR1', '_', tmp.season])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);  %% AR1 initialization
    
    %     data.([tmp.varname, '_assm'])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        
        %% read variables
        for iyear=min(cfg.iyears):max(cfg.iyears)
            tmp.iyear_str=num2str(iyear, '%04i');
            tmp.iind=iyear-min(cfg.iyears)+1;
            cfg.casename=['ensmean_', cfg.casename_m, '_i', tmp.iyear_str];
%             tmp.fy=iyear+lyear;
%             tmp.fy_str=num2str(tmp.fy, '%04i');
            %% monthly filename
            for mon=1:length(tmp.mons)
                tmp.mon=tmp.mons(mon);
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
                tmpvarid=netcdf.inqVarID(ncid,tmp.varname);
                [tmp.varname,tmp.xtype,tmp.dimids,tmp.natts]= netcdf.inqVar(ncid,tmpvarid);
                if length(tmp.dimids)>3
                     tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 0], [grid.nlon grid.nlat cfg.vlayer_cnt 1]));
%                      tmp.dd=tmp.dd.*grid.mask_ocn;
                     tmp.dd(abs(tmp.dd)>1e30)=NaN;
                     ddd=mean(tmp.dd,3,'omitnan'); % depth mean
                     tmp.ydata(1:grid.nlon,1:grid.nlat,mon) = ddd; %depth averaged value
                else
                     tmp.dd= netcdf.getVar(ncid,tmpvarid);
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
                    tmpvarid=netcdf.inqVarID(ncid,tmp.varname);
                    [tmp.varname,tmp.xtype,tmp.dimids,tmp.natts]= netcdf.inqVar(ncid,tmpvarid);
                    if length(tmp.dimids)>3
                         tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 0], [grid.nlon grid.nlat cfg.vlayer_cnt 1]));
    %                      tmp.dd=tmp.dd.*grid.mask_ocn;
                         tmp.dd(abs(tmp.dd)>1e30)=NaN;
                         ddd=mean(tmp.dd,3,'omitnan'); % depth mean
                         tmp.ydata(1:grid.nlon,1:grid.nlat,mon) = ddd; %depth averaged value
                    else
                         tmp.dd= netcdf.getVar(ncid,tmpvarid);
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
                    tmpvarid=netcdf.inqVarID(ncid,tmp.varname);
                    if length(tmp.dimids)>3
                         tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 0], [grid.nlon grid.nlat cfg.vlayer_cnt 1]));
%                          ddd=tmp.dd.*grid.dz_res; %% weight depth
%                          ddd=sum(ddd,3); % depth sum
%                          ddd=ddd./sum(grid.dz_res(:)); % divide total depth for avg
                         tmp.dd(abs(tmp.dd)>1e30)=NaN;
                         ddd=mean(tmp.dd,3,'omitnan'); % depth mean
                         tmp.ydata_lens2(1:grid.nlon,1:grid.nlat,mon) = ddd; %depth averaged value
                    else
                        tmp.ydata_lens2(1:grid.nlon,1:grid.nlat,mon) = netcdf.getVar(ncid,tmpvarid);
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
                        tmp.dd= netcdf.getVar(ncid,tmpvarid);
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
                    tmpvarid=netcdf.inqVarID(ncid,tmp.varname);
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
    %         data.([tmp.varname, '_model', '_', tmp.season])(:,:,iyear-min(cfg.iyears)+1+lyear)= tmp.ymean;
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
%             data.([tmp.varname, '_model', '_', tmp.season])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean;
%             data.([tmp.varname, '_obs'])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_obs;
%             data.([tmp.varname, '_lens2', '_', tmp.season])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_lens2;
%             data.([tmp.varname, '_assm'])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_assm;
            data.([tmp.varname, '_model', '_', tmp.season])(1:grid.nlon,1:grid.nlat,tmp.iind)= tmp.ymean;
            if isfield(tmp, 'ydata_mod_obs_masked')
                data.([tmp.varname, '_mod_obs_masked', '_', tmp.season])(1:grid.nlon,1:grid.nlat,tmp.iind)= tmp.ymean_mod_obs_masked;
                data.([tmp.varname, '_assm_obs_masked', '_', tmp.season])(1:grid.nlon,1:grid.nlat,tmp.iind)= tmp.ymean_assm_obs_masked;
                data.([tmp.varname, '_lens2_obs_masked', '_', tmp.season])(1:grid.nlon,1:grid.nlat,tmp.iind)= tmp.ymean_lens2_obs_masked;
            end
            data.([tmp.varname, '_obs', '_', tmp.season])(1:grid.nlon,1:grid.nlat,tmp.iind)= tmp.ymean_obs;
            data.([tmp.varname, '_assm', '_', tmp.season])(1:grid.nlon,1:grid.nlat,tmp.iind)= tmp.ymean_assm;
            data.([tmp.varname, '_lens2', '_', tmp.season])(1:grid.nlon,1:grid.nlat,tmp.iind)= tmp.ymean_lens2;
    %         %% AR1 (assign)
    %         tmp.ymean_AR1= mean(tmp.ydata_AR1, 3);
    %         data.([tmp.varname, '_AR1', '_', tmp.season])(:,:,iyear-min(cfg.iyears)+1+lyear)= tmp.ymean_AR1;
    %         tmp.ymean= mean(ncread(cfg.datafilename, ['assm_', tmp.varname], [1, 1, 1], [inf, inf, 12]), 3);
    %         data.([tmp.varname, '_assm'])(:,:,iyear-min(cfg.iyears)+1)= tmp.ymean;
        end
        
        %% get correlation coefficient
        for loni=1:grid.nlon
            for lati=1:grid.nlat
                 if (isnan(data.([tmp.varname, '_assm', '_', tmp.season])(loni,lati,1))~=1 & nansum(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,:))~=0)

                 %% corr assm
                     tmp.data = squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,:));
                     tmp.data_assm = squeeze(data.([tmp.varname, '_assm', '_', tmp.season])(loni,lati,:));
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_assm)), tmp.data_assm(isfinite(tmp.data_assm)));
                     
                 %% hcst-lens2 ~ assm
                     tmp.data_hcst_lens2 = squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,:)) - ...
                                    squeeze(data.([tmp.varname, '_lens2', '_', tmp.season])(loni,lati,:));
                     [tmp.corr_int, tmp.corr_int_p]=corrcoef(tmp.data_hcst_lens2(isfinite(tmp.data_assm)), tmp.data_assm(isfinite(tmp.data_assm)));
    
                     data.([tmp.varname, '_corr_assm', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_assm_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
                     data.([tmp.varname, '_corr_assm_int', '_', tmp.season])(loni,lati)=tmp.corr_int(1,2);
                     data.([tmp.varname, '_corr_assm_int_p', '_', tmp.season])(loni,lati)=tmp.corr_int_p(1,2);
    
                 %% corr lens2 ~ assm
                     tmp.lens2 = squeeze(data.([tmp.varname, '_lens2', '_', tmp.season])(loni,lati,:));
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.lens2(isfinite(tmp.data_assm)), tmp.data_assm(isfinite(tmp.data_assm)));
    
                     data.([tmp.varname, '_corr_assm_lens2', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_assm_lens2_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
    
                
                 %% corr obs ~ hcst
                     if isfield(tmp, 'ydata_mod_obs_masked')
                        tmp.data = squeeze(data.([tmp.varname, '_mod_obs_masked', '_', tmp.season])(loni,lati,:));
                     else
                        tmp.data = squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,:));
                     end
                     tmp.data_obs = squeeze(data.([tmp.varname, '_obs', '_', tmp.season])(loni,lati,:));
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
                     
                     if sum(squeeze(isfinite(tmp.data_obs)))<10
                         tmp.data_det = NaN(size(tmp.data));
                         tmp.data_obs_det = NaN(size(tmp.data));
                     else
                         tmp.data_det = Func_0028_detrend_linear_1d(tmp.data(isfinite(tmp.data_obs)));
                         tmp.data_obs_det = Func_0028_detrend_linear_1d(tmp.data_obs(isfinite(tmp.data_obs)));
                     end
                     [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data_det(isfinite(tmp.data_obs_det)), tmp.data_obs_det(isfinite(tmp.data_obs_det)));
    %                  [tmp.corr, tmp.corr_p]=corrcoef(squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,1+lyear:end-cfg.proj_year+1)), ...
    %                      squeeze(data.([tmp.varname, '_obs'])(loni,lati,1+lyear:end)));

                     if sum(squeeze(isfinite(tmp.data_obs)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                        tmp.corr_det=NaN(2,2);
                        tmp.corr_det_p=NaN(2,2);
                     end

                     data.([tmp.varname, '_corr_obs', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_obs_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
                     data.([tmp.varname, '_corr_obs_det', '_', tmp.season])(loni,lati)=tmp.corr_det(1,2);
                     data.([tmp.varname, '_corr_obs_det_p', '_', tmp.season])(loni,lati)=tmp.corr_det_p(1,2);

                 %% corr obs ~ assm
                     if isfield(tmp, 'ydata_mod_obs_masked')
                        tmp.data_assm = squeeze(data.([tmp.varname, '_assm_obs_masked', '_', tmp.season])(loni,lati,:));
                     else
                        tmp.data_assm = squeeze(data.([tmp.varname, '_assm', '_', tmp.season])(loni,lati,:));
                     end
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data_assm(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
                     
                     if sum(squeeze(isfinite(tmp.data_obs)))<10
                         tmp.data_assm_det = NaN(size(tmp.data));
                     else
                         tmp.data_assm_det = Func_0028_detrend_linear_1d(tmp.data_assm(isfinite(tmp.data_obs)));
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
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                     end

                     data.([tmp.varname, '_corr_obs_lens2', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_obs_lens2_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);

                 %% corr obs ~ hcst-lens2
                     tmp.data_hcst_lens2 = squeeze(tmp.data) - squeeze(tmp.data_lens2);
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data_hcst_lens2(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
                     if sum(squeeze(isfinite(tmp.data_obs)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                     end
                     data.([tmp.varname, '_corr_obs_int', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_obs_int_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);


    
    %                 %% AR1 (corr)
    %                  tmp.data_AR1 = squeeze(data.([tmp.varname, '_AR1', '_', tmp.season])(loni,lati,1+lyear:end-cfg.proj_year+1));
    %                  [tmp.corr_AR1, tmp.corr_AR1_p]=corrcoef(tmp.data_AR1, tmp.data_obs);
    %                  if (isnan(tmp.corr_AR1(1,2)))
    %                     tmp.corr_AR1(1,2)=0;
    %                  end
    %                  data.([tmp.varname, '_corr_AR1', '_', tmp.season])(loni,lati)=tmp.corr_AR1(1,2);
    %                  data.([tmp.varname, '_corr_AR1_p', '_', tmp.season])(loni,lati)=tmp.corr_AR1_p(1,2);
    %                  
    %                  tmp.data_AR1_det = Func_0028_detrend_linear_1d(data.([tmp.varname, '_AR1', '_', tmp.season])(loni,lati,1+lyear:end-cfg.proj_year+1));
    %                  [tmp.corr_AR1_det, tmp.corr_AR1_det_p]=corrcoef(tmp.data_AR1_det, tmp.data_obs_det);
    %                  if (isnan(tmp.corr_AR1_det(1,2)))
    %                     tmp.corr_AR1_det(1,2)=0;
    %                  end
    %                  data.([tmp.varname, '_corr_AR1_det', '_', tmp.season])(loni,lati)=tmp.corr_AR1_det(1,2);
    %                  data.([tmp.varname, '_corr_AR1_det_p', '_', tmp.season])(loni,lati)=tmp.corr_AR1_det_p(1,2);
    
                 else
                     data.([tmp.varname, '_corr_obs', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_p', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_det', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_det_p', '_', tmp.season])(loni,lati)=NaN;
    %                  data.([tmp.varname, '_corr_assm', '_', tmp.season])(loni,lati)=NaN;
    %                  data.([tmp.varname, '_corr_assm_p', '_', tmp.season])(loni,lati)=NaN;
    
    %                 %% AR1 (NaN)
    %                  data.([tmp.varname, '_corr_AR1', '_', tmp.season])(loni,lati)=NaN;
    %                  data.([tmp.varname, '_corr_AR1_p', '_', tmp.season])(loni,lati)=NaN;
    %                  data.([tmp.varname, '_corr_AR1_det', '_', tmp.season])(loni,lati)=NaN;
    %                  data.([tmp.varname, '_corr_AR1_det_p', '_', tmp.season])(loni,lati)=NaN;
    
                    %% ASSM (NaN)
                     data.([tmp.varname, '_corr_assm', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_p', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_int', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_int_p', '_', tmp.season])(loni,lati)=NaN;
    
                     %% LENS2 (NaN)
                     data.([tmp.varname, '_corr_assm_lens2', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_lens2_p', '_', tmp.season])(loni,lati)=NaN;

                     %% obs ~ ASSM (NaN)
                     data.([tmp.varname, '_corr_obs_assm', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_assm_p', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_assm_det', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_assm_det_p', '_', tmp.season])(loni,lati)=NaN;
                     
                     %% obs ~ LENS2 (NaN)
                     data.([tmp.varname, '_corr_obs_lens2', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_lens2_p', '_', tmp.season])(loni,lati)=NaN;

                     %% obs ~ HCST-LENS2 (NaN)
                     data.([tmp.varname, '_corr_obs_int', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_int_p', '_', tmp.season])(loni,lati)=NaN;
                    
                 end
            end
        end
        disp('abc')
        if isfield(tmp,'ydata_mod_obs_masked')
            tmp=rmfield(tmp, 'ydata_mod_obs_masked');
        end

        
        if ~exist(dirs.matroot,'dir'), mkdir(dirs.matroot); end
        fig_cfg.mat_name=[dirs.matroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_', tmp.season, '.mat'];
        save(fig_cfg.mat_name, 'data')
        clear data
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