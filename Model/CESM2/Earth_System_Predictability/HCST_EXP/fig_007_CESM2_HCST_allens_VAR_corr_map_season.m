    % %  Created 30-Jan-2023 by Yong-Yub Kim
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

% cfg.vars= {'TS', 'SST', 'PRECT', 'PSL'};
% cfg.vars= {'SSH', 'SST', 'PRECT', 'PSL', 'TS'};
% 
% cfg.vars= {'SSH', 'sumChl', 'photoC_TOT_zint', 'SST', 'PRECT', 'PSL'};
% cfg.vars= {'PRECT', 'PSL'};
% cfg.vars= {'photoC_TOT_zint'};
% cfg.vars= {'sumChl'};
% cfg.vars = { 'sumChl', 'photoC_TOT_zint',  'SSH'};
cfg.vars = {  'SST', 'PRECT', 'PSL'};

for vari=1:length(cfg.vars)
    cfg.var=cfg.vars{vari};
    % cfg.var='SST'; %SST PRECT PSL TS SSH sumChl
    cfg.obs_name=f_obs_name(cfg.var);
    cfg.obs_fname_mid=f_obs_name_mid(cfg.var);
    cfg.obs_varname=f_obs_varname(cfg.var);
    cfg.comp=Func_0025_CESM2_cmpname_var(cfg.var);
    cfg.obs_fname_module=f_obs_fname_module(cfg.comp);
    cfg.obs_iyears=1970:2020;
%     cfg.obs_iyears=f_obs_iyears(cfg.var);
    
    
    dirs.hcstroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/HCST_EXP/archive/', cfg.comp, '/', cfg.var];
    dirs.obsroot=['/Volumes/kyy_raid/kimyy/Observation/', cfg.obs_name, '/monthly_reg_', cfg.obs_fname_module(2:4)];
    dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', cfg.var];
    dirs.lens2root=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/LENS2/archive_analysis/', cfg.comp, '/', cfg.var];
    dirs.assmroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/ASSM_EXP/archive_analysis/', cfg.comp, '/', cfg.var];
    
    cfg.iyears=cfg.obs_iyears;
    % cfg.months=1:12;
    % cfg.scnm='HIST';
    cfg.gnm='f09_g17';
    % cfg.assm_factor='10';
    % cfg.ens_member='1';
    % cfg.proj_year=5;
%     cfg.season = {'MAM', 'JJA', 'SON', 'DJF', 'AMJ', 'JAS', 'OND', 'JFM', 'INI'};
    cfg.season = {'AMJ', 'JAS', 'OND', 'JFM', 'INI'};
%     cfg.season = {'INI'};
%     cfg.season = {'FMA'};

%     tmp.mon_cover=[1:60];
%     for ii=tmp.mon_cover
%         cfg.season{ii-min(tmp.mon_cover)+1} = num2str(ii,'%02i'); %% all month
%     end
    
    % cfg.component='ocn';
    % cfg.varnames={'temp', 'salt'};
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
    tmp.maskname = '/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/grids/ocn/RECCAP2_region_masks_all_v20210412_POP2_grid.nc';

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
    
    S = shaperead('landareas.shp');
    
    %% read & plot data
    tmp.varname=cfg.var;
    
%     %% AR1 (data read)
%     for iyear=min(cfg.iyears):max(cfg.iyears)
%         tmp.iyear_str=num2str(iyear, '%04i');
%         tmp.iind=iyear-min(cfg.iyears)+1;
%         cfg.casename_m=['ens_all'];
%         cfg.casename=[cfg.casename_m, '_i', tmp.iyear_str];
%         for mon=1:12
%             tmp.mon_str=num2str(mon, '%02i');
%             cfg.obs_fnm=[dirs.obsroot, tmp.fs, cfg.obs_fname_mid,tmp.iyear_str,tmp.mon_str, '.nc'];
%             tmp.all_data_obs(1:grid.nlon,1:grid.nlat,(tmp.iind-1)*12+mon)=ncread(cfg.obs_fnm, cfg.obs_varname);
%         end
%     end
%     
%     tmp.all_data_obs_des_det = Func_0029_deseasonalize_detrend_linear_3d(tmp.all_data_obs, 'omitnan');
    
    % pcolor(tmp.all_data_obs_des_det(1:grid.nlon,1:grid.nlat,1)'); shading flat; colorbar;
    % plot(squeeze(tmp.all_data_obs_des_det(146,135,:)))
    
%     %% AR1 (corr)
%     for loni=1:grid.nlon
%         for lati=1:grid.nlat
%              if (isnan(tmp.all_data_obs(loni,lati,1))~=1 & nansum(tmp.all_data_obs(loni,lati,:))~=0)
%                  [tmp.corr_AC_lag1, tmp.corr_AC_lag1_p] = ...
%                      corrcoef(squeeze(tmp.all_data_obs(loni,lati,1:end-1)), squeeze(tmp.all_data_obs(loni,lati,2:end)));
%                  data.([tmp.varname, '_AC_lag1'])(loni,lati)=tmp.corr_AC_lag1(1,2);
%                  data.([tmp.varname, '_AC_lag1_p'])(loni,lati)=tmp.corr_AC_lag1_p(1,2);
%                  
%                  [tmp.corr_AC_lag1_des_det, tmp.corr_AC_lag1_des_det_p] = ...
%                      corrcoef(squeeze(tmp.all_data_obs_des_det(loni,lati,1:end-1)), squeeze(tmp.all_data_obs_des_det(loni,lati,2:end)));
%                  data.([tmp.varname, '_AC_lag1_des_det'])(loni,lati)=tmp.corr_AC_lag1_des_det(1,2);
%                  data.([tmp.varname, '_AC_lag1_des_det_p'])(loni,lati)=tmp.corr_AC_lag1_des_det_p(1,2);
%              else
%                  data.([tmp.varname, '_AC_lag1'])(loni,lati)=NaN;
%                  data.([tmp.varname, '_AC_lag1_p'])(loni,lati)=NaN;
%                  data.([tmp.varname, '_AC_lag1_des_det'])(loni,lati)=NaN;
%                  data.([tmp.varname, '_AC_lag1_des_det_p'])(loni,lati)=NaN;
%              end
%         end
%     end
    
    % pcolor(data.SST_AC_lag1'); shading flat; colorbar; %(much better) 
    % pcolor(data.SST_AC_lag1_des_det'); shading flat; colorbar; % normal skill (deseasonlized+detrended)
    
    for lss=1:length(cfg.season)
        tmp.season=cfg.season{lss};
        tmp.mons = f_season_mons(tmp.season);
        cfg.casename_m=['ens_all'];
    
        dirs.datadir= [dirs.hcstroot, filesep, cfg.casename_m, filesep];
        dirs.assmdir= [dirs.assmroot, filesep, cfg.casename_m, filesep];
        dirs.lens2dir= [dirs.lens2root, filesep, 'ens_smbb', filesep];
        fprintf('%s_%s_%s  ',tmp.season,cfg.casename_m, tmp.varname); lap_time = tic;
        
        %% variables initialization
        data.([tmp.varname, '_model', '_', tmp.season])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data.([tmp.varname, '_assm', '_', tmp.season])=NaN(grid.nlon, grid.nlat, cfg.len_t_y); 
        data.([tmp.varname, '_lens2', '_', tmp.season])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);         
        data.([tmp.varname, '_bias', '_', tmp.season])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data.([tmp.varname, '_obs', '_', tmp.season])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data.([tmp.varname, '_AR1', '_', tmp.season])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);  %% AR1 initialization
    
    %     data.([tmp.varname, '_assm'])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        
        %% read variables
        for iyear=min(cfg.iyears):max(cfg.iyears)
            tmp.iyear_str=num2str(iyear, '%04i');
            tmp.iind=iyear-min(cfg.iyears)+1;
            cfg.casename=[cfg.casename_m, '_i', tmp.iyear_str];
    %         tmp.fy=iyear+lss;
    %         tmp.fy_str=num2str(tmp.fy, '%04i');
            %% monthly filename
            for mon=1:length(tmp.mons)
                tmp.mon=tmp.mons(mon);
                if tmp.mon>=12 && mod(tmp.mon,12)~=0
                    tmp.mon_str=num2str(tmp.mon-12*floor(tmp.mon/12), '%02i');
                    tmp.fy_str=num2str(iyear+floor(tmp.mon/12), '%04i');
                elseif tmp.mon>=12 && mod(tmp.mon,12) ==0
                    tmp.mon_str='12';
                    tmp.fy_str=num2str(iyear+floor(tmp.mon/12)-1, '%04i');
                else
                    tmp.mon_str=num2str(tmp.mon, '%02i');
                    tmp.fy_str=num2str(iyear, '%04i');
                end

                %% HCST
                cfg.mod_fnm=[dirs.datadir, tmp.fs, cfg.casename, tmp.fs, ...
                tmp.varname, '_', cfg.gnm, '.hcst.', cfg.casename, cfg.obs_fname_module, tmp.fy_str, '-', tmp.mon_str, '.nc'];
%                 tmp.sdata(1:grid.nlon,1:grid.nlat,mon)=ncread(cfg.mod_fnm, tmp.varname);
                ncid=netcdf.open(cfg.mod_fnm, 'NOWRITE');
                tmpvarid=netcdf.inqVarID(ncid,tmp.varname);
                tmp.sdata(1:grid.nlon,1:grid.nlat,mon) = netcdf.getVar(ncid,tmpvarid);
                netcdf.close(ncid);

                %% LENS2
                cfg.lens2_fnm=[dirs.lens2dir, tmp.fs, ...
                tmp.varname, '_ensmean_', tmp.fy_str, '-', tmp.mon_str, '.nc'];
%                 tmp.sdata(1:grid.nlon,1:grid.nlat,mon)=ncread(cfg.mod_fnm, tmp.varname);
                ncid=netcdf.open(cfg.lens2_fnm, 'NOWRITE');
                tmpvarid=netcdf.inqVarID(ncid,tmp.varname);
                tmp.sdata_lens2(1:grid.nlon,1:grid.nlat,mon) = netcdf.getVar(ncid,tmpvarid);
                netcdf.close(ncid);
                
                %% OBS
    %             cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'ersst_reg_cesm2.v5.',tmp.iyear_str,tmp.mon_str, '.nc'];
                cfg.obs_fnm=[dirs.obsroot, tmp.fs, cfg.obs_fname_mid,tmp.fy_str,tmp.mon_str, '.nc'];
%                 tmp.sdata_obs(1:grid.nlon,1:grid.nlat,mon)=ncread(cfg.obs_fnm, cfg.obs_varname);
                if exist(cfg.obs_fnm)~=0
                    ncid=netcdf.open(cfg.obs_fnm, 'NOWRITE');
                    tmpvarid=netcdf.inqVarID(ncid,cfg.obs_varname);
                    tmp.sdata_obs(1:grid.nlon,1:grid.nlat,mon) = netcdf.getVar(ncid,tmpvarid);
                    netcdf.close(ncid);
                else
                    tmp.sdata_obs(1:grid.nlon,1:grid.nlat,mon) = NaN;                    
                end

                %% ASSM
                cfg.assm_fnm=[dirs.assmdir, tmp.fs, tmp.varname, '_ensmean_', tmp.fy_str,'-', tmp.mon_str, '.nc'];
                if exist(cfg.assm_fnm)~=0
                    ncid=netcdf.open(cfg.assm_fnm, 'NOWRITE');
                    tmpvarid=netcdf.inqVarID(ncid,tmp.varname);
                    tmp.sdata_assm(1:grid.nlon,1:grid.nlat,mon) = netcdf.getVar(ncid,tmpvarid);
                    netcdf.close(ncid);
                else
                    tmp.sdata_assm(1:grid.nlon,1:grid.nlat,mon) = NaN;                    
                end
                
%                 if (strcmp(cfg.var, 'sumChl') ||  strcmp(cfg.var, 'TS')) %for Chls, it uses all available data
%                     tmp.obs_mask(:,:)=tmp.sdata_obs(1:grid.nlon,1:grid.nlat,mon)./tmp.sdata_obs(1:grid.nlon,1:grid.nlat,mon);
%                     tmp.sdata(1:grid.nlon,1:grid.nlat,mon) = tmp.sdata(1:grid.nlon,1:grid.nlat,mon) .* tmp.obs_mask(:,:); % masking with available observation
%                 end

%                 %% AR1 (integration)
%                 tmp.sdata_AR1(1:grid.nlon,1:grid.nlat,mon) = ...
%                     Func_0030_AR1_prog(tmp.all_data_obs(1:grid.nlon,1:grid.nlat,(tmp.iind-1)*12+1), data.([tmp.varname, '_AC_lag1']), tmp.mon-1, 0);
                
            end
    %         tmp.sdata(tmp.sdata==0)=NaN;
    %         data.time=ncread(cfg.datafilename, 'time');

%             if (strcmp(cfg.var, 'sumChl') ||  strcmp(cfg.var, 'TS')) %for Chls, it uses all available data
%                 tmp.smean= mean(tmp.sdata,3, 'omitnan');
%                 tmp.smean_obs= mean(tmp.sdata_obs,3, 'omitnan');
%                 tmp.smean_assm= mean(tmp.sdata_assm,3, 'omitnan');
%                 tmp.smean_lens2= mean(tmp.sdata_lens2,3, 'omitnan');
%             else
                tmp.smean= mean(tmp.sdata,3);
                tmp.smean_obs= mean(tmp.sdata_obs,3);
                tmp.smean_assm= mean(tmp.sdata_assm,3);
                tmp.smean_lens2= mean(tmp.sdata_lens2,3);
%                 tmp.smean_assm(tmp.smean_assm>10e30)=NaN;
%                 tmp.smean_assm(tmp.smean_assm==0)=NaN;
                
                 switch cfg.comp
                    case 'ocn'
                        tmp.smean = tmp.smean .* grid.mask_ocn;
                        tmp.smean_obs = tmp.smean_obs .* grid.mask_ocn;
                        tmp.smean_assm = tmp.smean_assm .* grid.mask_ocn;
                        tmp.smean_lens2 = tmp.smean_lens2 .* grid.mask_ocn;                        
                 end


%             end
            data.([tmp.varname, '_model', '_', tmp.season])(1:grid.nlon,1:grid.nlat,tmp.iind)= tmp.smean;
            data.([tmp.varname, '_obs', '_', tmp.season])(1:grid.nlon,1:grid.nlat,tmp.iind)= tmp.smean_obs;
            data.([tmp.varname, '_assm', '_', tmp.season])(1:grid.nlon,1:grid.nlat,tmp.iind)= tmp.smean_assm;
            data.([tmp.varname, '_lens2', '_', tmp.season])(1:grid.nlon,1:grid.nlat,tmp.iind)= tmp.smean_lens2;


%             %% AR1 (assign)
%             tmp.smean_AR1= mean(tmp.sdata_AR1, 3);
%             data.([tmp.varname, '_AR1', '_', tmp.season])(1:grid.nlon,1:grid.nlat,tmp.iind)= tmp.smean_AR1;
    
    %         tmp.ymean= mean(ncread(cfg.datafilename, ['assm_', tmp.varname], [1, 1, 1], [inf, inf, 12]), 3);
    %         data.([tmp.varname, '_assm'])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean;
    
    
            
    
        end
    
        tmp.ano_obs=data.([tmp.varname, '_obs', '_', tmp.season]) - mean(data.([tmp.varname, '_obs', '_', tmp.season]), 3, 'omitnan');
        tmp.ano_mod=data.([tmp.varname, '_model', '_', tmp.season]) - mean(data.([tmp.varname, '_model', '_', tmp.season]), 3, 'omitnan');
        tmp.ano_assm=data.([tmp.varname, '_assm', '_', tmp.season]) - mean(data.([tmp.varname, '_assm', '_', tmp.season]), 3, 'omitnan');
        tmp.ano_lens2=data.([tmp.varname, '_lens2', '_', tmp.season]) - mean(data.([tmp.varname, '_lens2', '_', tmp.season]), 3, 'omitnan');
    
        %% get nRMSE
        tmp.sig_0=sqrt( sum( tmp.ano_obs.^2, 3 ) ./ cfg.len_t_y ); %denominator (sigma_0)
        tmp.RMSE=sqrt( sum( ( tmp.ano_mod-tmp.ano_obs ).^2, 3 ) ./ cfg.len_t_y);
        data.([tmp.varname, '_nRMSE', '_', tmp.season])=tmp.RMSE ./ tmp.sig_0;
    %     pcolor(data.([tmp.varname, '_nRMSE', '_', tmp.season])'); shading flat; colorbar; colormap(jet); caxis([0,1]);
    
        %% get correlation coefficient
        for loni=1:grid.nlon
            for lati=1:grid.nlat
                 if (isnan(data.([tmp.varname, '_assm', '_', tmp.season])(loni,lati,1))~=1 & nansum(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,:))~=0)
                     
                     
                 %% corr OBS
                     tmp.data = squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,1:end));
                     tmp.data_obs = data.([tmp.varname, '_obs', '_', tmp.season])(loni,lati,1:end);
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data, tmp.data_obs);
                     
                     tmp.data_det = Func_0028_detrend_linear_1d(squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,1:end)));
                     tmp.data_obs_det = Func_0028_detrend_linear_1d(data.([tmp.varname, '_obs', '_', tmp.season])(loni,lati,1:end));
                     [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data_det, tmp.data_obs_det);
    
    %                  [tmp.corr, tmp.corr_p]=corrcoef(squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,1+lss:end-cfg.proj_year+1)), ...
    %                      squeeze(data.([tmp.varname, '_obs', '_', tmp.season])(loni,lati,1+lss:end)));
                     data.([tmp.varname, '_corr_obs', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_obs_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
                     data.([tmp.varname, '_corr_obs_det', '_', tmp.season])(loni,lati)=tmp.corr_det(1,2);
                     data.([tmp.varname, '_corr_obs_det_p', '_', tmp.season])(loni,lati)=tmp.corr_det_p(1,2);
    %                  [tmp.corr, tmp.corr_p]=corrcoef(squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,1+lss:end-cfg.proj_year+1)), ...
    %                      squeeze(data.([tmp.varname, '_assm'])(loni,lati,1+lss:end)));
    %                  data.([tmp.varname, '_corr_assm', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
    %                  data.([tmp.varname, '_corr_assm_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
                    
                %% corr assm
                     tmp.data = squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,1:end));
                     tmp.data_assm = squeeze(data.([tmp.varname, '_assm', '_', tmp.season])(loni,lati,1:end));
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_assm)), tmp.data_assm(isfinite(tmp.data_assm)));
                     
                     tmp.data_hcst_lens2 = squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,1:end)) - ...
                                    squeeze(data.([tmp.varname, '_lens2', '_', tmp.season])(loni,lati,1:end));
                     [tmp.corr_int, tmp.corr_int_p]=corrcoef(tmp.data_hcst_lens2(isfinite(tmp.data_assm)), tmp.data_assm(isfinite(tmp.data_assm)));
    
                     data.([tmp.varname, '_corr_assm', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_assm_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
                     data.([tmp.varname, '_corr_assm_int', '_', tmp.season])(loni,lati)=tmp.corr_int(1,2);
                     data.([tmp.varname, '_corr_assm_int_p', '_', tmp.season])(loni,lati)=tmp.corr_int_p(1,2);
 
                %% corr lens2-assm
                     tmp.lens2 = squeeze(data.([tmp.varname, '_lens2', '_', tmp.season])(loni,lati,1:end));
                     tmp.data_assm = data.([tmp.varname, '_assm', '_', tmp.season])(loni,lati,1:end);
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.lens2(isfinite(tmp.data_assm)), tmp.data_assm(isfinite(tmp.data_assm)));

                     data.([tmp.varname, '_corr_assm_lens2', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_assm_lens2_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);


%                     %% AR1 (corr)
%                      tmp.data_AR1 = squeeze(data.([tmp.varname, '_AR1', '_', tmp.season])(loni,lati,1:end));
%                      [tmp.corr_AR1, tmp.corr_AR1_p]=corrcoef(tmp.data_AR1, tmp.data_obs);
%                      data.([tmp.varname, '_corr_AR1', '_', tmp.season])(loni,lati)=tmp.corr_AR1(1,2);
%                      data.([tmp.varname, '_corr_AR1_p', '_', tmp.season])(loni,lati)=tmp.corr_AR1_p(1,2);
%                     
%                      tmp.data_AR1_det = Func_0028_detrend_linear_1d(data.([tmp.varname, '_AR1', '_', tmp.season])(loni,lati,1:end));
%                      [tmp.corr_AR1_det, tmp.corr_AR1_det_p]=corrcoef(tmp.data_AR1_det, tmp.data_obs_det);
%                      data.([tmp.varname, '_corr_AR1_det', '_', tmp.season])(loni,lati)=tmp.corr_AR1_det(1,2);
%                      data.([tmp.varname, '_corr_AR1_det_p', '_', tmp.season])(loni,lati)=tmp.corr_AR1_det_p(1,2);
    
    
                 else
                     data.([tmp.varname, '_corr_obs', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_p', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_det', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_det_p', '_', tmp.season])(loni,lati)=NaN;
    %                  data.([tmp.varname, '_corr_assm', '_', tmp.season])(loni,lati)=NaN;
    %                  data.([tmp.varname, '_corr_assm_p', '_', tmp.season])(loni,lati)=NaN;  
    
% %                      %% AR1 (NaN)
% %                      data.([tmp.varname, '_corr_AR1', '_', tmp.season])(loni,lati)=NaN;
% %                      data.([tmp.varname, '_corr_AR1_p', '_', tmp.season])(loni,lati)=NaN;
% %                      data.([tmp.varname, '_corr_AR1_det', '_', tmp.season])(loni,lati)=NaN;
% %                      data.([tmp.varname, '_corr_AR1_det_p', '_', tmp.season])(loni,lati)=NaN;
 
                     %% ASSM (NaN)
                     data.([tmp.varname, '_corr_assm', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_p', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_int', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_int_p', '_', tmp.season])(loni,lati)=NaN;

                     %% LENS2 (NaN)
                     data.([tmp.varname, '_corr_assm_lens2', '_', tmp.season])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_assm_lens2_p', '_', tmp.season])(loni,lati)=NaN;

                 end
            end
        end
        disp('abc')
    end
    % pcolor(data.([tmp.varname, '_corr_AR1', '_', 'MAM'])'); shading flat; colorbar;
    % pcolor(data.([tmp.varname, '_corr_AR1', '_', 'JJA'])'); shading flat; colorbar;
    % pcolor(data.([tmp.varname, '_corr_AR1', '_', 'SON'])'); shading flat; colorbar;
    % pcolor(data.([tmp.varname, '_corr_AR1', '_', 'DJF'])'); shading flat; colorbar;
    
    for lss_ind=1:length(cfg.season)
    %     tmp.lss_str=tmp.yearset{lss_ind};
        tmp.season=cfg.season{lss_ind};
        
% %         %% model & obs corr map --------------------------------------
% %         fig_cfg.name_rgn = 'Glob';
% %         fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
% %         fig_cfg.x_lim = [-180 180];
% %         fig_cfg.y_lim = [-80 89];
% %         fig_cfg.fig_size = [0,0,6,3.5];
% %         fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
% %         fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
% %         fig_cfg.title_pos = [0.5,0.93];
% %         fig_cfg.p_lim =0.1;
% %         fig_cfg.c_lim = [-1 1];
% %         [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
% %     
% %         tmp.X=grid.tlong([end, 1:end],:);
% %         tmp.Y=grid.tlat([end, 1:end],:);
% %         tmp.C=data.([tmp.varname, '_corr_obs', '_', tmp.season]);
% %         tmp.C=tmp.C([end, 1:end],:);
% % %         tmp.C_H=data.([tmp.varname, '_corr_AR1', '_', tmp.season]);
% % %         tmp.C_H=tmp.C_H([end, 1:end],:);
% % %         tmp.C_2=tmp.C;
% % %         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 
% %     
% %     
% %         [tmp.mean_corr, tmp.err] = ...
% %             Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_obs', '_', tmp.season]), grid.tlong, grid.tlat);
% %         fig_cfg.fig_name=[tmp.season, ', corr_ob_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
% %         fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
% %             'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
% %         %% map setting
% %         ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
% %             'fontname','freeserif'); 
% %     
% %         axis off; 
% %         hold on;
% %         setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
% %         set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
% %         text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
% %         'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
% %         'fontsize',14,'fontname','freeserif','interpreter','none')
% %     
% %         %% caxis & colorbar
% %         caxis(ax_m, fig_cfg.c_lim); 
% %         colormap(fig_cfg.c_map);
% %         cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
% %         set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
% %         title(cb,'R','fontsize',12);
% %     
% %         %% draw on ax_m
% %         h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
% %         shading flat;
% %         geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
% %     
% % %         if (strcmp(tmp.season, 'INI')~=1 & strcmp(cfg.var, 'sumChl')~=1) %for Chls, it uses all available data
% % %             %% <AR1 area -> hatch
% % %             pp2 = pcolorm(tmp.Y,tmp.X,tmp.C_2, 'parent', ax_m);
% % %             set(pp2,'linestyle','none','Tag','HatchingRegion');
% % %             hp = findobj(pp2,'Tag','HatchingRegion');
% % %             hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);
% % %         end
% %     
% %     %     %% >AR1 area -> dot
% %     %     tmp.Y2=tmp.Y(isfinite(tmp.C_2));
% %     %     tmp.X2=tmp.X(isfinite(tmp.C_2));
% %     %     tmp.dens=15;
% %     %     tmp.Y2=tmp.Y2(1:tmp.dens:end);
% %     %     tmp.X2=tmp.X2(1:tmp.dens:end);
% %     %     scatterm(tmp.Y2,tmp.X2, 5, 'k', 'o', 'parent',ax_m); 
% %     
% %         %% frame and label setting
% %         setm(ax_m,'frame','on','FLineWidth',1);
% %     
% %         label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
% %         label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
% %         mlabel; plabel;
% %         label_y=plabel; label_x=mlabel;
% %         for lxi=1:length(label_x)
% %             tmp.tmppos=label_x(lxi,1).Position;
% %             tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
% %             label_x(lxi,1).Position=tmp.tmppos;
% %             label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
% %         end
% %         for lyi=1:length(label_y)
% %             label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
% %         end
% %     
% %         %% save
% %         dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_map', filesep, 'model'];
% %         if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
% %         cfg.figname=[dirs.figdir, filesep, 'corr_obs_map_', tmp.varname, '_', tmp.season, '.tif'];
% %         print(fig_h, cfg.figname, '-dpng');
% %         RemoveWhiteSpace([], 'file', cfg.figname);
% %         fprintf('%7.1f sec\n', toc(lap_time) );
% %         close all;
    
% %     
% %     %% model & obs corr map (det) --------------------------------------
% %         fig_cfg.name_rgn = 'Glob';
% %         fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
% %         fig_cfg.x_lim = [-180 180];
% %         fig_cfg.y_lim = [-80 89];
% %         fig_cfg.fig_size = [0,0,6,3.5];
% %         fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
% %         fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
% %         fig_cfg.title_pos = [0.5,0.93];
% %         fig_cfg.p_lim =0.1;
% %         fig_cfg.c_lim = [-1 1];
% %         [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
% %     
% %         tmp.X=grid.tlong([end, 1:end],:);
% %         tmp.Y=grid.tlat([end, 1:end],:);
% %         tmp.C=data.([tmp.varname, '_corr_obs_det', '_', tmp.season]);
% %         tmp.C=tmp.C([end, 1:end],:);
% %         tmp.C_H=data.([tmp.varname, '_corr_AR1_det', '_', tmp.season]);
% %         tmp.C_H=tmp.C_H([end, 1:end],:);
% %         tmp.C_2=tmp.C;
% %         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 
% %     
% %         [tmp.mean_corr, tmp.err] = ...
% %             Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_obs_det', '_', tmp.season]), grid.tlong, grid.tlat);
% %         fig_cfg.fig_name=[tmp.season, ', corr_ob_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
% %         fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
% %             'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
% %         %% map setting
% %         ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
% %             'fontname','freeserif'); 
% %     
% %         axis off; 
% %         hold on;
% %         setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
% %         set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
% %         text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
% %         'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
% %         'fontsize',14,'fontname','freeserif','interpreter','none')
% %     
% %         %% caxis & colorbar
% %         caxis(ax_m, fig_cfg.c_lim); 
% %         colormap(fig_cfg.c_map);
% %         cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
% %         set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
% %         title(cb,'R','fontsize',12);
% %     
% %         %% draw on ax_m
% %         h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
% %         shading flat;
% %         geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
% %     
% %         if (strcmp(tmp.season, 'INI')~=1 & strcmp(cfg.var, 'sumChl')~=1) %for Chls, it uses all available data
% %             %% <AR1 area -> hatch
% %             pp2 = pcolorm(tmp.Y,tmp.X,tmp.C_2, 'parent', ax_m);
% %             set(pp2,'linestyle','none','Tag','HatchingRegion');
% %             hp = findobj(pp2,'Tag','HatchingRegion');
% %             hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','k','HatchLineWidth',0.5);
% %         end
% %         %% frame and label setting
% %         setm(ax_m,'frame','on','FLineWidth',1);
% %     
% %         label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
% %         label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
% %         mlabel; plabel;
% %         label_y=plabel; label_x=mlabel;
% %         for lxi=1:length(label_x)
% %             tmp.tmppos=label_x(lxi,1).Position;
% %             tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
% %             label_x(lxi,1).Position=tmp.tmppos;
% %             label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
% %         end
% %         for lyi=1:length(label_y)
% %             label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
% %         end
% %     
% %         %% save
% %         dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_map', filesep, 'model_det'];
% %         if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
% %         cfg.figname=[dirs.figdir, filesep, 'corr_obs_det_map_', tmp.varname, '_', tmp.season, '.tif'];
% %         print(fig_h, cfg.figname, '-dpng');
% %         RemoveWhiteSpace([], 'file', cfg.figname);
% %         fprintf('%7.1f sec\n', toc(lap_time) );
% %         close all;
    
    
% %     %% AR1 & obs corr map --------------------------------------
% %         fig_cfg.name_rgn = 'Glob';
% %         fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
% %         fig_cfg.x_lim = [-180 180];
% %         fig_cfg.y_lim = [-80 89];
% %         fig_cfg.fig_size = [0,0,6,3.5];
% %         fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
% %         fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
% %         fig_cfg.title_pos = [0.5,0.93];
% %         fig_cfg.p_lim =0.1;
% %         fig_cfg.c_lim = [-1 1];
% %         [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
% %     
% %         tmp.X=grid.tlong([end, 1:end],:);
% %         tmp.Y=grid.tlat([end, 1:end],:);
% %         tmp.C=data.([tmp.varname, '_corr_AR1', '_', tmp.season]);
% %         tmp.C=tmp.C([end, 1:end],:);
% %     
% %     
% %         [tmp.mean_corr, tmp.err] = ...
% %             Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_AR1', '_', tmp.season]), grid.tlong, grid.tlat);
% %         fig_cfg.fig_name=[tmp.season, ', corr_AR1_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
% %         fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
% %             'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
% %         %% map setting
% %         ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
% %             'fontname','freeserif'); 
% %     
% %         axis off; 
% %         hold on;
% %         setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
% %         set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
% %         text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
% %         'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
% %         'fontsize',14,'fontname','freeserif','interpreter','none')
% %     
% %         %% caxis & colorbar
% %         caxis(ax_m, fig_cfg.c_lim); 
% %         colormap(fig_cfg.c_map);
% %         cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
% %         set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
% %         title(cb,'R','fontsize',12);
% %     
% %         %% draw on ax_m
% %         h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
% %         shading flat;
% %         geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
% %     
% %     %             pp2 = pcolorm(tmp.Y,tmp.X,tmp.C_2, 'parent', ax_m);
% %     %             set(pp2,'linestyle','none','Tag','HatchingRegion');
% %     %             hp = findobj(pp2,'Tag','HatchingRegion');
% %     %             hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','k','HatchLineWidth',0.5);
% %     
% %         %% frame and label setting
% %         setm(ax_m,'frame','on','FLineWidth',1);
% %     
% %         label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
% %         label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
% %         mlabel; plabel;
% %         label_y=plabel; label_x=mlabel;
% %         for lxi=1:length(label_x)
% %             tmp.tmppos=label_x(lxi,1).Position;
% %             tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
% %             label_x(lxi,1).Position=tmp.tmppos;
% %             label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
% %         end
% %         for lyi=1:length(label_y)
% %             label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
% %         end
% %     
% %         %% save
% %         dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_map', filesep, 'AR1'];
% %         if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
% %         cfg.figname=[dirs.figdir, filesep, 'corr_AR1_map_', tmp.varname, '_', tmp.season, '.tif'];
% %         print(fig_h, cfg.figname, '-dpng');
% %         RemoveWhiteSpace([], 'file', cfg.figname);
% %         fprintf('%7.1f sec\n', toc(lap_time) );
% %         close all;
% %     
% %     %% AR1_det & obs_det corr map --------------------------------------
% %         fig_cfg.name_rgn = 'Glob';
% %         fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
% %         fig_cfg.x_lim = [-180 180];
% %         fig_cfg.y_lim = [-80 89];
% %         fig_cfg.fig_size = [0,0,6,3.5];
% %         fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
% %         fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
% %         fig_cfg.title_pos = [0.5,0.93];
% %         fig_cfg.p_lim =0.1;
% %         fig_cfg.c_lim = [-1 1];
% %         [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
% %     
% %         tmp.X=grid.tlong([end, 1:end],:);
% %         tmp.Y=grid.tlat([end, 1:end],:);
% %         tmp.C=data.([tmp.varname, '_corr_AR1_det', '_', tmp.season]);
% %         tmp.C=tmp.C([end, 1:end],:);
% %     
% %     
% %         [tmp.mean_corr, tmp.err] = ...
% %             Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_AR1_det', '_', tmp.season]), grid.tlong, grid.tlat);
% %         fig_cfg.fig_name=[tmp.season, ', corr_AR1_det_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
% %         fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
% %             'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
% %         %% map setting
% %         ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
% %             'fontname','freeserif'); 
% %     
% %         axis off; 
% %         hold on;
% %         setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
% %         set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
% %         text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
% %         'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
% %         'fontsize',14,'fontname','freeserif','interpreter','none')
% %     
% %         %% caxis & colorbar
% %         caxis(ax_m, fig_cfg.c_lim); 
% %         colormap(fig_cfg.c_map);
% %         cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
% %         set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
% %         title(cb,'R','fontsize',12);
% %     
% %         %% draw on ax_m
% %         h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
% %         shading flat;
% %         geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
% %     
% %     %             pp2 = pcolorm(tmp.Y,tmp.X,tmp.C_2, 'parent', ax_m);
% %     %             set(pp2,'linestyle','none','Tag','HatchingRegion');
% %     %             hp = findobj(pp2,'Tag','HatchingRegion');
% %     %             hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','k','HatchLineWidth',0.5);
% %     
% %         %% frame and label setting
% %         setm(ax_m,'frame','on','FLineWidth',1);
% %     
% %         label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
% %         label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
% %         mlabel; plabel;
% %         label_y=plabel; label_x=mlabel;
% %         for lxi=1:length(label_x)
% %             tmp.tmppos=label_x(lxi,1).Position;
% %             tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
% %             label_x(lxi,1).Position=tmp.tmppos;
% %             label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
% %         end
% %         for lyi=1:length(label_y)
% %             label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
% %         end
% %     
% %         %% save
% %         dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_map', filesep, 'AR1_det'];
% %         if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
% %         cfg.figname=[dirs.figdir, filesep, 'corr_AR1_det_map_', tmp.varname, '_', tmp.season, '.tif'];
% %         print(fig_h, cfg.figname, '-dpng');
% %         RemoveWhiteSpace([], 'file', cfg.figname);
% %         fprintf('%7.1f sec\n', toc(lap_time) );
% %         close all;
% %     
% %     
% %     
% %     %% model & obs nRMSE map --------------------------------------
% %         fig_cfg.name_rgn = 'Glob';
% %         fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
% %         fig_cfg.x_lim = [-180 180];
% %         fig_cfg.y_lim = [-80 89];
% %         fig_cfg.fig_size = [0,0,6,3.5];
% %         fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
% %         fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
% %         fig_cfg.title_pos = [0.5,0.93];
% %         fig_cfg.p_lim =0.1;
% %         fig_cfg.c_lim = [0 1];
% %         [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr', tmp.dropboxpath);
% %         fig_cfg.c_map= flip(fig_cfg.c_map); %wr -> rw
% %     
% %         tmp.X=grid.tlong([end, 1:end],:);
% %         tmp.Y=grid.tlat([end, 1:end],:);
% %         tmp.C=data.([tmp.varname, '_nRMSE', '_', tmp.season]);
% %         tmp.C=tmp.C([end, 1:end],:);
% %     
% %         [tmp.mean_corr, tmp.err] = ...
% %             Func_0011_get_area_weighted_mean(data.([tmp.varname, '_nRMSE', '_', tmp.season]), grid.tlong, grid.tlat);
% %         fig_cfg.fig_name=[tmp.season, ', nRMSE, ', tmp. varname];
% %         fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
% %             'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
% %         %% map setting
% %         ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
% %             'fontname','freeserif'); 
% %     
% %         axis off; 
% %         hold on;
% %         setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
% %         set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
% %         text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
% %         'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
% %         'fontsize',14,'fontname','freeserif','interpreter','none')
% %     
% %         %% caxis & colorbar
% %         caxis(ax_m, fig_cfg.c_lim); 
% %         colormap(fig_cfg.c_map);
% %         cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
% %         set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
% %         title(cb,'R','fontsize',12);
% %     
% %         %% draw on ax_m
% %         h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
% %         shading flat;
% %         geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
% %     
% %         %% frame and label setting
% %         setm(ax_m,'frame','on','FLineWidth',1);
% %     
% %         label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
% %         label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
% %         mlabel; plabel;
% %         label_y=plabel; label_x=mlabel;
% %         for lxi=1:length(label_x)
% %             tmp.tmppos=label_x(lxi,1).Position;
% %             tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
% %             label_x(lxi,1).Position=tmp.tmppos;
% %             label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
% %         end
% %         for lyi=1:length(label_y)
% %             label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
% %         end
% %     
% %         %% save
% %         dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_NRMSE_obs_map'];
% %         if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
% %         cfg.figname=[dirs.figdir, filesep, 'nRMSE_obs_map_', tmp.varname, '_', tmp.season, '.tif'];
% %         print(fig_h, cfg.figname, '-dpng');
% %         RemoveWhiteSpace([], 'file', cfg.figname);
% %         fprintf('%7.1f sec\n', toc(lap_time) );
% %         close all;
    
    %% model & assm corr map --------------------------------------
        fig_cfg.name_rgn = 'Glob';
        fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
        fig_cfg.x_lim = [-180 180];
        fig_cfg.y_lim = [-80 89];
        fig_cfg.fig_size = [0,0,6,3.5];
        fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
        fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
        fig_cfg.title_pos = [0.5,0.93];
        fig_cfg.p_lim =0.1;
        fig_cfg.c_lim = [-1 1];
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=data.([tmp.varname, '_corr_assm', '_', tmp.season]);
        tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 
    
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_assm', '_', tmp.season]), grid.tlong, grid.tlat);
        fig_cfg.fig_name=[tmp.season, ', corr_ob_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
        fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
            'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
        %% map setting
        ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
            'fontname','freeserif'); 
    
        axis off; 
        hold on;
        setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
        set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',14,'fontname','freeserif','interpreter','none')
    
        %% caxis & colorbar
        caxis(ax_m, fig_cfg.c_lim); 
        colormap(fig_cfg.c_map);
        cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
        set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
        title(cb,'R','fontsize',12);
    
        %% draw on ax_m
        h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
        shading flat;
        geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
    
%         if (strcmp(tmp.season, 'INI')~=1 & strcmp(cfg.var, 'sumChl')~=1) %for Chls, it uses all available data
%             %% <AR1 area -> hatch
%             pp2 = pcolorm(tmp.Y,tmp.X,tmp.C_2, 'parent', ax_m);
%             set(pp2,'linestyle','none','Tag','HatchingRegion');
%             hp = findobj(pp2,'Tag','HatchingRegion');
%             hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);
%         end
   
    
        %% frame and label setting
        setm(ax_m,'frame','on','FLineWidth',1);
    
        label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
        label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
        mlabel; plabel;
        label_y=plabel; label_x=mlabel;
        for lxi=1:length(label_x)
            tmp.tmppos=label_x(lxi,1).Position;
            tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
            label_x(lxi,1).Position=tmp.tmppos;
            label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
        end
        for lyi=1:length(label_y)
            label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
        end
    
        %% save
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_assm_map', filesep, 'model'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_assm_map_', tmp.varname, '_', tmp.season, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        fprintf('%7.1f sec\n', toc(lap_time) );
        close all;


     %% model-lens2 & assm corr map --------------------------------------
        fig_cfg.name_rgn = 'Glob';
        fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
        fig_cfg.x_lim = [-180 180];
        fig_cfg.y_lim = [-80 89];
        fig_cfg.fig_size = [0,0,6,3.5];
        fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
        fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
        fig_cfg.title_pos = [0.5,0.93];
        fig_cfg.p_lim =0.1;
        fig_cfg.c_lim = [-1 1];
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=data.([tmp.varname, '_corr_assm_int', '_', tmp.season]);
        tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 
    
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_assm_int', '_', tmp.season]), grid.tlong, grid.tlat);
        fig_cfg.fig_name=[tmp.season, ', corr_ob_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
        fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
            'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
        %% map setting
        ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
            'fontname','freeserif'); 
    
        axis off; 
        hold on;
        setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
        set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',14,'fontname','freeserif','interpreter','none')
    
        %% caxis & colorbar
        caxis(ax_m, fig_cfg.c_lim); 
        colormap(fig_cfg.c_map);
        cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
        set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
        title(cb,'R','fontsize',12);
    
        %% draw on ax_m
        h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
        shading flat;
        geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
    
%         if (strcmp(tmp.season, 'INI')~=1 & strcmp(cfg.var, 'sumChl')~=1) %for Chls, it uses all available data
%             %% <AR1 area -> hatch
%             pp2 = pcolorm(tmp.Y,tmp.X,tmp.C_2, 'parent', ax_m);
%             set(pp2,'linestyle','none','Tag','HatchingRegion');
%             hp = findobj(pp2,'Tag','HatchingRegion');
%             hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);
%         end
   
    
        %% frame and label setting
        setm(ax_m,'frame','on','FLineWidth',1);
    
        label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
        label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
        mlabel; plabel;
        label_y=plabel; label_x=mlabel;
        for lxi=1:length(label_x)
            tmp.tmppos=label_x(lxi,1).Position;
            tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
            label_x(lxi,1).Position=tmp.tmppos;
            label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
        end
        for lyi=1:length(label_y)
            label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
        end
    
        %% save
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_assm_int_map', filesep, 'model'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_assm_int_map_', tmp.varname, '_', tmp.season, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        fprintf('%7.1f sec\n', toc(lap_time) );
        close all;


    %% lens2 & assm corr map --------------------------------------
        fig_cfg.name_rgn = 'Glob';
        fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
        fig_cfg.x_lim = [-180 180];
        fig_cfg.y_lim = [-80 89];
        fig_cfg.fig_size = [0,0,6,3.5];
        fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
        fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
        fig_cfg.title_pos = [0.5,0.93];
        fig_cfg.p_lim =0.1;
        fig_cfg.c_lim = [-1 1];
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
    
        tmp.X=grid.tlong([end, 1:end],:);
        tmp.Y=grid.tlat([end, 1:end],:);
        tmp.C=data.([tmp.varname, '_corr_assm_lens2', '_', tmp.season]);
        tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data.([tmp.varname, '_corr_AR1', '_', tmp.season]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 
    
    
        [tmp.mean_corr, tmp.err] = ...
            Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_assm_lens2', '_', tmp.season]), grid.tlong, grid.tlat);
        fig_cfg.fig_name=[tmp.season, ', corr_ob_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
        fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
            'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
        %% map setting
        ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
            'fontname','freeserif'); 
    
        axis off; 
        hold on;
        setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
        set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',14,'fontname','freeserif','interpreter','none')
    
        %% caxis & colorbar
        caxis(ax_m, fig_cfg.c_lim); 
        colormap(fig_cfg.c_map);
        cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
        set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
        title(cb,'R','fontsize',12);
    
        %% draw on ax_m
        h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
        shading flat;
        geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
    
%         if (strcmp(tmp.season, 'INI')~=1 & strcmp(cfg.var, 'sumChl')~=1) %for Chls, it uses all available data
%             %% <AR1 area -> hatch
%             pp2 = pcolorm(tmp.Y,tmp.X,tmp.C_2, 'parent', ax_m);
%             set(pp2,'linestyle','none','Tag','HatchingRegion');
%             hp = findobj(pp2,'Tag','HatchingRegion');
%             hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);
%         end
   
    
        %% frame and label setting
        setm(ax_m,'frame','on','FLineWidth',1);
    
        label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
        label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
        mlabel; plabel;
        label_y=plabel; label_x=mlabel;
        for lxi=1:length(label_x)
            tmp.tmppos=label_x(lxi,1).Position;
            tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
            label_x(lxi,1).Position=tmp.tmppos;
            label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
        end
        for lyi=1:length(label_y)
            label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
        end
    
        %% save
        dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_assm_ext_map', filesep, 'model'];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'corr_assm_ext_map_', tmp.varname, '_', tmp.season, '.tif'];
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        fprintf('%7.1f sec\n', toc(lap_time) );
        close all;

    end

end


% tmp.assm=squeeze(data.([tmp.varname, '_assm', '_', tmp.season])(200,180,:));
% tmp.int=squeeze(data.([tmp.varname, '_assm', '_', tmp.season])(200,180,:)) - squeeze(data.([tmp.varname, '_lens2', '_', tmp.season])(200,180,:));
% tmp.lens2=squeeze(data.([tmp.varname, '_lens2', '_', tmp.season])(200,180,:));
% plot(tmp.assm-mean(tmp.assm), 'linewidth', 2)
% hold on
% plot(tmp.int-mean(tmp.int), 'linewidth', 2)
% plot(tmp.lens2-mean(tmp.lens2), 'linewidth', 2)
% legend({'assm', 'int', 'ext'})
% hold off
% set(gca, 'fontsize', 20)
% corr(tmp.int, tmp.assm)


function obsname_simple = f_obs_simple(obsname)
    switch obsname
        case 'en4.2_ba'
            obsname_simple='en4';
        case 'projdv7.3'
            obsname_simple='projd';
    end
end

function gmsst = f_gm_var(var_2d, area)
    mask_var=NaN(size(var_2d));
    mask_var(isfinite(var_2d))=1;
    area=area.*mask_var;
    var_2d=squeeze(var_2d);
    var_2d_sum=var_2d.*area;
    var_2d_sum=sum(var_2d_sum(:), 'omitnan');
    gmsst=var_2d_sum ./ sum(area(:), 'omitnan');
end

function varname_C = f_varname_C(varname)
    switch varname
        case 'temp'
            varname_C='TEMP';
        case 'salt'
            varname_C='SALT';
    end
end

function varname_unit = f_varname_unit(varname)
    switch varname
        case 'temp'
            varname_unit='\circC';
        case 'salt'
            varname_unit='g/kg';
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



