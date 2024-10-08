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

%% model configuration
% cfg.var='TS'; %SST PRECT PSL TS SSH sumChl
% cfg.vars = {'SSH', 'sumChl'};
% % cfg.vars = {'SST', 'PRECT', 'PSL', 'TS'};
% cfg.vars = {'SSH'};
% cfg.vars = {'SST', 'PRECT', 'PSL', 'TS', 'sumChl', 'photoC_TOT_zint'};
cfg.vars = { 'sumChl', 'photoC_TOT_zint',  'SSH'};
% cfg.vars = {  'SST', 'PRECT', 'PSL'};

for vari=1:length(cfg.vars)
    cfg.var=cfg.vars{vari};
cfg.obs_name=f_obs_name(cfg.var);
cfg.obs_fname_mid=f_obs_name_mid(cfg.var);
cfg.obs_varname=f_obs_varname(cfg.var);
cfg.comp=Func_0025_CESM2_cmpname_var(cfg.var);
cfg.obs_fname_module=f_obs_fname_module(cfg.comp);
cfg.obs_iyears=f_obs_iyears(cfg.var);

dirs.hcstroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/HCST_EXP/archive/', cfg.comp, '/', cfg.var];
dirs.obsroot=['/Volumes/kyy_raid/kimyy/Observation/', cfg.obs_name, '/monthly_reg_', cfg.obs_fname_module(2:4)];
dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', cfg.var];
dirs.lens2root=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/LENS2/archive_analysis/', cfg.comp, '/', cfg.var];
dirs.assmroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/ASSM_EXP/archive_analysis/', cfg.comp, '/', cfg.var];

cfg.iyears=cfg.obs_iyears;
cfg.gnm='f09_g17';
cfg.proj_year=5;

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

% % %% AR1 (data read)
% % for iyear=min(cfg.iyears):max(cfg.iyears)
% %     tmp.iyear_str=num2str(iyear, '%04i');
% %     tmp.iind=iyear-min(cfg.iyears)+1;
% %     cfg.casename_m=['ens_all'];
% %     cfg.casename=[cfg.casename_m, '_i', tmp.iyear_str];
% %     for mon=1:12
% %         tmp.mon_str=num2str(mon, '%02i');
% %         cfg.obs_fnm=[dirs.obsroot, tmp.fs, cfg.obs_fname_mid,tmp.iyear_str,tmp.mon_str, '.nc'];
% %         tmp.all_data_obs(:,:,(tmp.iind-1)*12+mon)=ncread(cfg.obs_fnm, cfg.obs_varname);
% %     end
% % end
% % 
% % tmp.all_data_obs_des_det = Func_0029_deseasonalize_detrend_linear_3d(tmp.all_data_obs, 'omitnan');

% pcolor(tmp.all_data_obs_des_det(:,:,1)'); shading flat; colorbar;
% plot(squeeze(tmp.all_data_obs_des_det(146,135,:)))

% % %% AR1 (corr)
% % for loni=1:grid.nlon
% %     for lati=1:grid.nlat
% %          if (isnan(tmp.all_data_obs(loni,lati,1))~=1 & nansum(tmp.all_data_obs(loni,lati,:))~=0)
% %              [tmp.corr_AC_lag1, tmp.corr_AC_lag1_p] = ...
% %                  corrcoef(squeeze(tmp.all_data_obs(loni,lati,1:end-1)), squeeze(tmp.all_data_obs(loni,lati,2:end)));
% %              data.([tmp.varname, '_AC_lag1'])(loni,lati)=tmp.corr_AC_lag1(1,2);
% %              data.([tmp.varname, '_AC_lag1_p'])(loni,lati)=tmp.corr_AC_lag1_p(1,2);
% %              
% %              [tmp.corr_AC_lag1_des_det, tmp.corr_AC_lag1_des_det_p] = ...
% %                  corrcoef(squeeze(tmp.all_data_obs_des_det(loni,lati,1:end-1)), squeeze(tmp.all_data_obs_des_det(loni,lati,2:end)));
% %              data.([tmp.varname, '_AC_lag1_des_det'])(loni,lati)=tmp.corr_AC_lag1_des_det(1,2);
% %              data.([tmp.varname, '_AC_lag1_des_det_p'])(loni,lati)=tmp.corr_AC_lag1_des_det_p(1,2);
% %          else
% %              data.([tmp.varname, '_AC_lag1'])(loni,lati)=NaN;
% %              data.([tmp.varname, '_AC_lag1_p'])(loni,lati)=NaN;
% %              data.([tmp.varname, '_AC_lag1_des_det'])(loni,lati)=NaN;
% %              data.([tmp.varname, '_AC_lag1_des_det_p'])(loni,lati)=NaN;
% %          end
% %     end
% % end

clear tmp.ydata tmp.ydata_lens2 tmp.ydata_obs tmp.ydata_assm

for lyear=1:cfg.proj_year-1
    tmp.lyear_str=num2str(lyear, '%02i');
    cfg.casename_m=['ens_all'];

    dirs.datadir= [dirs.hcstroot, filesep, cfg.casename_m, filesep];
    dirs.assmdir= [dirs.assmroot, filesep, cfg.casename_m, filesep];
    dirs.lens2dir= [dirs.lens2root, filesep, 'ens_smbb', filesep];
    fprintf('%02d_%s_%s  ',lyear,cfg.casename_m, tmp.varname); lap_time = tic;

    
    %% variables initialization
    data.([tmp.varname, '_model', '_l', tmp.lyear_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
    data.([tmp.varname, '_lens2', '_l', tmp.lyear_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
    data.([tmp.varname, '_bias', '_l', tmp.lyear_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
    data.([tmp.varname, '_obs'])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
    data.([tmp.varname, '_assm'])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
    data.([tmp.varname, '_AR1', '_l', tmp.lyear_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);  %% AR1 initialization

%     data.([tmp.varname, '_assm'])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
    
    %% read variables
    for iyear=min(cfg.iyears):max(cfg.iyears)
        tmp.iyear_str=num2str(iyear, '%04i');
        cfg.casename=[cfg.casename_m, '_i', tmp.iyear_str];
        tmp.fy=iyear+lyear;
        tmp.fy_str=num2str(tmp.fy, '%04i');
        %% monthly filename
        for mon=1:12
            tmp.mon_str=num2str(mon, '%02i');

            %% HCST
            cfg.mod_fnm=[dirs.datadir, tmp.fs, cfg.casename, tmp.fs, ...
            tmp.varname, '_', cfg.gnm, '.hcst.', cfg.casename, cfg.obs_fname_module, tmp.fy_str, '-', tmp.mon_str, '.nc'];
            ncid=netcdf.open(cfg.mod_fnm, 'NOWRITE');
            tmpvarid=netcdf.inqVarID(ncid,tmp.varname);
            tmp.ydata(1:grid.nlon,1:grid.nlat,mon) = netcdf.getVar(ncid,tmpvarid);
            netcdf.close(ncid);
%             tmp.ydata(:,:,mon)=ncread(cfg.mod_fnm, tmp.varname);
            
            %% LENS2
            cfg.lens2_fnm=[dirs.lens2dir, tmp.fs, ...
                tmp.varname, '_ensmean_', tmp.fy_str, '-', tmp.mon_str, '.nc'];
            ncid=netcdf.open(cfg.lens2_fnm, 'NOWRITE');
            tmpvarid=netcdf.inqVarID(ncid,tmp.varname);
            tmp.ydata_lens2(1:grid.nlon,1:grid.nlat,mon) = netcdf.getVar(ncid,tmpvarid);
            netcdf.close(ncid);

            %% OBS
            cfg.obs_fnm=[dirs.obsroot, tmp.fs, cfg.obs_fname_mid,tmp.fy_str,tmp.mon_str, '.nc'];
            if exist(cfg.obs_fnm)~=0
                ncid=netcdf.open(cfg.obs_fnm, 'NOWRITE');
                tmpvarid=netcdf.inqVarID(ncid,cfg.obs_varname);
                tmp.ydata_obs(1:grid.nlon,1:grid.nlat,mon) = netcdf.getVar(ncid,tmpvarid);
                netcdf.close(ncid);
            else
                tmp.ydata_obs(1:grid.nlon,1:grid.nlat,mon) = NaN;                    
            end
%             tmp.ydata_obs(:,:,mon)=ncread(cfg.obs_fnm, cfg.obs_varname);

            %% ASSM
            cfg.assm_fnm=[dirs.assmdir, tmp.fs, tmp.varname, '_ensmean_', tmp.fy_str,'-', tmp.mon_str, '.nc'];
            if exist(cfg.assm_fnm)~=0
                ncid=netcdf.open(cfg.assm_fnm, 'NOWRITE');
                tmpvarid=netcdf.inqVarID(ncid,tmp.varname);
                tmp.ydata_assm(1:grid.nlon,1:grid.nlat,mon) = netcdf.getVar(ncid,tmpvarid);

                netcdf.close(ncid);
            else
                tmp.ydata_assm(1:grid.nlon,1:grid.nlat,mon) = NaN;                    
            end

            
%             if (strcmp(cfg.var, 'sumChl') || strcmp(cfg.var, 'TS')) %for Chls, it uses all available data
%                 tmp.obs_mask(:,:)=tmp.ydata_obs(:,:,mon)./tmp.ydata_obs(:,:,mon);
%                 tmp.ydata(:,:,mon) = tmp.ydata(:,:,mon) .* tmp.obs_mask(:,:); % masking with available observation
%             end

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
            tmp.ymean_obs= mean(tmp.ydata_obs,3);
            tmp.ymean_assm= mean(tmp.ydata_assm,3);
            tmp.ymean_assm(tmp.ymean_assm>10e30)=NaN;
            tmp.ymean_lens2= mean(tmp.ydata_lens2,3);   

            switch cfg.comp
                case 'ocn'
                    tmp.ymean = tmp.ymean .* grid.mask_ocn;
                    tmp.ymean_obs = tmp.ymean_obs .* grid.mask_ocn;
                    tmp.ymean_assm = tmp.ymean_assm .* grid.mask_ocn;
                    tmp.ymean_lens2 = tmp.ymean_lens2 .* grid.mask_ocn;                        
             end

%         end
        data.([tmp.varname, '_model', '_l', tmp.lyear_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean;
        data.([tmp.varname, '_obs'])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_obs;
        data.([tmp.varname, '_lens2', '_l', tmp.lyear_str])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_lens2;
        data.([tmp.varname, '_assm'])(1:grid.nlon,1:grid.nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_assm;

%         %% AR1 (assign)
%         tmp.ymean_AR1= mean(tmp.ydata_AR1, 3);
%         data.([tmp.varname, '_AR1', '_l', tmp.lyear_str])(:,:,iyear-min(cfg.iyears)+1+lyear)= tmp.ymean_AR1;
%         tmp.ymean= mean(ncread(cfg.datafilename, ['assm_', tmp.varname], [1, 1, 1], [inf, inf, 12]), 3);
%         data.([tmp.varname, '_assm'])(:,:,iyear-min(cfg.iyears)+1)= tmp.ymean;
    end
    
    %% get correlation coefficient
    for loni=1:grid.nlon
        for lati=1:grid.nlat
             if (isnan(data.([tmp.varname, '_assm'])(loni,lati,1))~=1 & nansum(data.([tmp.varname, '_model_l', tmp.lyear_str])(loni,lati,:))~=0)
                 tmp.data = squeeze(data.([tmp.varname, '_model', '_l', tmp.lyear_str])(loni,lati,:));
                 tmp.data_obs = data.([tmp.varname, '_obs'])(loni,lati,:);
                 [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));

                 tmp.data_det = Func_0028_detrend_linear_1d(squeeze(data.([tmp.varname, '_model', '_l', tmp.lyear_str])(loni,lati,:)), 'omitnan');
                 tmp.data_obs_det = Func_0028_detrend_linear_1d(data.([tmp.varname, '_obs'])(loni,lati,:));
                 [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data_det(isfinite(tmp.data_obs_det)), tmp.data_obs_det(isfinite(tmp.data_obs_det)));
%                  [tmp.corr, tmp.corr_p]=corrcoef(squeeze(data.([tmp.varname, '_model', '_l', tmp.lyear_str])(loni,lati,1+lyear:end-cfg.proj_year+1)), ...
%                      squeeze(data.([tmp.varname, '_obs'])(loni,lati,1+lyear:end)));
                 data.([tmp.varname, '_corr_obs', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
                 data.([tmp.varname, '_corr_obs_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);
                 data.([tmp.varname, '_corr_obs_det', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det(1,2);
                 data.([tmp.varname, '_corr_obs_det_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_det_p(1,2);
%                  [tmp.corr, tmp.corr_p]=corrcoef(squeeze(data.([tmp.varname, '_model', '_l', tmp.lyear_str])(loni,lati,1+lyear:end-cfg.proj_year+1)), ...
%                      squeeze(data.([tmp.varname, '_assm'])(loni,lati,1+lyear:end)));
%                  data.([tmp.varname, '_corr_assm', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
%                  data.([tmp.varname, '_corr_assm_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);


             %% corr assm
                 tmp.data_assm = data.([tmp.varname, '_assm'])(loni,lati,:);
                 [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_assm)), tmp.data_assm(isfinite(tmp.data_assm)));
                 
                 tmp.data_hcst_lens2 = squeeze(data.([tmp.varname, '_model', '_l', tmp.lyear_str])(loni,lati,:)) - ...
                                squeeze(data.([tmp.varname, '_lens2', '_l', tmp.lyear_str])(loni,lati,:));
                 [tmp.corr_int, tmp.corr_int_p]=corrcoef(tmp.data_hcst_lens2(isfinite(tmp.data_assm)), tmp.data_assm(isfinite(tmp.data_assm)));

                 data.([tmp.varname, '_corr_assm', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
                 data.([tmp.varname, '_corr_assm_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);
                 data.([tmp.varname, '_corr_assm_int', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_int(1,2);
                 data.([tmp.varname, '_corr_assm_int_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_int_p(1,2);

             %% corr lens2-assm
                 tmp.lens2 = squeeze(data.([tmp.varname, '_lens2', '_l', tmp.lyear_str])(loni,lati,:));
                 [tmp.corr, tmp.corr_p]=corrcoef(tmp.lens2(isfinite(tmp.data_assm)), tmp.data_assm(isfinite(tmp.data_assm)));

                 data.([tmp.varname, '_corr_assm_lens2', '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
                 data.([tmp.varname, '_corr_assm_lens2_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);



%                 %% AR1 (corr)
%                  tmp.data_AR1 = squeeze(data.([tmp.varname, '_AR1', '_l', tmp.lyear_str])(loni,lati,1+lyear:end-cfg.proj_year+1));
%                  [tmp.corr_AR1, tmp.corr_AR1_p]=corrcoef(tmp.data_AR1, tmp.data_obs);
%                  if (isnan(tmp.corr_AR1(1,2)))
%                     tmp.corr_AR1(1,2)=0;
%                  end
%                  data.([tmp.varname, '_corr_AR1', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_AR1(1,2);
%                  data.([tmp.varname, '_corr_AR1_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_AR1_p(1,2);
%                  
%                  tmp.data_AR1_det = Func_0028_detrend_linear_1d(data.([tmp.varname, '_AR1', '_l', tmp.lyear_str])(loni,lati,1+lyear:end-cfg.proj_year+1));
%                  [tmp.corr_AR1_det, tmp.corr_AR1_det_p]=corrcoef(tmp.data_AR1_det, tmp.data_obs_det);
%                  if (isnan(tmp.corr_AR1_det(1,2)))
%                     tmp.corr_AR1_det(1,2)=0;
%                  end
%                  data.([tmp.varname, '_corr_AR1_det', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_AR1_det(1,2);
%                  data.([tmp.varname, '_corr_AR1_det_p', '_l', tmp.lyear_str])(loni,lati)=tmp.corr_AR1_det_p(1,2);

             else
                 data.([tmp.varname, '_corr_obs', '_l', tmp.lyear_str])(loni,lati)=NaN;
                 data.([tmp.varname, '_corr_obs_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                 data.([tmp.varname, '_corr_obs_det', '_l', tmp.lyear_str])(loni,lati)=NaN;
                 data.([tmp.varname, '_corr_obs_det_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
%                  data.([tmp.varname, '_corr_assm', '_l', tmp.lyear_str])(loni,lati)=NaN;
%                  data.([tmp.varname, '_corr_assm_p', '_l', tmp.lyear_str])(loni,lati)=NaN;

%                 %% AR1 (NaN)
%                  data.([tmp.varname, '_corr_AR1', '_l', tmp.lyear_str])(loni,lati)=NaN;
%                  data.([tmp.varname, '_corr_AR1_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
%                  data.([tmp.varname, '_corr_AR1_det', '_l', tmp.lyear_str])(loni,lati)=NaN;
%                  data.([tmp.varname, '_corr_AR1_det_p', '_l', tmp.lyear_str])(loni,lati)=NaN;

                %% ASSM (NaN)
                 data.([tmp.varname, '_corr_assm', '_l', tmp.lyear_str])(loni,lati)=NaN;
                 data.([tmp.varname, '_corr_assm_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
                 data.([tmp.varname, '_corr_assm_int', '_l', tmp.lyear_str])(loni,lati)=NaN;
                 data.([tmp.varname, '_corr_assm_int_p', '_l', tmp.lyear_str])(loni,lati)=NaN;

                 %% LENS2 (NaN)
                 data.([tmp.varname, '_corr_assm_lens2', '_l', tmp.lyear_str])(loni,lati)=NaN;
                 data.([tmp.varname, '_corr_assm_lens2_p', '_l', tmp.lyear_str])(loni,lati)=NaN;
             end
        end
    end
    disp('abc')
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


% tmp.yearset={'01', '02', '3_4', '5_9'};
tmp.yearset={'01', '02', '03', '04'};
for lyear_ind=1:length(tmp.yearset)
    tmp.lyear_str=tmp.yearset{lyear_ind};

% % %     %% model & obs corr map  ----------------------
% % %     fig_cfg.name_rgn = 'Glob';
% % %     fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
% % %     fig_cfg.x_lim = [-180 180];
% % %     fig_cfg.y_lim = [-80 89];
% % %     fig_cfg.fig_size = [0,0,6,3.5];
% % %     fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
% % %     fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
% % %     fig_cfg.title_pos = [0.5,0.93];
% % %     fig_cfg.p_lim =0.1;
% % %     fig_cfg.c_lim = [-1 1];
% % %     [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
% % % 
% % %     tmp.X=grid.tlong([end, 1:end],:);
% % %     tmp.Y=grid.tlat([end, 1:end],:);
% % %     tmp.C=data.([tmp.varname, '_corr_obs', '_l', tmp.lyear_str]);
% % %     tmp.C=tmp.C([end, 1:end],:);
% % %     tmp.C_H=data.([tmp.varname, '_corr_AR1', '_l', tmp.lyear_str]);
% % %     tmp.C_H=tmp.C_H([end, 1:end],:);
% % %     tmp.C_2=tmp.C;
% % %     tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 
% % %     
% % %     [tmp.mean_corr, tmp.err] = ...
% % %         Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_obs', '_l', tmp.lyear_str]), grid.tlong, grid.tlat);
% % %     fig_cfg.fig_name=['l',tmp.lyear_str, ', corr_ob_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
% % %     fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
% % %         'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
% % %     %% map setting
% % %     ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
% % %         'fontname','freeserif'); 
% % % 
% % %     axis off; 
% % %     hold on;
% % %     setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
% % %     set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
% % %     text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
% % %     'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
% % %     'fontsize',14,'fontname','freeserif','interpreter','none')
% % % 
% % %     %% caxis & colorbar
% % %     caxis(ax_m, fig_cfg.c_lim); 
% % %     colormap(fig_cfg.c_map);
% % %     cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
% % %     set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
% % %     title(cb,'R','fontsize',12);
% % % 
% % %     %% draw on ax_m
% % %     h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
% % %     shading flat;
% % %     geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
% % % 
% % %     %% <AR1 area -> hatch
% % %     if (strcmp(cfg.var, 'sumChl')~=1) %for Chls, it uses all available data
% % %         pp2 = pcolorm(tmp.Y,tmp.X,tmp.C_2, 'parent', ax_m);
% % %         set(pp2,'linestyle','none','Tag','HatchingRegion');
% % %         hp = findobj(pp2,'Tag','HatchingRegion');
% % %         hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);
% % %     end
% % %     %% frame and label setting
% % %     setm(ax_m,'frame','on','FLineWidth',1);
% % % 
% % %     label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
% % %     label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
% % %     mlabel; plabel;
% % %     label_y=plabel; label_x=mlabel;
% % %     for lxi=1:length(label_x)
% % %         tmp.tmppos=label_x(lxi,1).Position;
% % %         tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
% % %         label_x(lxi,1).Position=tmp.tmppos;
% % %         label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
% % %     end
% % %     for lyi=1:length(label_y)
% % %         label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
% % %     end
% % % 
% % %     %% save
% % %     dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_map', filesep, 'model'];
% % %     if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
% % %     cfg.figname=[dirs.figdir, filesep, 'corr_obs_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
% % %     print(fig_h, cfg.figname, '-dpng');
% % %     RemoveWhiteSpace([], 'file', cfg.figname);
% % %     fprintf('%7.1f sec\n', toc(lap_time) );
% % %     close all;
% % % 
% % % %% model & obs corr map (detrended)  ----------------------
% % %     fig_cfg.name_rgn = 'Glob';
% % %     fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
% % %     fig_cfg.x_lim = [-180 180];
% % %     fig_cfg.y_lim = [-80 89];
% % %     fig_cfg.fig_size = [0,0,6,3.5];
% % %     fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
% % %     fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
% % %     fig_cfg.title_pos = [0.5,0.93];
% % %     fig_cfg.p_lim =0.1;
% % %     fig_cfg.c_lim = [-1 1];
% % %     [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
% % % 
% % %     tmp.X=grid.tlong([end, 1:end],:);
% % %     tmp.Y=grid.tlat([end, 1:end],:);
% % %     tmp.C=data.([tmp.varname, '_corr_obs_det', '_l', tmp.lyear_str]);
% % %     tmp.C=tmp.C([end, 1:end],:);
% % %     tmp.C_H=data.([tmp.varname, '_corr_AR1_det', '_l', tmp.lyear_str]);
% % %     tmp.C_H=tmp.C_H([end, 1:end],:);
% % %     tmp.C_2=tmp.C;
% % %     tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 
% % % 
% % %     [tmp.mean_corr, tmp.err] = ...
% % %         Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_obs_det', '_l', tmp.lyear_str]), grid.tlong, grid.tlat);
% % %     fig_cfg.fig_name=['l',tmp.lyear_str, ', corr_ob_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
% % %     fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
% % %         'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
% % %     %% map setting
% % %     ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
% % %         'fontname','freeserif'); 
% % % 
% % %     axis off; 
% % %     hold on;
% % %     setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
% % %     set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
% % %     text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
% % %     'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
% % %     'fontsize',14,'fontname','freeserif','interpreter','none')
% % % 
% % %     %% caxis & colorbar
% % %     caxis(ax_m, fig_cfg.c_lim); 
% % %     colormap(fig_cfg.c_map);
% % %     cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
% % %     set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
% % %     title(cb,'R','fontsize',12);
% % % 
% % %     %% draw on ax_m
% % %     h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
% % %     shading flat;
% % %     geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
% % % 
% % %     %% <AR1 area -> hatch
% % %     if (strcmp(cfg.var, 'sumChl')~=1) %for Chls, it uses all available data
% % %         pp2 = pcolorm(tmp.Y,tmp.X,tmp.C_2, 'parent', ax_m);
% % %         set(pp2,'linestyle','none','Tag','HatchingRegion');
% % %         hp = findobj(pp2,'Tag','HatchingRegion');
% % %         hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','k','HatchLineWidth',0.5);
% % %     end
% % %     %% frame and label setting
% % %     setm(ax_m,'frame','on','FLineWidth',1);
% % % 
% % %     label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
% % %     label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
% % %     mlabel; plabel;
% % %     label_y=plabel; label_x=mlabel;
% % %     for lxi=1:length(label_x)
% % %         tmp.tmppos=label_x(lxi,1).Position;
% % %         tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
% % %         label_x(lxi,1).Position=tmp.tmppos;
% % %         label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
% % %     end
% % %     for lyi=1:length(label_y)
% % %         label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
% % %     end
% % % 
% % %     %% save
% % %     dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_map', filesep, 'model_det'];
% % %     if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
% % %     cfg.figname=[dirs.figdir, filesep, 'corr_obs_det_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
% % %     print(fig_h, cfg.figname, '-dpng');
% % %     RemoveWhiteSpace([], 'file', cfg.figname);
% % %     fprintf('%7.1f sec\n', toc(lap_time) );
% % %     close all;


% % %     %% AR1 & obs corr map  ----------------------
% % %     fig_cfg.name_rgn = 'Glob';
% % %     fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
% % %     fig_cfg.x_lim = [-180 180];
% % %     fig_cfg.y_lim = [-80 89];
% % %     fig_cfg.fig_size = [0,0,6,3.5];
% % %     fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
% % %     fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
% % %     fig_cfg.title_pos = [0.5,0.93];
% % %     fig_cfg.p_lim =0.1;
% % %     fig_cfg.c_lim = [-1 1];
% % %     [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
% % % 
% % %     tmp.X=grid.tlong([end, 1:end],:);
% % %     tmp.Y=grid.tlat([end, 1:end],:);
% % %     tmp.C=data.([tmp.varname, '_corr_AR1', '_l', tmp.lyear_str]);
% % %     tmp.C=tmp.C([end, 1:end],:);
% % %     
% % %     [tmp.mean_corr, tmp.err] = ...
% % %         Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_AR1', '_l', tmp.lyear_str]), grid.tlong, grid.tlat);
% % %     fig_cfg.fig_name=['l',tmp.lyear_str, ', corr_ob_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
% % %     fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
% % %         'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
% % %     %% map setting
% % %     ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
% % %         'fontname','freeserif'); 
% % % 
% % %     axis off; 
% % %     hold on;
% % %     setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
% % %     set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
% % %     text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
% % %     'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
% % %     'fontsize',14,'fontname','freeserif','interpreter','none')
% % % 
% % %     %% caxis & colorbar
% % %     caxis(ax_m, fig_cfg.c_lim); 
% % %     colormap(fig_cfg.c_map);
% % %     cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
% % %     set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
% % %     title(cb,'R','fontsize',12);
% % % 
% % %     %% draw on ax_m
% % %     h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
% % %     shading flat;
% % %     geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
% % % 
% % %     %% frame and label setting
% % %     setm(ax_m,'frame','on','FLineWidth',1);
% % % 
% % %     label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
% % %     label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
% % %     mlabel; plabel;
% % %     label_y=plabel; label_x=mlabel;
% % %     for lxi=1:length(label_x)
% % %         tmp.tmppos=label_x(lxi,1).Position;
% % %         tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
% % %         label_x(lxi,1).Position=tmp.tmppos;
% % %         label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
% % %     end
% % %     for lyi=1:length(label_y)
% % %         label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
% % %     end
% % % 
% % %     %% save
% % %     dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_map', filesep, 'AR1'];
% % %     if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
% % %     cfg.figname=[dirs.figdir, filesep, 'corr_AR1_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
% % %     print(fig_h, cfg.figname, '-dpng');
% % %     RemoveWhiteSpace([], 'file', cfg.figname);
% % %     fprintf('%7.1f sec\n', toc(lap_time) );
% % %     close all;
% % % 
% % % 
% % %     %% AR1 & obs corr map (detrended)  ----------------------
% % %     fig_cfg.name_rgn = 'Glob';
% % %     fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
% % %     fig_cfg.x_lim = [-180 180];
% % %     fig_cfg.y_lim = [-80 89];
% % %     fig_cfg.fig_size = [0,0,6,3.5];
% % %     fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
% % %     fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
% % %     fig_cfg.title_pos = [0.5,0.93];
% % %     fig_cfg.p_lim =0.1;
% % %     fig_cfg.c_lim = [-1 1];
% % %     [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
% % % 
% % %     tmp.X=grid.tlong([end, 1:end],:);
% % %     tmp.Y=grid.tlat([end, 1:end],:);
% % %     tmp.C=data.([tmp.varname, '_corr_AR1_det', '_l', tmp.lyear_str]);
% % %     tmp.C=tmp.C([end, 1:end],:);
% % % 
% % % 
% % %     [tmp.mean_corr, tmp.err] = ...
% % %         Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_AR1_det', '_l', tmp.lyear_str]), grid.tlong, grid.tlat);
% % %     fig_cfg.fig_name=['l',tmp.lyear_str, ', corr_ob_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
% % %     fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
% % %         'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
% % %     %% map setting
% % %     ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
% % %         'fontname','freeserif'); 
% % % 
% % %     axis off; 
% % %     hold on;
% % %     setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
% % %     set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
% % %     text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
% % %     'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
% % %     'fontsize',14,'fontname','freeserif','interpreter','none')
% % % 
% % %     %% caxis & colorbar
% % %     caxis(ax_m, fig_cfg.c_lim); 
% % %     colormap(fig_cfg.c_map);
% % %     cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
% % %     set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
% % %     title(cb,'R','fontsize',12);
% % % 
% % %     %% draw on ax_m
% % %     h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
% % %     shading flat;
% % %     geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
% % % 
% % %     %% frame and label setting
% % %     setm(ax_m,'frame','on','FLineWidth',1);
% % % 
% % %     label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
% % %     label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
% % %     mlabel; plabel;
% % %     label_y=plabel; label_x=mlabel;
% % %     for lxi=1:length(label_x)
% % %         tmp.tmppos=label_x(lxi,1).Position;
% % %         tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
% % %         label_x(lxi,1).Position=tmp.tmppos;
% % %         label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
% % %     end
% % %     for lyi=1:length(label_y)
% % %         label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
% % %     end
% % % 
% % %     %% save
% % %     dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_map', filesep, 'AR1_det'];
% % %     if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
% % %     cfg.figname=[dirs.figdir, filesep, 'corr_AR1_det_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
% % %     print(fig_h, cfg.figname, '-dpng');
% % %     RemoveWhiteSpace([], 'file', cfg.figname);
% % %     fprintf('%7.1f sec\n', toc(lap_time) );
% % %     close all;


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
    tmp.C=data.([tmp.varname, '_corr_assm', '_l', tmp.lyear_str]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data.([tmp.varname, '_corr_AR1', '_l', tmp.lyear_str]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_assm', '_l', tmp.lyear_str]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=['l',tmp.lyear_str, ', corr_assm_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    cfg.figname=[dirs.figdir, filesep, 'corr_assm_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
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
    tmp.C=data.([tmp.varname, '_corr_assm_int', '_l', tmp.lyear_str]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data.([tmp.varname, '_corr_AR1', '_l', tmp.lyear_str]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_assm_int', '_l', tmp.lyear_str]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=['l',tmp.lyear_str, ', corr_assm _, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    cfg.figname=[dirs.figdir, filesep, 'corr_assm_int_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];

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
    tmp.C=data.([tmp.varname, '_corr_assm_lens2', '_l', tmp.lyear_str]);
    tmp.C=tmp.C([end, 1:end],:);
%         tmp.C_H=data.([tmp.varname, '_corr_AR1', '_l', tmp.lyear_str]);
%         tmp.C_H=tmp.C_H([end, 1:end],:);
%         tmp.C_2=tmp.C;
%         tmp.C_2(tmp.C>tmp.C_H)=NaN;  % if AR1 corr > model = NaN 


    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_assm_lens2', '_l', tmp.lyear_str]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=['l',tmp.lyear_str, ', corr_assm_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
    cfg.figname=[dirs.figdir, filesep, 'corr_assm_ext_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    fprintf('%7.1f sec\n', toc(lap_time) );
    close all;
    
    



end
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