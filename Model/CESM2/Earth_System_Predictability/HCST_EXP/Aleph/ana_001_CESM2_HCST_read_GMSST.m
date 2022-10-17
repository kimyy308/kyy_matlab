% %  Created 05-Oct-2022 by Yong-Yub Kim
clc; clear all; close all;

%% set path
tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

%% model configuration
dirs.root='/mnt/lustre/proj/earth.system.predictability/HCST_EXP';
dirs.yoshi_root='/proj/yoshi/DATA/CESM2_ODA';
dirs.archive=[dirs.root, filesep, 'archive'];
dirs.saveroot='/mnt/lustre/proj/kimyy/Model/CESM2/ESP/HCST_EXP';

config.iyears=1992:2021;
config.months=1:12;
config.scenname='HIST';
config.gridname='f09_g17';
config.assm_factor='10';
config.ens_member='1';
config.proj_year=10;

% config.obsnames={'oras4', 'projdv7.3', 'en4.2'};
config.obsnames={'en4.2_ba'};
% config.obsnames_simple={'en4'};
config.ensnames={'ba-10p1'};

% config.ensnames={'ba-10p1', 'ba-10p2', 'ba-10p3', 'ba-10p4', 'ba-10p5', ...
%     'ba-20p1', 'ba-20p2', 'ba-20p3', 'ba-20p4', 'ba-20p5'};

config.component='ocn';
config.varnames={'TEMP','SALT'};

config.len_t_y = length(config.iyears);
config.len_t_m = length(config.months);
config.len_obs= length(config.obsnames);
config.len_ens= length(config.ensnames);

%% observation configuration
dirs.obs_en4root = '/mnt/lustre/proj/earth.system.predictability/OBS_DATA/POP2_g17_BA/EN4.2';
dirs.obs_projdroot = '/mnt/lustre/proj/earth.system.predictability/OBS_DATA/POP2_g17_BA/ProjD_v7.3';


%% grid set(lon, lat from obs)
% tosoerr.1957-11.nc
% time (day from 0000)
% temp, salt, tlong, tlat
config.obs_en4_filename=[dirs.obs_en4root, filesep, ...
                'tosoerr.', num2str(min(config.iyears)), '-', num2str(01, '%02i'), '.nc'];
grid.tlong=ncread(config.obs_en4_filename, 'TLONG');
grid.tlat=ncread(config.obs_en4_filename, 'TLAT');
grid.nlon=size(grid.tlong,1);
grid.nlat=size(grid.tlong,2);
grid.ntime=config.proj_year.*12;


%% grid set(mask from model)
tmp.obsname=config.obsnames{1};
iyear=min(config.iyears);
config.casename_m=[config.gridname, '.hcst.', tmp.obsname, '-', config.assm_factor, 'p', config.ens_member];
config.casename=[config.casename_m, '_i', num2str(iyear)];
dirs.datadir= [dirs.archive, filesep, config.casename_m, filesep, config.casename, ...
    filesep, 'ocn/hist'];

[tmp.error_status, tmp.value]=system(['ls ', dirs.datadir, '/*once*']);  % b.e21.BHISTsmbb.f09_g17.assm.oras4_ba-10p1.pop.h.once.nc
tmp.gridname = [tmp.value(1:end-1)];
grid.region_mask=ncread(tmp.gridname, 'REGION_MASK'); 
grid.ocean_mask=NaN(size(grid.region_mask));
grid.ocean_mask(grid.region_mask>0)=1;
grid.tarea = ncread(tmp.gridname, 'TAREA');

% % model filename example
% /mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive/
% f09_g17.hcst.en4.2_ba-10p1/f09_g17.hcst.en4.2_ba-10p1_i1993/ocn/hist/
% f09_g17.hcst.en4.2_ba-10p1_i1993.pop.h.1993-01.nc

%% read & save data
for obsind=1:length(config.obsnames)
    tmp.obsname=config.obsnames{obsind};
    tmp.obsname_simple= f_obs_simple(tmp.obsname);
    for iyear=min(config.iyears):max(config.iyears)
        tmp.iyear_str=num2str(iyear, '%04i');
        config.casename_m=[config.gridname, '.hcst.', tmp.obsname, '-', config.assm_factor, 'p', config.ens_member];
        config.casename=[config.casename_m, '_i', num2str(iyear)];
        dirs.datadir= [dirs.archive, filesep, config.casename_m, filesep, config.casename, ...
            filesep, 'ocn/hist'];
        
        for varind=1:length(config.varnames)
            tmp.varname=config.varnames{varind};
            grid.(['time_i', tmp.iyear_str])=NaN(1,grid.ntime);
            %% variables initialization
            data.([tmp.varname, '_model_', tmp.obsname_simple, '_i', tmp.iyear_str])=NaN(grid.nlon, grid.nlat, grid.ntime);
            data.([tmp.varname, '_bias_', tmp.obsname_simple, '_i', tmp.iyear_str])=NaN(grid.nlon, grid.nlat, grid.ntime);
            data.([tmp.varname, '_obs_', tmp.obsname_simple])=NaN(grid.nlon, grid.nlat, grid.ntime);
            data.([tmp.varname, '_sq_err_', tmp.obsname_simple, '_i', tmp.iyear_str])=NaN(grid.nlon, grid.nlat, grid.ntime);
            data.([tmp.varname, '_rmse_', tmp.obsname_simple, '_i', tmp.iyear_str])=NaN(grid.nlon, grid.nlat);
            data.(['gm_', tmp.varname, '_model_', tmp.obsname_simple, '_i', tmp.iyear_str])=NaN(1,grid.ntime);
            data.(['gm_', tmp.varname, '_bias_', tmp.obsname_simple, '_i', tmp.iyear_str])=NaN(1,grid.ntime);
            data.(['gm_', tmp.varname, '_obs_', tmp.obsname_simple, '_i', tmp.iyear_str])=NaN(1,grid.ntime);
            data.(['gm_', tmp.varname, '_sq_err_', tmp.obsname_simple, '_i', tmp.iyear_str])=NaN(1,grid.ntime);
            data.time=NaN(1,grid.ntime);
    
            for myear=iyear:iyear+config.proj_year-1
                for month=1:12
                    config.(['obs_',tmp.obsname_simple,'_filename'])=[dirs.(['obs_',tmp.obsname_simple,'root']), filesep, ...
                        'tosoerr.', num2str(myear), '-', num2str(month, '%02i'), '.nc'];
                    config.(['cesm2_',tmp.obsname_simple,'_filename']) = [dirs.datadir, filesep, ...
                        config.casename, '.pop.h.', num2str(myear), '-', num2str(month, '%02i'), '.nc'];
                    tmp.ind_var=(myear-iyear)*12+month;
                    data.([tmp.varname, '_model_', tmp.obsname_simple, '_i', tmp.iyear_str])(:,:,tmp.ind_var) = ...
                        ncread(config.(['cesm2_',tmp.obsname_simple,'_filename']), tmp.varname, [1 1 1 1], [inf inf 1 1]) ...
                        .* grid.ocean_mask;
                    if (exist(config.(['obs_',tmp.obsname_simple,'_filename']))==2)
                        data.([tmp.varname, '_obs_', tmp.obsname_simple])(:,:,tmp.ind_var) = ...
                            ncread(config.(['obs_',tmp.obsname_simple,'_filename']), tmp.varname, [1 1 1 1], [inf inf 1 1])...
                            .* grid.ocean_mask;
                    end
                    data.time(tmp.ind_var) = ncread(config.(['cesm2_',tmp.obsname_simple,'_filename']), 'time');
%                     disp('abc')
%                     pcolor(data.([tmp.varname, '_model_', tmp.obsname_simple, '_i', tmp.iyear_str])(:,:,tmp.ind_var)');
%                     shading flat; colorbar;
%                     pcolor(data.([tmp.varname, '_obs_', tmp.obsname_simple])(:,:,tmp.ind_var)');
%                     shading flat; colorbar;
                    data.([tmp.varname, '_bias_', tmp.obsname_simple, '_i', tmp.iyear_str])(:,:,tmp.ind_var) = ...
                        data.([tmp.varname, '_model_', tmp.obsname_simple, '_i', tmp.iyear_str])(:,:,tmp.ind_var) - ...
                        data.([tmp.varname, '_obs_', tmp.obsname_simple])(:,:,tmp.ind_var);
                    data.([tmp.varname, '_sq_err_', tmp.obsname_simple, '_i', tmp.iyear_str])(:,:,tmp.ind_var) = ...
                        data.([tmp.varname, '_bias_', tmp.obsname_simple, '_i', tmp.iyear_str])(:,:,tmp.ind_var).^2;
                    data.(['gm_', tmp.varname, '_model_', tmp.obsname_simple, '_i', tmp.iyear_str])(tmp.ind_var)= ...
                        f_gm_var(data.([tmp.varname, '_model_', tmp.obsname_simple, '_i', tmp.iyear_str])(:,:,tmp.ind_var), ...
                        grid.tarea);
                    data.(['gm_', tmp.varname, '_obs_', tmp.obsname_simple])(tmp.ind_var)= ...
                        f_gm_var(data.([tmp.varname, '_obs_', tmp.obsname_simple])(:,:,tmp.ind_var), ...
                        grid.tarea);
                    data.(['gm_', tmp.varname, '_bias_', tmp.obsname_simple, '_i', tmp.iyear_str])(tmp.ind_var)= ...
                        f_gm_var(data.([tmp.varname, '_bias_', tmp.obsname_simple, '_i', tmp.iyear_str])(:,:,tmp.ind_var), ...
                        grid.tarea);
                    data.(['gm_', tmp.varname, '_sq_err_', tmp.obsname_simple, '_i', tmp.iyear_str])(tmp.ind_var)= ...
                        f_gm_var(data.([tmp.varname, '_sq_err_', tmp.obsname_simple, '_i', tmp.iyear_str])(:,:,tmp.ind_var), ...
                        grid.tarea);
                end
            end
            data.([tmp.varname, '_rmse_', tmp.obsname_simple, '_i', tmp.iyear_str]) = ...
                sqrt(sum(data.([tmp.varname, '_sq_err_', tmp.obsname_simple, '_i', tmp.iyear_str]), 3) ./ grid.ntime);
%             pcolor(data.([tmp.varname, '_rmse_', tmp.obsname_simple, '_i', tmp.iyear_str])'); shading flat; colorbar;
            data.(['gm_', tmp.varname, '_rmse_', tmp.obsname_simple, '_i', tmp.iyear_str])= ...
                        f_gm_var(data.([tmp.varname, '_rmse_', tmp.obsname_simple, '_i', tmp.iyear_str])(:,:), ...
                        grid.tarea);
        end
%% save ncfile
        dirs.saveroot='/mnt/lustre/proj/kimyy/Model/CESM2/ESP/HCST_EXP';
        dirs.savedir=[dirs.saveroot, filesep, config.casename_m, filesep, 'GMSV'];
        mkdir(dirs.savedir);
        config.savefilename=[dirs.savedir, filesep, 'GMSV_', config.casename, '.nc'];
        
        if (exist(config.savefilename)==0)
            ncid = netcdf.create(config.savefilename,'NETCDF4');
    
            lon_dimid = netcdf.defDim(ncid, 'nlon', grid.nlon);
            lat_dimid = netcdf.defDim(ncid,'nlat', grid.nlat);
            time_dimid = netcdf.defDim(ncid, 'time', 0);
            one_dimid = netcdf.defDim(ncid, 'one', 1);
    
            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                'title', ['CESM2 Hindcast ', 'surface variables']);
            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                'source', ['CESM2, ',config.casename]);
            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                'initialized year', tmp.iyear_str);
            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                'author', 'Created by Y.Y. Kim');
            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                'date', date);
            
            timevarid=netcdf.defVar(ncid, 'time', 'NC_DOUBLE', time_dimid);
                netcdf.putAtt(ncid,timevarid,'long_name','time');
                netcdf.putAtt(ncid,timevarid,'units','days since 0000-01-01 00:00:00');
    %             netcdf.putAtt(ncid,timevarid,'bounds','time_bound');           
                netcdf.putAtt(ncid,timevarid,'calendar','noleap');           
            
            tlongvarid=netcdf.defVar(ncid, 'TLONG', 'NC_DOUBLE', [lon_dimid lat_dimid]);
                netcdf.putAtt(ncid,tlongvarid,'long_name','array of t-grid longitudes');
                netcdf.putAtt(ncid,tlongvarid,'units','degree_east');
                netcdf.defVarChunking(ncid, tlongvarid, 'CHUNKED', [grid.nlon/10, grid.nlat/10]);
                netcdf.defVarDeflate(ncid, tlongvarid, true, true, 1);
            tlatvarid=netcdf.defVar(ncid, 'TLAT', 'NC_DOUBLE', [lon_dimid lat_dimid]);
                netcdf.putAtt(ncid,tlatvarid,'long_name','array of t-grid latitudes');
                netcdf.putAtt(ncid,tlatvarid,'units','degree_north');
                netcdf.defVarChunking(ncid, tlatvarid, 'CHUNKED', [grid.nlon/10, grid.nlat/10]);
                netcdf.defVarDeflate(ncid, tlatvarid, true, true, 1);
    
            tempvarid=netcdf.defVar(ncid, 'TEMP', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                netcdf.putAtt(ncid,tempvarid,'long_name','Potential temperature (surface)');
                netcdf.putAtt(ncid,tempvarid,'units','degC');  
                netcdf.defVarChunking(ncid, tempvarid, 'CHUNKED', [grid.nlon/10, grid.nlat/10 1]);
                netcdf.defVarDeflate(ncid, tempvarid, true, true, 1);
            obs_tempvarid=netcdf.defVar(ncid, 'obs_temp', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                netcdf.putAtt(ncid,obs_tempvarid,'long_name','obs_temperature (surface)');
                netcdf.putAtt(ncid,obs_tempvarid,'units','degC');
                netcdf.putAtt(ncid,obs_tempvarid,'source', tmp.obsname);  
                netcdf.defVarChunking(ncid, obs_tempvarid, 'CHUNKED', [grid.nlon/10, grid.nlat/10 1]);
                netcdf.defVarDeflate(ncid, obs_tempvarid, true, true, 1);
            bias_tempvarid=netcdf.defVar(ncid, 'bias_temp', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                netcdf.putAtt(ncid,bias_tempvarid,'long_name','bias_temperature');
                netcdf.putAtt(ncid,bias_tempvarid,'units','degC');  
                netcdf.defVarChunking(ncid, bias_tempvarid, 'CHUNKED', [grid.nlon/10, grid.nlat/10 1]);
                netcdf.defVarDeflate(ncid, bias_tempvarid, true, true, 1);
            sq_err_tempvarid=netcdf.defVar(ncid, 'sq_err_temp', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                netcdf.putAtt(ncid,sq_err_tempvarid,'long_name','sq_err_temperature');
                netcdf.putAtt(ncid,sq_err_tempvarid,'units','degC^2');  
                netcdf.defVarChunking(ncid, sq_err_tempvarid, 'CHUNKED', [grid.nlon/10, grid.nlat/10 1]);
                netcdf.defVarDeflate(ncid, sq_err_tempvarid, true, true, 1);
            rmse_tempvarid=netcdf.defVar(ncid, 'rmse_temp', 'NC_FLOAT', [lon_dimid lat_dimid]);
                netcdf.putAtt(ncid,rmse_tempvarid,'long_name','rmse_temperature');
                netcdf.putAtt(ncid,rmse_tempvarid,'units','degC');  
                netcdf.defVarChunking(ncid, rmse_tempvarid, 'CHUNKED', [grid.nlon/10, grid.nlat/10]);
                netcdf.defVarDeflate(ncid, rmse_tempvarid, true, true, 1);
            gm_tempvarid=netcdf.defVar(ncid, 'gm_temp', 'NC_FLOAT', [time_dimid]);
                netcdf.putAtt(ncid,gm_tempvarid,'long_name','global mean Potential temperature (surface)');
                netcdf.putAtt(ncid,gm_tempvarid,'units','degC');
            gm_obs_tempvarid=netcdf.defVar(ncid, 'gm_obs_temp', 'NC_FLOAT', [time_dimid]);
                netcdf.putAtt(ncid,gm_obs_tempvarid,'long_name','global mean obs_temperature (surface)');
                netcdf.putAtt(ncid,gm_obs_tempvarid,'units','degC');
            gm_bias_tempvarid=netcdf.defVar(ncid, 'gm_bias_temp', 'NC_FLOAT', [time_dimid]);
                netcdf.putAtt(ncid,gm_bias_tempvarid,'long_name','global mean bias_temperature (surface)');
                netcdf.putAtt(ncid,gm_bias_tempvarid,'units','degC');
            gm_sq_err_tempvarid=netcdf.defVar(ncid, 'gm_sq_err_temp', 'NC_FLOAT', [time_dimid]);
                netcdf.putAtt(ncid,gm_sq_err_tempvarid,'long_name','global mean sq_err_temperature (surface)');
                netcdf.putAtt(ncid,gm_sq_err_tempvarid,'units','degC^2');
            gm_rmse_tempvarid=netcdf.defVar(ncid, 'gm_rmse_temp', 'NC_FLOAT', [one_dimid]);
                netcdf.putAtt(ncid,gm_rmse_tempvarid,'long_name','global mean rmse_temperature (surface)');
                netcdf.putAtt(ncid,gm_rmse_tempvarid,'units','degC');
    
            saltvarid=netcdf.defVar(ncid, 'SALT', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                netcdf.putAtt(ncid,saltvarid,'long_name','Salinity (surface)');
                netcdf.putAtt(ncid,saltvarid,'units','gram/kilogram');  
                netcdf.defVarChunking(ncid, saltvarid, 'CHUNKED', [grid.nlon/10, grid.nlat/10 1]);
                netcdf.defVarDeflate(ncid, saltvarid, true, true, 1);
            obs_saltvarid=netcdf.defVar(ncid, 'obs_salt', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                netcdf.putAtt(ncid,obs_saltvarid,'long_name','obs_salinity (surface)');
                netcdf.putAtt(ncid,obs_saltvarid,'units','gram/kilogram');
                netcdf.putAtt(ncid,obs_saltvarid,'source', tmp.obsname);  
                netcdf.defVarChunking(ncid, obs_saltvarid, 'CHUNKED', [grid.nlon/10, grid.nlat/10 1]);
                netcdf.defVarDeflate(ncid, obs_saltvarid, true, true, 1);
            bias_saltvarid=netcdf.defVar(ncid, 'bias_salt', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                netcdf.putAtt(ncid,bias_saltvarid,'long_name','bias_salinity');
                netcdf.putAtt(ncid,bias_saltvarid,'units','gram/kilogram');  
                netcdf.defVarChunking(ncid, bias_saltvarid, 'CHUNKED', [grid.nlon/10, grid.nlat/10 1]);
                netcdf.defVarDeflate(ncid, bias_saltvarid, true, true, 1);
            sq_err_saltvarid=netcdf.defVar(ncid, 'sq_err_salt', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                netcdf.putAtt(ncid,sq_err_saltvarid,'long_name','sq_err_salinity');
                netcdf.putAtt(ncid,sq_err_saltvarid,'units','(gram/kilogram)^2');  
                netcdf.defVarChunking(ncid, sq_err_saltvarid, 'CHUNKED', [grid.nlon/10, grid.nlat/10 1]);
                netcdf.defVarDeflate(ncid, sq_err_saltvarid, true, true, 1);
            rmse_saltvarid=netcdf.defVar(ncid, 'rmse_salt', 'NC_FLOAT', [lon_dimid lat_dimid]);
                netcdf.putAtt(ncid,rmse_saltvarid,'long_name','rmse_salterature');
                netcdf.putAtt(ncid,rmse_saltvarid,'units','gram/kilogram');  
                netcdf.defVarChunking(ncid, rmse_saltvarid, 'CHUNKED', [grid.nlon/10, grid.nlat/10]);
                netcdf.defVarDeflate(ncid, rmse_saltvarid, true, true, 1);
            gm_saltvarid=netcdf.defVar(ncid, 'gm_salt', 'NC_FLOAT', [time_dimid]);
                netcdf.putAtt(ncid,gm_saltvarid,'long_name','global mean Salinity (surface)');
                netcdf.putAtt(ncid,gm_saltvarid,'units','gram/kilogram');
            gm_obs_saltvarid=netcdf.defVar(ncid, 'gm_obs_salt', 'NC_FLOAT', [time_dimid]);
                netcdf.putAtt(ncid,gm_obs_saltvarid,'long_name','global mean obs_salinity (surface)');
                netcdf.putAtt(ncid,gm_obs_saltvarid,'units','gram/kilogram');
            gm_bias_saltvarid=netcdf.defVar(ncid, 'gm_bias_salt', 'NC_FLOAT', [time_dimid]);
                netcdf.putAtt(ncid,gm_bias_saltvarid,'long_name','global mean bias_salinity (surface)');
                netcdf.putAtt(ncid,gm_bias_saltvarid,'units','gram/kilogram');
            gm_sq_err_saltvarid=netcdf.defVar(ncid, 'gm_sq_err_salt', 'NC_FLOAT', [time_dimid]);
                netcdf.putAtt(ncid,gm_sq_err_saltvarid,'long_name','global mean sq_err_salinity (surface)');
                netcdf.putAtt(ncid,gm_sq_err_saltvarid,'units','(gram/kilogram)^2');
            gm_rmse_saltvarid=netcdf.defVar(ncid, 'gm_rmse_salt', 'NC_FLOAT', [one_dimid]);
                netcdf.putAtt(ncid,gm_rmse_saltvarid,'long_name','global mean rmse_salinity (surface)');
                netcdf.putAtt(ncid,gm_rmse_saltvarid,'units','gram/kilogram');    
                
            netcdf.endDef(ncid);
    
            netcdf.putVar(ncid, timevarid, 0, grid.ntime, data.time);
            netcdf.putVar(ncid, tlongvarid, [0 0], [grid.nlon grid.nlat], grid.tlong);
            netcdf.putVar(ncid, tlatvarid, [0 0], [grid.nlon grid.nlat], grid.tlat);
    
            netcdf.putVar(ncid, tempvarid, [0 0 0], [grid.nlon grid.nlat grid.ntime], data.(['TEMP', '_model_', tmp.obsname_simple, '_i', tmp.iyear_str]));
            netcdf.putVar(ncid, obs_tempvarid, [0 0 0], [grid.nlon grid.nlat grid.ntime], data.(['TEMP', '_obs_', tmp.obsname_simple]));
            netcdf.putVar(ncid, bias_tempvarid, [0 0 0], [grid.nlon grid.nlat grid.ntime], data.(['TEMP', '_bias_', tmp.obsname_simple, '_i', tmp.iyear_str]));
            netcdf.putVar(ncid, sq_err_tempvarid, [0 0 0], [grid.nlon grid.nlat grid.ntime], data.(['TEMP', '_sq_err_', tmp.obsname_simple, '_i', tmp.iyear_str]));
            netcdf.putVar(ncid, rmse_tempvarid, [0 0], [grid.nlon grid.nlat], data.(['TEMP', '_rmse_', tmp.obsname_simple, '_i', tmp.iyear_str]));
            netcdf.putVar(ncid, gm_tempvarid, [0], [grid.ntime], data.(['gm_', 'TEMP', '_model_', tmp.obsname_simple, '_i', tmp.iyear_str]));
            netcdf.putVar(ncid, gm_obs_tempvarid, [0], [grid.ntime], data.(['gm_', 'TEMP', '_obs_', tmp.obsname_simple]));
            netcdf.putVar(ncid, gm_bias_tempvarid, [0], [grid.ntime], data.(['gm_', 'TEMP', '_bias_', tmp.obsname_simple, '_i', tmp.iyear_str]));
            netcdf.putVar(ncid, gm_sq_err_tempvarid, [0], [grid.ntime], data.(['gm_', 'TEMP', '_sq_err_', tmp.obsname_simple, '_i', tmp.iyear_str]));
            netcdf.putVar(ncid, gm_rmse_tempvarid, [0], [1], data.(['gm_', 'TEMP', '_rmse_', tmp.obsname_simple, '_i', tmp.iyear_str]));
    
            netcdf.putVar(ncid, saltvarid, [0 0 0], [grid.nlon grid.nlat grid.ntime], data.(['SALT', '_model_', tmp.obsname_simple, '_i', tmp.iyear_str]));
            netcdf.putVar(ncid, obs_saltvarid, [0 0 0], [grid.nlon grid.nlat grid.ntime], data.(['SALT', '_obs_', tmp.obsname_simple]));
            netcdf.putVar(ncid, bias_saltvarid, [0 0 0], [grid.nlon grid.nlat grid.ntime], data.(['SALT', '_bias_', tmp.obsname_simple, '_i', tmp.iyear_str]));
            netcdf.putVar(ncid, sq_err_saltvarid, [0 0 0], [grid.nlon grid.nlat grid.ntime], data.(['SALT', '_sq_err_', tmp.obsname_simple, '_i', tmp.iyear_str]));
            netcdf.putVar(ncid, rmse_saltvarid, [0 0], [grid.nlon grid.nlat], data.(['SALT', '_rmse_', tmp.obsname_simple, '_i', tmp.iyear_str]));
            netcdf.putVar(ncid, gm_saltvarid, [0], [grid.ntime], data.(['gm_', 'SALT', '_model_', tmp.obsname_simple, '_i', tmp.iyear_str]));
            netcdf.putVar(ncid, gm_obs_saltvarid, [0], [grid.ntime], data.(['gm_', 'SALT', '_obs_', tmp.obsname_simple]));
            netcdf.putVar(ncid, gm_bias_saltvarid, [0], [grid.ntime], data.(['gm_', 'SALT', '_bias_', tmp.obsname_simple, '_i', tmp.iyear_str]));
            netcdf.putVar(ncid, gm_sq_err_saltvarid, [0], [grid.ntime], data.(['gm_', 'SALT', '_sq_err_', tmp.obsname_simple, '_i', tmp.iyear_str]));
            netcdf.putVar(ncid, gm_rmse_saltvarid, [0], [1], data.(['gm_', 'SALT', '_rmse_', tmp.obsname_simple, '_i', tmp.iyear_str]));
    
            netcdf.close(ncid);
            
            disp(['saved file is ', config.savefilename]);
        end
    end
end


function var_obs = f_varname_obs(var_model)
    switch var_model
        case 'SALT'
            var_obs='bsal';  % PSS-78       
        case 'TEMP'
            var_obs='theta'; % ITS-90
        case 'DIC' % mmol/m^3
            var_obs='dic'; %umol/kg
        case 'PO4' %mmol/m^3
            var_obs='phos';  % standard method with low accuracy, umol/kg
%             var_obs='llp';  %high precision, nmol/kg
        case 'ALK' % meq/m^3
            var_obs='alk'; % 'ueq/kg'
        case 'NO3' %mmol/m^3
%             var_obs='nit'; %(NO2 + NO3 in ALOHA data) umol/kg
            var_obs='lln'; %(NO2 + NO3 in ALOHA data) nmol/kg
        case 'SiO3' %mmol/m^3
            var_obs='sil';  %umol/kg
        case 'PD' %g/cm^3
            var_obs='sigma';  %kg/m3
        case 'NO2'
            var_obs='no2'; %nmol/kg
    end
end

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


