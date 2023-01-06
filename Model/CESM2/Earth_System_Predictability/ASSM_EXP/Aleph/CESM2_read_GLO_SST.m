% %  Created 15-Dec-2022 by Yong-Yub Kim
clc; clear all; close all;

%% set path
tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

%% model configuration
dirs.root='/mnt/lustre/proj/earth.system.predictability/ASSM_EXP';
dirs.yoshi_root='/proj/yoshi/DATA/CESM2_ODA';
dirs.archive=[dirs.root, filesep, 'archive'];
dirs.saveroot='/mnt/lustre/proj/kimyy/Model/CESM2/ESP/ASSM_EXP';

config.years=1960:2021;
config.months=1:12;
config.scenname='HIST';
config.gridname='f09_g17';
% config.obsnames={'oras4', 'projdv7.3', 'en4.2'};
config.obsnames={'projdv7.3', 'en4.2'};


% config.obsnames={'oras4'};

% config.ensnames={'ba-10p1'};

% config.ensnames={'ba-10p1', 'ba-10p2', 'ba-10p3', 'ba-10p4', 'ba-10p5'};

config.ensnames={'ba-10p1', 'ba-10p2', 'ba-10p3', 'ba-10p4', 'ba-10p5', ...
    'ba-20p1', 'ba-20p2', 'ba-20p3', 'ba-20p4', 'ba-20p5'};

config.component='ocn';
% config.varnames={'SALT', 'TEMP', 'DIC', 'PO4', 'ALK', 'NO3', 'SiO3', 'PD', 'NO2'};
config.varnames={'TEMP'};

config.len_t_y = length(config.years);
config.len_t_m = length(config.months);
config.len_t = config.len_t_y * config.len_t_m;
config.len_obs= length(config.obsnames);
config.len_ens= length(config.ensnames);



% system(['ls ', dirs.datadir])
% b.e21.BHISTsmbb.f09_g17.assm.oras4_ba-10p1.pop.h.once.nc

% RCM_info.name={'test2127', 'test2128', 'test2129', 'test2130', 'test2131', 'ens2201'};

%% observation configuration
% config_obs.staname = 'ALOHA';
% config_obs.sta_lon = 360.0-158.0;
% config_obs.sta_lat = 22.0+45.0/60.0;

config_obs.staname = 'GLO';
config_obs.sta_lon = [0, 360]; % degrees west to 0~360 degrees east
config_obs.sta_lat = [-100 100]; % degrees north

%bsal; bottle salinity
%theta; potential temperature
%dic; dissolved inorganic carbon, DIC
%phos; phosphate
%alk; alkaliniy
%nit; Nitrate + Nitrite
%sil; Silicate
%sigma; potential density


%% get model data
for varind=1:length(config.varnames)
    tmp.varname=config.varnames{varind};
  
    for obsind=1:length(config.obsnames)
        tmp.obsname=config.obsnames{obsind};
        mkdir([dirs.saveroot, filesep, config_obs.staname]);
%         config.savename=[dirs.saveroot, filesep, config_obs.staname, filesep, 'CESM2_', config_obs.staname, '_', ...
%             tmp.obsname, '_', tmp.varname, '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
%         config.savename_cfg=[dirs.saveroot, filesep, config_obs.staname, filesep, 'cfg_CESM2_', config_obs.staname, '_', ...
%             tmp.obsname, '_', tmp.varname, '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
%         if (exist(config.savename) ==0)
            for ensind=1:length(config.ensnames)
                tmp.ensname=config.ensnames{ensind};
                config.casename_m=[config.gridname, '.hcst.', tmp.obsname, '-', tmp.ensname];

                config.savename=[dirs.saveroot, filesep, config_obs.staname, filesep, 'CESM2_', config_obs.staname, '_', ...
                    tmp.obsname, '_', tmp.ensname, '_', tmp.varname, '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.nc'];

                dirs.datadir=[dirs.archive, filesep, 'b.e21.B', config.scenname, 'smbb.', ...
                    config.gridname, '.assm.', tmp.obsname, '_', tmp.ensname, filesep, ...
                    config.component, filesep, 'hist'];
                dirs.yoshi_datadir=[dirs.yoshi_root, filesep, ...
                    'assm.', tmp.obsname, '_', tmp.ensname];
    
                %% read grid
                if ensind==1
                    [tmp.error_status, tmp.value]=system(['ls ', dirs.datadir, '| grep once']);
                    tmp.gridname = [dirs.datadir, filesep, tmp.value(1:end-1)];
                    CESM2_grid.lon_t=ncread(tmp.gridname, 'TLONG'); 
                    CESM2_grid.lat_t=ncread(tmp.gridname, 'TLAT'); 
                    CESM2_grid.region_mask=ncread(tmp.gridname, 'REGION_MASK'); 
                    CESM2_grid.ocean_mask=NaN(size(CESM2_grid.region_mask));
                    CESM2_grid.ocean_mask(CESM2_grid.region_mask>0)=1;
                    CESM2_grid.z_t=ncread(tmp.gridname, 'z_t');  % (cm)
                    CESM2_grid.dz=ncread(tmp.gridname, 'dz');
                    CESM2_grid.tarea=ncread(tmp.gridname, 'TAREA'); 

                    CESM2_grid.len_z_t=length(CESM2_grid.z_t);
                    
%                     [CESM2_grid.obs_lon_ind_w, CESM2_grid.obs_lon_ind_e, CESM2_grid.obs_lat_ind_s, CESM2_grid.obs_lat_ind_n] = ...
%                         Func_0012_findind_Y(0.1, [config_obs.sta_lon, config_obs.sta_lat], ...
%                         CESM2_grid.lon_t.*CESM2_grid.ocean_mask, ...
%                         CESM2_grid.lat_t.*CESM2_grid.ocean_mask, 'CESM2'); % find valid lon, lat index near station
                    CESM2_grid.obs_lon_ind_w=1;
                    CESM2_grid.obs_lat_ind_s=1;
                    [CESM2_grid.obs_lon_ind_e, CESM2_grid.obs_lat_ind_n]=size(CESM2_grid.lon_t);
                    CESM2_grid.len_x= CESM2_grid.obs_lon_ind_e-CESM2_grid.obs_lon_ind_w+1;
                    CESM2_grid.len_y= CESM2_grid.obs_lat_ind_n-CESM2_grid.obs_lat_ind_s+1;
                    CESM2_grid.cut_lon_t = CESM2_grid.lon_t(CESM2_grid.obs_lon_ind_w:CESM2_grid.obs_lon_ind_e, ...
                        CESM2_grid.obs_lat_ind_s:CESM2_grid.obs_lat_ind_n);
                    CESM2_grid.cut_lat_t = CESM2_grid.lat_t(CESM2_grid.obs_lon_ind_w:CESM2_grid.obs_lon_ind_e, ...
                        CESM2_grid.obs_lat_ind_s:CESM2_grid.obs_lat_ind_n);
                    CESM2_grid.dz_3d= permute(repmat(CESM2_grid.dz',1,CESM2_grid.len_x,CESM2_grid.len_y),[2,3,1]);
                    CESM2_grid.z_t_3d= permute(repmat(CESM2_grid.z_t',1,CESM2_grid.len_x,CESM2_grid.len_y),[2,3,1]);
                %% initialize

%                     CESM2_data.(tmp.varname)(1:config.len_ens, CESM2_grid.len_x, CESM2_grid.len_y, 1:config.len_t_y,1:config.len_t_m)=NaN;
                    
                end
                CESM2_data.(tmp.varname)(CESM2_grid.len_x, CESM2_grid.len_y, config.len_t_y*config.len_t_m)=NaN;

                %% read data
                for yind=1:length(config.years)
    
                    tmp.year=config.years(yind); tmp.yearstr=num2str(tmp.year, '%04i');
                    disp([tmp.varname, ', ', tmp.obsname, ', ', tmp.ensname, ', ', tmp.yearstr]);
    
                    for mind=1:length(config.months)
                        tmp.month=config.months(mind); tmp.monthstr=num2str(tmp.month, '%02i');
                        tmp.tind=(yind-1)*length(config.months)+mind;
                        [tmp.error_status, tmp.value]=system(['ls ', dirs.yoshi_datadir, ...
                            '/*pop*', tmp.yearstr, '-', tmp.monthstr, '.* | grep ', 'h.', tmp.yearstr, '-', tmp.monthstr]);
                        if strcmp(tmp.value(end-3:end-1), '.nc')==1 && strcmp('NO2', tmp.varname)~=1

                            tmp.filename = tmp.value(1:end-1);
    
                            tmp.data= ...
                                squeeze(ncread(tmp.filename, tmp.varname, ...
                                [CESM2_grid.obs_lon_ind_w, CESM2_grid.obs_lat_ind_s, 1, 1], ...
                                [CESM2_grid.len_x, CESM2_grid.len_y, 1, inf])); % just read surface
%                             CESM2_data.(tmp.varname)(ensind,:,:,yind,mind)= tmp.data; % [ens, x, y, year, month, z]
                            
                            CESM2_data.(tmp.varname)(:,:,tmp.tind)= tmp.data; % [ens, x, y, year, month, z]
                            CESM2_data.time(tmp.tind)=ncread(tmp.filename, 'time');

                            CESM2_data.(['gm_',tmp.varname])(tmp.tind) = f_gm_var(tmp.data, CESM2_grid.tarea);
                        else
                            tmp.filename = NaN;
%                             CESM2_data.(tmp.varname)(ensind,:,:,yind,mind)= NaN;
                            CESM2_data.(tmp.varname)(:,:,tmp.tind)= NaN; % [ens, x, y, year, month, z]
                            CESM2_data.(['gm_',tmp.varname])(tmp.tind) = NaN;
                            CESM2_data.time(tmp.tind) = CESM2_data.time(tmp.tind-1) + 30;
                        end
                        clear functions  % to read faster
                    end
                end
                if (exist(config.savename)==0)
                    ncid = netcdf.create(config.savename,'NETCDF4');
            
                    lon_dimid = netcdf.defDim(ncid, 'nlon', CESM2_grid.len_x);
                    lat_dimid = netcdf.defDim(ncid,'nlat', CESM2_grid.len_y);
                    time_dimid = netcdf.defDim(ncid, 'time', 0);
                    one_dimid = netcdf.defDim(ncid, 'one', 1);
            
                    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                        'title', ['CESM2 Hindcast ', 'surface variables']);
                    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                        'source', ['CESM2, ',config.casename_m]);
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
                        netcdf.defVarChunking(ncid, tlongvarid, 'CHUNKED', [CESM2_grid.len_x/10, CESM2_grid.len_y/10]);
                        netcdf.defVarDeflate(ncid, tlongvarid, true, true, 1);
                    tlatvarid=netcdf.defVar(ncid, 'TLAT', 'NC_DOUBLE', [lon_dimid lat_dimid]);
                        netcdf.putAtt(ncid,tlatvarid,'long_name','array of t-grid latitudes');
                        netcdf.putAtt(ncid,tlatvarid,'units','degree_north');
                        netcdf.defVarChunking(ncid, tlatvarid, 'CHUNKED', [CESM2_grid.len_x/10, CESM2_grid.len_y/10]);
                        netcdf.defVarDeflate(ncid, tlatvarid, true, true, 1);
            
                    tempvarid=netcdf.defVar(ncid, 'TEMP', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                        netcdf.putAtt(ncid,tempvarid,'long_name','Potential temperature (surface)');
                        netcdf.putAtt(ncid,tempvarid,'units','degC');  
                        netcdf.defVarChunking(ncid, tempvarid, 'CHUNKED', [CESM2_grid.len_x/10, CESM2_grid.len_y/10 1]);
                        netcdf.defVarDeflate(ncid, tempvarid, true, true, 1);
                    gm_tempvarid=netcdf.defVar(ncid, 'gm_temp', 'NC_FLOAT', [time_dimid]);
                        netcdf.putAtt(ncid,gm_tempvarid,'long_name','global mean Potential temperature (surface)');
                        netcdf.putAtt(ncid,gm_tempvarid,'units','degC');
                    
                                                        
                    netcdf.endDef(ncid);
            
                    netcdf.putVar(ncid, timevarid, 0, config.len_t , CESM2_data.time);
                    netcdf.putVar(ncid, tlongvarid, [0 0], [CESM2_grid.len_x CESM2_grid.len_y], CESM2_grid.cut_lon_t);
                    netcdf.putVar(ncid, tlatvarid, [0 0], [CESM2_grid.len_x CESM2_grid.len_y], CESM2_grid.cut_lat_t);
            
                    netcdf.putVar(ncid, tempvarid, [0 0 0], [CESM2_grid.len_x CESM2_grid.len_y config.len_t], CESM2_data.(tmp.varname));
                    
                    netcdf.putVar(ncid, gm_tempvarid, [0], [config.len_t], CESM2_data.(['gm_',tmp.varname]));
                                                   
                    netcdf.close(ncid);
                    
                    disp(['saved file is ', config.savename]);
                end


            end
%             save(config.savename, '-struct', 'CESM2_data', '-v7.3');
%             save(config.savename_cfg, 'CESM2_data', 'CESM2_grid', 'config', 'config_obs', 'dirs', '-v7.3');
%             clear CESM2_data
%         else
%             load(config.savename, 'CESM2_data')
%         end
    end

end

% tmp.dd=squeeze(mean(CESM2_data.DIC_vm_200m,4));
% plot(tmp.dd)

%     
%% save ncfile
%                 dirs.saveroot='/mnt/lustre/proj/kimyy/Model/CESM2/ESP/HCST_EXP';
%                 dirs.savedir=[dirs.saveroot, filesep, config.casename_m, filesep, 'GMSV'];
%                 mkdir(dirs.savedir);
%                 config.savefilename=[dirs.savedir, filesep, 'GMSV_', config.casename, '.nc'];
                
                



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


function gmsst = f_gm_var(var_2d, area)
    mask_var=NaN(size(var_2d));
    mask_var(isfinite(var_2d))=1;
    area=area.*mask_var;
    var_2d=squeeze(var_2d);
    var_2d_sum=var_2d.*area;
    var_2d_sum=sum(var_2d_sum(:), 'omitnan');
    gmsst=var_2d_sum ./ sum(area(:), 'omitnan');
end
