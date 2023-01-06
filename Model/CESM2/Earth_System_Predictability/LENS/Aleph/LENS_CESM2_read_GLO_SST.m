% %  Created 16-Dec-2022 by Yong-Yub Kim
clc; clear all; close all;

%% set path
t.dbp = '/mnt/lustre/proj/kimyy/Dropbox'; % t.dropboxpath
fs=filesep; %filesep
addpath(genpath([t.dbp, fs, 'source', fs, 'matlab', fs, 'function']));
            [t.dbp, t.er_s] = Func_0008_set_dropbox_path(computer);

% t.er_s : t.error_status

%% model configuration
% dirs.r='/mnt/lustre/proj/earth.system.predictability/ASSM_EXP';
dirs.r='/mnt/lustre/proj/jedwards/archive'; %dirs.root
% dirs.yoshi_root='/proj/yoshi/DATA/CESM2_ODA';
dirs.arc=[dirs.r, fs, 'archive']; %dirs.archive
dirs.svr='/mnt/lustre/proj/kimyy/Model/CESM2/ESP/LENS/surface';  %dirs.saveroot

cfg.y=1960:2021; %config.years
cfg.m=1:12; %config.months
cfg.sce='HIST'; %cfg.scenname
cfg.gr='f09_g17'; %cfg.gridname


cfg.ens={'ba-10p1', 'ba-10p2', 'ba-10p3', 'ba-10p4', 'ba-10p5', ...
    'ba-20p1', 'ba-20p2', 'ba-20p3', 'ba-20p4', 'ba-20p5'}; %cfg.ensnames

t.cmp='ocn'; %t.component
% cfg.v={'SALT', 'TEMP', 'DIC', 'PO4', 'ALK', 'NO3', 'SiO3', 'PD', 'NO2'};
cfg.v={'TEMP'}; %cfg.varnames

cfg.l_y = length(cfg.y); %cfg.l_t_y
cfg.l_m = length(cfg.m); %cfg.l_m
cfg.l_t = cfg.l_y * cfg.l_m; %cfg.l_t
cfg.l_e= length(cfg.ens); %cfg.l_e

% system(['ls ', dirs.dat])
% b.e21.BHISTsmbb.f09_g17.assm.oras4_ba-10p1.pop.h.once.nc

%% observation configuration
cfg_o.st = 'GLO'; %config_obs.staname
cfg_o.st_lon = [0, 360]; % degrees west to 0~360 degrees east, cfg_o.st_lon
cfg_o.st_lat = [-100 100]; % degrees north, cfg_o.st_lat

cfg.ys=[num2str(min(cfg.y)), '_', num2str(max(cfg.y))]; %cfg.year_range_string


t.el=system(['ls ', dirs.r, '/*/ | grep -v ext | rev | cut -c 2-9 | rev']); %ensemble name list (nonuniq)
% ls -d */ | grep -v ext | rev | cut -c 2-9 | rev

%% get model data
for eind=1:length(cfg.ens)
    t.ens=cfg.ens{eind}; %t.ensname
    for vind=1:length(cfg.v)
    t.v=cfg.v{vind}; %t.varname
  
  
        mkdir([dirs.svr, fs, cfg_o.st]); %[saveroot, filesep, obs_station]
%         cfg.svnm=[dirs.svr, fs, cfg_o.st, fs, 'CESM2_', cfg_o.st, '_', ...
%             t.obsname, '_', t.v, '_', num2str(min(cfg.y)), '_', num2str(max(cfg.y)), '.mat'];
%         cfg.svnm_cfg=[dirs.svr, fs, cfg_o.st, fs, 'cfg_CESM2_', cfg_o.st, '_', ...
%             t.obsname, '_', t.v, '_', num2str(min(cfg.y)), '_', num2str(max(cfg.y)), '.mat'];
%         if (exist(cfg.svnm) ==0)
        cfg.casename_m=[cfg.gr, '.hcst.', t.obsname, '-', t.ens];

        cfg.svnm=strcat(dirs.svr, fs, cfg_o.st, fs, ...
            ['CESM2_', cfg_o.st, '_', t.ens, '_', t.v, '_', cfg.ys, '.nc']); %cfg.savefilename

        dirs.dat=[dirs.arc, fs, 'b.e21.B', cfg.sce, 'smbb.', ...
            cfg.gr, '.assm.', t.obsname, '_', t.ens, fs, ...
            t.cmp, fs, 'hist']; %dirs.datadir
        dirs.yoshi_datadir=[dirs.yoshi_root, fs, ...
            'assm.', t.obsname, '_', t.ens];

        %% read grid
        if eind==1
            [t.er_s, t.value]=system(['ls ', dirs.dat, '| grep once']);
            t.gridname = [dirs.dat, fs, t.value(1:end-1)];
            CESM2_grid.lon_t=ncread(t.gridname, 'TLONG'); 
            CESM2_grid.lat_t=ncread(t.gridname, 'TLAT'); 
            CESM2_grid.region_mask=ncread(t.gridname, 'REGION_MASK'); 
            CESM2_grid.ocean_mask=NaN(size(CESM2_grid.region_mask));
            CESM2_grid.ocean_mask(CESM2_grid.region_mask>0)=1;
            CESM2_grid.z_t=ncread(t.gridname, 'z_t');  % (cm)
            CESM2_grid.dz=ncread(t.gridname, 'dz');
            CESM2_grid.tarea=ncread(t.gridname, 'TAREA'); 

            CESM2_grid.len_z_t=length(CESM2_grid.z_t);
            
%                     [CESM2_grid.obs_lon_ind_w, CESM2_grid.obs_lon_ind_e, CESM2_grid.obs_lat_ind_s, CESM2_grid.obs_lat_ind_n] = ...
%                         Func_0012_findind_Y(0.1, [cfg_o.st_lon, cfg_o.st_lat], ...
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

%                     CESM2_data.(t.v)(1:cfg.l_e, CESM2_grid.len_x, CESM2_grid.len_y, 1:cfg.l_y,1:cfg.l_m)=NaN;
            
        end
        CESM2_data.(t.v)(CESM2_grid.len_x, CESM2_grid.len_y, cfg.l_y*cfg.l_m)=NaN;

        %% read data
        for yind=1:length(cfg.y)

            t.year=cfg.y(yind); t.yearstr=num2str(t.year, '%04i');
            disp([t.v, ', ', t.obsname, ', ', t.ens, ', ', t.yearstr]);

            for mind=1:length(cfg.m)
                t.month=cfg.m(mind); t.monthstr=num2str(t.month, '%02i');
                t.tind=(yind-1)*length(cfg.m)+mind;
                [t.er_s, t.value]=system(['ls ', dirs.yoshi_datadir, ...
                    '/*pop*', t.yearstr, '-', t.monthstr, '.* | grep ', 'h.', t.yearstr, '-', t.monthstr]);
                if strcmp(t.value(end-3:end-1), '.nc')==1 && strcmp('NO2', t.v)~=1

                    t.filename = t.value(1:end-1);

                    t.data= ...
                        squeeze(ncread(t.filename, t.v, ...
                        [CESM2_grid.obs_lon_ind_w, CESM2_grid.obs_lat_ind_s, 1, 1], ...
                        [CESM2_grid.len_x, CESM2_grid.len_y, 1, inf])); % just read surface
%                             CESM2_data.(t.v)(eind,:,:,yind,mind)= t.data; % [ens, x, y, year, month, z]
                    
                    CESM2_data.(t.v)(:,:,t.tind)= t.data; % [ens, x, y, year, month, z]
                    CESM2_data.time(t.tind)=ncread(t.filename, 'time');

                    CESM2_data.(['gm_',t.v])(t.tind) = f_gm_var(t.data, CESM2_grid.tarea);
                else
                    t.filename = NaN;
%                             CESM2_data.(t.v)(eind,:,:,yind,mind)= NaN;
                    CESM2_data.(t.v)(:,:,t.tind)= NaN; % [ens, x, y, year, month, z]
                    CESM2_data.(['gm_',t.v])(t.tind) = NaN;
                    CESM2_data.time(t.tind) = CESM2_data.time(t.tind-1) + 30;
                end
                clear functions  % to read faster
            end
        end
        if (exist(cfg.svnm)==0)
            ncid = netcdf.create(cfg.svnm,'NETCDF4');
    
            lon_dimid = netcdf.defDim(ncid, 'nlon', CESM2_grid.len_x);
            lat_dimid = netcdf.defDim(ncid,'nlat', CESM2_grid.len_y);
            time_dimid = netcdf.defDim(ncid, 'time', 0);
            one_dimid = netcdf.defDim(ncid, 'one', 1);
    
            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                'title', ['CESM2 Hindcast ', 'surface variables']);
            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                'source', ['CESM2, ',cfg.casename_m]);
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
    
            netcdf.putVar(ncid, timevarid, 0, cfg.l_t , CESM2_data.time);
            netcdf.putVar(ncid, tlongvarid, [0 0], [CESM2_grid.len_x CESM2_grid.len_y], CESM2_grid.cut_lon_t);
            netcdf.putVar(ncid, tlatvarid, [0 0], [CESM2_grid.len_x CESM2_grid.len_y], CESM2_grid.cut_lat_t);
    
            netcdf.putVar(ncid, tempvarid, [0 0 0], [CESM2_grid.len_x CESM2_grid.len_y cfg.l_t], CESM2_data.(t.v));
            
            netcdf.putVar(ncid, gm_tempvarid, [0], [cfg.l_t], CESM2_data.(['gm_',t.v]));
                                           
            netcdf.close(ncid);
            
            disp(['saved file is ', cfg.svnm]);
        end


    end
%             save(cfg.svnm, '-struct', 'CESM2_data', '-v7.3');
%             save(cfg.svnm_cfg, 'CESM2_data', 'CESM2_grid', 'config', 'config_obs', 'dirs', '-v7.3');
%             clear CESM2_data
%         else
%             load(cfg.svnm, 'CESM2_data')
%         end

end

% t.dd=squeeze(mean(CESM2_data.DIC_vm_200m,4));
% plot(t.dd)

%     
%% save ncfile
%                 dirs.svr='/mnt/lustre/proj/kimyy/Model/CESM2/ESP/HCST_EXP';
%                 dirs.savedir=[dirs.svr, fs, cfg.casename_m, fs, 'GMSV'];
%                 mkdir(dirs.savedir);
%                 cfg.savefilename=[dirs.savedir, fs, 'GMSV_', cfg.casename, '.nc'];
                
                



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
