% %  Created 06-Jan-2023 by Yong-Yub Kim
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

cfg.y=1960:2030; %config.years
cfg.m=1:12; %config.months
cfg.sce='HIST'; %cfg.scenname
cfg.gr='f09_g17'; %cfg.gridname


cfg.ens={'ba-10p1', 'ba-10p2', 'ba-10p3', 'ba-10p4', 'ba-10p5', ...
    'ba-20p1', 'ba-20p2', 'ba-20p3', 'ba-20p4', 'ba-20p5'}; %cfg.ensnames

% cfg.v={'SALT', 'TEMP', 'DIC', 'PO4', 'ALK', 'NO3', 'SiO3', 'PD', 'NO2'};
% cfg.v={'TEMP', 'SALT', 'diatChl', 'diazChl', 'spChl', 'dst_a1', 'dst_a2', 'dst_a3', 'U', 'V'}; %cfg.varnames
cfg.v={'TEMP', 'SALT',  'dst_a1', 'dst_a2', 'dst_a3', 'U', 'V'}; %cfg.varnames

% cfg.v={'dst_a1'};

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


% [t.er_s, t.el]=
system(['ls -d ', dirs.r, '/*/ | grep -v ext | rev | cut -c 2-9 | rev > ~/tmpn']); %save ensemble name list (nonuniq)
cfg.ens=unique(textread('~/tmpn', '%s')); %tmp.ensnamelist

% ls -d */ | grep -v ext | rev | cut -c 2-9 | rev

%% get model data

for y=min(cfg.y):max(cfg.y)
    for m=1:12
        t.i_t=m+(y-min(cfg.y))*12;
        cfg.t(t.i_t)=datenum(y,m,15);
    end
end


for eind=1:length(cfg.ens)
% for eind=1:1

    t.ens=cfg.ens{eind}; %t.ensname
    for vind=1:length(cfg.v)
        t.v=cfg.v{vind}; %t.varname
        switch t.v
            case {'TEMP', 'SALT', 'diatChl', 'diazChl', 'spChl'}
                t.cmp='ocn'; %t.component
            case {'dst_a1', 'AODVIS', 'AODDUST', 'U', 'V'}
                t.cmp='atm';
        end
    
    %    system(['ls ', dirs.r, '/*', t.ens, '*/', t.cmp, '/proc/tseries/month_1/* | grep ', '''', '\.', '''', t.v, '''', '\.', ''''])
        t.dl=dir([dirs.r, '/*', t.ens, '*/', t.cmp, '/proc/tseries/month_1/*.', t.v, '.*']); %dir information list
        t.n_f=length(t.dl); %t.number of files
    
    
    
        mkdir([dirs.svr, fs, cfg_o.st, fs, t.cmp, fs, t.v]); %[saveroot, filesep, obs_station, filesep, variable]
    %     cfg.casename_m=[cfg.gr, '.hcst.', t.obsname, '-', t.ens];
    
        cfg.svnm=strcat(dirs.svr, fs, cfg_o.st, fs, t.cmp, fs, t.v, fs, ...
            ['CESM2_', cfg_o.st, '_', t.ens, '_', t.v, '_', cfg.ys, '.nc']); %cfg.savefilename (will be saved)
    
    %     dirs.dat=[dirs.arc, fs, 'b.e21.B', cfg.sce, 'smbb.', ...
    %         cfg.gr, '.assm.', t.obsname, '_', t.ens, fs, ...
    %         t.cmp, fs, 'hist']; %dirs.datadir
    
    
    %     CESM2_data.(t.v)(CESM2_grid.l_x, CESM2_grid.l_y, cfg.l_y*cfg.l_m)=NaN;
    
        t.y=min(cfg.y);
        switch t.cmp
            case 'ocn'
                t.i_sf=1;
            case 'atm'
                t.i_sf=32;
        end
        while t.y<=max(cfg.y)
            for fi=1:t.n_f
                t.f_sp=split(t.dl(fi).name, '.'); %t.filename_splited
                t.f_yms=split(t.f_sp{end-1}, '-'); %t.filename_year month
                t.f_y_min=str2num(t.f_yms{1}(1:4));
                t.f_y_max=str2num(t.f_yms{2}(1:4));
                if (t.f_y_min <= t.y && t.y <= t.f_y_max) % search filename with right period
                    t.fn=[t.dl(fi).folder, fs, t.dl(fi).name];
                    if t.y==min(cfg.y)
                        t.gn = t.fn;
                        switch t.cmp
                            case 'ocn'
                                CESM2_grid.lon_t=ncread(t.gn, 'TLONG'); 
                                CESM2_grid.lat_t=ncread(t.gn, 'TLAT'); 
                                CESM2_grid.region_mask=ncread(t.gn, 'REGION_MASK'); 
                                CESM2_grid.ocn_mask=NaN(size(CESM2_grid.region_mask));
                                CESM2_grid.ocn_mask(CESM2_grid.region_mask>0)=1;
                                CESM2_grid.z_t=ncread(t.gn, 'z_t');  % (cm)
                                CESM2_grid.dz=ncread(t.gn, 'dz');
                                CESM2_grid.tarea=ncread(t.gn, 'TAREA'); 
                            case 'atm'
                                CESM2_grid.lon=ncread(t.gn, 'lon'); %1.5
                                CESM2_grid.lat=ncread(t.gn, 'lat'); %0.9424
                                [CESM2_grid.lat_t, CESM2_grid.lon_t] = meshgrid(CESM2_grid.lat, CESM2_grid.lon);
                                t.km_lat=m_lldist([0 0], [CESM2_grid.lat(1) CESM2_grid.lat(2)]);
                                for lati=1:length(CESM2_grid.lat)
                                    t.km_lon(1:length(CESM2_grid.lon),lati)=m_lldist([CESM2_grid.lon(1) CESM2_grid.lon(2)],...
                                        [CESM2_grid.lat(lati) CESM2_grid.lat(lati)]);
                                end
                                CESM2_grid.tarea=t.km_lon.*t.km_lat;
                        end
                        
            
%                         CESM2_grid.len_z_t=length(CESM2_grid.z_t);
                        
                        if strcmp(cfg_o.st, 'GLO')
                            CESM2_grid.id_w=1;
                            CESM2_grid.id_s=1;
                        else
                            switch t.cmp
                                case 'ocn'
                                [CESM2_grid.id_w, CESM2_grid.id_e, CESM2_grid.id_s, CESM2_grid.id_n] = ...
                                    Func_0012_findind_Y(0.1, [cfg_o.st_lon, cfg_o.st_lat], ...
                                    CESM2_grid.lon_t.*CESM2_grid.ocn_mask, ...
                                    CESM2_grid.lat_t.*CESM2_grid.ocn_mask, 'CESM2'); % find valid lon, lat index near station
                                case 'atm'
                                 [CESM2_grid.id_w, CESM2_grid.id_e, CESM2_grid.id_s, CESM2_grid.id_n] = ...
                                    Func_0012_findind_Y(0.1, [cfg_o.st_lon, cfg_o.st_lat], ...
                                    CESM2_grid.lon_t, ...
                                    CESM2_grid.lat_t, 'CESM2'); % find valid lon, lat index near station
                            end
                        end
                        [CESM2_grid.id_e, CESM2_grid.id_n]=size(CESM2_grid.lon_t);

                        CESM2_grid.l_x= CESM2_grid.id_e-CESM2_grid.id_w+1;
                        CESM2_grid.l_y= CESM2_grid.id_n-CESM2_grid.id_s+1;
                        CESM2_grid.cut_lon_t = CESM2_grid.lon_t(CESM2_grid.id_w:CESM2_grid.id_e, ...
                            CESM2_grid.id_s:CESM2_grid.id_n);
                        CESM2_grid.cut_lat_t = CESM2_grid.lat_t(CESM2_grid.id_w:CESM2_grid.id_e, ...
                            CESM2_grid.id_s:CESM2_grid.id_n);
%                         CESM2_grid.dz_3d= permute(repmat(CESM2_grid.dz',1,CESM2_grid.l_x,CESM2_grid.l_y),[2,3,1]);
%                         CESM2_grid.z_t_3d= permute(repmat(CESM2_grid.z_t',1,CESM2_grid.l_x,CESM2_grid.l_y),[2,3,1]);
    
                        CESM2_data.(t.v)(CESM2_grid.l_x, CESM2_grid.l_y, cfg.l_y*cfg.l_m)=NaN;
                        t.i_tmin=1+(t.f_y_min-min(cfg.y))*cfg.l_m; % minimum t index
                        t.i_tmax=((t.f_y_max - min(cfg.y))+1)*cfg.l_m; % length t to read
                        CESM2_data.(t.v)(CESM2_grid.id_w:CESM2_grid.id_e, CESM2_grid.id_s:CESM2_grid.id_n, 1:t.i_tmax) = ...
                            ncread(t.fn, t.v, [CESM2_grid.id_w, CESM2_grid.id_s, t.i_sf, t.i_tmin], ...
                            [CESM2_grid.l_x, CESM2_grid.l_y, 1, t.i_tmax]);
                        t.i_t=1+t.i_tmax;
                        t.y=t.f_y_max+1;
                    else
                        if t.f_y_max < max(cfg.y)
                            t.i_tmax=((t.f_y_max - t.f_y_min)+1)*cfg.l_m;
                            CESM2_data.(t.v)(CESM2_grid.id_w:CESM2_grid.id_e, CESM2_grid.id_s:CESM2_grid.id_n, t.i_t:t.i_t+t.i_tmax-1)= ...
                                ncread(t.fn, t.v, [CESM2_grid.id_w, CESM2_grid.id_s, t.i_sf, 1], ...
                                [CESM2_grid.l_x, CESM2_grid.l_y, 1, t.i_tmax]);
                            t.i_t=t.i_t + t.i_tmax;
                            t.y=t.f_y_max+1;
                        else
                            t.i_tmax=((max(cfg.y) - t.f_y_min)+1)*cfg.l_m;
                            CESM2_data.(t.v)(CESM2_grid.id_w:CESM2_grid.id_e, CESM2_grid.id_s:CESM2_grid.id_n, t.i_t:t.i_t+t.i_tmax-1)= ...
                                ncread(t.fn, t.v, [CESM2_grid.id_w, CESM2_grid.id_s, t.i_sf, 1], ...
                                [CESM2_grid.l_x, CESM2_grid.l_y, 1, t.i_tmax]);
                            t.y=max(cfg.y)+1;
                        end
                    end
                    disp(num2str(t.y))
                end
            end
            clear functions  % to read faster
        end
        CESM2_data.(['gm_',t.v]) = NaN(cfg.l_t,1);
        for tind=1:cfg.l_t
            CESM2_data.(['gm_',t.v])(tind) = f_gm_var(CESM2_data.(t.v)(:,:,tind), CESM2_grid.tarea);
        end

% pcolor(CESM2_data.TEMP(:,:,744)'); shading flat; colorbar;
       

        if (exist(cfg.svnm)==0)
            ncid = netcdf.create(cfg.svnm,'NETCDF4');
    
            lon_dimid = netcdf.defDim(ncid, 'nlon', CESM2_grid.l_x);
            lat_dimid = netcdf.defDim(ncid,'nlat', CESM2_grid.l_y);
            time_dimid = netcdf.defDim(ncid, 'time', 0);
            one_dimid = netcdf.defDim(ncid, 'one', 1);
    
            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                'title', ['CESM2 LENS2 ', 'surface variables']);
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
                netcdf.defVarChunking(ncid, tlongvarid, 'CHUNKED', [CESM2_grid.l_x/10, CESM2_grid.l_y/10]);
                netcdf.defVarDeflate(ncid, tlongvarid, true, true, 1);
            tlatvarid=netcdf.defVar(ncid, 'TLAT', 'NC_DOUBLE', [lon_dimid lat_dimid]);
                netcdf.putAtt(ncid,tlatvarid,'long_name','array of t-grid latitudes');
                netcdf.putAtt(ncid,tlatvarid,'units','degree_north');
                netcdf.defVarChunking(ncid, tlatvarid, 'CHUNKED', [CESM2_grid.l_x/10, CESM2_grid.l_y/10]);
                netcdf.defVarDeflate(ncid, tlatvarid, true, true, 1);
    
            tempvarid=netcdf.defVar(ncid, t.v, 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                netcdf.putAtt(ncid,tempvarid,'long_name', [t.v,' (surface)']);
                netcdf.defVarChunking(ncid, tempvarid, 'CHUNKED', [CESM2_grid.l_x/10, CESM2_grid.l_y/10 1]);
                netcdf.defVarDeflate(ncid, tempvarid, true, true, 1);
            gm_tempvarid=netcdf.defVar(ncid, ['gm_', t.v], 'NC_FLOAT', [time_dimid]);
                netcdf.putAtt(ncid,gm_tempvarid,'long_name','global mean');
                netcdf.putAtt(ncid,gm_tempvarid,'units','degC');
            
                                                
            netcdf.endDef(ncid);
    
            netcdf.putVar(ncid, timevarid, 0, cfg.l_t , cfg.t);
            netcdf.putVar(ncid, tlongvarid, [0 0], [CESM2_grid.l_x CESM2_grid.l_y], CESM2_grid.cut_lon_t);
            netcdf.putVar(ncid, tlatvarid, [0 0], [CESM2_grid.l_x CESM2_grid.l_y], CESM2_grid.cut_lat_t);
    
            netcdf.putVar(ncid, tempvarid, [0 0 0], [CESM2_grid.l_x CESM2_grid.l_y cfg.l_t], CESM2_data.(t.v));
            
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