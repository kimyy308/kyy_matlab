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
grid.regions = [130 190 15 50];


cfg.vars = {'TEMP', 'SSH', 'SALT', 'DIC', 'DIC_ALT_CO2', 'FG_CO2'};

% tmp.dimids= [1, 2, 4];
cfg.vlayer=1; % surf, vertical slice 
% cfg.vlayer=1:10; % 10layer. don't put more than 15
cfg.vlayer=10; % 100m, vertical slice 
cfg.vlayer=20; % 200m, vertical slice 
% cfg.vlayer=15; %150m
cfg.vlayer=27; %305m
cfg.vlayer=31; %408m
cfg.vlayer_1st=min(cfg.vlayer);
cfg.vlayer_cnt=max(cfg.vlayer)-cfg.vlayer_1st+1;

for vari=1:length(cfg.vars)
    tmp.fs=filesep;
    cfg.var=cfg.vars{vari};
    cfg.obs_name=f_obs_name(cfg.var);
    cfg.obs_fname_mid=f_obs_name_mid(cfg.var);
    cfg.obs_varname=f_obs_varname(cfg.var);
    cfg.comp=Func_0025_CESM2_cmpname_var(cfg.var);
    cfg.obs_fname_module=f_obs_fname_module(cfg.comp);
    
    cfg.obs_iyears=1960:2020;

    cfg.obs_iyears2=f_obs_iyears(cfg.var);
  
      
    dirs.obsroot=['/mnt/lustre/proj/kimyy/Observation/', cfg.obs_name];
    dirs.matroot=['/mnt/lustre/proj/kimyy/tr_sysong/mat/', cfg.comp,'/', cfg.var];
    dirs.lens2root=['/mnt/lustre/proj/kimyy/tr_sysong/LENS2/', cfg.comp, '/', cfg.var];
    dirs.assmroot=['/mnt/lustre/proj/kimyy/tr_sysong/ASSM_EXP/', cfg.comp, '/', cfg.var];    
    
    cfg.iyears=cfg.obs_iyears;
    cfg.gnm='f09_g17';
    cfg.proj_year=5;

    cfg.max_iy=max(cfg.iyears);
    cfg.min_iy=min(cfg.iyears);
    cfg.len_t_y = length(cfg.iyears);
       
    tmp.gridname = [dirs.assmroot, tmp.fs, '../grid.nc'];
    tmp.maskname = '/mnt/lustre/proj/kimyy/Model/CESM2/ESP/grids/ocn/RECCAP2_region_masks_all_v20210412_POP2_grid.nc';
    
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

    [grid.id_w, grid.id_e, grid.id_s, grid.id_n] = Func_0012_findind_Y(1.0, grid.regions, ...
                grid.tlong, grid.tlat, 'CESM2'); % find valid lon, lat index near station

    grid.cut_tlong=grid.tlong(grid.id_w:grid.id_e, grid.id_s:grid.id_n);
    grid.cut_tlat=grid.tlat(grid.id_w:grid.id_e, grid.id_s:grid.id_n);
    grid.cut_nlon=size(grid.cut_tlong,1);
    grid.cut_nlat=size(grid.cut_tlat,2);

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
        cfg.casename_m=['all'];
    
        dirs.assmdir= [dirs.assmroot,  filesep];
        dirs.lens2dir= [dirs.lens2root, filesep];
        lap_time = tic;
    
        
%% variables initialization
        data.([tmp.varname, '_lens2'])=NaN(grid.cut_nlon, grid.cut_nlat, cfg.len_t_y);
        data.([tmp.varname, '_obs'])=NaN(grid.cut_nlon, grid.cut_nlat, cfg.len_t_y);
        data.([tmp.varname, '_assm'])=NaN(grid.cut_nlon, grid.cut_nlat, cfg.len_t_y);
        

    %     data.([tmp.varname, '_assm'])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        
%% read variables
        for iyear=min(cfg.iyears):max(cfg.iyears)
            tmp.iyear_str=num2str(iyear, '%04i');
            cfg.casename=['ensmean_', cfg.casename_m, '_i', tmp.iyear_str];
            tmp.fy=iyear;
            tmp.fy_str=num2str(tmp.fy, '%04i');

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
                    elseif (strcmp(cfg.obs_name, 'EN4')==1)
                        cfg.obs_fnm=['/mnt/lustre/proj/kimyy/Observation/EN4/EN4.2', tmp.fs, 'toso-en4.', tmp.fy_str, '-', tmp.m_str, '.nc'];
                    elseif (strcmp(cfg.obs_name, 'projdv7.3')==1)
                        cfg.obs_fnm=['/mnt/lustre/proj/kimyy/Observation/projdv7.3/ProjD_v7.3', tmp.fs, 'toso-proj7.', tmp.fy_str, '-', tmp.m_str, '.nc'];
                    else
                        cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'monthly_reg_', cfg.obs_fname_module(2:4),tmp.fs,cfg.obs_fname_mid,tmp.fy_str,tmp.m_str,'.nc'];
                    end

                    if exist(cfg.obs_fnm)~=0
                        ncid=netcdf.open(cfg.obs_fnm, 'NOWRITE');
                        tmpvarid=netcdf.inqVarID(ncid,cfg.obs_varname);
                        tmp.dd =  netcdf.getVar(ncid,tmpvarid, [grid.id_w grid.id_s cfg.vlayer_1st-1 0], [grid.cut_nlon grid.cut_nlat cfg.vlayer_cnt 1]);
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
                        tmp.ydata_obs(1:grid.cut_nlon,1:grid.cut_nlat,mi) = tmp.dd;
                        netcdf.close(ncid);
                        
                    else
                        tmp.ydata_obs(1:grid.cut_nlon,1:grid.cut_nlat,mi) = NaN;                    
                    end
                end
                tmp.ymean_obs(1:grid.cut_nlon,1:grid.cut_nlat)=mean(tmp.ydata_obs,3);

            else
                tmp.ymean_obs(1:grid.cut_nlon,1:grid.cut_nlat) = NaN;          
            end
%             tmp.ydata_obs(:,:,mon)=ncread(cfg.obs_fnm, cfg.obs_varname);

%% ASSM
            for mi=1:12
                tmp.m_str=num2str(mi,'%02i');
                cfg.assm_fnm=[dirs.assmdir, tmp.fs, tmp.varname, '_ensmean_', tmp.fy_str, '-', tmp.m_str, '.nc'];
                ncid=netcdf.open(cfg.assm_fnm, 'NOWRITE');
                tmpvarid=netcdf.inqVarID(ncid,tmp.fvarname);
                [tmp.fvarname, tmp.xtype, tmp.dimids, tmp.natts]= netcdf.inqVar(ncid, tmpvarid);
                if length(tmp.dimids)>3
                     tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [grid.id_w grid.id_s cfg.vlayer_1st-1 0], [grid.cut_nlon grid.cut_nlat cfg.vlayer_cnt 1]));
                     tmp.dd(abs(tmp.dd)>1e30)=NaN;
%                      ddd=mean(tmp.dd,3,'omitnan'); % depth mean
                     tmp.ydata_assm(1:grid.cut_nlon,1:grid.cut_nlat,mi) = tmp.dd; %depth averaged value             
                else
                    tmp.ydata_assm(1:grid.cut_nlon,1:grid.cut_nlat,mi) = netcdf.getVar(ncid,tmpvarid, ...
                        [grid.id_w grid.id_s 0], [grid.cut_nlon grid.cut_nlat 1]);
                end
                if(strcmp(tmp.varname(1:2),'mu')~=1 && strcmp(tmp.varname(1:2),'su')~=1)
                    data.units=netcdf.getAtt(ncid,tmpvarid,'units');
                else
                    data.units=' ';
                end
                netcdf.close(ncid);
            end
            tmp.ymean_assm=mean(tmp.ydata_assm,3);
                 

%% save variables as structure

            data.([tmp.varname, '_obs'])(1:grid.cut_nlon,1:grid.cut_nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_obs;
            data.([tmp.varname, '_assm'])(1:grid.cut_nlon,1:grid.cut_nlat,iyear-min(cfg.iyears)+1)= tmp.ymean_assm;

        end  %% initialized year loop end

            %% LENS2 mean
            ty=1;
            
            while ty<=cfg.len_t_y
                tmp.year=cfg.iyears(ty);
                tmp.fdir=dirs.lens2dir;
                tmp.flist=dir( [tmp.fdir, '/',cfg.var, '*'] );
                for kk = 1: length (tmp.flist)
                    tmp.fname_in = tmp.flist(kk).name;
                    tmp.fname_split = strsplit(tmp.fname_in, {'.'} );
                    tmp.fname_period = tmp.fname_split{1}(end-12:end);
                    fyear_str   = strsplit( tmp.fname_period, '-' );
                    fyear_start = str2num( fyear_str{1}(1:4) );
                    fyear_end   = str2num( fyear_str{2}(1:4) );
                    if( tmp.year >= fyear_start && tmp.year <= fyear_end )
        %                 disp(tmp.fname_period);
                        flag_file_in = true;            break;
                    end
                end
                tmp.fname=[tmp.fdir, '/', cfg.var, '_ensmean_', tmp.fname_period, '.nc'];
                ind_start=(tmp.year - fyear_start)*12+1;
                if fyear_end<=max(cfg.iyears)
                    ind_count=(fyear_end-tmp.year+1)*12;
                else
                    ind_count=(max(cfg.iyears)-tmp.year+1)*12;
                end
                
                ncid=netcdf.open(tmp.fname, 'NOWRITE');
                tmpvarid=netcdf.inqVarID(ncid,tmp.fvarname);
                if length(tmp.dimids)>3
                     tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [grid.id_w grid.id_s cfg.vlayer_1st-1 ind_start-1], [grid.cut_nlon grid.cut_nlat cfg.vlayer_cnt ind_count]));
                     tmp.dd(abs(tmp.dd)>1e30)=NaN;
%                      ddd=mean(tmp.dd,3,'omitnan'); % depth mean
                     tmp.ydata_lens2(1:grid.cut_nlon,1:grid.cut_nlat, (ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd;
                else
                    tmp.dd=netcdf.getVar(ncid,tmpvarid, [grid.id_w grid.id_s ind_start-1], [grid.cut_nlon grid.cut_nlat ind_count]);
                     tmp.dd(abs(tmp.dd)>1e30)=NaN;                        
                    tmp.ydata_lens2(1:grid.cut_nlon,1:grid.cut_nlat,(ty-1)*12+1:(ty-1)*12+ind_count) =tmp.dd;
                end
                netcdf.close(ncid);

%                 if length(tmp.dimids)>3
%                      tmp.fname2='/proj/kimyy/tmp/test_mat.nc';
%                      system(['cdo -O -w ', cfg.vm_str, ' -sellevel,', cfg.vl_z1_str,'/',cfg.vl_z2_str, ' ', tmp.fname, ' ', tmp.fname2]);
%                      try
%                          ncid = netcdf.open(tmp.fname2, 'NOWRITE');
%                          tmpvarid = netcdf.inqVarID(ncid, cfg.var);
%                          tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 ind_start-1], [grid.nlon grid.nlat ind_count]));
%                          tmp.dd(abs(tmp.dd)>1e30)=NaN;
%                          data_lens2.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd; %depth averaged value
%                      catch
%                          pause(60);
%                          system(['rm -f ', tmp.fname2]);
%                          system(['cdo -O -w ', cfg.vm_str, ' -sellevel,', cfg.vl_z1_str,'/',cfg.vl_z2_str, ' ', tmp.fname, ' ', tmp.fname2]);
%                          ncid = netcdf.open(tmp.fname2, 'NOWRITE');
%                          tmpvarid = netcdf.inqVarID(ncid, cfg.var);
%                          tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 ind_start-1], [grid.nlon grid.nlat ind_count]));
%                          tmp.dd(abs(tmp.dd)>1e30)=NaN;
%                          data_lens2.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd; %depth averaged value
%                      end
%         
%                 else
%                     ncid = netcdf.open(tmp.fname, 'NOWRITE');
%                     tmpvarid = netcdf.inqVarID(ncid, cfg.var);
%                     tmp.dd=  netcdf.getVar(ncid,tmpvarid, [0 0 ind_start-1], [grid.nlon grid.nlat ind_count]);
%                     tmp.dd(abs(tmp.dd)>1e30)=NaN;
%                     data_lens2.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd;
%                 end
%                 netcdf.close(ncid);
        
                ty=ty+ind_count/12;
            end
            tmp.ydata2_lens2=reshape(tmp.ydata_lens2, [grid.cut_nlon, grid.cut_nlat, 12, cfg.len_t_y]);
            tmp.ymean_lens2=squeeze(mean(tmp.ydata2_lens2,3));
            data.([tmp.varname, '_lens2'])(1:grid.cut_nlon,1:grid.cut_nlat,:)= tmp.ymean_lens2;

%             for mi=1:12
%                 tmp.m_str=num2str(mi,'%02i');
%                 cfg.lens2_fnm=[dirs.lens2dir, tmp.fs, ...
%                     tmp.varname, '_ensmean_', tmp.fy_str, '-', tmp.m_str, '.nc'];
%                 ncid=netcdf.open(cfg.lens2_fnm, 'NOWRITE');
%                 tmpvarid=netcdf.inqVarID(ncid,tmp.fvarname);
%                 if length(tmp.dimids)>3
%                      tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [grid.id_w grid.id_s cfg.vlayer_1st-1 0], [grid.cut_nlon grid.cut_nlat cfg.vlayer_cnt 1]));
%                      tmp.dd(abs(tmp.dd)>1e30)=NaN;
%                      ddd=mean(tmp.dd,3,'omitnan'); % depth mean
%                      tmp.ydata_lens2(1:grid.nlon,1:grid.nlat, mi) = ddd; %depth averaged value
%                 else
%                     tmp.dd=netcdf.getVar(ncid,tmpvarid);
%                      tmp.dd(abs(tmp.dd)>1e30)=NaN;                        
%                     tmp.ydata_lens2(1:grid.nlon,1:grid.nlat,mi) =tmp.dd;
%                 end
%                 netcdf.close(ncid);
%             end
%             tmp.ymean_lens2=mean(tmp.ydata_lens2,3);


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
        data.([tmp.varname, '_lens2_stdt'])= std(data.([tmp.varname, '_lens2']),0,3,'omitnan');

   %% get model error(bias) with obs 
        data.([tmp.varname, '_lens2_err_obs']) = ...
            data.([tmp.varname, '_lens2']) - data.([tmp.varname, '_obs']);
        data.([tmp.varname, '_assm_err_obs']) = ...
            data.([tmp.varname, '_assm']) - data.([tmp.varname, '_obs']);

   %% get RMSE with obs (drift-removed)
        data.([tmp.varname, '_lens2_rmse_obs']) = ...
            sqrt( mean( ...
            (data.([tmp.varname, '_lens2_err_obs']) ...
            - mean(data.([tmp.varname, '_lens2_err_obs']), 3, 'omitnan')).^2, ...
            3, 'omitnan') );
        data.([tmp.varname, '_assm_rmse_obs']) = ...
            sqrt( mean( ...
            (data.([tmp.varname, '_assm_err_obs']) ...
            - mean(data.([tmp.varname, '_assm_err_obs']), 3, 'omitnan')).^2, ...
            3, 'omitnan') );
    
    %% get hindcast climatology as the function of lead year
    data.([tmp.varname, '_assm_clim']) = mean(data.([tmp.varname, '_assm'])(:,:,:),3, 'omitnan'); %1964 ~ 2020
    data.([tmp.varname, '_lens2_clim']) = mean(data.([tmp.varname, '_lens2'])(:,:,:),3, 'omitnan'); %1964 ~ 2020
    data.([tmp.varname, '_obs_clim']) = mean(data.([tmp.varname, '_obs'])(:,:,:),3, 'omitnan'); %1964 ~ 2020
        
    disp('statistics calculated')
    fprintf('%7.1f sec\n', toc(lap_time) );

% % %     %% get svd modes (first 10 modes)
% % %         tmp.svd_modes=5;
% % %         % lens2 <-> hcst
% % %         [data2.([tmp.varname, '_svd_lens2_model_lvmap_left']), ...
% % %             data2.([tmp.varname, '_svd_lens2_model_pcs_left']), ...
% % %             data2.([tmp.varname, '_svd_lens2_model_lvmap_right']), ...
% % %             data2.([tmp.varname, '_svd_lens2_model_pcs_right']), ...
% % %             data2.([tmp.varname, '_svd_lens2_model_lambda']), ...
% % %             data2.([tmp.varname, '_svd_lens2_model_scf'])] = ...
% % %             mca( data2.([tmp.varname, '_lens2_ano']), ... 
% % %                  data2.([tmp.varname, '_model_ano']), tmp.svd_modes);
% % %         for orderi=1:tmp.svd_modes
% % %             data2.([tmp.varname, '_svd_order_lens2_model']) = ...
% % %                 order(mean(abs(data2.([tmp.varname, '_svd_lens2_model_pcs_left'])(orderi,:)),'all'));
% % %             data2.([tmp.varname, '_svd_lens2_model_lvmap_left'])(:,:,orderi) = ...
% % %                 data2.([tmp.varname, '_svd_lens2_model_lvmap_left'])(:,:,orderi) ...
% % %                 .* 10.^data2.([tmp.varname, '_svd_order_lens2_model']);
% % %             data2.([tmp.varname, '_svd_lens2_model_pcs_left'])(orderi,:) = ...
% % %                 data2.([tmp.varname, '_svd_lens2_model_pcs_left'])(orderi,:) ...
% % %                 .* 10.^(-data2.([tmp.varname, '_svd_order_lens2_model']));
% % %             data2.([tmp.varname, '_svd_lens2_model_lvmap_right'])(:,:,orderi) = ...
% % %                 data2.([tmp.varname, '_svd_lens2_model_lvmap_right'])(:,:,orderi) ...
% % %                 .* 10.^data2.([tmp.varname, '_svd_order_lens2_model']);
% % %             data2.([tmp.varname, '_svd_lens2_model_pcs_right'])(orderi,:) = ...
% % %                 data2.([tmp.varname, '_svd_lens2_model_pcs_right'])(orderi,:) ...
% % %                 .* 10.^(-data2.([tmp.varname, '_svd_order_lens2_model']));
% % %         end
% % %         tmp.pcs=data2.([tmp.varname, '_svd_lens2_model_pcs_right'])(1,:);
% % %         tmp.pcs=reshape(tmp.pcs,[1,1,length(tmp.pcs)]);
% % %         data2.([tmp.varname, '_svd_lens2_model_right_1mode_recon']) = ...
% % %             data2.([tmp.varname, '_svd_lens2_model_lvmap_right'])(:,:,1).*tmp.pcs(1,1,:);
% % % 
% % %         % lens2 <-> assm
% % %         [data2.([tmp.varname, '_svd_lens2_assm_lvmap_left']), ...
% % %             data2.([tmp.varname, '_svd_lens2_assm_pcs_left']), ...
% % %             data2.([tmp.varname, '_svd_lens2_assm_lvmap_right']), ...
% % %             data2.([tmp.varname, '_svd_lens2_assm_pcs_right']), ...
% % %             data2.([tmp.varname, '_svd_lens2_assm_lambda']), ...
% % %             data2.([tmp.varname, '_svd_lens2_assm_scf'])] = ...
% % %             mca( data2.([tmp.varname, '_lens2_ano']), ... 
% % %                  data2.([tmp.varname, '_assm_ano']), tmp.svd_modes);
% % %         for orderi=1:tmp.svd_modes
% % %             data2.([tmp.varname, '_svd_order_lens2_assm']) = ...
% % %                 order(mean(abs(data2.([tmp.varname, '_svd_lens2_assm_pcs_left'])(orderi,:)),'all'));
% % %             data2.([tmp.varname, '_svd_lens2_assm_lvmap_left'])(:,:,orderi) = ...
% % %                 data2.([tmp.varname, '_svd_lens2_assm_lvmap_left'])(:,:,orderi) ...
% % %                 .* 10.^data2.([tmp.varname, '_svd_order_lens2_assm']);
% % %             data2.([tmp.varname, '_svd_lens2_assm_pcs_left'])(orderi,:) = ...
% % %                 data2.([tmp.varname, '_svd_lens2_assm_pcs_left'])(orderi,:) ...
% % %                 .* 10.^(-data2.([tmp.varname, '_svd_order_lens2_assm']));
% % %             data2.([tmp.varname, '_svd_lens2_assm_lvmap_right'])(:,:,orderi) = ...
% % %                 data2.([tmp.varname, '_svd_lens2_assm_lvmap_right'])(:,:,orderi) ...
% % %                 .* 10.^data2.([tmp.varname, '_svd_order_lens2_assm']);
% % %             data2.([tmp.varname, '_svd_lens2_assm_pcs_right'])(orderi,:)= ...
% % %                 data2.([tmp.varname, '_svd_lens2_assm_pcs_right'])(orderi,:) ...
% % %                 .* 10.^(-data2.([tmp.varname, '_svd_order_lens2_assm']));
% % %         end
% % %         tmp.pcs(1,1,:)=data2.([tmp.varname, '_svd_lens2_assm_pcs_right'])(1,:);
% % %         tmp.pcs=reshape(tmp.pcs,[1,1,length(tmp.pcs)]);
% % %         data2.([tmp.varname, '_svd_lens2_assm_right_1mode_recon']) = ...
% % %             data2.([tmp.varname, '_svd_lens2_assm_lvmap_right'])(:,:,1).*tmp.pcs(1,1,:);
% % % 
% % %         % lens2 <-> obs
% % %         if isnan(sum( data2.([tmp.varname, '_obs_ano'])(:)))
% % %             data2.([tmp.varname, '_svd_lens2_obs_lvmap_left']) = ...
% % %                 NaN(size(data2.([tmp.varname, '_svd_lens2_assm_lvmap_left'])));
% % %             data2.([tmp.varname, '_svd_lens2_obs_pcs_left']) = ...
% % %                 NaN(size(data2.([tmp.varname, '_svd_lens2_assm_pcs_left'])));
% % %             data2.([tmp.varname, '_svd_lens2_obs_lvmap_right']) = ...
% % %                 NaN(size(data2.([tmp.varname, '_svd_lens2_assm_lvmap_left'])));
% % %             data2.([tmp.varname, '_svd_lens2_obs_pcs_right']) = ...
% % %                 NaN(size(data2.([tmp.varname, '_svd_lens2_assm_pcs_left'])));
% % %             data2.([tmp.varname, '_svd_lens2_obs_lambda']) = ...
% % %                 NaN(size(data2.([tmp.varname, '_svd_lens2_model_lambda'])));
% % %             data2.([tmp.varname, '_svd_lens2_obs_scf']) = ...
% % %                 NaN(size(data2.([tmp.varname, '_svd_lens2_model_scf'])));
% % %         else
% % %             [data2.([tmp.varname, '_svd_lens2_obs_lvmap_left']), ...
% % %                 data2.([tmp.varname, '_svd_lens2_obs_pcs_left']), ...
% % %                 data2.([tmp.varname, '_svd_lens2_obs_lvmap_right']), ...
% % %                 data2.([tmp.varname, '_svd_lens2_obs_pcs_right']), ...
% % %                 data2.([tmp.varname, '_svd_lens2_obs_lambda']), ...
% % %                 data2.([tmp.varname, '_svd_lens2_obs_scf'])] = ...
% % %                 mca( data2.([tmp.varname, '_lens2_ano']), ... 
% % %                      data2.([tmp.varname, '_obs_ano']), tmp.svd_modes);
% % %             for orderi=1:tmp.svd_modes
% % %                 data2.([tmp.varname, '_svd_order_lens2_obs']) = ...
% % %                     order(mean(abs(data2.([tmp.varname, '_svd_lens2_obs_pcs_left'])(orderi,:)),'all'));
% % %                 data2.([tmp.varname, '_svd_lens2_obs_lvmap_left'])(:,:,orderi) = ...
% % %                     data2.([tmp.varname, '_svd_lens2_obs_lvmap_left'])(:,:,orderi) ...
% % %                     .* 10.^data2.([tmp.varname, '_svd_order_lens2_obs']);
% % %                 data2.([tmp.varname, '_svd_lens2_obs_pcs_left'])(orderi,:) = ...
% % %                     data2.([tmp.varname, '_svd_lens2_obs_pcs_left'])(orderi,:) ...
% % %                     .* 10.^(-data2.([tmp.varname, '_svd_order_lens2_obs']));
% % %                 data2.([tmp.varname, '_svd_lens2_obs_lvmap_right'])(:,:,orderi) = ...
% % %                     data2.([tmp.varname, '_svd_lens2_obs_lvmap_right'])(:,:,orderi) ...
% % %                     .* 10.^data2.([tmp.varname, '_svd_order_lens2_obs']);
% % %                 data2.([tmp.varname, '_svd_lens2_obs_pcs_right'])(orderi,:) = ...
% % %                     data2.([tmp.varname, '_svd_lens2_obs_pcs_right'])(orderi,:)...
% % %                     .* 10.^(-data2.([tmp.varname, '_svd_order_lens2_obs']));
% % %             end
% % %         end
% % %       
% % %         tmp.pcs(1,1,:)=data2.([tmp.varname, '_svd_lens2_obs_pcs_right'])(1,:);
% % %         tmp.pcs=reshape(tmp.pcs,[1,1,length(tmp.pcs)]);
% % %         data2.([tmp.varname, '_svd_lens2_obs_right_1mode_recon']) = ...
% % %             data2.([tmp.varname, '_svd_lens2_obs_lvmap_right'])(:,:,1).*tmp.pcs(1,1,:);
% % % 
% % %     disp('svd finished')
% % %     fprintf('%7.1f sec\n', toc(lap_time) );
 
        %% get correlation coefficient
        for loni=1:grid.cut_nlon
            for lati=1:grid.cut_nlat
                 if (isnan(data.([tmp.varname, '_assm'])(loni,lati,1))~=1 && sum(data.([tmp.varname, '_lens2'])(loni,lati,:), 'omitnan')~=0)
    
                 %% corr lens2 ~ assm
                     tmp.lens2 = squeeze(data.([tmp.varname, '_lens2'])(loni,lati,:));
                     tmp.data_assm = squeeze(data.([tmp.varname, '_assm'])(loni,lati,:));
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.lens2(isfinite(tmp.data_assm)), tmp.data_assm(isfinite(tmp.data_assm)));

                     data.([tmp.varname, '_corr_assm_lens2'])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_assm_lens2_p'])(loni,lati)=tmp.corr_p(1,2);

                 %% corr obs ~ assm
                     tmp.data_assm = squeeze(data.([tmp.varname, '_assm'])(loni,lati,:));
                     tmp.data_obs = squeeze(data.([tmp.varname, '_obs'])(loni, lati,:));
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data_assm(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));

                     if sum(squeeze(isfinite(tmp.data_obs)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                     end

                     data.([tmp.varname, '_corr_obs_assm'])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_obs_assm_p'])(loni,lati)=tmp.corr_p(1,2);

                 %% corr obs ~ lens2
                        tmp.data_lens2 = squeeze(data.([tmp.varname, '_lens2'])(loni,lati,:));
                     [tmp.corr, tmp.corr_p]=corrcoef(tmp.data_lens2(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));

                     if sum(squeeze(isfinite(tmp.data_obs)))<10
                        tmp.corr=NaN(2,2);
                        tmp.corr_p=NaN(2,2);
                        tmp.corr_det=NaN(2,2);
                        tmp.corr_det_p=NaN(2,2);
                        tmp.data_lens2_trend=NaN;
                     end

                     data.([tmp.varname, '_corr_obs_lens2'])(loni,lati)=tmp.corr(1,2);
                     data.([tmp.varname, '_corr_obs_lens2_p'])(loni,lati)=tmp.corr_p(1,2);          
    
                 else
                     %% obs ~ ASSM (NaN)
                     data.([tmp.varname, '_corr_obs_assm'])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_assm_p'])(loni,lati)=NaN;
                     
                     %% obs ~ LENS2 (NaN)
                     data.([tmp.varname, '_corr_obs_lens2'])(loni,lati)=NaN;
                     data.([tmp.varname, '_corr_obs_lens2_p'])(loni,lati)=NaN;
                 end
            end
        end
        disp('correlation1 finished')
        fprintf('%7.1f sec\n', toc(lap_time) );

        if isfield(tmp,'ymean_mod_obs_masked')
            tmp=rmfield(tmp, 'ymean_mod_obs_masked');
        end
        
        if ~exist(dirs.matroot,'dir'), mkdir(dirs.matroot); end
        fig_cfg.mat_name=[dirs.matroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_lon',num2str(grid.regions(1)), '_', num2str(grid.regions(2)), '_', ...
            '_lat',num2str(grid.regions(3)), '_', num2str(grid.regions(4)), ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), '_', ...
            'obs_', cfg.obs_name, '.mat'];

        save(fig_cfg.mat_name, 'data', 'grid', 'cfg')
        clear data tmp
        fprintf('%7.1f sec\n', toc(lap_time) );
    
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
        case 'TEMP'
            obsname_simple='EN4';
%             obsname_simple='projdv7.3';
        case 'SALT'
            obsname_simple='EN4';
%             obsname_simple='projdv7.3';
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
        case 'TEMP'
            obsname_simple='TEMP';
        case 'SALT'
            obsname_simple='SALT';
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

