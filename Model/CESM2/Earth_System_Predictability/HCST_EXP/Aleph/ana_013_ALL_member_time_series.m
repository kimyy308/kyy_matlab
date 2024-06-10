% %  Created 21-Feb-2024 by Yong-Yub Kim

clc; clear all; close all;

%% set path
tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

%% model configuration
dirs.root='/mnt/lustre/proj/earth.system.predictability';
dirs.assm_root=[dirs.root, '/ASSM_EXP_timeseries/archive'];
dirs.hcst_root=[dirs.root, '/HCST_EXP_timeseries/archive'];
dirs.lens2_root='/proj/jedwards/archive';

dirs.saveroot=[dirs.root, '/statistics'];

cfg.iyears=1960:2020;
cfg.months=1:12;
cfg.scenname='HIST';
cfg.gridname='f09_g17';
cfg.proj_year=5;

grid.regions = [150 210 75 90];
grid.regions = [285 310 35 40]; % NASPG


% cfg.vars={'TS', 'SST', 'PRECT', 'TWS', 'GPP', 'FAREA_BURNED', 'photoC_TOT_zint_100m', ...
%     'COL_FIRE_CLOSS', 'SOILWATER_10CM', 'TLAI', 'SSH', 'PSL', 'AEROD_v'};
cfg.vars={'SST'};
cfg.vars={'TEMP', 'SALT'};


for vari=1:length(cfg.vars)
    cfg.var=cfg.vars{vari};
    tmp.varname=cfg.var;
    cfg.vlayer=1:37; % surf, vertical slice 
% % cfg.vlayer=1:10; % 10layer. don't put more than 15
% cfg.vlayer=10; % 100m, vertical slice 
% cfg.vlayer=20; % 200m, vertical slice 
% % cfg.vlayer=15; %150m
% cfg.vlayer=27; %305m
% cfg.vlayer=31; %408m
% cfg.vlayer=37; %707m

    cfg.vlayer_1st=min(cfg.vlayer);
    cfg.vlayer_cnt=max(cfg.vlayer)-cfg.vlayer_1st+1;
    cfg.vlayer_ind=cfg.vlayer_1st:cfg.vlayer_1st+cfg.vlayer_cnt-1;

% cfg.assm_factors={'10', '20'};
% cfg.ens_members={'1', '2', '3', '4', '5'};
% cfg.obsnames={'en4.2', 'projdv7.3'};  %{'oras4', 'projdv7.3', 'en4.2'}

    cfg.len_t_y=length(cfg.iyears);
    
    cfg.comp=Func_0025_CESM2_cmpname_var(cfg.var);
    cfg.obs_fname_module=f_obs_fname_module(cfg.comp);
    cfg.obs_name=f_obs_name(cfg.var);
    cfg.obs_fname_mid=f_obs_name_mid(cfg.var);
    cfg.obs_varname=f_obs_varname(cfg.var);

    dirs.obsroot=['/mnt/lustre/proj/kimyy/Observation/', cfg.obs_name];
    dirs.grid_root=['/mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive_transfer/', cfg.comp];

%% set grid
    tmp.gridname = [dirs.grid_root, tmp.fs, 'grid.nc'];
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
            cfg.vl_z1=grid.z_t(cfg.vlayer_1st);
            cfg.vl_z2=grid.z_t(cfg.vlayer_1st+cfg.vlayer_cnt-1);
            cfg.vl_z1_str=num2str(cfg.vl_z1*100);
            cfg.vl_z2_str=num2str(cfg.vl_z2*100);
            if cfg.vl_z1==cfg.vl_z2
                cfg.vm_str=' ';
            else
                cfg.vm_str='-vertmean';
            end
    
        case {'atm', 'lnd'}
            grid.lon=ncread(tmp.gridname, 'lon');
            grid.lat=ncread(tmp.gridname, 'lat');
            [grid.tlat, grid.tlong]=meshgrid(grid.lat, grid.lon);
            grid.area=ncread(tmp.gridname, 'AREA');cfg.vlayer_1st=min(cfg.vlayer);
            cfg.vlayer_cnt=max(cfg.vlayer)-cfg.vlayer_1st+1;
            grid.lfrac=ncread(tmp.gridname, 'LANDFRAC');
            grid.area=grid.area.*grid.lfrac;
    end

    grid.nlon=size(grid.tlong,1);
    grid.nlat=size(grid.tlat,2);
    grid.dz_3d=repmat(grid.dz_res, [grid.nlon, grid.nlat, cfg.vlayer_cnt]);

    [grid.id_w, grid.id_e, grid.id_s, grid.id_n] = Func_0012_findind_Y(1.0, grid.regions, ...
                grid.tlong, grid.tlat, 'CESM2'); % find valid lon, lat index near station
    if length(grid.regions)==2
        tmp.idlon=round((grid.id_w+grid.id_e)/2);
        tmp.idlat=round((grid.id_s+grid.id_n)/2);
        grid.id_w=tmp.idlon;
        grid.id_e=tmp.idlon;
        grid.id_s=tmp.idlat;
        grid.id_n=tmp.idlat;
    end
    grid.cut_tlong=grid.tlong(grid.id_w:grid.id_e, grid.id_s:grid.id_n);
    grid.cut_tlat=grid.tlat(grid.id_w:grid.id_e, grid.id_s:grid.id_n);
    grid.cut_nlon=size(grid.cut_tlong,1);
    grid.cut_nlat=size(grid.cut_tlat,2);
    grid.cut_dz_3d=grid.dz_3d(grid.id_w:grid.id_e, grid.id_s:grid.id_n, cfg.vlayer);


%% read OBS
cfg.max_iy=max(cfg.iyears);
for ty=1:length(cfg.iyears)
    tmp.year=cfg.iyears(ty);
    tmp.fy=tmp.year;
    tmp.fy_str=num2str(tmp.fy);
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

% 
        if exist(cfg.obs_fnm)~=0
            ncid=netcdf.open(cfg.obs_fnm, 'NOWRITE');
            tmpvarid=netcdf.inqVarID(ncid,cfg.obs_varname);
            tmp.dd =  netcdf.getVar(ncid,tmpvarid, [grid.id_w-1 grid.id_s-1 cfg.vlayer_1st-1 0], [grid.cut_nlon grid.cut_nlat cfg.vlayer_cnt 1]);
            tmp.dd=double(tmp.dd);
            if (strcmp(cfg.obs_name, 'ERA5')==1 || strcmp(cfg.var, 'TLAI')==1)
                tmp.add_offset=netcdf.getAtt(ncid,tmpvarid,'add_offset');
                tmp.scale_factor=netcdf.getAtt(ncid,tmpvarid,'scale_factor');
                tmp.dd=tmp.dd.*tmp.scale_factor+tmp.add_offset;
            end
             tmp.dd(abs(tmp.dd)>1e30)=NaN;
            %% unit correction
            if strcmp(cfg.obs_name, 'GPCC')==1
                if strcmp(cfg.var, 'PRECT')
                    tmp.dd=tmp.dd./1000.0./86400/eomday(tmp.fy,mi);
                elseif strcmp(cfg.var, 'RAIN')
                    tmp.dd=tmp.dd./86400/eomday(tmp.fy,mi);
                end
            elseif (strcmp(cfg.obs_name, 'CMEMS')==1 && strcmp(cfg.var, 'SOILWATER_10CM')==1)
                tmp.dd=tmp.dd.*1000.*(10./3);
            elseif (strcmp(cfg.obs_name, 'CMEMS')==1 && strcmp(cfg.var, 'SSH')==1)
                tmp.dd=tmp.dd./100; % cm -> m
            elseif (strcmp(cfg.obs_name, 'GLEAM')==1 && strcmp(cfg.var, 'SOILWATER_10CM')==1)
                tmp.dd=tmp.dd.*(1000)./10; % 1. m^3 -> kg, 2. 100cm(1m) -> 10cm,  m3/m3 -> 10cm(surface) soil kg/m2
            elseif (strcmp(cfg.obs_name, 'NOAA')==1 && strcmp(cfg.var, 'TWS')==1)
                tmp.dd(tmp.dd<=0)=NaN;
                tmp.dd(tmp.dd>740)=NaN;
            elseif strcmp(cfg.var, 'FAREA_BURNED')
%                                 tmp.dd=tmp.dd./86400./eomday(tmp.fy,mi)./grid.area;
                    tmp.dd=tmp.dd./86400./grid.area;
            elseif strcmp(cfg.var, 'COL_FIRE_CLOSS')
                tmp.dd=tmp.dd./86400./eomday(tmp.fy,mi);
            end
%                 tmp.ydata_obs(1:grid.cut_nlon,1:grid.cut_nlat,mi) = tmp.dd;

            if cfg.vlayer_cnt==1 %% 2D or 3D
                data_obs.([cfg.var])(:,:,(ty-1)*12+mi) = tmp.dd;
            else
                data_obs.([cfg.var])(:,:,:,(ty-1)*12+mi) = tmp.dd;
            end
            netcdf.close(ncid);
            
        else
            if cfg.vlayer_cnt==1 %% 2D
                data_obs.([cfg.var])(:,:,(ty-1)*12+mi) = NaN(grid.cut_nlon, grid.cut_nlat);
            else %% 3D
                data_obs.([cfg.var])(:,:,:,(ty-1)*12+mi) = NaN(grid.cut_nlon, grid.cut_nlat, cfg.vlayer_cnt);
            end
        end
    end
end
% yearly mean
if cfg.vlayer_cnt==1 %% 2D or 3D
    tmp.reshp=reshape(data_obs.([cfg.var]), [grid.cut_nlon, grid.cut_nlat, 12, cfg.len_t_y]);
    data_obs.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,3));
else
    tmp.reshp=reshape(data_obs.([cfg.var]), [grid.cut_nlon, grid.cut_nlat, cfg.vlayer_cnt, 12, cfg.len_t_y]);
    data_obs.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,4));
end

%% spatial mean
if cfg.vlayer_cnt==1 %% 2D
    if length(grid.regions)==4
    [m_data_obs.([cfg.var]), tmp.err] = ...
        Func_0011_get_area_weighted_mean(data_obs.([cfg.var]), grid.cut_tlong, grid.cut_tlat);
    [m_data_obs.([cfg.var, '_ym']), tmp.err] = ...
        Func_0011_get_area_weighted_mean(data_obs.([cfg.var, '_ym']), grid.cut_tlong, grid.cut_tlat);
    elseif length(grid.regions)==2
        m_data_obs.([cfg.var]) = data_obs.([cfg.var]);
        m_data_obs.([cfg.var, '_ym']) = data_obs.([cfg.var, '_ym']);
    end
else %% 3D
    grid.cut_dz_3d_masked=grid.cut_dz_3d .* ...
        data_obs.([cfg.var, '_ym'])(:,:,:,1) ./ data_obs.([cfg.var, '_ym'])(:,:,:,1);
    grid.cut_dz_3d_norm=grid.cut_dz_3d_masked./ sum(grid.cut_dz_3d_masked,3,'omitnan');
    data_obs_zint.([cfg.var])=squeeze(sum(data_obs.([cfg.var]).*grid.cut_dz_3d_norm, 3, 'omitnan'));
    data_obs_zint.([cfg.var, '_ym'])=squeeze(sum(data_obs.([cfg.var, '_ym']).*grid.cut_dz_3d_norm, 3, 'omitnan'));

    if length(grid.regions)==4
    [m_data_obs.([cfg.var]), tmp.err] = ...
        Func_0011_get_area_weighted_mean(data_obs_zint.([cfg.var]), grid.cut_tlong, grid.cut_tlat);
    [m_data_obs.([cfg.var, '_ym']), tmp.err] = ...
        Func_0011_get_area_weighted_mean(data_obs_zint.([cfg.var, '_ym']), grid.cut_tlong, grid.cut_tlat);
    elseif length(grid.regions)==2
        m_data_obs.([cfg.var]) = data_obs_zint.([cfg.var]);
        m_data_obs.([cfg.var, '_ym']) = data_obs_zint.([cfg.var, '_ym']);
    end
end

%% get ASSM members
cfg_assm.list = dir( [dirs.assm_root, '/', '*BHISTsmbb*' ]); 
for kk = 1 : length( cfg_assm.list )
    tmp.fname_in    = cfg_assm.list(kk).name;
    tmp.fname_split = strsplit( tmp.fname_in, {'.'} );
    cfg_assm.members{kk}=[ tmp.fname_split{6}, '.', tmp.fname_split{7}];  
end      
cfg_assm.len_mem=length(cfg_assm.members);

%% get HCST members
cfg_hcst.list = dir( [dirs.hcst_root, '/', '*f09*' ]); 
for kk = 1 : length( cfg_hcst.list )
    tmp.fname_in    = cfg_hcst.list(kk).name;
    tmp.fname_split = strsplit( tmp.fname_in, {'.'} );
    cfg_hcst.members{kk}=[ tmp.fname_split{3}, '.', tmp.fname_split{4}];     
end
cfg_hcst.len_mem=length(cfg_hcst.members);

%% get LENS2 members
cfg_lens2.list = dir( [dirs.lens2_root, '/', '*BHISTsmbb*' ]); 
for kk = 1 : length( cfg_lens2.list )
    tmp.fname_in    = cfg_lens2.list(kk).name;
    tmp.fname_split = strsplit( tmp.fname_in, {'.'} );
    cfg_lens2.members{kk}=[ tmp.fname_split{5}, '.', tmp.fname_split{6}];     
end
cfg_lens2.len_mem=length(cfg_lens2.members);


%% read ASSM members
tic;
if cfg.vlayer_cnt==1 %% 2D or 3D
    data_assm.([cfg.var])=NaN(cfg_assm.len_mem, grid.cut_nlon, grid.cut_nlat, cfg.len_t_y*12); % initialization
else
    data_assm.([cfg.var])=NaN(cfg_assm.len_mem, grid.cut_nlon, grid.cut_nlat, cfg.vlayer_cnt, cfg.len_t_y*12); % initialization
end
if isfield(tmp, 'dimids')==1
    tmp=rmfield(tmp, 'dimids');
end
for mi= 1:cfg_assm.len_mem
    ty=1;
    tmp.member=cfg_assm.members{mi};
    disp(['assm, ', tmp.member]);
    while ty<=cfg.len_t_y
        tmp.year=cfg.iyears(ty);
        tmp.scen = f_scen(tmp.year);
        tmp.casenm = ['b.e21.', tmp.scen, '.f09_g17.assm.',tmp.member];
        tmp.fdir=[dirs.assm_root, '/', tmp.casenm, '/' cfg.comp, '/proc/tseries/month_1'];
        tmp.flist=dir( [tmp.fdir, '/','*.',cfg.var, '.*'] );
        for kk = 1: length (tmp.flist)
            tmp.fname_in = tmp.flist(kk).name;
            tmp.fname_split = strsplit(tmp.fname_in, {'.'} );
            tmp.fname_period = tmp.fname_split{11};
            fyear_str   = strsplit( tmp.fname_period, '-' );
            fyear_start = str2num( fyear_str{1}(1:4) );
            fyear_end   = str2num( fyear_str{2}(1:4) );
            if( tmp.year >= fyear_start && tmp.year <= fyear_end )
%                 disp(tmp.fname_period);
                flag_file_in = true;            break;
            end
        end
        tmp.fname=[tmp.fdir, '/', tmp.casenm,  cfg.obs_fname_module, cfg.var, '.', tmp.fname_period, '.nc'];
        ind_start=(tmp.year - fyear_start)*12+1;
        if fyear_end<=max(cfg.iyears)
            ind_count=(fyear_end-tmp.year+1)*12;
        else
            ind_count=(max(cfg.iyears)-tmp.year+1)*12;
        end
        
        ncid = netcdf.open(tmp.fname, 'NOWRITE');
        tmpvarid = netcdf.inqVarID(ncid, cfg.var);
        if ~isfield(tmp, 'dimids')
            [tmp.fvarname,tmp.xtype,tmp.dimids,tmp.natts]= netcdf.inqVar(ncid,tmpvarid);
        end
        if length(tmp.dimids)>3
            tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, ...
                     [grid.id_w-1 grid.id_s-1 cfg.vlayer_1st-1 ind_start-1], ...
                     [grid.cut_nlon grid.cut_nlat cfg.vlayer_cnt ind_count]));
            tmp.dd(abs(tmp.dd)>1e30)=NaN;
            if cfg.vlayer_cnt==1 %% 2D or 3D
                data_assm.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd;
            else
                data_assm.([cfg.var])(mi,:,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd; %depth averaged value                
            end
        else
            tmp.dd=  netcdf.getVar(ncid,tmpvarid, ...
                [grid.id_w-1 grid.id_s-1 ind_start-1], ...
                [grid.cut_nlon grid.cut_nlat ind_count]);
            tmp.dd(abs(tmp.dd)>1e30)=NaN;
            data_assm.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd;
        end
        netcdf.close(ncid);
        ty=ty+ind_count/12;
    end
end

%% yearly mean
if cfg.vlayer_cnt==1 %% 2D or 3D
    tmp.reshp=reshape(data_assm.([cfg.var]), [cfg_assm.len_mem, grid.cut_nlon, grid.cut_nlat, 12, cfg.len_t_y]);
    data_assm.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,4));
% data_assm_em.([cfg.var, '_ym'])=squeeze(mean(data_assm.([cfg.var, '_ym']),1));
else %% 3D
    tmp.reshp=reshape(data_assm.([cfg.var]), [cfg_assm.len_mem, grid.cut_nlon, grid.cut_nlat, cfg.vlayer_cnt, 12, cfg.len_t_y]);
    data_assm.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,5));
end

%% spatial mean
if cfg.vlayer_cnt==1 %% 2D
    if length(grid.regions)==4
        for mi=1:size(data_assm.([cfg.var]),1)
            [m_data_assm.([cfg.var])(mi,:), tmp.err] = ...
                Func_0011_get_area_weighted_mean(squeeze(data_assm.([cfg.var])(mi,:,:,:)), grid.cut_tlong, grid.cut_tlat);
            [m_data_assm.([cfg.var, '_ym'])(mi,:), tmp.err] = ...
                Func_0011_get_area_weighted_mean(squeeze(data_assm.([cfg.var, '_ym'])(mi,:,:,:)), grid.cut_tlong, grid.cut_tlat);
        end
    elseif length(grid.regions)==2
        m_data_assm.([cfg.var]) = data_assm.([cfg.var]);
        m_data_assm.([cfg.var, '_ym']) = data_assm.([cfg.var, '_ym']);
    end
else %% 3D
    grid.cut_dz_3d_masked=grid.cut_dz_3d .* ...
        data_assm.([cfg.var, '_ym'])(1,:,:,:,1) ./ data_assm.([cfg.var, '_ym'])(1,:,:,:,1);
    grid.cut_dz_3d_norm=grid.cut_dz_3d_masked./ sum(grid.cut_dz_3d_masked,4,'omitnan');
    grid.cut_dz_4d_norm=repmat(grid.cut_dz_3d_norm, [size(data_assm.([cfg.var, '_ym']),1). size(grid.cut_dz_3d_norm)]);
    data_assm_zint.([cfg.var])=squeeze(sum(data_assm.([cfg.var]).*grid.cut_dz_4d_norm, 4, 'omitnan'));
    data_assm_zint.([cfg.var, '_ym'])=squeeze(sum(data_assm.([cfg.var, '_ym']).*grid.cut_dz_4d_norm, 3, 'omitnan'));
    
    if length(grid.regions)==4
        for mi=1:size(data_assm.([cfg.var]),1)
            [m_data_assm.([cfg.var])(mi,:), tmp.err] = ...
                Func_0011_get_area_weighted_mean(data_assm_zint.([cfg.var])(mi,:,:,:), grid.cut_tlong, grid.cut_tlat);
            [m_data_assm.([cfg.var, '_ym'])(mi,:), tmp.err] = ...
                Func_0011_get_area_weighted_mean(data_assm_zint.([cfg.var, '_ym'])(mi,:,:,:), grid.cut_tlong, grid.cut_tlat);
        end
    elseif length(grid.regions)==2
        m_data_assm.([cfg.var]) = data_assm_zint.([cfg.var]);
        m_data_assm.([cfg.var, '_ym']) = data_assm_zint.([cfg.var, '_ym']);
    end
end
toc;


%% read LENS2 members
tic;
if cfg.vlayer_cnt==1 %% 2D or 3D
    data_lens2.([cfg.var])=NaN(cfg_lens2.len_mem, grid.cut_nlon, grid.cut_nlat, cfg.len_t_y*12); % initialization
else
    data_lens2.([cfg.var])=NaN(cfg_lens2.len_mem, grid.cut_nlon, grid.cut_nlat, cfg.vlayer_cnt, cfg.len_t_y*12); % initialization
end
for mi= 1:cfg_lens2.len_mem
    ty=1;
    tmp.member=cfg_lens2.members{mi};
    disp(['lens2, ', tmp.member]);
    while ty<=cfg.len_t_y
        tmp.year=cfg.iyears(ty);
        tmp.scen = f_scen(tmp.year);
        tmp.casenm = ['b.e21.', tmp.scen, '.f09_g17.',tmp.member];
        tmp.fdir=[dirs.lens2_root, '/', tmp.casenm, '/' cfg.comp, '/proc/tseries/month_1'];
        tmp.flist=dir( [tmp.fdir, '/','*.',cfg.var, '.*'] );
        for kk = 1: length (tmp.flist)
            tmp.fname_in = tmp.flist(kk).name;
            tmp.fname_split = strsplit(tmp.fname_in, {'.'} );
            tmp.fname_period = tmp.fname_split{10};
            fyear_str   = strsplit( tmp.fname_period, '-' );
            fyear_start = str2num( fyear_str{1}(1:4) );
            fyear_end   = str2num( fyear_str{2}(1:4) );
            if( tmp.year >= fyear_start && tmp.year <= fyear_end )
%                 disp(tmp.fname_period);
                flag_file_in = true;            break;
            end
        end
        tmp.fname=[tmp.fdir, '/', tmp.casenm, cfg.obs_fname_module, cfg.var, '.', tmp.fname_period, '.nc'];
        ind_start=(tmp.year - fyear_start)*12+1;
        if fyear_end<=max(cfg.iyears)
            ind_count=(fyear_end-tmp.year+1)*12;
        else
            ind_count=(max(cfg.iyears)-tmp.year+1)*12;
        end
        
        ncid = netcdf.open(tmp.fname, 'NOWRITE');
        tmpvarid = netcdf.inqVarID(ncid, cfg.var);
        if length(tmp.dimids)>3
            tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, ...
                     [grid.id_w-1 grid.id_s-1 cfg.vlayer_1st-1 ind_start-1], ...
                     [grid.cut_nlon grid.cut_nlat cfg.vlayer_cnt ind_count]));
            tmp.dd(abs(tmp.dd)>1e30)=NaN;
            if cfg.vlayer_cnt==1 %% 2D or 3D
                data_lens2.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd;
            else
                data_lens2.([cfg.var])(mi,:,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd; %depth averaged value                
            end
        else
            tmp.dd=  netcdf.getVar(ncid,tmpvarid, ...
                [grid.id_w-1 grid.id_s-1 ind_start-1], ...
                [grid.cut_nlon grid.cut_nlat ind_count]);
            tmp.dd(abs(tmp.dd)>1e30)=NaN;
            data_lens2.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd;
        end
        netcdf.close(ncid);
        ty=ty+ind_count/12;
    end
end
%% yearly mean
if cfg.vlayer_cnt==1 %% 2D
    tmp.reshp=reshape(data_lens2.([cfg.var]), [cfg_lens2.len_mem, grid.cut_nlon, grid.cut_nlat, 12, cfg.len_t_y]);
    data_lens2.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,4));
% data_lens2_em.([cfg.var, '_ym'])=squeeze(mean(data_lens2.([cfg.var, '_ym']),1));
else %% 3D
    tmp.reshp=reshape(data_lens2.([cfg.var]), [cfg_lens2.len_mem, grid.cut_nlon, grid.cut_nlat, cfg.vlayer_cnt, 12, cfg.len_t_y]);
    data_lens2.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,5));
end

% %% spatial mean
% if length(grid.regions)==4
%     for mi=1:size(data_lens2.([cfg.var]),1)
%         [m_data_lens2.([cfg.var])(mi,:), tmp.err] = ...
%             Func_0011_get_area_weighted_mean(squeeze(data_lens2.([cfg.var])(mi,:,:,:)), grid.cut_tlong, grid.cut_tlat);
%         [m_data_lens2.([cfg.var, '_ym'])(mi,:), tmp.err] = ...
%             Func_0011_get_area_weighted_mean(squeeze(data_lens2.([cfg.var, '_ym'])(mi,:,:,:)), grid.cut_tlong, grid.cut_tlat);
%     end
% elseif length(grid.regions)==2
%     m_data_lens2.([cfg.var]) = data_lens2.([cfg.var]);
%     m_data_lens2.([cfg.var, '_ym']) = data_lens2.([cfg.var, '_ym']);
% end

%% spatial mean
if cfg.vlayer_cnt==1 %% 2D
    if length(grid.regions)==4
        for mi=1:size(data_lens2.([cfg.var]),1)
            [m_data_lens2.([cfg.var])(mi,:), tmp.err] = ...
                Func_0011_get_area_weighted_mean(squeeze(data_lens2.([cfg.var])(mi,:,:,:)), grid.cut_tlong, grid.cut_tlat);
            [m_data_lens2.([cfg.var, '_ym'])(mi,:), tmp.err] = ...
                Func_0011_get_area_weighted_mean(squeeze(data_lens2.([cfg.var, '_ym'])(mi,:,:,:)), grid.cut_tlong, grid.cut_tlat);
        end
    elseif length(grid.regions)==2
        m_data_lens2.([cfg.var]) = data_lens2.([cfg.var]);
        m_data_lens2.([cfg.var, '_ym']) = data_lens2.([cfg.var, '_ym']);
    end
else %% 3D
    grid.cut_dz_3d_masked=grid.cut_dz_3d .* ...
        data_lens2.([cfg.var, '_ym'])(1,:,:,:,1) ./ data_lens2.([cfg.var, '_ym'])(1,:,:,:,1);
    grid.cut_dz_3d_norm=grid.cut_dz_3d_masked./ sum(grid.cut_dz_3d_masked,4,'omitnan');
    grid.cut_dz_4d_norm=repmat(grid.cut_dz_3d_norm, [size(data_lens2.([cfg.var, '_ym']),1). size(grid.cut_dz_3d_norm)]);
    data_lens2_zint.([cfg.var])=squeeze(sum(data_lens2.([cfg.var]).*grid.cut_dz_4d_norm, 4, 'omitnan'));
    data_lens2_zint.([cfg.var, '_ym'])=squeeze(sum(data_lens2.([cfg.var, '_ym']).*grid.cut_dz_4d_norm, 3, 'omitnan'));
    
    if length(grid.regions)==4
        for mi=1:size(data_lens2.([cfg.var]),1)
            [m_data_lens2.([cfg.var])(mi,:), tmp.err] = ...
                Func_0011_get_area_weighted_mean(data_lens2_zint.([cfg.var])(mi,:,:,:), grid.cut_tlong, grid.cut_tlat);
            [m_data_lens2.([cfg.var, '_ym'])(mi,:), tmp.err] = ...
                Func_0011_get_area_weighted_mean(data_lens2_zint.([cfg.var, '_ym'])(mi,:,:,:), grid.cut_tlong, grid.cut_tlat);
        end
    elseif length(grid.regions)==2
        m_data_lens2.([cfg.var]) = data_lens2_zint.([cfg.var]);
        m_data_lens2.([cfg.var, '_ym']) = data_lens2_zint.([cfg.var, '_ym']);
    end
end
toc;



%% read HCST members

for ly=1:5 %1:5
    tmp.ly_str=['ly',num2str(ly)];
    if cfg.vlayer_cnt==1 %% 2D or 3D
        data_hcst.([cfg.var])=NaN(cfg_hcst.len_mem, grid.cut_nlon, grid.cut_nlat, cfg.len_t_y*12); % initialization
    else
        data_hcst.([cfg.var])=NaN(cfg_hcst.len_mem, grid.cut_nlon, grid.cut_nlat, cfg.vlayer_cnt, cfg.len_t_y*12); % initialization        
    end
    for mi= 1:cfg_hcst.len_mem
        ty=1;
        tmp.member=cfg_hcst.members{mi};
        disp(['hcst, ', tmp.member]);
        while ty<=cfg.len_t_y
            tmp.year=cfg.iyears(ty);
            tmp.yearstr=num2str(tmp.year);
            tmp.casenm = ['f09_g17.hcst.',tmp.member];
            tmp.casenm_i = [tmp.casenm, '_i', tmp.yearstr];
            tmp.fdir=[dirs.hcst_root, '/', tmp.casenm, '/' tmp.casenm_i, '/', cfg.comp, '/proc/tseries/month_1'];
            tmp.flist=dir( [tmp.fdir, '/','*.',cfg.var, '.*'] );
            tmp.fname=[tmp.fdir, '/', tmp.flist(1).name];
           
            ind_start=(ly - 1)*12+1;
            ind_count=12;
            
            ncid = netcdf.open(tmp.fname, 'NOWRITE');
            tmpvarid = netcdf.inqVarID(ncid, cfg.var);
            if length(tmp.dimids)>3
                tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, ...
                         [grid.id_w-1 grid.id_s-1 cfg.vlayer_1st-1 ind_start-1], ...
                         [grid.cut_nlon grid.cut_nlat cfg.vlayer_cnt ind_count]));
                tmp.dd(abs(tmp.dd)>1e30)=NaN;
                if cfg.vlayer_cnt==1 %% 2D or 3D
                    data_hcst.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd;
                else
                    data_hcst.([cfg.var])(mi,:,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd; %depth averaged value                
                end
            else
                tmp.dd=  netcdf.getVar(ncid,tmpvarid, ...
                    [grid.id_w-1 grid.id_s-1 ind_start-1], ...
                    [grid.cut_nlon grid.cut_nlat ind_count]);
                tmp.dd(abs(tmp.dd)>1e30)=NaN;
                data_hcst.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd;
            end
            netcdf.close(ncid);
    
            ty=ty+ind_count/12;
        end
    end
    % yearly mean
    if cfg.vlayer_cnt==1 %% 2D
        tmp.reshp=reshape(data_hcst.([cfg.var]), [cfg_hcst.len_mem, grid.cut_nlon, grid.cut_nlat, 12, cfg.len_t_y]);
        data_hcst.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,4));
%     data_hcst_em.([cfg.var, '_ym'])=squeeze(mean(data_hcst.([cfg.var, '_ym']),1));
    else
        tmp.reshp=reshape(data_hcst.([cfg.var]), [cfg_hcst.len_mem, grid.cut_nlon, grid.cut_nlat, cfg.vlayer_cnt, 12, cfg.len_t_y]);
        data_hcst.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,5));
    end
    
    %% spatial mean
    if cfg.vlayer_cnt==1 %% 2D
        if length(grid.regions)==4
            for mi=1:size(data_hcst.([cfg.var]),1)
                [m_data_hcst.([cfg.var]).(tmp.ly_str)(mi,:), tmp.err] = ...
                    Func_0011_get_area_weighted_mean(squeeze(data_hcst.([cfg.var])(mi,:,:,:)), grid.cut_tlong, grid.cut_tlat);
                [m_data_hcst.([cfg.var, '_ym']).(tmp.ly_str)(mi,:), tmp.err] = ...
                    Func_0011_get_area_weighted_mean(squeeze(data_hcst.([cfg.var, '_ym'])(mi,:,:,:)), grid.cut_tlong, grid.cut_tlat);
            end
        elseif length(grid.regions)==2
            m_data_hcst.([cfg.var]).(tmp.ly_str) = data_hcst.([cfg.var]);
            m_data_hcst.([cfg.var, '_ym']).(tmp.ly_str) = data_hcst.([cfg.var, '_ym']);
        end
    else %% 3D
        grid.cut_dz_3d_masked=grid.cut_dz_3d .* ...
            data_hcst.([cfg.var, '_ym'])(1,:,:,:,1) ./ data_hcst.([cfg.var, '_ym'])(1,:,:,:,1);
        grid.cut_dz_3d_norm=grid.cut_dz_3d_masked./ sum(grid.cut_dz_3d_masked,4,'omitnan');
        grid.cut_dz_4d_norm=repmat(grid.cut_dz_3d_norm, [size(data_hcst.([cfg.var, '_ym']),1). size(grid.cut_dz_3d_norm)]);
        data_hcst_zint.([cfg.var])=squeeze(sum(data_hcst.([cfg.var]).*grid.cut_dz_4d_norm, 4, 'omitnan'));
        data_hcst_zint.([cfg.var, '_ym'])=squeeze(sum(data_hcst.([cfg.var, '_ym']).*grid.cut_dz_4d_norm, 3, 'omitnan'));
        
        if length(grid.regions)==4
            for mi=1:size(data_hcst.([cfg.var]),1)
                [m_data_hcst.([cfg.var]).(tmp.ly_str)(mi,:), tmp.err] = ...
                    Func_0011_get_area_weighted_mean(data_hcst_zint.([cfg.var])(mi,:,:,:), grid.cut_tlong, grid.cut_tlat);
                [m_data_hcst.([cfg.var, '_ym']).(tmp.ly_str)(mi,:), tmp.err] = ...
                    Func_0011_get_area_weighted_mean(data_hcst_zint.([cfg.var, '_ym'])(mi,:,:,:), grid.cut_tlong, grid.cut_tlat);
            end
        elseif length(grid.regions)==2
            m_data_hcst.([cfg.var]).(tmp.ly_str) = data_hcst_zint.([cfg.var]);
            m_data_hcst.([cfg.var, '_ym']).(tmp.ly_str) = data_hcst_zint.([cfg.var, '_ym']);
        end
    end
    toc;
    
end




%% save matfile
if length(grid.regions)==2
    grid.regions2(1)=grid.regions(1); grid.regions2(2)=grid.regions(1);
    grid.regions2(3)=grid.regions(2); grid.regions2(4)=grid.regions(2);
else
    grid.regions2=grid.regions;
end

mkdir([dirs.saveroot, '/ts_data']);
matfilename=[dirs.saveroot, '/ts_data/', 'ts_data_all_',cfg.var, ...
            '_lon',num2str(grid.regions2(1)), '_', num2str(grid.regions2(2)), '_', ...
            '_lat',num2str(grid.regions2(3)), '_', num2str(grid.regions2(4)), ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), '_', ...
            'obs_', cfg.obs_name, '.mat'];
save(matfilename, 'cfg_assm', 'cfg_hcst', 'cfg_lens2', ...
    'data_obs', 'm_data_obs', 'm_data_assm', ...
    'm_data_hcst', 'm_data_lens2', 'grid');

% save(matfilename, 'cfg_assm', 'cfg_hcst', 'cfg_lens2', 'corrval_hcst', 'grid');
% matfilename=[dirs.saveroot, '/corr_raw/', 'corr_assm_',cfg.var,'.mat'];
% save(matfilename, 'cfg_assm', 'cfg_hcst', 'cfg_lens2', 'corrval_assm', 'grid');
% matfilename=[dirs.saveroot, '/corr_raw/', 'corr_lens2_',cfg.var,'.mat'];
% save(matfilename, 'cfg_assm', 'cfg_hcst', 'cfg_lens2', 'corrval_lens2', 'grid');

end

function scenname = f_scen(year)
    if year<=2014
        scenname='BHISTsmbb';
    else
        scenname='BSSP370smbb';
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
            obsname_simple='ersst_reg_5deg.v5.';
        case 'PRECT'
            obsname_simple='GPCC_reg_5deg.v5.';
        case 'RAIN'
            obsname_simple='GPCC_reg_5deg.v5.';
        case 'PSL'
            obsname_simple='ERA5_msl_reg_5deg.';
        case 'SOILWATER_10CM'
%             obsname_simple='SM_reg_5deg.';
            obsname_simple='GLEAM_reg_5deg.v5.';
        case 'TWS'
            obsname_simple='TSW_reg_5deg.';
        case 'SSH'
            obsname_simple='CMEMS_reg_5deg.';
        case 'TS'
%             obsname_simple='HadCRUT5_reg_5deg.';
            obsname_simple='ERA5_t2m_reg_5deg.';
        case 'sumChl'
            obsname_simple='OC_CCI_reg_5deg.';
        case 'TLAI'
            obsname_simple='LAI_reg_5deg.'; % NOAA LAI
        case 'FAREA_BURNED'
            obsname_simple='burned_area_reg_5deg.'; % MODIS Fire_cci v5.1
            obsname_simple='AVHRR-LTDR_reg_5deg.'; % AVHRR-LTDR
        case 'COL_FIRE_CLOSS'
            obsname_simple='FIRE_CLOSS_reg_5deg.'; % GFED
        case 'photoC_TOT_zint'
%             obsname_simple='s_vgpm_reg_5deg.'; % VGPM
%             obsname_simple='ensmean_reg_5deg.'; % VGPM   
            obsname_simple='CMEMS_reg_5deg.'; %Globcolour;
        case 'photoC_TOT_zint_100m'
%             obsname_simple='s_vgpm_reg_5deg.'; % VGPM
%             obsname_simple='ensmean_reg_5deg.'; % VGPM
            obsname_simple='CMEMS_reg_5deg.'; %Globcolour;
        case 'GPP'
%             obsname_simple='ORNL_DAAC_reg_5deg.'; % ORNL_DAAC
            obsname_simple='VODCA2GPP_reg_5deg.'; % VODCA2GPP
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

