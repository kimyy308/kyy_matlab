% %  Created 21-Feb-2024 by Yong-Yub Kim

clc; clear all; close all;

%% set maximum thread
LASTN = maxNumCompThreads(10);

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

% grid.regions = [0 360 -60 60];
% grid.regions = [190 240 -5 5];

grid.regions = [0 360 -90 90];
% grid.regions = [0 360 -67 -38];


% cfg.vars={'TS', 'SST', 'PRECT', 'TWS', 'GPP', 'FAREA_BURNED', 'photoC_TOT_zint_100m', ...
%     'COL_FIRE_CLOSS', 'SOILWATER_10CM', 'TLAI', 'SSH', 'PSL', 'AEROD_v'};
% cfg.vars={'SST'};
% cfg.vars={'PSL'};
cfg.vars={'SSH'};
% cfg.vars={'TEMP', 'SALT'};
% cfg.vars={'SSH'};
% cfg.vars={'photoC_TOT_zint_100m'};
% cfg.vars={'NO3'};

flags.AMO=1;
flags.ENSO=1;
flags.PDO=1;
flags.SAM=0;



for vari=1:length(cfg.vars)
    cfg.var=cfg.vars{vari};
    tmp.varname=cfg.var;
    cfg.vlayer=1;


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
            grid.nlon=size(grid.tlong,1);
            grid.nlat=size(grid.tlat,2);
            grid.dz_3d=repmat(grid.dz_res, [grid.nlon, grid.nlat, cfg.vlayer_cnt]);
        case {'atm', 'lnd'}
            grid.lon=ncread(tmp.gridname, 'lon');
            grid.lat=ncread(tmp.gridname, 'lat');
            [grid.tlat, grid.tlong]=meshgrid(grid.lat, grid.lon);
            grid.area=ncread(tmp.gridname, 'AREA');cfg.vlayer_1st=min(cfg.vlayer);
            cfg.vlayer_cnt=max(cfg.vlayer)-cfg.vlayer_1st+1;
            grid.lfrac=ncread(tmp.gridname, 'LANDFRAC');
            grid.area=grid.area.*grid.lfrac;
            grid.nlon=size(grid.tlong,1);
            grid.nlat=size(grid.tlat,2);
    end

    

    [grid.id_w, grid.id_e, grid.id_s, grid.id_n] = Func_0012_findind_Y(0.1, grid.regions, ...
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
    switch cfg.comp
        case {'ocn', 'ice'}
            grid.cut_dz_3d=grid.dz_3d(grid.id_w:grid.id_e, grid.id_s:grid.id_n, cfg.vlayer);
            grid.cut_dz_3d2=reshape(grid.cut_dz_3d, [1 size(grid.cut_dz_3d)]);
    end


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
        elseif (strcmp(cfg.obs_name, 'CMEMS')==1 && strcmp(tmp.varname, 'SSH')==1)
            cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'sla/monthly_reg_pop/', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
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
            [aa,bb,dimids,dd]=netcdf.inqVar(ncid,tmpvarid);
            if length(dimids)==3
                tmp.dd =  netcdf.getVar(ncid,tmpvarid, [grid.id_w-1 grid.id_s-1 0], [grid.cut_nlon grid.cut_nlat 1]);                
            else
                tmp.dd =  netcdf.getVar(ncid,tmpvarid, [grid.id_w-1 grid.id_s-1 cfg.vlayer_1st-1 0], [grid.cut_nlon grid.cut_nlat cfg.vlayer_cnt 1]);
            end
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



% % yearly mean
% if cfg.vlayer_cnt==1 %% 2D or 3D
%     tmp.reshp=reshape(data_obs.([cfg.var]), [grid.cut_nlon, grid.cut_nlat, 12, cfg.len_t_y]);
%     data_obs.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,3));
% else
%     tmp.reshp=reshape(data_obs.([cfg.var]), [grid.cut_nlon, grid.cut_nlat, cfg.vlayer_cnt, 12, cfg.len_t_y]);
%     data_obs.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,4));
% end



% % %% spatial mean
% % if cfg.vlayer_cnt==1 %% 2D
% %     if length(grid.regions)==4
% %     [m_data_obs.([cfg.var]), tmp.err] = ...
% %         Func_0011_get_area_weighted_mean(data_obs.([cfg.var]), grid.cut_tlong, grid.cut_tlat);
% %     [m_data_obs.([cfg.var, '_ym']), tmp.err] = ...
% %         Func_0011_get_area_weighted_mean(data_obs.([cfg.var, '_ym']), grid.cut_tlong, grid.cut_tlat);
% %     elseif length(grid.regions)==2
% %         m_data_obs.([cfg.var]) = data_obs.([cfg.var]);
% %         m_data_obs.([cfg.var, '_ym']) = data_obs.([cfg.var, '_ym']);
% %     end
% % else %% 3D
% %     grid.cut_dz_3d_masked=grid.cut_dz_3d .* ...
% %         data_obs.([cfg.var, '_ym'])(:,:,:,1) ./ data_obs.([cfg.var, '_ym'])(:,:,:,1);
% %     grid.cut_dz_3d_norm=grid.cut_dz_3d_masked./ sum(grid.cut_dz_3d_masked,3,'omitnan');
% %     data_obs_zint.([cfg.var])=squeeze(sum(data_obs.([cfg.var]).*grid.cut_dz_3d_norm, 3, 'omitnan'));
% %     data_obs_zint.([cfg.var, '_ym'])=squeeze(sum(data_obs.([cfg.var, '_ym']).*grid.cut_dz_3d_norm, 3, 'omitnan'));
% % 
% %     if length(grid.regions)==4
% %     [m_data_obs.([cfg.var]), tmp.err] = ...
% %         Func_0011_get_area_weighted_mean(data_obs_zint.([cfg.var]), grid.cut_tlong, grid.cut_tlat);
% %     [m_data_obs.([cfg.var, '_ym']), tmp.err] = ...
% %         Func_0011_get_area_weighted_mean(data_obs_zint.([cfg.var, '_ym']), grid.cut_tlong, grid.cut_tlat);
% %     elseif length(grid.regions)==2
% %         m_data_obs.([cfg.var]) = data_obs_zint.([cfg.var]);
% %         m_data_obs.([cfg.var, '_ym']) = data_obs_zint.([cfg.var, '_ym']);
% %     end
% % end

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

% %% yearly mean
% if cfg.vlayer_cnt==1 %% 2D or 3D
%     tmp.reshp=reshape(data_assm.([cfg.var]), [cfg_assm.len_mem, grid.cut_nlon, grid.cut_nlat, 12, cfg.len_t_y]);
%     data_assm.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,4));
% % data_assm_em.([cfg.var, '_ym'])=squeeze(mean(data_assm.([cfg.var, '_ym']),1));
% else %% 3D
%     tmp.reshp=reshape(data_assm.([cfg.var]), [cfg_assm.len_mem, grid.cut_nlon, grid.cut_nlat, cfg.vlayer_cnt, 12, cfg.len_t_y]);
%     data_assm.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,5));
% end

%% spatial mean
% % if cfg.vlayer_cnt==1 %% 2D
% %     if length(grid.regions)==4
% %         for mi=1:size(data_assm.([cfg.var]),1)
% %             [m_data_assm.([cfg.var])(mi,:), tmp.err] = ...
% %                 Func_0011_get_area_weighted_mean(squeeze(data_assm.([cfg.var])(mi,:,:,:)), grid.cut_tlong, grid.cut_tlat);
% %             [m_data_assm.([cfg.var, '_ym'])(mi,:), tmp.err] = ...
% %                 Func_0011_get_area_weighted_mean(squeeze(data_assm.([cfg.var, '_ym'])(mi,:,:,:)), grid.cut_tlong, grid.cut_tlat);
% %         end
% %     elseif length(grid.regions)==2
% %         m_data_assm.([cfg.var]) = data_assm.([cfg.var]);
% %         m_data_assm.([cfg.var, '_ym']) = data_assm.([cfg.var, '_ym']);
% %     end
% % else %% 3D
% %     grid.cut_dz_3d_masked=grid.cut_dz_3d2 .* ...
% %         data_assm.([cfg.var, '_ym'])(1,:,:,:,1) ./ data_assm.([cfg.var, '_ym'])(1,:,:,:,1) ;
% %     grid.cut_dz_3d_norm=grid.cut_dz_3d_masked./ sum(grid.cut_dz_3d_masked,4,'omitnan');
% %     grid.cut_dz_4d_norm=repmat(grid.cut_dz_3d_norm, [size(data_assm.([cfg.var, '_ym']),1), 1]);
% %     data_assm_zint.([cfg.var])=squeeze(sum(data_assm.([cfg.var]).*grid.cut_dz_4d_norm, 4, 'omitnan'));
% %     data_assm_zint.([cfg.var, '_ym'])=squeeze(sum(data_assm.([cfg.var, '_ym']).*grid.cut_dz_4d_norm, 4, 'omitnan'));
% %     data_assm_zint.([cfg.var])(data_assm_zint.([cfg.var])==0)=NaN;
% %     data_assm_zint.([cfg.var, '_ym'])(data_assm_zint.([cfg.var, '_ym'])==0)=NaN;
% %     if length(grid.regions)==4
% %         for mi=1:size(data_assm.([cfg.var]),1)
% %             [m_data_assm.([cfg.var])(mi,:), tmp.err] = ...
% %                 Func_0011_get_area_weighted_mean(squeeze(data_assm_zint.([cfg.var])(mi,:,:,:)), grid.cut_tlong, grid.cut_tlat);
% %             [m_data_assm.([cfg.var, '_ym'])(mi,:), tmp.err] = ...
% %                 Func_0011_get_area_weighted_mean(squeeze(data_assm_zint.([cfg.var, '_ym'])(mi,:,:,:)), grid.cut_tlong, grid.cut_tlat);
% %         end
% %     elseif length(grid.regions)==2
% %         m_data_assm.([cfg.var]) = data_assm_zint.([cfg.var]);
% %         m_data_assm.([cfg.var, '_ym']) = data_assm_zint.([cfg.var, '_ym']);
% %     end
% % end
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
% %% yearly mean
% if cfg.vlayer_cnt==1 %% 2D
%     tmp.reshp=reshape(data_lens2.([cfg.var]), [cfg_lens2.len_mem, grid.cut_nlon, grid.cut_nlat, 12, cfg.len_t_y]);
%     data_lens2.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,4));
% % data_lens2_em.([cfg.var, '_ym'])=squeeze(mean(data_lens2.([cfg.var, '_ym']),1));
% else %% 3D
%     tmp.reshp=reshape(data_lens2.([cfg.var]), [cfg_lens2.len_mem, grid.cut_nlon, grid.cut_nlat, cfg.vlayer_cnt, 12, cfg.len_t_y]);
%     data_lens2.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,5));
% end


% %% spatial mean
% if cfg.vlayer_cnt==1 %% 2D
%     if length(grid.regions)==4
%         for mi=1:size(data_lens2.([cfg.var]),1)
%             [m_data_lens2.([cfg.var])(mi,:), tmp.err] = ...
%                 Func_0011_get_area_weighted_mean(squeeze(data_lens2.([cfg.var])(mi,:,:,:)), grid.cut_tlong, grid.cut_tlat);
%             [m_data_lens2.([cfg.var, '_ym'])(mi,:), tmp.err] = ...
%                 Func_0011_get_area_weighted_mean(squeeze(data_lens2.([cfg.var, '_ym'])(mi,:,:,:)), grid.cut_tlong, grid.cut_tlat);
%         end
%     elseif length(grid.regions)==2
%         m_data_lens2.([cfg.var]) = data_lens2.([cfg.var]);
%         m_data_lens2.([cfg.var, '_ym']) = data_lens2.([cfg.var, '_ym']);
%     end
% else %% 3D
%     grid.cut_dz_3d_masked=grid.cut_dz_3d2 .* ...
%         data_lens2.([cfg.var, '_ym'])(1,:,:,:,1) ./ data_lens2.([cfg.var, '_ym'])(1,:,:,:,1);
%     grid.cut_dz_3d_norm=grid.cut_dz_3d_masked./ sum(grid.cut_dz_3d_masked,4,'omitnan');
%     grid.cut_dz_4d_norm=repmat(grid.cut_dz_3d_norm, [size(data_lens2.([cfg.var, '_ym']),1), 1]);
%     data_lens2_zint.([cfg.var])=squeeze(sum(data_lens2.([cfg.var]).*grid.cut_dz_4d_norm, 4, 'omitnan'));
%     data_lens2_zint.([cfg.var, '_ym'])=squeeze(sum(data_lens2.([cfg.var, '_ym']).*grid.cut_dz_4d_norm, 4, 'omitnan'));
%     data_lens2_zint.([cfg.var])(data_lens2_zint.([cfg.var])==0)=NaN;
%     data_lens2_zint.([cfg.var, '_ym'])(data_lens2_zint.([cfg.var, '_ym'])==0)=NaN;
%     if length(grid.regions)==4
%         for mi=1:size(data_lens2.([cfg.var]),1)
%             [m_data_lens2.([cfg.var])(mi,:), tmp.err] = ...
%                 Func_0011_get_area_weighted_mean(squeeze(data_lens2_zint.([cfg.var])(mi,:,:,:)), grid.cut_tlong, grid.cut_tlat);
%             [m_data_lens2.([cfg.var, '_ym'])(mi,:), tmp.err] = ...
%                 Func_0011_get_area_weighted_mean(squeeze(data_lens2_zint.([cfg.var, '_ym'])(mi,:,:,:)), grid.cut_tlong, grid.cut_tlat);
%         end
%     elseif length(grid.regions)==2
%         m_data_lens2.([cfg.var]) = data_lens2_zint.([cfg.var]);
%         m_data_lens2.([cfg.var, '_ym']) = data_lens2_zint.([cfg.var, '_ym']);
%     end
% end
toc;



%% read HCST members

for ly=1:5 %1:5
    tmp.ly_str=['ly',num2str(ly)];
    if cfg.vlayer_cnt==1 %% 2D or 3D
        data_hcst.([cfg.var]).(tmp.ly_str)=NaN(cfg_hcst.len_mem, grid.cut_nlon, grid.cut_nlat, cfg.len_t_y*12); % initialization
    else
        data_hcst.([cfg.var]).(tmp.ly_str)=NaN(cfg_hcst.len_mem, grid.cut_nlon, grid.cut_nlat, cfg.vlayer_cnt, cfg.len_t_y*12); % initialization        
    end
    for mi= 1:cfg_hcst.len_mem
        ty=1;
        tmp.member=cfg_hcst.members{mi};
        disp(['hcst, ', tmp.member, ', ', tmp.ly_str]);
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
                    data_hcst.([cfg.var]).(tmp.ly_str)(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd;
                else
                    data_hcst.([cfg.var]).(tmp.ly_str)(mi,:,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd; %depth averaged value                
                end
            else
                tmp.dd=  netcdf.getVar(ncid,tmpvarid, ...
                    [grid.id_w-1 grid.id_s-1 ind_start-1], ...
                    [grid.cut_nlon grid.cut_nlat ind_count]);
                tmp.dd(abs(tmp.dd)>1e30)=NaN;
                data_hcst.([cfg.var]).(tmp.ly_str)(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd;
            end
            netcdf.close(ncid);
    
            ty=ty+ind_count/12;
        end
    end
%     % yearly mean
%     if cfg.vlayer_cnt==1 %% 2D
%         tmp.reshp=reshape(data_hcst.([cfg.var]), [cfg_hcst.len_mem, grid.cut_nlon, grid.cut_nlat, 12, cfg.len_t_y]);
%         data_hcst.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,4));
% %     data_hcst_em.([cfg.var, '_ym'])=squeeze(mean(data_hcst.([cfg.var, '_ym']),1));
%     else
%         tmp.reshp=reshape(data_hcst.([cfg.var]), [cfg_hcst.len_mem, grid.cut_nlon, grid.cut_nlat, cfg.vlayer_cnt, 12, cfg.len_t_y]);
%         data_hcst.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,5));
%     end
    
%     %% spatial mean
%     if cfg.vlayer_cnt==1 %% 2D
%         if length(grid.regions)==4
%             for mi=1:size(data_hcst.([cfg.var]),1)
%                 [m_data_hcst.([cfg.var]).(tmp.ly_str)(mi,:), tmp.err] = ...
%                     Func_0011_get_area_weighted_mean(squeeze(data_hcst.([cfg.var])(mi,:,:,:)), grid.cut_tlong, grid.cut_tlat);
%                 [m_data_hcst.([cfg.var, '_ym']).(tmp.ly_str)(mi,:), tmp.err] = ...
%                     Func_0011_get_area_weighted_mean(squeeze(data_hcst.([cfg.var, '_ym'])(mi,:,:,:)), grid.cut_tlong, grid.cut_tlat);
%             end
%         elseif length(grid.regions)==2
%             m_data_hcst.([cfg.var]).(tmp.ly_str) = data_hcst.([cfg.var]);
%             m_data_hcst.([cfg.var, '_ym']).(tmp.ly_str) = data_hcst.([cfg.var, '_ym']);
%         end
%     else %% 3D
%         grid.cut_dz_3d_masked=grid.cut_dz_3d2 .* ...
%             data_hcst.([cfg.var, '_ym'])(1,:,:,:,1) ./ data_hcst.([cfg.var, '_ym'])(1,:,:,:,1);
%         grid.cut_dz_3d_norm=grid.cut_dz_3d_masked./ sum(grid.cut_dz_3d_masked,4,'omitnan');
%         grid.cut_dz_4d_norm=repmat(grid.cut_dz_3d_norm, [size(data_hcst.([cfg.var, '_ym']),1), 1]);
%         data_hcst_zint.([cfg.var])=squeeze(sum(data_hcst.([cfg.var]).*grid.cut_dz_4d_norm, 4, 'omitnan'));
%         data_hcst_zint.([cfg.var, '_ym'])=squeeze(sum(data_hcst.([cfg.var, '_ym']).*grid.cut_dz_4d_norm, 3, 'omitnan'));
%         data_hcst_zint.([cfg.var])(data_hcst_zint.([cfg.var])==0)=NaN;
%         data_hcst_zint.([cfg.var, '_ym'])(data_hcst_zint.([cfg.var, '_ym'])==0)=NaN;
%         if length(grid.regions)==4
%             for mi=1:size(data_hcst.([cfg.var]),1)
%                 [m_data_hcst.([cfg.var]).(tmp.ly_str)(mi,:), tmp.err] = ...
%                     Func_0011_get_area_weighted_mean(squeeze(data_hcst_zint.([cfg.var])(mi,:,:,:)), grid.cut_tlong, grid.cut_tlat);
%                 [m_data_hcst.([cfg.var, '_ym']).(tmp.ly_str)(mi,:), tmp.err] = ...
%                     Func_0011_get_area_weighted_mean(squeeze(data_hcst_zint.([cfg.var, '_ym'])(mi,:,:,:)), grid.cut_tlong, grid.cut_tlat);
%             end
%         elseif length(grid.regions)==2
%             m_data_hcst.([cfg.var]).(tmp.ly_str) = data_hcst_zint.([cfg.var]);
%             m_data_hcst.([cfg.var, '_ym']).(tmp.ly_str) = data_hcst_zint.([cfg.var, '_ym']);
%         end
%     end
    toc;
    
end



if flags.ENSO==1


    %% ENSO
    %% ENSO - grid information
    data_ENSO.regions = [190 240 -5 5];
    [data_ENSO.id_w, data_ENSO.id_e, data_ENSO.id_s, data_ENSO.id_n] = ...
        Func_0012_findind_Y(0.1, data_ENSO.regions, ...
                grid.tlong, grid.tlat, 'CESM2'); % find valid lon, lat index near station
    data_ENSO.cut_tlong=grid.tlong(data_ENSO.id_w:data_ENSO.id_e, data_ENSO.id_s:data_ENSO.id_n);
    data_ENSO.cut_tlat=grid.tlat(data_ENSO.id_w:data_ENSO.id_e, data_ENSO.id_s:data_ENSO.id_n);
    data_ENSO.cut_nlon=size(data_ENSO.cut_tlong,1);
    data_ENSO.cut_nlat=size(data_ENSO.cut_tlat,2);
    
    
    %% ENSO - time set
    for ti=1:length(cfg.iyears)
        for mi=1:12
            data_ENSO.time((ti-1)*12+mi)=cfg.iyears(ti)+mi*1/12-1/24;
        end
    end
    
    %% ENSO - OBS
    tmp.data=data_obs.(cfg.var)(data_ENSO.id_w:data_ENSO.id_e, data_ENSO.id_s:data_ENSO.id_n, :);
    tmp.data(tmp.data<-900)=NaN;
    [data_ENSO.obs, tmp.err] = ...
        Func_0011_get_area_weighted_mean(squeeze(tmp.data), data_ENSO.cut_tlong, data_ENSO.cut_tlat);
    
    %% ENSO - ASSM
    data_ENSO.assm_members=cfg_assm.members;
    for mi=1:size(data_assm.(cfg.var),1)
        tmp.data=squeeze(data_assm.(cfg.var)(mi,data_ENSO.id_w:data_ENSO.id_e, data_ENSO.id_s:data_ENSO.id_n,:));
        [data_ENSO.assm(mi,:), tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(tmp.data), data_ENSO.cut_tlong, data_ENSO.cut_tlat);
    end
    
    % hold on
    % for mi=1:size(data_assm.(cfg.var),1)
    %     plot(data_ENSO.time,data_ENSO.lens2(mi,:));
    % end
    % hold off
    
    %% ENSO - LENS2
    data_ENSO.lens2_members=cfg_lens2.members;
    for mi=1:size(data_lens2.(cfg.var),1)
        tmp.data=squeeze(data_lens2.(cfg.var)(mi,data_ENSO.id_w:data_ENSO.id_e, data_ENSO.id_s:data_ENSO.id_n,:));
        [data_ENSO.lens2(mi,:), tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(tmp.data), data_ENSO.cut_tlong, data_ENSO.cut_tlat);
    end
    
    %% ENSO - HCST
    data_ENSO.hcst_members=cfg_hcst.members;
    for ly=1:5 %1:5
        tmp.ly_str=['ly',num2str(ly)];
        for mi=1:size(data_hcst.(cfg.var).(tmp.ly_str),1)
            tmp.data=squeeze(data_hcst.(cfg.var).(tmp.ly_str)(mi,data_ENSO.id_w:data_ENSO.id_e, data_ENSO.id_s:data_ENSO.id_n,:));
            [data_ENSO.hcst(ly,mi,:), tmp.err] = ...
                Func_0011_get_area_weighted_mean(squeeze(tmp.data), data_ENSO.cut_tlong, data_ENSO.cut_tlat);
        end
    end
    
    % hold on
    % for mi=1:size(data_hcst.(cfg.var).ly1,1)
    %     plot(data_ENSO.time,squeeze(data_ENSO.hcst(4,mi,:)));
    % end
    % hold off
    
    %% ENSO - save matfile 
    mkdir([dirs.saveroot, '/clim_indices']);
    matfilename=[dirs.saveroot, '/clim_indices/', 'clim_indices_', cfg.var, '_all_','ENSO_', ...
                'obs_', cfg.obs_name, '.mat'];
    save(matfilename, 'cfg', 'cfg_assm', 'cfg_hcst', 'cfg_lens2', 'data_ENSO');
end


if flags.AMO==1

    %% AMO index
    %% AMO - grid information
    data_AMO.ATL_regions = [280 360 0 60];
    [data_AMO.ATL_id_w, data_AMO.ATL_id_e, data_AMO.ATL_id_s, data_AMO.ATL_id_n] = ...
        Func_0012_findind_Y(0.1, data_AMO.ATL_regions, ...
                grid.tlong, grid.tlat, 'CESM2'); % find valid lon, lat index near station
    data_AMO.ATL_cut_tlong=grid.tlong(data_AMO.ATL_id_w:data_AMO.ATL_id_e, data_AMO.ATL_id_s:data_AMO.ATL_id_n);
    data_AMO.ATL_cut_tlat=grid.tlat(data_AMO.ATL_id_w:data_AMO.ATL_id_e, data_AMO.ATL_id_s:data_AMO.ATL_id_n);
    data_AMO.ATL_cut_nlon=size(data_AMO.ATL_cut_tlong,1);
    data_AMO.ATL_cut_nlat=size(data_AMO.ATL_cut_tlat,2);
    
    data_AMO.GLO_regions = [0 360 -60 60];
    [data_AMO.GLO_id_w, data_AMO.GLO_id_e, data_AMO.GLO_id_s, data_AMO.GLO_id_n] = ...
        Func_0012_findind_Y(0.1, data_AMO.GLO_regions, ...
                grid.tlong, grid.tlat, 'CESM2'); % find valid lon, lat index near station
    data_AMO.GLO_cut_tlong=grid.tlong(data_AMO.GLO_id_w:data_AMO.GLO_id_e, data_AMO.GLO_id_s:data_AMO.GLO_id_n);
    data_AMO.GLO_cut_tlat=grid.tlat(data_AMO.GLO_id_w:data_AMO.GLO_id_e, data_AMO.GLO_id_s:data_AMO.GLO_id_n);
    data_AMO.GLO_cut_nlon=size(data_AMO.GLO_cut_tlong,1);
    data_AMO.GLO_cut_nlat=size(data_AMO.GLO_cut_tlat,2);
    
    
    %% AMO - time set
    for ti=1:length(cfg.iyears)
        for mi=1:12
            data_AMO.time((ti-1)*12+mi)=cfg.iyears(ti)+mi*1/12-1/24;
        end
    end
    
    %% AMO - OBS 
    tmp.data=data_obs.(cfg.var)(data_AMO.GLO_id_w:data_AMO.GLO_id_e, data_AMO.GLO_id_s:data_AMO.GLO_id_n, :);
    tmp.data(tmp.data<-900)=NaN;
    [data_AMO.GLO_obs, tmp.err] = ...
        Func_0011_get_area_weighted_mean(squeeze(tmp.data), data_AMO.GLO_cut_tlong, data_AMO.GLO_cut_tlat);
    
    tmp.data_dseason= ...
        reshape(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]) - ...
        mean(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]),4, 'omitnan'), [size(tmp.data)]);
    [data_AMO.GLO_obs_dseason, tmp.err] = ...
        Func_0011_get_area_weighted_mean(squeeze(tmp.data_dseason), data_AMO.GLO_cut_tlong, data_AMO.GLO_cut_tlat);
    
    tmp.data=data_obs.(cfg.var)(data_AMO.ATL_id_w:data_AMO.ATL_id_e, data_AMO.ATL_id_s:data_AMO.ATL_id_n, :);
    tmp.data(tmp.data<-900)=NaN;
    [data_AMO.ATL_obs, tmp.err] = ...
        Func_0011_get_area_weighted_mean(squeeze(tmp.data), data_AMO.ATL_cut_tlong, data_AMO.ATL_cut_tlat);
    
    tmp.data_dseason= ...
        reshape(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]) - ...
        mean(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]),4, 'omitnan'), [size(tmp.data)]);
    [data_AMO.ATL_obs_dseason, tmp.err] = ...
        Func_0011_get_area_weighted_mean(squeeze(tmp.data_dseason), data_AMO.ATL_cut_tlong, data_AMO.ATL_cut_tlat);
    
    
    % %% normal lowpass
    % data_AMO.obs=data_AMO.ATL_obs-data_AMO.GLO_obs;
    % data_AMO.lp_period_month=120;
    % data_AMO.obs_lp=lowpass(data_AMO.obs, 1/data_AMO.lp_period_month,  'steepness', 0.90);
    % % data_AMO.obs_lp(1:data_AMO.lp_period_month/2)=NaN;
    % % data_AMO.obs_lp(end-(data_AMO.lp_period_month/2-1):end)=NaN;
    % plot(data_AMO.time,data_AMO.obs_lp)
    
    % %% normal lowpass_dseason
    data_AMO.obs_dseason=data_AMO.ATL_obs_dseason-data_AMO.GLO_obs_dseason;
    data_AMO.lp_period_month=120;
    data_AMO.obs_dseason_lp=lowpass(data_AMO.obs_dseason, 1/data_AMO.lp_period_month,  'steepness', 0.9999);
    % data_AMO.obs_lp(1:data_AMO.lp_period_month/2)=NaN;
    % data_AMO.obs_lp(end-(data_AMO.lp_period_month/2-1):end)=NaN;
    plot(data_AMO.time,data_AMO.obs_dseason_lp)
    hold on
    plot(data_AMO.time, data_AMO.obs_dseason)
    hold off
    
    % %% movmean
    % data_AMO.obs_movm=movmean(data_AMO.obs, 120, 'Endpoints', 'discard');
    % plot(data_AMO.time(60:end-60),data_AMO.obs_movm)
    
    % %% yearly - lowpass
    % data_AMO.obs_y=mean(reshape(data_AMO.obs, [12, length(data_AMO.obs)/12]),1);
    % data_AMO.lp_period_y=10;
    % data_AMO.obs_y_lp=lowpass(data_AMO.obs_y, 1/data_AMO.lp_period_y,  'steepness', 0.80);
    % plot(cfg.iyears,data_AMO.obs_y_lp)
    % hold on
    % plot(cfg.iyears,data_AMO.obs_y)
    
    %% AMO - ASSM
    data_AMO.assm_members=cfg_assm.members;
    for mi=1:size(data_assm.(cfg.var),1)
        tmp.data=squeeze(data_assm.(cfg.var)(mi,data_AMO.GLO_id_w:data_AMO.GLO_id_e, data_AMO.GLO_id_s:data_AMO.GLO_id_n,:));
        tmp.data(tmp.data==0)=NaN;
        [data_AMO.GLO_assm(mi,:), tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(tmp.data), data_AMO.GLO_cut_tlong, data_AMO.GLO_cut_tlat);
        
        tmp.data_dseason= ...
        reshape(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]) - ...
        mean(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]),4), [size(tmp.data)]);
        [data_AMO.GLO_assm_dseason(mi,:), tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(tmp.data_dseason), data_AMO.GLO_cut_tlong, data_AMO.GLO_cut_tlat);
    
        data_AMO.assm_members{mi}
    end
    
    for mi=1:size(data_assm.(cfg.var),1)
        tmp.data=squeeze(data_assm.(cfg.var)(mi,data_AMO.ATL_id_w:data_AMO.ATL_id_e, data_AMO.ATL_id_s:data_AMO.ATL_id_n,:));
        tmp.data(tmp.data==0)=NaN;
        [data_AMO.ATL_assm(mi,:), tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(tmp.data), data_AMO.ATL_cut_tlong, data_AMO.ATL_cut_tlat);
    
        tmp.data_dseason= ...
        reshape(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]) - ...
        mean(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]),4), [size(tmp.data)]);
        [data_AMO.ATL_assm_dseason(mi,:), tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(tmp.data_dseason), data_AMO.ATL_cut_tlong, data_AMO.ATL_cut_tlat);
    
        data_AMO.assm_members{mi}
    end
    data_AMO.assm=data_AMO.ATL_assm-data_AMO.GLO_assm;
    data_AMO.assm_dseason=data_AMO.ATL_assm_dseason-data_AMO.GLO_assm_dseason;
    
    
    for mi=1:size(data_assm.(cfg.var),1)
        data_AMO.assm_dseason_lp(mi,:)=lowpass(data_AMO.assm_dseason(mi,:), 1/data_AMO.lp_period_month,  'steepness', 0.9999);
        data_AMO.assm_y=mean(reshape(data_AMO.assm(mi,:), [12, size(data_AMO.assm,2)/12]),1);
    %     data_AMO.assm_y_lp(mi,:)=lowpass(data_AMO.assm_y, 1/data_AMO.lp_period_y, 'steepness', 0.80);
    end
    
    % hold on
    % for mi=1:size(data_assm.(cfg.var),1)
    %     plot(data_AMO.time,data_AMO.assm_dseason_lp(mi,:), 'r');
    % %     plot(cfg.iyears,data_AMO.assm_y_lp(mi,:), 'r');
    % end
    % hold off
    
    %% AMO - LENS2
    data_AMO.lens2_members=cfg_lens2.members;
    for mi=1:size(data_lens2.(cfg.var),1)
        tmp.data=squeeze(data_lens2.(cfg.var)(mi,data_AMO.GLO_id_w:data_AMO.GLO_id_e, data_AMO.GLO_id_s:data_AMO.GLO_id_n,:));
        tmp.data(tmp.data==0)=NaN;
        [data_AMO.GLO_lens2(mi,:), tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(tmp.data), data_AMO.GLO_cut_tlong, data_AMO.GLO_cut_tlat);
    
        tmp.data_dseason= ...
        reshape(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]) - ...
        mean(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]),4), [size(tmp.data)]);
        [data_AMO.GLO_lens2_dseason(mi,:), tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(tmp.data_dseason), data_AMO.GLO_cut_tlong, data_AMO.GLO_cut_tlat);
    
        data_AMO.lens2_members{mi}
    end
    
    for mi=1:size(data_lens2.(cfg.var),1)
        tmp.data=squeeze(data_lens2.(cfg.var)(mi,data_AMO.ATL_id_w:data_AMO.ATL_id_e, data_AMO.ATL_id_s:data_AMO.ATL_id_n,:));
        tmp.data(tmp.data==0)=NaN;
        [data_AMO.ATL_lens2(mi,:), tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(tmp.data), data_AMO.ATL_cut_tlong, data_AMO.ATL_cut_tlat);
    
        tmp.data_dseason= ...
        reshape(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]) - ...
        mean(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]),4), [size(tmp.data)]);
        [data_AMO.ATL_lens2_dseason(mi,:), tmp.err] = ...
            Func_0011_get_area_weighted_mean(squeeze(tmp.data_dseason), data_AMO.ATL_cut_tlong, data_AMO.ATL_cut_tlat);
    
        data_AMO.lens2_members{mi}
    end
    data_AMO.lens2=data_AMO.ATL_lens2-data_AMO.GLO_lens2;
    data_AMO.lens2_dseason=data_AMO.ATL_lens2_dseason-data_AMO.GLO_lens2_dseason;
    
    
    for mi=1:size(data_lens2.(cfg.var),1)
        data_AMO.lens2_dseason_lp(mi,:)=lowpass(data_AMO.lens2_dseason(mi,:), 1/data_AMO.lp_period_month,  'steepness', 0.9999);
    
        data_AMO.lens2_y=mean(reshape(data_AMO.lens2(mi,:), [12, size(data_AMO.lens2,2)/12]),1);
    %     data_AMO.lens2_y_lp(mi,:)=lowpass(data_AMO.lens2_y, 1/data_AMO.lp_period_y, 'steepness', 0.80);
    end
    
    % hold on
    % for mi=1:size(data_lens2.(cfg.var),1)
    %     plot(data_AMO.time,data_AMO.lens2_dseason_lp(mi,:), 'g');    
    % %     plot(cfg.iyears,data_AMO.lens2_y_lp(mi,:), 'g');
    % end
    % hold off
    
    
    %% AMO - HCST
    data_AMO.hcst_members=cfg_hcst.members;
    for ly=1:5 %1:5
        tmp.ly_str=['ly',num2str(ly)];
        for mi=1:size(data_hcst.(cfg.var).(tmp.ly_str),1)
            tmp.data=squeeze(data_hcst.(cfg.var).(tmp.ly_str)(mi,data_AMO.GLO_id_w:data_AMO.GLO_id_e, data_AMO.GLO_id_s:data_AMO.GLO_id_n,:));
            tmp.data(tmp.data==0)=NaN;
            [data_AMO.GLO_hcst(ly,mi,:), tmp.err] = ...
                Func_0011_get_area_weighted_mean(squeeze(tmp.data), data_AMO.GLO_cut_tlong, data_AMO.GLO_cut_tlat);
    
            tmp.data_dseason= ...
            reshape(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]) - ...
            mean(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]),4), [size(tmp.data)]);
            [data_AMO.GLO_hcst_dseason(ly,mi,:), tmp.err] = ...
                Func_0011_get_area_weighted_mean(squeeze(tmp.data_dseason), data_AMO.GLO_cut_tlong, data_AMO.GLO_cut_tlat);
    
            disp(['ly',num2str(ly), data_AMO.hcst_members{mi}]);
        end
    end
    
    for ly=1:5 %1:5
        tmp.ly_str=['ly',num2str(ly)];
        for mi=1:size(data_hcst.(cfg.var).(tmp.ly_str),1)
            tmp.data=squeeze(data_hcst.(cfg.var).(tmp.ly_str)(mi,data_AMO.ATL_id_w:data_AMO.ATL_id_e, data_AMO.ATL_id_s:data_AMO.ATL_id_n,:));
            tmp.data(tmp.data==0)=NaN;
            [data_AMO.ATL_hcst(ly,mi,:), tmp.err] = ...
                Func_0011_get_area_weighted_mean(squeeze(tmp.data), data_AMO.ATL_cut_tlong, data_AMO.ATL_cut_tlat);
    
            tmp.data_dseason= ...
            reshape(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]) - ...
            mean(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]),4), [size(tmp.data)]);
            [data_AMO.ATL_hcst_dseason(ly,mi,:), tmp.err] = ...
                Func_0011_get_area_weighted_mean(squeeze(tmp.data_dseason), data_AMO.ATL_cut_tlong, data_AMO.ATL_cut_tlat);
    
            ly*100+mi
        end
    end
    
    data_AMO.hcst=data_AMO.ATL_hcst-data_AMO.GLO_hcst;
    data_AMO.hcst_dseason=data_AMO.ATL_hcst_dseason-data_AMO.GLO_hcst_dseason;
    
    for ly=1:5 %1:5
        for mi=1:size(data_hcst.(cfg.var).ly1,1)
            data_AMO.hcst_dseason_lp(ly,mi,:)=lowpass(squeeze(data_AMO.hcst_dseason(ly,mi,:)), 1/data_AMO.lp_period_month,  'steepness', 0.9999);
            data_AMO.hcst_y=mean(reshape(data_AMO.hcst(ly,mi,:), [12, size(data_AMO.hcst,3)/12]),1);
    %         data_AMO.hcst_y_lp(ly,mi,:)=lowpass(data_AMO.hcst_y, 1/data_AMO.lp_period_y, 'steepness', 0.80);
        end
    end
    
    % hold on
    % ly=5;
    % for mi=1:size(data_hcst.(cfg.var).ly1,1)
    %     plot(data_AMO.time+ly-1,squeeze(data_AMO.hcst_dseason_lp(ly,mi,:)), 'b');
    % %     plot(cfg.iyears,squeeze(data_AMO.hcst_y_lp(ly,mi,:)), 'k');
    % end
    % hold off
    
    %% AMO - save matfile 
    mkdir([dirs.saveroot, '/clim_indices']);
    matfilename=[dirs.saveroot, '/clim_indices/', 'clim_indices_', cfg.var, '_all_','AMO_', ...
                'obs_', cfg.obs_name, '.mat'];
    % nccreate('mync.nc', 'AMO_hcst', 'Dimensions', {'ly', 5, 'mem', 20, 't', 732});
    % ncwrite('mync.nc', 'AMO_hcst', data_AMO.hcst_dseason_lp);
    save(matfilename, 'cfg', 'cfg_assm', 'cfg_hcst', 'cfg_lens2', 'data_AMO');
end



if flags.PDO==1
    
    AMO_matfilename=[dirs.saveroot, '/clim_indices/', 'clim_indices_', cfg.var, '_all_','AMO_', ...
                'obs_', cfg.obs_name, '.mat'];
  
    load(AMO_matfilename, 'cfg', 'cfg_assm', 'cfg_hcst', 'cfg_lens2', 'data_AMO');

    %% PDO index
    %% PDO - grid information
    data_PDO.regions = [120 260 20 60];
    [data_PDO.id_w, data_PDO.id_e, data_PDO.id_s, data_PDO.id_n] = ...
        Func_0012_findind_Y(0.1, data_PDO.regions, ...
                grid.tlong, grid.tlat, 'CESM2'); % find valid lon, lat index near station
    data_PDO.cut_tlong=grid.tlong(data_PDO.id_w:data_PDO.id_e, data_PDO.id_s:data_PDO.id_n);
    data_PDO.cut_tlat=grid.tlat(data_PDO.id_w:data_PDO.id_e, data_PDO.id_s:data_PDO.id_n);
    data_PDO.cut_nlon=size(data_PDO.cut_tlong,1);
    data_PDO.cut_nlat=size(data_PDO.cut_tlat,2);
    
    % [lv, pc, var_exp] = Func_0024_EOF_3d(data,X);
    
    %% PDO - time set
    for ti=1:length(cfg.iyears)
        for mi=1:12
            data_PDO.time((ti-1)*12+mi)=cfg.iyears(ti)+mi*1/12-1/24;
        end
    end
    
    %% PDO & NPGO - OBS 
    tmp.data=data_obs.(cfg.var)(data_PDO.id_w:data_PDO.id_e, data_PDO.id_s:data_PDO.id_n, :);
    tmp.data(tmp.data<-900)=NaN;
    
    tmp.data_dseason= ...
        reshape(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]) - ...
        mean(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]),4, 'omitnan'), [size(tmp.data)]);
    
    tmp.GLO_obs_3d(1,1,:)=data_AMO.GLO_obs_dseason;
    tmp.data_dseason=tmp.data_dseason-tmp.GLO_obs_3d;
    tmp.data_dseason=tmp.data_dseason-mean(tmp.data_dseason,3,'omitnan'); % make anomaly before PDO calculation
    
    maxnind=find(data_AMO.GLO_obs_dseason==0,1,'last');

    tmp.num_modes=3;
    [data_PDO.lv_obs, data_PDO.pct_obs, data_PDO.var_exp_obs] = Func_0024_EOF_3d(tmp.data_dseason(:,:,maxnind+1:end),tmp.num_modes, data_PDO.cut_tlat);
    % plot(data_PDO.time,data_PDO.pct_obs(:,1))
    
    % [data_PDO.GLO_obs, tmp.err] = ...
    %     Func_0011_get_area_weighted_mean(squeeze(tmp.data), data_PDO.GLO_cut_tlong, data_PDO.GLO_cut_tlat);
    
    
    %% normal lowpass_dseason
    data_PDO.lp_period_month=120;
    for modei=1:length(tmp.num_modes)
        data_PDO.pct_obs_lp(:,modei)=lowpass(data_PDO.pct_obs(:,modei), 1/data_PDO.lp_period_month,  'steepness', 0.9999);
    end
    
    
    %% PDO & NPGO - ASSM
    data_PDO.assm_members=cfg_assm.members;
    
    for mi=1:size(data_assm.(cfg.var),1)
        tmp.data=squeeze(data_assm.(cfg.var)(mi,data_PDO.id_w:data_PDO.id_e, data_PDO.id_s:data_PDO.id_n, :));
        tmp.data(tmp.data==0)=NaN;
        tmp.data_mean=mean(tmp.data,3);
        tmp.data(isnan(repmat(tmp.data_mean, [1 1 size(tmp.data,3)])))=NaN;
    
        tmp.data_dseason= ...
            reshape(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]) - ...
            mean(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]),4), [size(tmp.data)]);
        tmp.GLO_assm_3d(1,1,:)=data_AMO.GLO_assm_dseason(mi,:);
        tmp.data_dseason=tmp.data_dseason-tmp.GLO_assm_3d; %% remove global increase
        tmp.data_dseason=tmp.data_dseason-mean(tmp.data_dseason,3); % make anomaly before PDO calculation
    
        [data_PDO.lv_assm(mi,:,:,:), data_PDO.pct_assm(mi,:,:), data_PDO.var_exp_assm(mi,:)] = Func_0024_EOF_3d(tmp.data_dseason,tmp.num_modes, data_PDO.cut_tlat);
    
        data_PDO.assm_members{mi};
    end
    
    for mi=1:size(data_assm.(cfg.var),1)
        for modei=1:length(tmp.num_modes)
            data_PDO.pct_assm_lp(mi,:,modei)=lowpass(data_PDO.pct_assm(mi,:,modei), 1/data_AMO.lp_period_month,  'steepness', 0.9999);
        end
    end
    
    %% PDO & NPGO - LENS2
    data_PDO.lens2_members=cfg_lens2.members;
    
    for mi=1:size(data_lens2.(cfg.var),1)
        tmp.data=squeeze(data_lens2.(cfg.var)(mi,data_PDO.id_w:data_PDO.id_e, data_PDO.id_s:data_PDO.id_n, :));
        tmp.data(tmp.data==0)=NaN;
        tmp.data_mean=mean(tmp.data,3);
        tmp.data(isnan(repmat(tmp.data_mean, [1 1 size(tmp.data,3)])))=NaN;
    
        tmp.data_dseason= ...
            reshape(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]) - ...
            mean(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]),4), [size(tmp.data)]);
        tmp.GLO_lens2_3d(1,1,:)=data_AMO.GLO_lens2_dseason(mi,:);
        tmp.data_dseason=tmp.data_dseason-tmp.GLO_lens2_3d; %% remove global increase
        tmp.data_dseason=tmp.data_dseason-mean(tmp.data_dseason,3); % make anomaly before PDO calculation
    
        [data_PDO.lv_lens2(mi,:,:,:), data_PDO.pct_lens2(mi,:,:), data_PDO.var_exp_lens2(mi,:)] = ...
            Func_0024_EOF_3d(tmp.data_dseason,tmp.num_modes, data_PDO.cut_tlat);
    
        data_PDO.lens2_members{mi};
    end
    
    for mi=1:size(data_lens2.(cfg.var),1)
        for modei=1:length(tmp.num_modes)
            data_PDO.pct_lens2_lp(mi,:,modei)=lowpass(data_PDO.pct_lens2(mi,:,modei), 1/data_AMO.lp_period_month,  'steepness', 0.9999);
        end
    end
    
    
    %% PDO & NPGO - HCST
    data_PDO.hcst_members=cfg_hcst.members;
    
    for ly=1:5 %1:5
        tmp.ly_str=['ly',num2str(ly)];
        for mi=1:size(data_hcst.(cfg.var).(tmp.ly_str),1)
            tmp.data=squeeze(data_hcst.(cfg.var).(tmp.ly_str)(mi,data_PDO.id_w:data_PDO.id_e, data_PDO.id_s:data_PDO.id_n, :));
            tmp.data(tmp.data==0)=NaN;
            tmp.data_mean=mean(tmp.data,3);
            tmp.data(isnan(repmat(tmp.data_mean, [1 1 size(tmp.data,3)])))=NaN;
    
            tmp.data_dseason= ...
            reshape(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]) - ...
            mean(reshape(tmp.data, [size(tmp.data,1), size(tmp.data,2), 12, size(tmp.data,3)/12]),4), [size(tmp.data)]);
    
            tmp.GLO_hcst_3d(1,1,:)=data_AMO.GLO_hcst_dseason(ly,mi,:);
            tmp.data_dseason=tmp.data_dseason-tmp.GLO_hcst_3d; %% remove global increase
            tmp.data_dseason=tmp.data_dseason-mean(tmp.data_dseason,3); % make anomaly before PDO calculation
            
            [data_PDO.lv_hcst(ly,mi,:,:,:), data_PDO.pct_hcst(ly,mi,:,:), data_PDO.var_exp_hcst(ly,mi,:)] = ...
                Func_0024_EOF_3d(tmp.data_dseason,tmp.num_modes, data_PDO.cut_tlat);
    
    %         [data_AMO.ATL_hcst_dseason(ly,mi,:), tmp.err] = ...
    %             Func_0011_get_area_weighted_mean(squeeze(tmp.data_dseason), data_AMO.ATL_cut_tlong, data_AMO.ATL_cut_tlat);
    
            ly*100+mi
        end
    end
    
    for ly=1:5 %1:5
        tmp.ly_str=['ly',num2str(ly)];
        for mi=1:size(data_hcst.(cfg.var).(tmp.ly_str),1)
            for modei=1:length(tmp.num_modes)
                data_PDO.pct_hcst_lp(ly,mi,:,modei)=lowpass(squeeze(data_PDO.pct_hcst(ly,mi,:,modei)), 1/data_AMO.lp_period_month,  'steepness', 0.9999);
            end
        end
    end
    
    %% PDO - save matfile 
    mkdir([dirs.saveroot, '/clim_indices']);
    matfilename=[dirs.saveroot, '/clim_indices/', 'clim_indices_', cfg.var, '_all_','PDO_', ...
                'obs_', cfg.obs_name, '.mat'];
    save(matfilename, 'cfg', 'cfg_assm', 'cfg_hcst', 'cfg_lens2', 'data_PDO');

%     hold on
%     for pi=1:20
%         if mean(data_PDO.pct_assm(pi,end,1))<0
%             plot(data_PDO.time,data_PDO.pct_assm(pi,:,1), 'r')
%         else
%             plot(data_PDO.time,-data_PDO.pct_assm(pi,:,1), 'r')            
%         end
%     end
%     plot(data_PDO.time,data_PDO.pct_obs(:,1), 'color', 'k', 'linewidth', 2)
% 
% 
%     hold on
%     for pi=1:20
%             plot(data_AMO.time,data_AMO.assm_dseason_lp(pi,:), 'r')
%     end
%     plot(data_PDO.time,data_AMO.obs_dseason_lp(:), 'color', 'k', 'linewidth', 2)

end

if flags.SAM==1
    grid.regions_40= [0 360 -40 -40];
    [grid.id_w_40, grid.id_e_40, grid.id_s_40, grid.id_n_40] = Func_0012_findind_Y(0.1, grid.regions_40, ...
                grid.cut_tlong, grid.cut_tlat, 'CESM2'); % find valid lon, lat index near station
    tmp.idlat=round((grid.id_s_40+grid.id_n_40)/2);
    grid.id_s_40=tmp.idlat;
    grid.id_n_40=tmp.idlat;

    grid.regions_65= [0 360 -65 -65];
    [grid.id_w_65, grid.id_e_65, grid.id_s_65, grid.id_n_65] = Func_0012_findind_Y(0.5, grid.regions_65, ...
                grid.cut_tlong, grid.cut_tlat, 'CESM2'); % find valid lon, lat index near station
    tmp.idlat=round((grid.id_s_65+grid.id_n_65)/2);
    grid.id_s_65=tmp.idlat;
    grid.id_n_65=tmp.idlat;

    grid.cut_tlong_40=grid.cut_tlong(grid.id_w_40:grid.id_e_40, grid.id_s_40:grid.id_n_40);
    grid.cut_tlat_40=grid.cut_tlat(grid.id_w_40:grid.id_e_40, grid.id_s_40:grid.id_n_40);
    grid.cut_nlon_40=size(grid.cut_tlong,1);

    grid.cut_tlong_65=grid.cut_tlong(grid.id_w_65:grid.id_e_65, grid.id_s_65:grid.id_n_65);
    grid.cut_tlat_65=grid.cut_tlat(grid.id_w_65:grid.id_e_65, grid.id_s_65:grid.id_n_65);
    grid.cut_nlon_65=size(grid.cut_tlong,1);

    %% SAM - OBS
    tmp.data_40=data_obs.(cfg.var)(grid.id_w_40:grid.id_e_40, grid.id_s_40:grid.id_n_40, :);
    tmp.data_40(tmp.data_40<-900)=NaN;
    tmp.data_40(tmp.data_40==0)=NaN;
    data_SAM.obs_40=squeeze(mean(tmp.data_40,1));
    data_SAM.obs_40=data_SAM.obs_40-mean(data_SAM.obs_40);
    data_SAM.obs_40=data_SAM.obs_40./std(data_SAM.obs_40);

    tmp.data_65=data_obs.(cfg.var)(grid.id_w_65:grid.id_e_65, grid.id_s_65:grid.id_n_65, :);
    tmp.data_65(tmp.data_65<-900)=NaN;
    tmp.data_65(tmp.data_65==0)=NaN;
    data_SAM.obs_65=squeeze(mean(tmp.data_65,1));
    data_SAM.obs_65=data_SAM.obs_65-mean(data_SAM.obs_65);
    data_SAM.obs_65=data_SAM.obs_65./std(data_SAM.obs_65);

    data_SAM.obs=data_SAM.obs_40-data_SAM.obs_65;
    
     %% SAM - time set
    for ti=1:length(cfg.iyears)
        for mi=1:12
            data_SAM.time((ti-1)*12+mi)=cfg.iyears(ti)+mi*1/12-1/24;
        end
    end
    data_SAM.lp_period_month=120;
    data_SAM.obs_lp=lowpass(data_SAM.obs, 1/data_SAM.lp_period_month,  'steepness', 0.989);

    plot(data_SAM.time,data_SAM.obs_lp)


    %% SAM-ASSM

    data_SAM.assm_members=cfg_assm.members;
    
    for mi=1:size(data_assm.(cfg.var),1)
        tmp.data_40=squeeze(data_assm.(cfg.var)(mi,grid.id_w_40:grid.id_e_40, grid.id_s_40:grid.id_n_40, :));
        tmp.data_40(tmp.data_40<-900)=NaN;
        tmp.data_40(tmp.data_40==0)=NaN;
        data_SAM.assm_40(mi,:)=squeeze(mean(tmp.data_40,1));
        data_SAM.assm_40(mi,:)=data_SAM.assm_40(mi,:)-mean(data_SAM.assm_40(mi,:));
        data_SAM.assm_40(mi,:)=data_SAM.assm_40(mi,:)./std(data_SAM.assm_40(mi,:));
    
        tmp.data_65=squeeze(data_assm.(cfg.var)(mi,grid.id_w_65:grid.id_e_65, grid.id_s_65:grid.id_n_65, :));
        tmp.data_65(tmp.data_65<-900)=NaN;
        tmp.data_65(tmp.data_65==0)=NaN;
        data_SAM.assm_65(mi,:)=squeeze(mean(tmp.data_65,1));
        data_SAM.assm_65(mi,:)=data_SAM.assm_65(mi,:)-mean(data_SAM.assm_65(mi,:));
        data_SAM.assm_65(mi,:)=data_SAM.assm_65(mi,:)./std(data_SAM.assm_65(mi,:));
    
        data_SAM.assm(mi,:)=data_SAM.assm_40(mi,:)-data_SAM.assm_65(mi,:);
        data_SAM.assm_lp(mi,:)=lowpass(squeeze(data_SAM.assm(mi,:)), 1/data_SAM.lp_period_month,  'steepness', 0.989);

    end

%     hold on
%     for pi=1:20
%             plot(data_SAM.time,data_SAM.assm_lp(pi,:), 'r', 'linewidth', 0.5)
%     end
%     plot(data_SAM.time,mean(data_SAM.assm_lp,1), 'color', 'r', 'linewidth', 3)
%     plot(data_SAM.time,data_SAM.obs_lp(:), 'color', 'k', 'linewidth', 3)
    

    %% SAM-LENS2

    data_SAM.lens2_members=cfg_lens2.members;
    
    for mi=1:size(data_lens2.(cfg.var),1)
        tmp.data_40=squeeze(data_lens2.(cfg.var)(mi,grid.id_w_40:grid.id_e_40, grid.id_s_40:grid.id_n_40, :));
        tmp.data_40(tmp.data_40<-900)=NaN;
        tmp.data_40(tmp.data_40==0)=NaN;
        data_SAM.lens2_40(mi,:)=squeeze(mean(tmp.data_40,1));
        data_SAM.lens2_40(mi,:)=data_SAM.lens2_40(mi,:)-mean(data_SAM.lens2_40(mi,:));
        data_SAM.lens2_40(mi,:)=data_SAM.lens2_40(mi,:)./std(data_SAM.lens2_40(mi,:));
    
        tmp.data_65=squeeze(data_lens2.(cfg.var)(mi,grid.id_w_65:grid.id_e_65, grid.id_s_65:grid.id_n_65, :));
        tmp.data_65(tmp.data_65<-900)=NaN;
        tmp.data_65(tmp.data_65==0)=NaN;
        data_SAM.lens2_65(mi,:)=squeeze(mean(tmp.data_65,1));
        data_SAM.lens2_65(mi,:)=data_SAM.lens2_65(mi,:)-mean(data_SAM.lens2_65(mi,:));
        data_SAM.lens2_65(mi,:)=data_SAM.lens2_65(mi,:)./std(data_SAM.lens2_65(mi,:));
    
        data_SAM.lens2(mi,:)=data_SAM.lens2_40(mi,:)-data_SAM.lens2_65(mi,:);
        data_SAM.lens2_lp(mi,:)=lowpass(squeeze(data_SAM.lens2(mi,:)), 1/data_SAM.lp_period_month,  'steepness', 0.989);

    end
  


    %% SAM-HCST
    data_SAM.hcst_members=cfg_hcst.members;
    
    for ly=1:5
        tmp.ly_str=['ly',num2str(ly)];
        for mi=1:size(data_hcst.(cfg.var).(tmp.ly_str),1)
            tmp.data_40=squeeze(data_hcst.(cfg.var).(tmp.ly_str)(mi,grid.id_w_40:grid.id_e_40, grid.id_s_40:grid.id_n_40, :));
            tmp.data_40(tmp.data_40<-900)=NaN;
            tmp.data_40(tmp.data_40==0)=NaN;
            data_SAM.hcst_40(ly,mi,:)=squeeze(mean(tmp.data_40,1));
            data_SAM.hcst_40(ly,mi,:)=data_SAM.hcst_40(ly,mi,:)-mean(data_SAM.hcst_40(ly,mi,:));
            data_SAM.hcst_40(ly,mi,:)=data_SAM.hcst_40(ly,mi,:)./std(data_SAM.hcst_40(ly,mi,:));
        
            tmp.data_65=squeeze(data_hcst.(cfg.var).(tmp.ly_str)(mi,grid.id_w_65:grid.id_e_65, grid.id_s_65:grid.id_n_65, :));
            tmp.data_65(tmp.data_65<-900)=NaN;
            tmp.data_65(tmp.data_65==0)=NaN;
            data_SAM.hcst_65(ly,mi,:)=squeeze(mean(tmp.data_65,1));
            data_SAM.hcst_65(ly,mi,:)=data_SAM.hcst_65(ly,mi,:)-mean(data_SAM.hcst_65(ly,mi,:));
            data_SAM.hcst_65(ly,mi,:)=data_SAM.hcst_65(ly,mi,:)./std(data_SAM.hcst_65(ly,mi,:));
        
            data_SAM.hcst(ly,mi,:)=data_SAM.hcst_40(ly,mi,:)-data_SAM.hcst_65(ly,mi,:);
            data_SAM.hcst_lp(ly,mi,:)=lowpass(squeeze(data_SAM.hcst(ly,mi,:)), 1/data_SAM.lp_period_month,  'steepness', 0.989);
        end
    end
    
%     hold on
%     for pi=1:50
%             plot(data_SAM.time,data_SAM.lens2_lp(pi,:), 'g', 'linewidth', 0.5)
%     end
%     plot(data_SAM.time,mean(data_SAM.lens2_lp,1), 'color', 'g', 'linewidth', 3)
% 
%     for pi=1:20
%             plot(data_SAM.time,data_SAM.assm_lp(pi,:), 'r', 'linewidth', 0.5)
%     end
%     plot(data_SAM.time,mean(data_SAM.assm_lp,1), 'color', 'r', 'linewidth', 3)
% 
%     ly=5;
%     for pi=1:20
%             plot(data_SAM.time,squeeze(data_SAM.hcst_lp(ly,pi,:)), 'b', 'linewidth', 0.5)
%     end
%     plot(data_SAM.time+ly-1,squeeze(mean(data_SAM.hcst_lp(ly,:,:),2)), 'color', 'b', 'linewidth', 3)
% 
%     plot(data_SAM.time,data_SAM.obs_lp(:), 'color', 'k', 'linewidth', 3)


    %% SAM - save matfile 
    mkdir([dirs.saveroot, '/clim_indices']);
    matfilename=[dirs.saveroot, '/clim_indices/', 'clim_indices_', cfg.var, '_all_','SAM_', ...
                'obs_', cfg.obs_name, '.mat'];
    save(matfilename, 'cfg', 'cfg_assm', 'cfg_hcst', 'cfg_lens2', 'data_SAM');

end



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

