% %  Created 21-Feb-2024 by Yong-Yub Kim
% % Updated 08-Mar-2024 by Yong-Yub Kim, observation paths

clc; clear all; close all;
warning off;

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

cfg.iyears=1964:2020;
cfg.months=1:12;
cfg.scenname='HIST';
cfg.gridname='f09_g17';
cfg.proj_year=5;


% cfg.vars={'TS', 'SST', 'PRECT', 'TWS', 'GPP', 'FAREA_BURNED', 'photoC_TOT_zint_100m', ...
%     'COL_FIRE_CLOSS', 'SOILWATER_10CM', 'TLAI', 'SSH', 'PSL', 'AEROD_v'};
cfg.vars={'TREFHT', 'TS', 'SST', 'PRECT', 'TWS', 'GPP', 'FAREA_BURNED', 'photoC_TOT_zint_100m', ...
     'SSH', 'PSL', 'NO3', 'COL_FIRE_CLOSS', 'SOILWATER_10CM', 'TLAI', 'Fe', 'SiO3', 'PO4'};
% cfg.vars={'photoC_TOT_zint_100m'};
% cfg.vars={'COL_FIRE_CLOSS', 'SOILWATER_10CM', 'TLAI', 'SSH', 'PSL', 'AEROD_v'};
% cfg.vars={'SSH', 'PSL', 'AEROD_v'};
% cfg.vars={'TS'};
% cfg.vars={'TLAI'};
cfg.vars={'ABCDEFG'};
% cfg.vars={'SSH'};

% cfg.vars={'sithick'};
for vari=1:length(cfg.vars)
    cfg.var=cfg.vars{vari};
    disp(cfg.var);
    tmp.varname=cfg.var;
    cfg.vlayer=1; % surf, vertical slice 
% cfg.vlayer=1:10; % 10layer. don't put more than 15

    cfg.vlayer_1st=min(cfg.vlayer);
    cfg.vlayer_cnt=max(cfg.vlayer)-cfg.vlayer_1st+1;


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


tmp.dimids=[1,1,1];

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
            cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'GlobColour/monthly_reg_pop/PP/', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
        elseif (strcmp(cfg.obs_name, 'CMEMS')==1 && strcmp(tmp.varname, 'photoC_TOT_zint')==1)
            cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'GlobColour/monthly_reg_pop/PP/', tmp.fs, cfg.obs_fname_mid, tmp.fy_str,tmp.m_str, '.nc'];
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
        elseif (strcmp(cfg.obs_name, 'VODCA2GPP')==1)
            cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'monthly_reg_cam',tmp.fs,cfg.obs_fname_mid,tmp.fy_str,tmp.m_str,'.nc'];
        elseif (strcmp(cfg.obs_name, 'SMYLE')==1)
            cfg.obs_fnm='na';
            cfg.obs_fnm_ts=['/mnt/lustre/proj/kimyy/Observation/CRU_TS/monthly_reg_cam/CRU_TS_reg_cesm2.',tmp.fy_str,tmp.m_str,'.nc'];
            cfg.obs_fnm_sst=['/mnt/lustre/proj/kimyy/Observation/HadISST1/monthly_reg_cam/HadISST1_reg_cesm2.',tmp.fy_str,tmp.m_str,'.nc'];
        else
            cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'monthly_reg_', cfg.obs_fname_module(2:4),tmp.fs,cfg.obs_fname_mid,tmp.fy_str,tmp.m_str,'.nc'];
        end
% 
        if exist(cfg.obs_fnm)~=0
            ncid=netcdf.open(cfg.obs_fnm, 'NOWRITE');
            tmpvarid=netcdf.inqVarID(ncid,cfg.obs_varname);
            tmp.dd =  netcdf.getVar(ncid,tmpvarid);
            tmp.dd=double(tmp.dd);
            if (strcmp(cfg.obs_name, 'ERA5')==1 || strcmp(cfg.var, 'TLAI')==1)
                tmp.add_offset=netcdf.getAtt(ncid,tmpvarid,'add_offset');
                tmp.scale_factor=netcdf.getAtt(ncid,tmpvarid,'scale_factor');
                tmp.dd=tmp.dd.*tmp.scale_factor+tmp.add_offset;
            end
             tmp.dd(abs(tmp.dd)>1e30)=NaN;
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
%                 tmp.ydata_obs(1:grid.nlon,1:grid.nlat,mi) = tmp.dd;
            data_obs.([cfg.var])(:,:,(ty-1)*12+mi) = tmp.dd;
            netcdf.close(ncid);
            
        else
%             tmp.ydata_obs(1:grid.nlon,1:grid.nlat,mi) = NaN(grid.nlon, grid.nlat);  
            data_obs.([cfg.var])(:,:,(ty-1)*12+mi) = NaN(grid.nlon, grid.nlat);
            if (strcmp(cfg.obs_name, 'SMYLE')==1)
                %TS
                if exist(cfg.obs_fnm_ts)~=0
                    ncid=netcdf.open(cfg.obs_fnm_ts, 'NOWRITE');
                    tmpvarid=netcdf.inqVarID(ncid,'tmp');
                    tmp.dd =  netcdf.getVar(ncid,tmpvarid);
                    tmp.dd=double(tmp.dd);
                     tmp.dd(abs(tmp.dd)>1e30)=NaN;
                     tmp.dd(abs(tmp.dd)>100)=NaN;
                    tmp.ydata_obs(1:grid.nlon,1:grid.nlat) = tmp.dd;
                    netcdf.close(ncid);
                else
                    tmp.ydata_obs(1:grid.nlon,1:grid.nlat) = NaN;                    
                end
                %SST
                if exist(cfg.obs_fnm_sst)~=0
                    ncid=netcdf.open(cfg.obs_fnm_sst, 'NOWRITE');
                    tmpvarid=netcdf.inqVarID(ncid,'sst');
                    tmp.dd =  netcdf.getVar(ncid,tmpvarid);
                    tmp.dd=double(tmp.dd);
                     tmp.dd(abs(tmp.dd)>1e30)=NaN;
                     tmp.dd(abs(tmp.dd)>100)=NaN;
                    tmp.ydata_obs2(1:grid.nlon,1:grid.nlat) = tmp.dd;
                    netcdf.close(ncid);
                else
                    tmp.ydata_obs2(1:grid.nlon,1:grid.nlat) = NaN;                    
                end
                tmp.ydata_obs(isnan(tmp.ydata_obs))=tmp.ydata_obs2(isnan(tmp.ydata_obs));
                data_obs.([cfg.var])(:,:,(ty-1)*12+mi)=tmp.ydata_obs;
            end
        end
    end
end
% yearly mean
tmp.reshp=reshape(data_obs.([cfg.var]), [grid.nlon, grid.nlat, 12, cfg.len_t_y]);
data_obs.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,3));
data_obs.([cfg.var, '_2ym'])=movmean(data_obs.([cfg.var, '_ym']),2,3, 'Endpoints', 'discard');
data_obs.([cfg.var, '_3ym'])=movmean(data_obs.([cfg.var, '_ym']),3,3, 'Endpoints', 'discard');
data_obs.([cfg.var, '_4ym'])=movmean(data_obs.([cfg.var, '_ym']),4,3, 'Endpoints', 'discard');


% 63~64 : 19~20
% movmean(tmp.ratio_for_EOF(:,85:105,:),12,3,'Endpoints', 'discard')


%% get ASSM members
cfg_assm.list = dir( [dirs.assm_root, '/', '*BHISTsmbb*' ]); 
for kk = 1 : length( cfg_assm.list )
    tmp.fname_in    = cfg_assm.list(kk).name;
    tmp.fname_split = strsplit( tmp.fname_in, {'.'} );
    cfg_assm.members{kk}=[ tmp.fname_split{6}, '.', tmp.fname_split{7}];
%     fyear_str   = strsplit( fname_split{end-1}, '-' );
%     fyear_start = str2num( fyear_str{1}(1:4) );
%     fyear_end   = str2num( fyear_str{2}(1:4) );
%     if( tempyear >= fyear_start && tempyear <= fyear_end &a&     ...                 
%             strcmp( fname_split{2}, 'Omon' ) &&         ...
%             strcmp( fname_split{3}, testname ) &&      ...                 
%             strcmp( fname_split{4}, 'historical' ) )
%         flag_file_in = true;            break;
%     end         
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
data_assm.([cfg.var])=NaN(cfg_assm.len_mem, grid.nlon, grid.nlat, cfg.len_t_y*12); % initialization
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
        
        if ~isfield(tmp, 'dimids')
            ncid = netcdf.open(tmp.fname, 'NOWRITE');
            tmpvarid = netcdf.inqVarID(ncid, cfg.var);
            [tmp.fvarname,tmp.xtype,tmp.dimids,tmp.natts]= netcdf.inqVar(ncid,tmpvarid);
            netcdf.close(ncid);
        end
        if length(tmp.dimids)>3
             tmp.fname2=['/proj/kimyy/tmp/test_mat_',tmp.varname, '.nc'];
             system(['cdo -O -w ', cfg.vm_str, ' -sellevel,', cfg.vl_z1_str,'/',cfg.vl_z2_str, ' ', tmp.fname, ' ', tmp.fname2]);
             try
                 ncid = netcdf.open(tmp.fname2, 'NOWRITE');
                 tmpvarid = netcdf.inqVarID(ncid, cfg.var);
                 tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 0 ind_start-1], [grid.nlon grid.nlat 1 ind_count]));
                 tmp.dd(abs(tmp.dd)>1e30)=NaN;
                 data_assm.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd; %depth averaged value
             catch 
                 pause(60)
                 system(['rm -f ', tmp.fname2]);
                 system(['cdo -O -w ', cfg.vm_str, ' -sellevel,', cfg.vl_z1_str,'/',cfg.vl_z2_str, ' ', tmp.fname, ' ', tmp.fname2]);
                 ncid = netcdf.open(tmp.fname2, 'NOWRITE');
                 tmpvarid = netcdf.inqVarID(ncid, cfg.var);
                 tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 0 ind_start-1], [grid.nlon grid.nlat 1 ind_count]));
                 tmp.dd(abs(tmp.dd)>1e30)=NaN;
                 data_assm.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd; %depth averaged value
             end
        else
            ncid = netcdf.open(tmp.fname, 'NOWRITE');
            tmpvarid = netcdf.inqVarID(ncid, cfg.var);
            tmp.dd=  netcdf.getVar(ncid,tmpvarid, [0 0 ind_start-1], [grid.nlon grid.nlat ind_count]);
            tmp.dd(abs(tmp.dd)>1e30)=NaN;
            data_assm.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd;
        end
        netcdf.close(ncid);
        ty=ty+ind_count/12;
    end
end
% yearly mean
tmp.reshp=reshape(data_assm.([cfg.var]), [cfg_assm.len_mem, grid.nlon, grid.nlat, 12, cfg.len_t_y]);
data_assm.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,4));
data_assm= rmfield(data_assm, cfg.var);
tmp= rmfield(tmp, 'reshp');
data_assm.([cfg.var, '_2ym'])=movmean(data_assm.([cfg.var, '_ym']),2,4, 'Endpoints', 'discard');
data_assm.([cfg.var, '_3ym'])=movmean(data_assm.([cfg.var, '_ym']),3,4, 'Endpoints', 'discard');
data_assm.([cfg.var, '_4ym'])=movmean(data_assm.([cfg.var, '_ym']),4,4, 'Endpoints', 'discard');
data_assm_em.([cfg.var, '_ym'])=squeeze(mean(data_assm.([cfg.var, '_ym']),1));
data_assm_em.([cfg.var, '_2ym'])=squeeze(mean(data_assm.([cfg.var, '_2ym']),1));
data_assm_em.([cfg.var, '_3ym'])=squeeze(mean(data_assm.([cfg.var, '_3ym']),1));
data_assm_em.([cfg.var, '_4ym'])=squeeze(mean(data_assm.([cfg.var, '_4ym']),1));
data_assm_em.spr=squeeze(std(data_assm.([cfg.var, '_ym']),0,1));
data_assm_em.q_17=squeeze(quantile(data_assm.([cfg.var, '_ym']),0.17,1));
data_assm_em.q_83=squeeze(quantile(data_assm.([cfg.var, '_ym']),0.83,1));
data_assm_em.ensmin=squeeze(quantile(data_assm.([cfg.var, '_ym']),0,1));
data_assm_em.ensmax=squeeze(quantile(data_assm.([cfg.var, '_ym']),1,1));
toc;

%% corr, OBS <-> ASSM
tic;
for mi2= 1:cfg_assm.len_mem
    tmp.member2=cfg_assm.members{mi2};
    corrval_assm.obs_assm.corr_member{1, mi2}=[cfg.obs_name, ' x ', tmp.member2];
    disp(['corr, ', corrval_assm.obs_assm.corr_member{1, mi2}]);
    for loni=1:grid.nlon
        for lati=1:grid.nlat
            [tmp.corr, tmp.corr_p]=corrcoef(data_obs.([cfg.var,'_ym'])(loni,lati,:), ...
                data_assm.([cfg.var,'_ym'])(mi2,loni,lati,:), 'Rows', 'complete');
            corrval_assm.obs_assm.val(mi2, loni,lati)=single(tmp.corr(1,2));
%             corrval_assm.obs_assm.p(mi2, loni,lati)=single(tmp.corr_p(1,2));
        end
    end
end
corrval_assm.obs_assm.val_mean=squeeze(mean(corrval_assm.obs_assm.val,1));
corrval_assm.obs_assm.val_median=squeeze(median(corrval_assm.obs_assm.val,1));

if strcmp(tmp.varname, 'SSH')
    for mi2= 1:cfg_assm.len_mem
        tmp.member2=cfg_assm.members{mi2};
        disp(['det(obs) corr, ', corrval_assm.obs_assm.corr_member{1, mi2}]);
        for loni=1:grid.nlon
            for lati=1:grid.nlat
                [tmp.det_data, tmp.trend] = ...
                    Func_0028_detrend_linear_1d(data_obs.([cfg.var,'_ym'])(loni,lati,:), 'omitnan');
		tmp.data2=data_assm.([cfg.var,'_ym'])(mi2,loni,lati,:);
		tmp.data2(isnan(tmp.det_data))=NaN;
                [tmp.det_data2, tmp.trend2] = ...
                    Func_0028_detrend_linear_1d(tmp.data2, 'omitnan');
		
                [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
                    tmp.det_data2, 'Rows', 'complete');
%                [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
%                    data_assm.([cfg.var,'_ym'])(mi2,loni,lati,:), 'Rows', 'complete');
                corrval_assm.obsdet_assm.val(mi2, loni,lati)=single(tmp.corr(1,2));
    %             corrval_assm.obs_assm.p(mi2, loni,lati)=single(tmp.corr_p(1,2));
            end
        end
    end
    corrval_assm.obsdet_assm.val_mean=squeeze(mean(corrval_assm.obsdet_assm.val,1));
    corrval_assm.obsdet_assm.val_median=squeeze(median(corrval_assm.obsdet_assm.val,1));
end

% movm(2~4)
for mvi=2:4
    str_mvi=num2str(mvi);
    disp(['ly', str_mvi]);
    for mi2= 1:cfg_assm.len_mem
        tmp.member2=cfg_assm.members{mi2};
        corrval_assm.obs_assm.corr_member{1, mi2}=[cfg.obs_name, ' x ', tmp.member2];
        disp(['corr, ly', str_mvi, ', ', corrval_assm.obs_assm.corr_member{1, mi2}]);
        for loni=1:grid.nlon
            for lati=1:grid.nlat
                [tmp.corr, tmp.corr_p]=corrcoef(data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:), ...
                    data_assm.([cfg.var,'_', str_mvi, 'ym'])(mi2,loni,lati,:), 'Rows', 'complete');
                corrval_assm.(['obs_assm_', str_mvi, 'ym']).val(mi2, loni,lati)=single(tmp.corr(1,2));
    %             corrval_assm.obs_assm.p(mi2, loni,lati)=single(tmp.corr_p(1,2));
            end
        end
    end
    corrval_assm.(['obs_assm_', str_mvi, 'ym']).val_mean=squeeze(mean(corrval_assm.(['obs_assm_', str_mvi, 'ym']).val,1));
    corrval_assm.(['obs_assm_', str_mvi, 'ym']).val_median=squeeze(median(corrval_assm.(['obs_assm_', str_mvi, 'ym']).val,1));
end

if strcmp(tmp.varname, 'SSH')
    for mvi=2:4
        str_mvi=num2str(mvi);
        disp(['ly', str_mvi]);
        for mi2= 1:cfg_assm.len_mem
            tmp.member2=cfg_assm.members{mi2};
            disp(['corr_det, ly', str_mvi, ', ', corrval_assm.obs_assm.corr_member{1, mi2}]);
            for loni=1:grid.nlon
                for lati=1:grid.nlat
                    [tmp.det_data, tmp.trend] = ...
                    Func_0028_detrend_linear_1d(data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:), 'omitnan');
   	       	    tmp.data2=data_assm.([cfg.var,'_', str_mvi, 'ym'])(mi2,loni,lati,:);
		    tmp.data2(isnan(tmp.det_data))=NaN;
                    [tmp.det_data2, tmp.trend2] = ...
                        Func_0028_detrend_linear_1d(tmp.data2, 'omitnan');

                    [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
                        tmp.det_data2, 'Rows', 'complete');
%                    [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
%                        data_assm.([cfg.var,'_', str_mvi, 'ym'])(mi2,loni,lati,:), 'Rows', 'complete');
                    corrval_assm.(['obsdet_assm_', str_mvi, 'ym']).val(mi2, loni,lati)=single(tmp.corr(1,2));
        %             corrval_assm.obs_assm.p(mi2, loni,lati)=single(tmp.corr_p(1,2));
                end
            end
        end
        corrval_assm.(['obsdet_assm_', str_mvi, 'ym']).val_mean=squeeze(mean(corrval_assm.(['obsdet_assm_', str_mvi, 'ym']).val,1));
        corrval_assm.(['obsdet_assm_', str_mvi, 'ym']).val_median=squeeze(median(corrval_assm.(['obsdet_assm_', str_mvi, 'ym']).val,1));
    end
end


% data_assm.ensmean=squeeze(mean(data_assm.([cfg.var]),1));
% data_assm.ensmean_ym=squeeze(mean(data_assm.([cfg.var, '_ym']),1));

% ensmean based corr
tmp.bootstrpn=100;
for loni=1:grid.nlon
    for lati=1:grid.nlat
        [tmp.corr, tmp.corr_p]=corrcoef(data_obs.([cfg.var,'_ym'])(loni,lati,:), ...
            data_assm_em.([cfg.var, '_ym'])(loni,lati,:), 'Rows', 'complete');
        corrval_assm.obs_assm_em.val(loni,lati)=single(tmp.corr(1,2));
        tmp.a=data_obs.([cfg.var,'_ym'])(loni,lati,:);
        tmp.b=data_assm_em.([cfg.var, '_ym'])(loni,lati,:);
%         if isfinite(tmp.corr(1,2))
%             tmp.btst=bootstrp(tmp.bootstrpn, 'corr', tmp.a(logical(isfinite(tmp.a).*isfinite(tmp.b))), tmp.b(logical(isfinite(tmp.b).*isfinite(tmp.a))));
%             corrval_assm.obs_assm_em.val_btst(loni,lati,1:tmp.bootstrpn)=tmp.btst;
% %             corrval_assm.obs_assm_em.val_005(loni,lati)=quantile(tmp.btst,0.005);
% %             corrval_assm.obs_assm_em.val_025(loni,lati)=quantile(tmp.btst,0.025);
% %             corrval_assm.obs_assm_em.val_050(loni,lati)=quantile(tmp.btst,0.05);
% %             corrval_assm.obs_assm_em.val_170(loni,lati)=quantile(tmp.btst,0.17);
% %             corrval_assm.obs_assm_em.val_995(loni,lati)=quantile(tmp.btst,0.995);
% %             corrval_assm.obs_assm_em.val_975(loni,lati)=quantile(tmp.btst,0.975);
% %             corrval_assm.obs_assm_em.val_950(loni,lati)=quantile(tmp.btst,0.95);
% %             corrval_assm.obs_assm_em.val_830(loni,lati)=quantile(tmp.btst,0.83);
%         else
%             corrval_assm.obs_assm_em.val_btst(loni,lati,1:tmp.bootstrpn)=NaN;
%         end
    end
end

if strcmp(tmp.varname, 'SSH')
    for loni=1:grid.nlon
        for lati=1:grid.nlat
            [tmp.det_data, tmp.trend] = ...
                    Func_0028_detrend_linear_1d(data_obs.([cfg.var,'_ym'])(loni,lati,:), 'omitnan');
            tmp.data2=data_assm_em.([cfg.var, '_ym'])(loni,lati,:);
            tmp.data2(isnan(tmp.det_data))=NaN;
            [tmp.det_data2, tmp.trend] = ...
                    Func_0028_detrend_linear_1d(tmp.data2, 'omitnan');

            [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
                tmp.det_data2, 'Rows', 'complete');
           % [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
           %     data_assm_em.([cfg.var, '_ym'])(loni,lati,:), 'Rows', 'complete');
            corrval_assm.obsdet_assm_em.val(loni,lati)=single(tmp.corr(1,2));
            tmp.a=data_obs.([cfg.var,'_ym'])(loni,lati,:);
            tmp.b=data_assm_em.([cfg.var, '_ym'])(loni,lati,:);
        end
    end
end


% mv
for mvi=2:4
    str_mvi=num2str(mvi);
    disp(['ly', str_mvi]);
    for loni=1:grid.nlon
        for lati=1:grid.nlat
            [tmp.corr, tmp.corr_p]=corrcoef(data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:), ...
                data_assm_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:), 'Rows', 'complete');
            corrval_assm.(['obs_assm_em_', str_mvi, 'ym']).val(loni,lati)=single(tmp.corr(1,2));

            tmp.a=data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:);
            tmp.b=data_assm_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:);
%             if isfinite(tmp.corr(1,2))
%                 tmp.btst=bootstrp(tmp.bootstrpn, 'corr', tmp.a(logical(isfinite(tmp.a).*isfinite(tmp.b))), tmp.b(logical(isfinite(tmp.b).*isfinite(tmp.a))));
%                 corrval_assm.(['obs_assm_em_', str_mvi, 'ym']).val_btst(loni,lati,1:tmp.bootstrpn)=tmp.btst;
% %                 corrval_assm.(['obs_assm_em_', str_mvi, 'ym']).val_005(loni,lati)=quantile(tmp.btst,0.005);
% %                 corrval_assm.(['obs_assm_em_', str_mvi, 'ym']).val_025(loni,lati)=quantile(tmp.btst,0.025);
% %                 corrval_assm.(['obs_assm_em_', str_mvi, 'ym']).val_050(loni,lati)=quantile(tmp.btst,0.05);
% %                 corrval_assm.(['obs_assm_em_', str_mvi, 'ym']).val_170(loni,lati)=quantile(tmp.btst,0.17);
% %                 corrval_assm.(['obs_assm_em_', str_mvi, 'ym']).val_995(loni,lati)=quantile(tmp.btst,0.995);
% %                 corrval_assm.(['obs_assm_em_', str_mvi, 'ym']).val_975(loni,lati)=quantile(tmp.btst,0.975);
% %                 corrval_assm.(['obs_assm_em_', str_mvi, 'ym']).val_950(loni,lati)=quantile(tmp.btst,0.95);
% %                 corrval_assm.(['obs_assm_em_', str_mvi, 'ym']).val_830(loni,lati)=quantile(tmp.btst,0.83);
%             else
%                 corrval_assm.(['obs_assm_em_', str_mvi, 'ym']).val_btst(loni,lati,1:tmp.bootstrpn)=NaN;
%             end
        end
    end
end
toc;
% pcolor(data_obs.TREFHT_ym(:,:,1)'); shading flat; colorbar;
% pcolor(corrval_assm.(['obs_assm_', str_mvi, 'ym']).val_median'); shading flat; colorbar;
% pcolor(corrval_assm.('obs_assm').val_median'); shading flat; colorbar;
% pcolor(corrval_assm.(['obs_assm_em_', str_mvi, 'ym']).val'); shading flat; colorbar;

if strcmp(tmp.varname, 'SSH')
    for mvi=2:4
        str_mvi=num2str(mvi);
        disp(['ly', str_mvi]);
        for loni=1:grid.nlon
            for lati=1:grid.nlat
                [tmp.det_data, tmp.trend] = ...
                    Func_0028_detrend_linear_1d(data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:), 'omitnan');
                tmp.data2=data_assm_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:);
                tmp.data2(isnan(tmp.det_data))=NaN;
                [tmp.det_data2, tmp.trend] = ...
                    Func_0028_detrend_linear_1d(tmp.data2, 'omitnan');
                   
                [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
                    tmp.det_data2, 'Rows', 'complete');
               % [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
               %     data_assm_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:), 'Rows', 'complete');
                corrval_assm.(['obsdet_assm_em_', str_mvi, 'ym']).val(loni,lati)=single(tmp.corr(1,2));
    
                tmp.a=data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:);
                tmp.b=data_assm_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:);
            end
        end
    end
end



%% read LENS2 members
tic;
data_lens2.([cfg.var])=NaN(cfg_lens2.len_mem, grid.nlon, grid.nlat, cfg.len_t_y*12); % initialization

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
        
        if length(tmp.dimids)>3
             tmp.fname2=['/proj/kimyy/tmp/test_mat_',tmp.varname, '.nc'];
             system(['cdo -O -w ', cfg.vm_str, ' -sellevel,', cfg.vl_z1_str,'/',cfg.vl_z2_str, ' ', tmp.fname, ' ', tmp.fname2]);
             try
                 ncid = netcdf.open(tmp.fname2, 'NOWRITE');
                 tmpvarid = netcdf.inqVarID(ncid, cfg.var);
                 tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 0 ind_start-1], [grid.nlon grid.nlat 1 ind_count]));
                 tmp.dd(abs(tmp.dd)>1e30)=NaN;
                 data_lens2.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd; %depth averaged value
             catch
                 pause(60);
                 system(['rm -f ', tmp.fname2]);
                 system(['cdo -O -w ', cfg.vm_str, ' -sellevel,', cfg.vl_z1_str,'/',cfg.vl_z2_str, ' ', tmp.fname, ' ', tmp.fname2]);
                 ncid = netcdf.open(tmp.fname2, 'NOWRITE');
                 tmpvarid = netcdf.inqVarID(ncid, cfg.var);
                 tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 0 ind_start-1], [grid.nlon grid.nlat 1 ind_count]));
                 tmp.dd(abs(tmp.dd)>1e30)=NaN;
                 data_lens2.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd; %depth averaged value
             end

        else
            ncid = netcdf.open(tmp.fname, 'NOWRITE');
            tmpvarid = netcdf.inqVarID(ncid, cfg.var);
            tmp.dd=  netcdf.getVar(ncid,tmpvarid, [0 0 ind_start-1], [grid.nlon grid.nlat ind_count]);
            tmp.dd(abs(tmp.dd)>1e30)=NaN;
            data_lens2.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd;
        end
        netcdf.close(ncid);

        ty=ty+ind_count/12;
    end
end
tmp.reshp=reshape(data_lens2.([cfg.var]), [cfg_lens2.len_mem, grid.nlon, grid.nlat, 12, cfg.len_t_y]);
data_lens2.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,4));
data_lens2= rmfield(data_lens2, cfg.var);
tmp= rmfield(tmp, 'reshp');
data_lens2.([cfg.var, '_2ym'])=movmean(data_lens2.([cfg.var, '_ym']),2,4, 'Endpoints', 'discard');
data_lens2.([cfg.var, '_3ym'])=movmean(data_lens2.([cfg.var, '_ym']),3,4, 'Endpoints', 'discard');
data_lens2.([cfg.var, '_4ym'])=movmean(data_lens2.([cfg.var, '_ym']),4,4, 'Endpoints', 'discard');
data_lens2_em.([cfg.var, '_ym'])=squeeze(mean(data_lens2.([cfg.var, '_ym']),1));
data_lens2_em.([cfg.var, '_2ym'])=squeeze(mean(data_lens2.([cfg.var, '_2ym']),1));
data_lens2_em.([cfg.var, '_3ym'])=squeeze(mean(data_lens2.([cfg.var, '_3ym']),1));
data_lens2_em.([cfg.var, '_4ym'])=squeeze(mean(data_lens2.([cfg.var, '_4ym']),1));
data_lens2_em.spr=squeeze(std(data_lens2.([cfg.var, '_ym']),0,1));
data_lens2_em.q_17=squeeze(quantile(data_lens2.([cfg.var, '_ym']),0.17,1));
data_lens2_em.q_83=squeeze(quantile(data_lens2.([cfg.var, '_ym']),0.83,1));
data_lens2_em.ensmin=squeeze(quantile(data_lens2.([cfg.var, '_ym']),0,1));
data_lens2_em.ensmax=squeeze(quantile(data_lens2.([cfg.var, '_ym']),1,1));
toc;


%% corr, OBS <-> LENS2
tic;
for mi2= 1:cfg_lens2.len_mem
    tmp.member2=cfg_lens2.members{mi2};
    corrval_lens2.obs_lens2.corr_member{1, mi2}=[cfg.obs_name, ' x ', tmp.member2];
    disp(['corr, ', corrval_lens2.obs_lens2.corr_member{1, mi2}]);
    for loni=1:grid.nlon
        for lati=1:grid.nlat
            [tmp.corr, tmp.corr_p]=corrcoef(data_obs.([cfg.var,'_ym'])(loni,lati,:), ...
                data_lens2.([cfg.var,'_ym'])(mi2,loni,lati,:), 'Rows', 'complete');
            corrval_lens2.obs_lens2.val(mi2, loni,lati)=single(tmp.corr(1,2));
            
        end
    end
end
corrval_lens2.obs_lens2.val_mean=squeeze(mean(corrval_lens2.obs_lens2.val,1));
corrval_lens2.obs_lens2.val_median=squeeze(median(corrval_lens2.obs_lens2.val,1));

if strcmp(tmp.varname, 'SSH')
    for mi2= 1:cfg_lens2.len_mem
        tmp.member2=cfg_lens2.members{mi2};
        disp(['corr, ', corrval_lens2.obs_lens2.corr_member{1, mi2}]);
        for loni=1:grid.nlon
            for lati=1:grid.nlat
                [tmp.det_data, tmp.trend] = ...
                    Func_0028_detrend_linear_1d(data_obs.([cfg.var,'_ym'])(loni,lati,:), 'omitnan');
                tmp.data2=data_lens2.([cfg.var,'_ym'])(mi2,loni,lati,:);
                tmp.data2(isnan(tmp.det_data))=NaN;
                [tmp.det_data2, tmp.trend] = ...
                    Func_0028_detrend_linear_1d(tmp.data2, 'omitnan');
                
                [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
                    tmp.det_data2, 'Rows', 'complete');
               % [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
               %     data_lens2.([cfg.var,'_ym'])(mi2,loni,lati,:), 'Rows', 'complete');
                corrval_lens2.obsdet_lens2.val(mi2, loni,lati)=single(tmp.corr(1,2));
            end
        end
    end
    corrval_lens2.obsdet_lens2.val_mean=squeeze(mean(corrval_lens2.obsdet_lens2.val,1));
    corrval_lens2.obsdet_lens2.val_median=squeeze(median(corrval_lens2.obsdet_lens2.val,1));
end


% movm(2~4)
for mvi=2:4
    str_mvi=num2str(mvi);
    disp(['ly', str_mvi]);
    for mi2= 1:cfg_lens2.len_mem
        tmp.member2=cfg_lens2.members{mi2};
        corrval_lens2.obs_lens2.corr_member{1, mi2}=[cfg.obs_name, ' x ', tmp.member2];
        disp(['corr, ly', str_mvi, ', ', corrval_lens2.obs_lens2.corr_member{1, mi2}]);
        for loni=1:grid.nlon
            for lati=1:grid.nlat
                [tmp.corr, tmp.corr_p]=corrcoef(data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:), ...
                    data_lens2.([cfg.var,'_', str_mvi, 'ym'])(mi2,loni,lati,:), 'Rows', 'complete');
                corrval_lens2.(['obs_lens2_', str_mvi, 'ym']).val(mi2, loni,lati)=single(tmp.corr(1,2));
    %             corrval_lens2.obs_lens2.p(mi2, loni,lati)=single(tmp.corr_p(1,2));
            end
        end
    end
    corrval_lens2.(['obs_lens2_', str_mvi, 'ym']).val_mean=squeeze(mean(corrval_lens2.(['obs_lens2_', str_mvi, 'ym']).val,1));
    corrval_lens2.(['obs_lens2_', str_mvi, 'ym']).val_median=squeeze(median(corrval_lens2.(['obs_lens2_', str_mvi, 'ym']).val,1));
end

if strcmp(tmp.varname, 'SSH')
    for mvi=2:4
        str_mvi=num2str(mvi);
        disp(['ly', str_mvi]);
        for mi2= 1:cfg_lens2.len_mem
            tmp.member2=cfg_lens2.members{mi2};
            disp(['corr, ly', str_mvi, ', ', corrval_lens2.obs_lens2.corr_member{1, mi2}]);
            for loni=1:grid.nlon
                for lati=1:grid.nlat
                    [tmp.det_data, tmp.trend] = ...
                    Func_0028_detrend_linear_1d(data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:), 'omitnan');
                    tmp.data2=data_lens2.([cfg.var,'_', str_mvi, 'ym'])(mi2,loni,lati,:);
                    tmp.data2(isnan(tmp.det_data))=NaN;
                    [tmp.det_data2, tmp.trend] = ...
                    Func_0028_detrend_linear_1d(tmp.data2, 'omitnan');
 
                   % [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
                   %     data_lens2.([cfg.var,'_', str_mvi, 'ym'])(mi2,loni,lati,:), 'Rows', 'complete');
                    [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
                        tmp.det_data2, 'Rows', 'complete');
                    corrval_lens2.(['obsdet_lens2_', str_mvi, 'ym']).val(mi2, loni,lati)=single(tmp.corr(1,2));
        %             corrval_lens2.obs_lens2.p(mi2, loni,lati)=single(tmp.corr_p(1,2));
                end
            end
        end
        corrval_lens2.(['obsdet_lens2_', str_mvi, 'ym']).val_mean=squeeze(mean(corrval_lens2.(['obsdet_lens2_', str_mvi, 'ym']).val,1));
        corrval_lens2.(['obsdet_lens2_', str_mvi, 'ym']).val_median=squeeze(median(corrval_lens2.(['obsdet_lens2_', str_mvi, 'ym']).val,1));
    end
end


% ensmean based corr
for loni=1:grid.nlon
    for lati=1:grid.nlat
        [tmp.corr, tmp.corr_p]=corrcoef(data_obs.([cfg.var,'_ym'])(loni,lati,:), ...
            data_lens2_em.([cfg.var, '_ym'])(loni,lati,:), 'Rows', 'complete');
        corrval_lens2.obs_lens2_em.val(loni,lati)=single(tmp.corr(1,2));

        tmp.a=data_obs.([cfg.var,'_ym'])(loni,lati,:);
        tmp.b=data_lens2_em.([cfg.var, '_ym'])(loni,lati,:);
%         if isfinite(tmp.corr(1,2))
%             tmp.btst=bootstrp(tmp.bootstrpn, 'corr', tmp.a(logical(isfinite(tmp.a).*isfinite(tmp.b))), tmp.b(logical(isfinite(tmp.b).*isfinite(tmp.a))));
%             corrval_lens2.obs_lens2_em.val_btst(loni,lati,1:tmp.bootstrpn)=tmp.btst;
% %             corrval_lens2.obs_lens2_em.val_005(loni,lati)=quantile(tmp.btst,0.005);
% %             corrval_lens2.obs_lens2_em.val_025(loni,lati)=quantile(tmp.btst,0.025);
% %             corrval_lens2.obs_lens2_em.val_050(loni,lati)=quantile(tmp.btst,0.05);
% %             corrval_lens2.obs_lens2_em.val_170(loni,lati)=quantile(tmp.btst,0.17);
% %             corrval_lens2.obs_lens2_em.val_995(loni,lati)=quantile(tmp.btst,0.995);
% %             corrval_lens2.obs_lens2_em.val_975(loni,lati)=quantile(tmp.btst,0.975);
% %             corrval_lens2.obs_lens2_em.val_950(loni,lati)=quantile(tmp.btst,0.95);
% %             corrval_lens2.obs_lens2_em.val_830(loni,lati)=quantile(tmp.btst,0.83);
%         else
%             corrval_lens2.obs_lens2_em.val_btst(loni,lati,1:tmp.bootstrpn)=NaN;
%         end
    end
end
if strcmp(tmp.varname, 'SSH')
    for loni=1:grid.nlon
        for lati=1:grid.nlat
             [tmp.det_data, tmp.trend] = ...
                    Func_0028_detrend_linear_1d(data_obs.([cfg.var,'_ym'])(loni,lati,:), 'omitnan');
            tmp.data2=data_lens2_em.([cfg.var, '_ym'])(loni,lati,:);
            tmp.data2(isnan(tmp.det_data))=NaN;
             [tmp.det_data2, tmp.trend] = ...
                    Func_0028_detrend_linear_1d(tmp.data2, 'omitnan');
          
            [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
                tmp.det_data2, 'Rows', 'complete');
           % [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
           %     data_lens2_em.([cfg.var, '_ym'])(loni,lati,:), 'Rows', 'complete');
            corrval_lens2.obsdet_lens2_em.val(loni,lati)=single(tmp.corr(1,2));
    
            tmp.a=data_obs.([cfg.var,'_ym'])(loni,lati,:);
            tmp.b=data_lens2_em.([cfg.var, '_ym'])(loni,lati,:);
        end
    end
end


% movm(2~4)
for mvi=2:4
    str_mvi=num2str(mvi);
    disp(['ly', str_mvi]);
    for loni=1:grid.nlon
        for lati=1:grid.nlat
            [tmp.corr, tmp.corr_p]=corrcoef(data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:), ...
                data_lens2_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:), 'Rows', 'complete');
            corrval_lens2.(['obs_lens2_em_', str_mvi, 'ym']).val(loni,lati)=single(tmp.corr(1,2));

            tmp.a=data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:);
            tmp.b=data_lens2_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:);
%             if isfinite(tmp.corr(1,2))
%                 tmp.btst=bootstrp(tmp.bootstrpn, 'corr', tmp.a(logical(isfinite(tmp.a).*isfinite(tmp.b))), tmp.b(logical(isfinite(tmp.b).*isfinite(tmp.a))));
%                 corrval_lens2.(['obs_lens2_em_', str_mvi, 'ym']).val_btst(loni,lati,1:tmp.bootstrpn)=tmp.btst;
% %                 corrval_lens2.(['obs_lens2_em_', str_mvi, 'ym']).val_005(loni,lati)=quantile(tmp.btst,0.005);
% %                 corrval_lens2.(['obs_lens2_em_', str_mvi, 'ym']).val_025(loni,lati)=quantile(tmp.btst,0.025);
% %                 corrval_lens2.(['obs_lens2_em_', str_mvi, 'ym']).val_050(loni,lati)=quantile(tmp.btst,0.05);
% %                 corrval_lens2.(['obs_lens2_em_', str_mvi, 'ym']).val_170(loni,lati)=quantile(tmp.btst,0.17);
% %                 corrval_lens2.(['obs_lens2_em_', str_mvi, 'ym']).val_995(loni,lati)=quantile(tmp.btst,0.995);
% %                 corrval_lens2.(['obs_lens2_em_', str_mvi, 'ym']).val_975(loni,lati)=quantile(tmp.btst,0.975);
% %                 corrval_lens2.(['obs_lens2_em_', str_mvi, 'ym']).val_950(loni,lati)=quantile(tmp.btst,0.95);
% %                 corrval_lens2.(['obs_lens2_em_', str_mvi, 'ym']).val_830(loni,lati)=quantile(tmp.btst,0.83);
%             else
%                 corrval_lens2.(['obs_lens2_em_', str_mvi, 'ym']).val_btst(loni,lati,1:tmp.bootstrpn)=NaN;
%             end
        end
    end
end

if strcmp(tmp.varname, 'SSH')
    for mvi=2:4
        str_mvi=num2str(mvi);
        disp(['ly', str_mvi]);
        for loni=1:grid.nlon
            for lati=1:grid.nlat
		[tmp.det_data, tmp.trend] = ...
                    Func_0028_detrend_linear_1d(data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:), 'omitnan');
                tmp.data2=data_lens2_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:);
                tmp.data2(isnan(tmp.det_data))=NaN;
		[tmp.det_data2, tmp.trend] = ...
                    Func_0028_detrend_linear_1d(tmp.data2, 'omitnan');
                
                [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
                    tmp.det_data2, 'Rows', 'complete');
               % [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
               %     data_lens2_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:), 'Rows', 'complete');
                corrval_lens2.(['obsdet_lens2_em_', str_mvi, 'ym']).val(loni,lati)=single(tmp.corr(1,2));
    
                tmp.a=data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:);
                tmp.b=data_lens2_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:);
            end
        end
    end
end


%% set valid grid indice for data compression
tmp.data=data_assm_em.([cfg.var, '_ym'])(:,:,1);
if strcmp(tmp.varname, 'SST')
    tmp.data(tmp.data==0)=NaN;
end
tmp.vl=~isnan(tmp.data);
grid.valid_ind=find(tmp.vl(:));
[grid.valid_ind_i, grid.valid_ind_j]=find(tmp.vl);

%% corr, ASSM <-> LENS2
tic;
for mi1= 1:cfg_assm.len_mem
    tmp.member1=cfg_assm.members{mi1};
    for mi2= 1:cfg_lens2.len_mem
        tmp.member2=cfg_lens2.members{mi2};
        corrval_lens2.assm_lens2.corr_member{1, (mi1-1)*cfg_lens2.len_mem+mi2}=[tmp.member1, ' x ', tmp.member2];
        disp([tmp.varname, 'corr, ', corrval_lens2.assm_lens2.corr_member{1, (mi1-1)*cfg_lens2.len_mem+mi2}]);
%         for loni=1:grid.nlon
%             for lati=1:grid.nlat
%                 [tmp.corr, tmp.corr_p]=corrcoef(data_assm.([cfg.var,'_ym'])(mi1,loni,lati,:), ...
%                     data_lens2.([cfg.var,'_ym'])(mi2,loni,lati,:), 'Rows', 'complete');
%                 corrval_lens2.assm_lens2.val((mi1-1)*cfg_lens2.len_mem+mi2, loni,lati)=single(tmp.corr(1,2));
%                 corrval_lens2.assm_lens2.p((mi1-1)*cfg_lens2.len_mem+mi2, loni,lati)=single(tmp.corr_p(1,2));
%             end
%         end
        for vali=1:length(grid.valid_ind)
            loni=grid.valid_ind_i(vali);
            lati=grid.valid_ind_j(vali);
                [tmp.corr, tmp.corr_p]=corrcoef(data_assm.([cfg.var,'_ym'])(mi1,loni,lati,:), ...
                    data_lens2.([cfg.var,'_ym'])(mi2,loni,lati,:), 'Rows', 'complete');
                corrval_lens2.assm_lens2.val((mi1-1)*cfg_lens2.len_mem+mi2, vali)=single(tmp.corr(1,2));
        end
    end
end
corrval_lens2.assm_lens2.val_mean=squeeze(mean(corrval_lens2.assm_lens2.val,1));
corrval_lens2.assm_lens2.val_median=squeeze(median(corrval_lens2.assm_lens2.val,1));

for mvi=2:4
    str_mvi=num2str(mvi);
    disp(['ly', str_mvi]);
    for mi1= 1:cfg_assm.len_mem
        tmp.member1=cfg_assm.members{mi1};
        for mi2= 1:cfg_lens2.len_mem
            tmp.member2=cfg_lens2.members{mi2};
            corrval_lens2.assm_lens2.corr_member{1, (mi1-1)*cfg_lens2.len_mem+mi2}=[tmp.member1, ' x ', tmp.member2];
            disp([tmp.varname, 'corr, ly', str_mvi, ', ', corrval_lens2.assm_lens2.corr_member{1, (mi1-1)*cfg_lens2.len_mem+mi2}]);
            for vali=1:length(grid.valid_ind)
                loni=grid.valid_ind_i(vali);
                lati=grid.valid_ind_j(vali);
                    [tmp.corr, tmp.corr_p]=corrcoef(data_assm.([cfg.var,'_', str_mvi, 'ym'])(mi1,loni,lati,:), ...
                        data_lens2.([cfg.var,'_', str_mvi, 'ym'])(mi2,loni,lati,:), 'Rows', 'complete');
                    corrval_lens2.(['assm_lens2_', str_mvi, 'ym']).val((mi1-1)*cfg_lens2.len_mem+mi2, vali)=single(tmp.corr(1,2));
            end
        end
    end
    corrval_lens2.(['assm_lens2_', str_mvi, 'ym']).val_mean=squeeze(mean(corrval_lens2.(['assm_lens2_', str_mvi, 'ym']).val,1));
    corrval_lens2.(['assm_lens2_', str_mvi, 'ym']).val_median=squeeze(median(corrval_lens2.(['assm_lens2_', str_mvi, 'ym']).val,1));
end

% ensmean based corr
% data_lens2.ensmean=squeeze(mean(data_lens2.([cfg.var]),1));
% data_lens2_em.([cfg.var, '_ym'])=squeeze(mean(data_lens2.([cfg.var, '_ym']),1));
for loni=1:grid.nlon
    for lati=1:grid.nlat
        [tmp.corr, tmp.corr_p]=corrcoef(data_assm_em.([cfg.var, '_ym'])(loni,lati,:), ...
            data_lens2_em.([cfg.var, '_ym'])(loni,lati,:));
        corrval_lens2.assm_lens2_em.val(loni,lati)=single(tmp.corr(1,2));

        tmp.a=data_assm_em.([cfg.var, '_ym'])(loni,lati,:);
        tmp.b=data_lens2_em.([cfg.var, '_ym'])(loni,lati,:);
%         if isfinite(tmp.corr(1,2))
%             tmp.btst=bootstrp(tmp.bootstrpn, 'corr', tmp.a(logical(isfinite(tmp.a).*isfinite(tmp.b))), tmp.b(logical(isfinite(tmp.b).*isfinite(tmp.a))));
%             corrval_lens2.assm_lens2_em.val_btst(loni,lati,1:tmp.bootstrpn)=tmp.btst;
% %             corrval_lens2.assm_lens2_em.val_005(loni,lati)=quantile(tmp.btst,0.005);
% %             corrval_lens2.assm_lens2_em.val_025(loni,lati)=quantile(tmp.btst,0.025);
% %             corrval_lens2.assm_lens2_em.val_050(loni,lati)=quantile(tmp.btst,0.05);
% %             corrval_lens2.assm_lens2_em.val_170(loni,lati)=quantile(tmp.btst,0.17);
% %             corrval_lens2.assm_lens2_em.val_995(loni,lati)=quantile(tmp.btst,0.995);
% %             corrval_lens2.assm_lens2_em.val_975(loni,lati)=quantile(tmp.btst,0.975);
% %             corrval_lens2.assm_lens2_em.val_950(loni,lati)=quantile(tmp.btst,0.95);
% %             corrval_lens2.assm_lens2_em.val_830(loni,lati)=quantile(tmp.btst,0.83);
%         else
%             corrval_lens2.assm_lens2_em.val_btst(loni,lati,1:tmp.bootstrpn)=NaN;
%         end
    end
end

for mvi=2:4
    str_mvi=num2str(mvi);
    disp(['ly', str_mvi]);
    for loni=1:grid.nlon
        for lati=1:grid.nlat
            [tmp.corr, tmp.corr_p]=corrcoef(data_assm_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:), ...
                data_lens2_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:));
            corrval_lens2.(['assm_lens2_em_', str_mvi, 'ym']).val(loni,lati)=single(tmp.corr(1,2));

            tmp.a=data_assm_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:);
            tmp.b=data_lens2_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:);
%             if isfinite(tmp.corr(1,2))
%                 tmp.btst=bootstrp(tmp.bootstrpn, 'corr', tmp.a(logical(isfinite(tmp.a).*isfinite(tmp.b))), tmp.b(logical(isfinite(tmp.b).*isfinite(tmp.a))));
%                 corrval_lens2.(['assm_lens2_em_', str_mvi, 'ym']).val_btst(loni,lati,1:tmp.bootstrpn)=tmp.btst;
% %                 corrval_lens2.(['assm_lens2_em_', str_mvi, 'ym']).val_005(loni,lati)=quantile(tmp.btst,0.005);
% %                 corrval_lens2.(['assm_lens2_em_', str_mvi, 'ym']).val_025(loni,lati)=quantile(tmp.btst,0.025);
% %                 corrval_lens2.(['assm_lens2_em_', str_mvi, 'ym']).val_050(loni,lati)=quantile(tmp.btst,0.05);
% %                 corrval_lens2.(['assm_lens2_em_', str_mvi, 'ym']).val_170(loni,lati)=quantile(tmp.btst,0.17);
% %                 corrval_lens2.(['assm_lens2_em_', str_mvi, 'ym']).val_995(loni,lati)=quantile(tmp.btst,0.995);
% %                 corrval_lens2.(['assm_lens2_em_', str_mvi, 'ym']).val_975(loni,lati)=quantile(tmp.btst,0.975);
% %                 corrval_lens2.(['assm_lens2_em_', str_mvi, 'ym']).val_950(loni,lati)=quantile(tmp.btst,0.95);
% %                 corrval_lens2.(['assm_lens2_em_', str_mvi, 'ym']).val_830(loni,lati)=quantile(tmp.btst,0.83);
%             else
%                 corrval_lens2.(['assm_lens2_em_', str_mvi, 'ym']).val_btst(loni,lati,1:tmp.bootstrpn)=NaN;
%             end
        end
    end
end

toc;



% pcolor(grid.lon, grid.lat,squeeze(corrval.assm_lens2_em.val-corrval.assm_lens2.val_mean)'); shading flat; colorbar;
% pcolor(grid.lon, grid.lat,squeeze(corrval.assm_lens2_em.val-corrval.assm_lens2.val_median)'); shading flat; colorbar;
% caxis([-0.3 0.7])
% caxis([0 0.4])
% figure; pcolor(grid.lon, grid.lat,squeeze(corrval.assm_lens2.val_mean)'); shading flat; colorbar;
% pcolor(grid.lon, grid.lat,squeeze(corrval.assm_lens2_em.val)'); shading flat; colorbar;

% tic
% for i=1:100000
% corrcoef(rand(61,1), rand(61,1));
% end
% toc
% 
% tic
% for i=1:100000
% corr(rand(61,1), rand(61,1));
% end
% toc






%% save matfile (OBS, ASSM)

matfilename=[dirs.saveroot, '/corr_raw/', 'corr_assm_',cfg.var,'_v', num2str(cfg.vlayer_1st), '_v', num2str(max(cfg.vlayer)), '.mat'];
corrval_assm.grid=grid;
corrval_assm.data=data_assm_em;
corrval_assm.data_obs=data_obs;
save(matfilename, '-struct', 'corrval_assm', '-v7.3');

matfilename=[dirs.saveroot, '/corr_raw/', 'corr_lens2_',cfg.var,'_v', num2str(cfg.vlayer_1st), '_v', num2str(max(cfg.vlayer)), '.mat'];
corrval_lens2.grid=grid;
corrval_lens2.data=data_lens2_em;
save(matfilename, '-struct', 'corrval_lens2', '-v7.3');

clear corrval_assm corrval_lens2 
data_lens2= rmfield(data_lens2, [cfg.var,'_ym']);
data_lens2= rmfield(data_lens2, [cfg.var,'_2ym']);
data_lens2= rmfield(data_lens2, [cfg.var,'_3ym']);
data_lens2= rmfield(data_lens2, [cfg.var,'_4ym']);


%% read HCST members

for ly=1:5 %1:5
    tmp.ly_str=['ly',num2str(ly)];
    data_hcst.([cfg.var])=NaN(cfg_hcst.len_mem, grid.nlon, grid.nlat, cfg.len_t_y*12); % initialization
    
    for mi= 1:cfg_hcst.len_mem
        ty=1;
        tmp.member=cfg_hcst.members{mi};
        disp(['hcst, ', tmp.member]);
        while ty<=cfg.len_t_y
            tmp.year=cfg.iyears(ty)-ly+1;
            tmp.yearstr=num2str(tmp.year);
            tmp.casenm = ['f09_g17.hcst.',tmp.member];
            tmp.casenm_i = [tmp.casenm, '_i', tmp.yearstr];
            tmp.fdir=[dirs.hcst_root, '/', tmp.casenm, '/' tmp.casenm_i, '/', cfg.comp, '/proc/tseries/month_1'];
            tmp.flist=dir( [tmp.fdir, '/','*.',cfg.var, '.*'] );
            tmp.fname=[tmp.fdir, '/', tmp.flist(1).name];
           
            ind_start=(ly - 1)*12+1;
            ind_count=12;
            
            if length(tmp.dimids)>3
                tmp.fname2=['/proj/kimyy/tmp/test_mat_',tmp.varname, '.nc'];
                system(['cdo -O -w ', cfg.vm_str, ' -sellevel,', cfg.vl_z1_str,'/',cfg.vl_z2_str, ' ', tmp.fname, ' ', tmp.fname2]);
                try
                    ncid = netcdf.open(tmp.fname2, 'NOWRITE');
                    tmpvarid = netcdf.inqVarID(ncid, cfg.var);
                    tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 0 ind_start-1], [grid.nlon grid.nlat 1 ind_count]));
                    tmp.dd(abs(tmp.dd)>1e30)=NaN;
                    data_hcst.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd; %depth averaged value
                catch 
                    pause(60)
                    system(['rm -f ', tmp.fname2]);
                    system(['cdo -O -w ', cfg.vm_str, ' -sellevel,', cfg.vl_z1_str,'/',cfg.vl_z2_str, ' ', tmp.fname, ' ', tmp.fname2]);
                    ncid = netcdf.open(tmp.fname2, 'NOWRITE');
                    tmpvarid = netcdf.inqVarID(ncid, cfg.var);
                    tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 0 ind_start-1], [grid.nlon grid.nlat 1 ind_count]));
                    tmp.dd(abs(tmp.dd)>1e30)=NaN;
                    data_hcst.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd; %depth averaged value
                end
            else
                ncid = netcdf.open(tmp.fname, 'NOWRITE');
                tmpvarid = netcdf.inqVarID(ncid, cfg.var);
                tmp.dd=  netcdf.getVar(ncid,tmpvarid, [0 0 ind_start-1], [grid.nlon grid.nlat ind_count]);
                tmp.dd(abs(tmp.dd)>1e30)=NaN;
                data_hcst.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd;
            end
            netcdf.close(ncid);
            ty=ty+ind_count/12;
        end
    end
    % yearly mean
    tmp.reshp=reshape(data_hcst.([cfg.var]), [cfg_hcst.len_mem, grid.nlon, grid.nlat, 12, cfg.len_t_y]);
    data_hcst.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,4));
    tmp= rmfield(tmp, 'reshp');
    data_hcst= rmfield(data_hcst, cfg.var);
    data_hcst_em.(tmp.ly_str).([cfg.var, '_ym'])=squeeze(mean(data_hcst.([cfg.var, '_ym']),1));
    data_hcst_em.(tmp.ly_str).spr=squeeze(std(data_hcst.([cfg.var, '_ym']),0,1));
    data_hcst_em.(tmp.ly_str).q_17=squeeze(quantile(data_hcst.([cfg.var, '_ym']),0.17,1));
    data_hcst_em.(tmp.ly_str).q_83=squeeze(quantile(data_hcst.([cfg.var, '_ym']),0.83,1));
    data_hcst_em.(tmp.ly_str).ensmin=squeeze(quantile(data_hcst.([cfg.var, '_ym']),0,1));
    data_hcst_em.(tmp.ly_str).ensmax=squeeze(quantile(data_hcst.([cfg.var, '_ym']),1,1));

    if ly==2
        disp('ly2sum')
        data_hcst.([cfg.var, '_2ym'])=data_hcst.([cfg.var, '_ym'])(:,:,:,1:end-1);
        data_hcst.([cfg.var, '_3ym'])=data_hcst.([cfg.var, '_ym'])(:,:,:,1:end-2);
        data_hcst.([cfg.var, '_4ym'])=data_hcst.([cfg.var, '_ym'])(:,:,:,1:end-3);
    end
    if ly==3
        disp('ly3sum')
        data_hcst.([cfg.var, '_2ym'])=data_hcst.([cfg.var, '_2ym'])+data_hcst.([cfg.var, '_ym'])(:,:,:,1:end-1);
        data_hcst.([cfg.var, '_3ym'])=data_hcst.([cfg.var, '_3ym'])+data_hcst.([cfg.var, '_ym'])(:,:,:,1:end-2);
        data_hcst.([cfg.var, '_4ym'])=data_hcst.([cfg.var, '_4ym'])+data_hcst.([cfg.var, '_ym'])(:,:,:,1:end-3);
    end
    if ly==4
        disp('ly4sum')
        data_hcst.([cfg.var, '_3ym'])=data_hcst.([cfg.var, '_3ym'])+data_hcst.([cfg.var, '_ym'])(:,:,:,1:end-2);
        data_hcst.([cfg.var, '_4ym'])=data_hcst.([cfg.var, '_4ym'])+data_hcst.([cfg.var, '_ym'])(:,:,:,1:end-3);
    end
    if ly==5
        disp('ly5sum')
        data_hcst.([cfg.var, '_4ym'])=data_hcst.([cfg.var, '_4ym'])+data_hcst.([cfg.var, '_ym'])(:,:,:,1:end-3);
    end

    %% corr, ASSM <-> HCST
    tic;
    for mi1= 1:cfg_assm.len_mem
        tmp.member1=cfg_assm.members{mi1};
        for mi2= 1:cfg_hcst.len_mem
            tmp.member2=cfg_hcst.members{mi2};
            corrval_hcst.assm_hcst.(tmp.ly_str).corr_member{1, (mi1-1)*cfg_hcst.len_mem+mi2}=[tmp.member1, ' x ', tmp.member2];
            disp([tmp.varname,', ly:', num2str(ly), ', corr, ', corrval_hcst.assm_hcst.(tmp.ly_str).corr_member{1, (mi1-1)*cfg_hcst.len_mem+mi2}]);
%             for loni=1:grid.nlon
%                 for lati=1:grid.nlat
%                     [tmp.corr, tmp.corr_p]=corrcoef(data_assm.([cfg.var,'_ym'])(mi1,loni,lati,:), ...
%                         data_hcst.([cfg.var,'_ym'])(mi2,loni,lati,:), 'Rows', 'complete');
%                     corrval_hcst.assm_hcst.(tmp.ly_str).val((mi1-1)*cfg_hcst.len_mem+mi2, loni,lati)=single(tmp.corr(1,2));
%                     corrval_hcst.assm_hcst.(tmp.ly_str).p((mi1-1)*cfg_hcst.len_mem+mi2, loni,lati)=single(tmp.corr_p(1,2));
% 
%                     %% internal signal
%                     [tmp.corr, tmp.corr_p]=corrcoef(data_assm.([cfg.var,'_ym'])(mi1,loni,lati,:), ...
%                         squeeze(data_hcst.([cfg.var,'_ym'])(mi2,loni,lati,:)) - ...
%                         squeeze(data_lens2_em.([cfg.var, '_ym'])(loni,lati,:)), 'Rows', 'complete');
%                     corrval_hcst.assm_hcst_int.(tmp.ly_str).val((mi1-1)*cfg_hcst.len_mem+mi2, loni,lati)=single(tmp.corr(1,2));
%                     corrval_hcst.assm_hcst_int.(tmp.ly_str).p((mi1-1)*cfg_hcst.len_mem+mi2, loni,lati)=single(tmp.corr_p(1,2));
% 
%                 end
%             end
            for vali=1:length(grid.valid_ind)
                loni=grid.valid_ind_i(vali);
                lati=grid.valid_ind_j(vali);
                 [tmp.corr, tmp.corr_p]=corrcoef(data_assm.([cfg.var,'_ym'])(mi1,loni,lati,:), ...
                    data_hcst.([cfg.var,'_ym'])(mi2,loni,lati,:), 'Rows', 'complete');
                corrval_hcst.assm_hcst.(tmp.ly_str).val((mi1-1)*cfg_hcst.len_mem+mi2, vali)=single(tmp.corr(1,2));
%                 corrval_hcst.assm_hcst.(tmp.ly_str).p((mi1-1)*cfg_hcst.len_mem+mi2, vali)=single(tmp.corr_p(1,2));

                %% internal signal
                [tmp.corr, tmp.corr_p]=corrcoef(data_assm.([cfg.var,'_ym'])(mi1,loni,lati,:), ...
                    squeeze(data_hcst.([cfg.var,'_ym'])(mi2,loni,lati,:)) - ...
                    squeeze(data_lens2_em.([cfg.var, '_ym'])(loni,lati,:)), 'Rows', 'complete');
                corrval_hcst.assm_hcst_int.(tmp.ly_str).val((mi1-1)*cfg_hcst.len_mem+mi2, vali)=single(tmp.corr(1,2));
%                 corrval_hcst.assm_hcst_int.(tmp.ly_str).p((mi1-1)*cfg_hcst.len_mem+mi2, vali)=single(tmp.corr_p(1,2));
            end
        end
    end
    corrval_hcst.assm_hcst.(tmp.ly_str).val_mean=squeeze(mean(corrval_hcst.assm_hcst.(tmp.ly_str).val,1));
    corrval_hcst.assm_hcst.(tmp.ly_str).val_median=squeeze(median(corrval_hcst.assm_hcst.(tmp.ly_str).val,1));
    
    % ensmean based corr
%     data_assm.ensmean=squeeze(mean(data_assm.([cfg.var]),1));
%     data_hcst.ensmean=squeeze(mean(data_hcst.([cfg.var]),1));
%     data_assm.ensmean_ym=squeeze(mean(data_assm.([cfg.var, '_ym']),1));
%     data_hcst.ensmean_ym=squeeze(mean(data_hcst.([cfg.var, '_ym']),1));
    data_assm.ensmean_ym=data_assm_em.([cfg.var, '_ym']);
    data_hcst.ensmean_ym=data_hcst_em.(tmp.ly_str).([cfg.var, '_ym']);
    for loni=1:grid.nlon
        for lati=1:grid.nlat
            [tmp.corr, tmp.corr_p]=corrcoef(data_assm.ensmean_ym(loni,lati,:), ...
                data_hcst.ensmean_ym(loni,lati,:));
            corrval_hcst.assm_hcst_em.(tmp.ly_str).val(loni,lati)=single(tmp.corr(1,2));
            
            tmp.a=data_assm.ensmean_ym(loni,lati,:);
            tmp.b=data_hcst.ensmean_ym(loni,lati,:);
%             if isfinite(tmp.corr(1,2))
%                 tmp.btst=bootstrp(tmp.bootstrpn, 'corr', tmp.a(logical(isfinite(tmp.a).*isfinite(tmp.b))), tmp.b(logical(isfinite(tmp.b).*isfinite(tmp.a))));
%                 corrval_hcst.assm_hcst_em.val_btst(loni,lati,1:tmp.bootstrpn)=tmp.btst;
% %                 corrval_hcst.assm_hcst_em.val_005(loni,lati)=quantile(tmp.btst,0.005);
% %                 corrval_hcst.assm_hcst_em.val_025(loni,lati)=quantile(tmp.btst,0.025);
% %                 corrval_hcst.assm_hcst_em.val_050(loni,lati)=quantile(tmp.btst,0.05);
% %                 corrval_hcst.assm_hcst_em.val_170(loni,lati)=quantile(tmp.btst,0.17);
% %                 corrval_hcst.assm_hcst_em.val_995(loni,lati)=quantile(tmp.btst,0.995);
% %                 corrval_hcst.assm_hcst_em.val_975(loni,lati)=quantile(tmp.btst,0.975);
% %                 corrval_hcst.assm_hcst_em.val_950(loni,lati)=quantile(tmp.btst,0.95);
% %                 corrval_hcst.assm_hcst_em.val_830(loni,lati)=quantile(tmp.btst,0.83);
%             else
%                 corrval_hcst.assm_hcst_em.val_btst(loni,lati,1:tmp.bootstrpn)=NaN;
%             end

            %% internal signal
            [tmp.corr, tmp.corr_p]=corrcoef(data_assm.ensmean_ym(loni,lati,:), ...
                squeeze(data_hcst.ensmean_ym(loni,lati,:)) - ...
                squeeze(data_lens2_em.([cfg.var, '_ym'])(loni,lati,:)));
            corrval_hcst.assm_hcst_em_int.(tmp.ly_str).val(loni,lati)=single(tmp.corr(1,2));
%             corrval_hcst.assm_hcst_em_int.(tmp.ly_str).p(loni,lati)=single(tmp.corr_p(1,2));
        end
    end
    toc;
    
%     pcolor(grid.lon, grid.lat, squeeze(corrval_hcst.assm_hcst_em.val-corrval_hcst.assm_hcst.(tmp.ly_str).val_mean)'); shading flat; colorbar;
%     pcolor(grid.lon, grid.lat,squeeze(corrval_hcst.assm_hcst_em.val-corrval_hcst.assm_hcst.(tmp.ly_str).val_median)'); shading flat; colorbar;
%     caxis([0 0.4])
%     pcolor(grid.lon, grid.lat,squeeze(corrval_hcst.assm_hcst.(tmp.ly_str).val_mean)'); shading flat; colorbar;
%     pcolor(grid.lon, grid.lat,squeeze(corrval_hcst.assm_hcst_em.val)'); shading flat; colorbar;


    %% corr, OBS <-> HCST
    tic;
    for mi2= 1:cfg_hcst.len_mem
        tmp.member2=cfg_hcst.members{mi2};
        corrval_hcst.obs_hcst.(tmp.ly_str).corr_member{1, mi2}=[cfg.obs_name, ' x ', tmp.member2];
        disp(['corr, ', corrval_hcst.obs_hcst.(tmp.ly_str).corr_member{1, mi2}]);
        for loni=1:grid.nlon
            for lati=1:grid.nlat
                [tmp.corr, tmp.corr_p]=corrcoef(data_obs.([cfg.var,'_ym'])(loni,lati,:), ...
                    data_hcst.([cfg.var,'_ym'])(mi2,loni,lati,:), 'Rows', 'complete');
                corrval_hcst.obs_hcst.(tmp.ly_str).val(mi2, loni,lati)=single(tmp.corr(1,2));

                %% internal signal
                [tmp.corr, tmp.corr_p]=corrcoef(data_obs.([cfg.var,'_ym'])(loni,lati,:), ...
                    squeeze(data_hcst.([cfg.var,'_ym'])(mi2,loni,lati,:)) - ...
                    squeeze(data_lens2_em.([cfg.var, '_ym'])(loni,lati,:)), 'Rows', 'complete');
                corrval_hcst.obs_hcst_int.(tmp.ly_str).val(mi2, loni,lati)=single(tmp.corr(1,2));
%                 corrval_hcst.obs_hcst_int.(tmp.ly_str).p(mi2, loni,lati)=single(tmp.corr_p(1,2));
            end
        end
    end
    corrval_hcst.obs_hcst.(tmp.ly_str).val_mean=squeeze(mean(corrval_hcst.obs_hcst.(tmp.ly_str).val,1));
    corrval_hcst.obs_hcst.(tmp.ly_str).val_median=squeeze(median(corrval_hcst.obs_hcst.(tmp.ly_str).val,1));
 
   if strcmp(tmp.varname, 'SSH')
    for mi2= 1:cfg_hcst.len_mem
        tmp.member2=cfg_hcst.members{mi2};
        corrval_hcst.obs_hcst.(tmp.ly_str).corr_member{1, mi2}=[cfg.obs_name, ' x ', tmp.member2];
        disp(['corr, ', corrval_hcst.obs_hcst.(tmp.ly_str).corr_member{1, mi2}]);
        for loni=1:grid.nlon
            for lati=1:grid.nlat
		[tmp.det_data, tmp.trend] = ...
                    Func_0028_detrend_linear_1d(data_obs.([cfg.var,'_ym'])(loni,lati,:), 'omitnan');
                tmp.data2=data_hcst.([cfg.var,'_ym'])(mi2,loni,lati,:);
                tmp.data2(isnan(tmp.det_data));
		[tmp.det_data2, tmp.trend] = ...
                    Func_0028_detrend_linear_1d(tmp.data2, 'omitnan');

                [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
                    tmp.det_data2, 'Rows', 'complete');
               % [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
               %     data_hcst.([cfg.var,'_ym'])(mi2,loni,lati,:), 'Rows', 'complete');
                corrval_hcst.obsdet_hcst.(tmp.ly_str).val(mi2, loni,lati)=single(tmp.corr(1,2));
            end
        end
    end
    corrval_hcst.obsdet_hcst.(tmp.ly_str).val_mean=squeeze(mean(corrval_hcst.obsdet_hcst.(tmp.ly_str).val,1));
    corrval_hcst.obsdet_hcst.(tmp.ly_str).val_median=squeeze(median(corrval_hcst.obsdet_hcst.(tmp.ly_str).val,1));
    end
    
    % ensmean based corr
    for loni=1:grid.nlon
        for lati=1:grid.nlat
            [tmp.corr, tmp.corr_p]=corrcoef(data_obs.([cfg.var,'_ym'])(loni,lati,:), ...
                data_hcst.ensmean_ym(loni,lati,:), 'Rows', 'complete');
            corrval_hcst.obs_hcst_em.(tmp.ly_str).val(loni,lati)=single(tmp.corr(1,2));

            tmp.a=data_obs.([cfg.var,'_ym'])(loni,lati,:);
            tmp.b=data_hcst.ensmean_ym(loni,lati,:);
%             if isfinite(tmp.corr(1,2))
%                 tmp.btst=bootstrp(tmp.bootstrpn, 'corr', tmp.a(logical(isfinite(tmp.a).*isfinite(tmp.b))), tmp.b(logical(isfinite(tmp.b).*isfinite(tmp.a))));
%                 corrval_hcst.obs_hcst_em.val_btst(loni,lati,1:tmp.bootstrpn)=tmp.btst;
% %                 corrval_hcst.obs_hcst_em.val_005(loni,lati)=quantile(tmp.btst,0.005);
% %                 corrval_hcst.obs_hcst_em.val_025(loni,lati)=quantile(tmp.btst,0.025);
% %                 corrval_hcst.obs_hcst_em.val_050(loni,lati)=quantile(tmp.btst,0.05);
% %                 corrval_hcst.obs_hcst_em.val_170(loni,lati)=quantile(tmp.btst,0.17);
% %                 corrval_hcst.obs_hcst_em.val_995(loni,lati)=quantile(tmp.btst,0.995);
% %                 corrval_hcst.obs_hcst_em.val_975(loni,lati)=quantile(tmp.btst,0.975);
% %                 corrval_hcst.obs_hcst_em.val_950(loni,lati)=quantile(tmp.btst,0.95);
% %                 corrval_hcst.obs_hcst_em.val_830(loni,lati)=quantile(tmp.btst,0.83);
%             else
%                 corrval_hcst.obs_hcst_em.val_btst(loni,lati,1:tmp.bootstrpn)=NaN;
%             end

            %% internal signal
            [tmp.corr, tmp.corr_p]=corrcoef(data_obs.([cfg.var,'_ym'])(loni,lati,:), ...
                squeeze(data_hcst.ensmean_ym(loni,lati,:)) - ...
                squeeze(data_lens2_em.([cfg.var, '_ym'])(loni,lati,:)), 'Rows', 'complete');
            corrval_hcst.obs_hcst_em_int.(tmp.ly_str).val(loni,lati)=single(tmp.corr(1,2));
%             corrval_hcst.obs_hcst_em_int.(tmp.ly_str).p(loni,lati)=single(tmp.corr_p(1,2));

        end
    end

    if strcmp(tmp.varname, 'SSH')
    for loni=1:grid.nlon
        for lati=1:grid.nlat
            [tmp.det_data, tmp.trend] = ...
                    Func_0028_detrend_linear_1d(data_obs.([cfg.var,'_ym'])(loni,lati,:), 'omitnan');
            tmp.data2=data_hcst.ensmean_ym(loni,lati,:);
            tmp.data2(isnan(tmp.det_data))=NaN;
            [tmp.det_data2, tmp.trend] = ...
                    Func_0028_detrend_linear_1d(tmp.data2, 'omitnan');
        
            [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
                tmp.det_data2, 'Rows', 'complete');
           % [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
           %     data_hcst.ensmean_ym(loni,lati,:), 'Rows', 'complete');
            corrval_hcst.obsdet_hcst_em.(tmp.ly_str).val(loni,lati)=single(tmp.corr(1,2));

            tmp.a=data_obs.([cfg.var,'_ym'])(loni,lati,:);
            tmp.b=data_hcst.ensmean_ym(loni,lati,:);

        end
    end
    end
    toc;
    
end

data_hcst.([cfg.var, '_2ym']) = data_hcst.([cfg.var, '_2ym'])./2;
data_hcst.([cfg.var, '_3ym']) = data_hcst.([cfg.var, '_3ym'])./3;
data_hcst.([cfg.var, '_4ym']) = data_hcst.([cfg.var, '_4ym'])./4;
data_hcst_em.([cfg.var, '_2ym'])=squeeze(mean(data_hcst.([cfg.var, '_2ym']),1));
data_hcst_em.([cfg.var, '_3ym'])=squeeze(mean(data_hcst.([cfg.var, '_3ym']),1));
data_hcst_em.([cfg.var, '_4ym'])=squeeze(mean(data_hcst.([cfg.var, '_4ym']),1));

%% corr, ASSM <-> HCST
for mvi=2:4
    str_mvi=num2str(mvi);
    disp(['mvy', str_mvi]);
    for mi1= 1:cfg_assm.len_mem
        tmp.member1=cfg_assm.members{mi1};
        for mi2= 1:cfg_hcst.len_mem
            tmp.member2=cfg_hcst.members{mi2};
            corrval_hcst.(['assm_hcst_',str_mvi, 'ym']).corr_member{1, (mi1-1)*cfg_hcst.len_mem+mi2}=[tmp.member1, ' x ', tmp.member2];
            disp([tmp.varname,', mvy:', str_mvi, ', corr, ', corrval_hcst.(['assm_hcst_',str_mvi, 'ym']).corr_member{1, (mi1-1)*cfg_hcst.len_mem+mi2}]);
            for vali=1:length(grid.valid_ind)
                loni=grid.valid_ind_i(vali);
                lati=grid.valid_ind_j(vali);
                 [tmp.corr, tmp.corr_p]=corrcoef(data_assm.([cfg.var,'_', str_mvi, 'ym'])(mi1,loni,lati,:), ...
                    data_hcst.([cfg.var,'_', str_mvi, 'ym'])(mi2,loni,lati,:), 'Rows', 'complete');
                corrval_hcst.(['assm_hcst_',str_mvi, 'ym']).val((mi1-1)*cfg_hcst.len_mem+mi2, vali)=single(tmp.corr(1,2));
    
                %% internal signal
                [tmp.corr, tmp.corr_p]=corrcoef(data_assm.([cfg.var,'_', str_mvi, 'ym'])(mi1,loni,lati,:), ...
                    squeeze(data_hcst.([cfg.var,'_', str_mvi, 'ym'])(mi2,loni,lati,:)) - ...
                    squeeze(data_lens2_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:)), 'Rows', 'complete');
                corrval_hcst.(['assm_hcst_int_', str_mvi, 'ym']).val((mi1-1)*cfg_hcst.len_mem+mi2, vali)=single(tmp.corr(1,2));
            end
        end
    end
    corrval_hcst.(['assm_hcst_',str_mvi, 'ym']).val_mean=squeeze(mean(corrval_hcst.(['assm_hcst_',str_mvi, 'ym']).val,1));
    corrval_hcst.(['assm_hcst_',str_mvi, 'ym']).val_median=squeeze(median(corrval_hcst.(['assm_hcst_',str_mvi, 'ym']).val,1));
end

% ensmean based corr (mvi)
for mvi=2:4
    str_mvi=num2str(mvi);
    disp(['mvy', str_mvi]);
    for loni=1:grid.nlon
        for lati=1:grid.nlat
            [tmp.corr, tmp.corr_p]=corrcoef(data_assm_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:), ...
                data_hcst_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:));
            corrval_hcst.(['assm_hcst_em_',str_mvi,'ym']).val(loni,lati)=single(tmp.corr(1,2));

            tmp.a=data_assm_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:);
            tmp.b=data_hcst_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:);
%             if isfinite(tmp.corr(1,2))
%                 tmp.btst=bootstrp(tmp.bootstrpn, 'corr', tmp.a(logical(isfinite(tmp.a).*isfinite(tmp.b))), tmp.b(logical(isfinite(tmp.b).*isfinite(tmp.a))));
%                 corrval_hcst.(['assm_hcst_em_',str_mvi,'ym']).val_btst(loni,lati,1:tmp.bootstrpn)=tmp.btst;
% %                 corrval_hcst.(['assm_hcst_em_',str_mvi,'ym']).val_005(loni,lati)=quantile(tmp.btst,0.005);
% %                 corrval_hcst.(['assm_hcst_em_',str_mvi,'ym']).val_025(loni,lati)=quantile(tmp.btst,0.025);
% %                 corrval_hcst.(['assm_hcst_em_',str_mvi,'ym']).val_050(loni,lati)=quantile(tmp.btst,0.05);
% %                 corrval_hcst.(['assm_hcst_em_',str_mvi,'ym']).val_170(loni,lati)=quantile(tmp.btst,0.17);
% %                 corrval_hcst.(['assm_hcst_em_',str_mvi,'ym']).val_995(loni,lati)=quantile(tmp.btst,0.995);
% %                 corrval_hcst.(['assm_hcst_em_',str_mvi,'ym']).val_975(loni,lati)=quantile(tmp.btst,0.975);
% %                 corrval_hcst.(['assm_hcst_em_',str_mvi,'ym']).val_950(loni,lati)=quantile(tmp.btst,0.95);
% %                 corrval_hcst.(['assm_hcst_em_',str_mvi,'ym']).val_830(loni,lati)=quantile(tmp.btst,0.83);
%             else
%                 corrval_hcst.(['assm_hcst_em_',str_mvi,'ym']).val_btst(loni,lati,1:tmp.bootstrpn)=NaN;
%             end

            %% internal signal
            [tmp.corr, tmp.corr_p]=corrcoef(data_assm_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:), ...
                squeeze(data_hcst_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:)) - ...
                squeeze(data_lens2_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:)));
            corrval_hcst.(['assm_hcst_em_int_',str_mvi,'ym']).val(loni,lati)=single(tmp.corr(1,2));
        end
    end
end

%% corr, OBS <-> HCST
for mvi=2:4
    str_mvi=num2str(mvi);
    disp(['mvy', str_mvi]);
    for mi2= 1:cfg_hcst.len_mem
        tmp.member2=cfg_hcst.members{mi2};
        corrval_hcst.(['obs_hcst_',str_mvi, 'ym']).corr_member{1, mi2}=[cfg.obs_name, ' x ', tmp.member2];
        disp(['corr, mvi: ', str_mvi, ', ', corrval_hcst.(['obs_hcst_',str_mvi, 'ym']).corr_member{1, mi2}]);
        for loni=1:grid.nlon
            for lati=1:grid.nlat
                [tmp.corr, tmp.corr_p]=corrcoef(data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:), ...
                    data_hcst.([cfg.var,'_', str_mvi, 'ym'])(mi2,loni,lati,:), 'Rows', 'complete');
                corrval_hcst.(['obs_hcst_',str_mvi, 'ym']).val(mi2, loni,lati)=single(tmp.corr(1,2));
                %% internal signal
                [tmp.corr, tmp.corr_p]=corrcoef(data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:), ...
                    squeeze(data_hcst.([cfg.var,'_', str_mvi, 'ym'])(mi2,loni,lati,:)) - ...
                    squeeze(data_lens2_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:)), 'Rows', 'complete');
                corrval_hcst.(['obs_hcst_int_',str_mvi, 'ym']).val(mi2, loni,lati)=single(tmp.corr(1,2));
            end
        end
    end
    corrval_hcst.(['obs_hcst_',str_mvi, 'ym']).val_mean=squeeze(mean(corrval_hcst.(['obs_hcst_',str_mvi, 'ym']).val,1));
    corrval_hcst.(['obs_hcst_',str_mvi, 'ym']).val_median=squeeze(median(corrval_hcst.(['obs_hcst_',str_mvi, 'ym']).val,1));
   
    
    if strcmp(tmp.varname, 'SSH')
    for mi2= 1:cfg_hcst.len_mem
        tmp.member2=cfg_hcst.members{mi2};
        disp(['corr, mvi: ', str_mvi, ', ', corrval_hcst.(['obs_hcst_',str_mvi, 'ym']).corr_member{1, mi2}]);
        for loni=1:grid.nlon
            for lati=1:grid.nlat
                [tmp.det_data, tmp.trend] = ...
                    Func_0028_detrend_linear_1d(data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:), 'omitnan');
                tmp.data2=data_hcst.([cfg.var,'_', str_mvi, 'ym'])(mi2,loni,lati,:);
                tmp.data2(isnan(tmp.det_data))=NaN;
                [tmp.det_data2, tmp.trend] = ...
                    Func_0028_detrend_linear_1d(tmp.data2, 'omitnan');
               % [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
               %     data_hcst.([cfg.var,'_', str_mvi, 'ym'])(mi2,loni,lati,:), 'Rows', 'complete');
                [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
                    tmp.det_data2, 'Rows', 'complete');
                corrval_hcst.(['obsdet_hcst_',str_mvi, 'ym']).val(mi2, loni,lati)=single(tmp.corr(1,2));
            end
        end
    end
    end
 
    % ensmean based corr
    for loni=1:grid.nlon
        for lati=1:grid.nlat
            
            [tmp.corr, tmp.corr_p]=corrcoef(data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:), ...
                data_hcst_em.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:), 'Rows', 'complete');
            corrval_hcst.(['obs_hcst_em_',str_mvi, 'ym']).val(loni,lati)=single(tmp.corr(1,2));
            
           
            tmp.a=data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:);
            tmp.b=data_hcst_em.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:);
%             if isfinite(tmp.corr(1,2))
%                 tmp.btst=bootstrp(tmp.bootstrpn, 'corr', tmp.a(logical(isfinite(tmp.a).*isfinite(tmp.b))), tmp.b(logical(isfinite(tmp.b).*isfinite(tmp.a))));
%                 corrval_hcst.(['obs_hcst_em_',str_mvi, 'ym']).val_btst(loni,lati,1:tmp.bootstrpn)=tmp.btst;
% %                 corrval_hcst.(['obs_hcst_em_',str_mvi, 'ym']).val_005(loni,lati)=quantile(tmp.btst,0.005);
% %                 corrval_hcst.(['obs_hcst_em_',str_mvi, 'ym']).val_025(loni,lati)=quantile(tmp.btst,0.025);
% %                 corrval_hcst.(['obs_hcst_em_',str_mvi, 'ym']).val_050(loni,lati)=quantile(tmp.btst,0.05);
% %                 corrval_hcst.(['obs_hcst_em_',str_mvi, 'ym']).val_170(loni,lati)=quantile(tmp.btst,0.17);
% %                 corrval_hcst.(['obs_hcst_em_',str_mvi, 'ym']).val_995(loni,lati)=quantile(tmp.btst,0.995);
% %                 corrval_hcst.(['obs_hcst_em_',str_mvi, 'ym']).val_975(loni,lati)=quantile(tmp.btst,0.975);
% %                 corrval_hcst.(['obs_hcst_em_',str_mvi, 'ym']).val_950(loni,lati)=quantile(tmp.btst,0.95);
% %                 corrval_hcst.(['obs_hcst_em_',str_mvi, 'ym']).val_830(loni,lati)=quantile(tmp.btst,0.83);
%             else
%                 corrval_hcst.(['obs_hcst_em_',str_mvi, 'ym']).val_btst(loni,lati,1:tmp.bootstrpn)=NaN;
%             end

            %% internal signal
            [tmp.corr, tmp.corr_p]=corrcoef(data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:), ...
                squeeze(data_hcst_em.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:)) - ...
                squeeze(data_lens2_em.([cfg.var, '_', str_mvi, 'ym'])(loni,lati,:)), 'Rows', 'complete');
            corrval_hcst.(['obs_hcst_em_int_',str_mvi, 'ym']).val(loni,lati)=single(tmp.corr(1,2));
%             corrval_hcst.obs_hcst_em_int.(tmp.ly_str).p(loni,lati)=single(tmp.corr_p(1,2));

        end
    end
    if strcmp(tmp.varname, 'SSH')
    for loni=1:grid.nlon
        for lati=1:grid.nlat
            [tmp.det_data, tmp.trend] = ...
                Func_0028_detrend_linear_1d(data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:), 'omitnan');
            tmp.data2=data_hcst_em.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:);
            tmp.data2(isnan(tmp.det_data))=NaN;
            [tmp.det_data2, tmp.trend] = ...
                Func_0028_detrend_linear_1d(tmp.data2, 'omitnan');
             
            [tmp.corr, tmp.corr_p]=corrcoef(tmp.det_data, ...
                tmp.det_data2, 'Rows', 'complete');
           % [tmp.corr, tmp.corr_p]=corrcoef(data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:), ...
           %     data_hcst_em.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:), 'Rows', 'complete');
            corrval_hcst.(['obsdet_hcst_em_',str_mvi, 'ym']).val(loni,lati)=single(tmp.corr(1,2));
            tmp.a=data_obs.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:);
            tmp.b=data_hcst_em.([cfg.var,'_', str_mvi, 'ym'])(loni,lati,:);
        end
    end

    end
end
toc;

% pcolor(grid.lon, grid.lat, squeeze(corrval_lens2.obs_lens2_em.val-corrval_lens2.obs_lens2.val_mean)'); shading flat; colorbar;
% pcolor(grid.lon, grid.lat,squeeze(corrval_lens2.obs_lens2_em.val-corrval_lens2.obs_lens2.val_median)'); shading flat; colorbar;
% % caxis([-0.3 0.7])
% pcolor(grid.lon, grid.lat,squeeze(corrval_lens2.obs_lens2.val_mean)'); shading flat; colorbar;
% pcolor(grid.lon, grid.lat,squeeze(corrval_lens2.obs_lens2_em.val)'); shading flat; colorbar;


%% save matfile
matfilename=[dirs.saveroot, '/corr_raw/', 'corr_hcst_',cfg.var,'_v', num2str(cfg.vlayer_1st), '_v', num2str(max(cfg.vlayer)), '.mat'];
mkdir([dirs.saveroot, '/corr_raw']);
% save(matfilename, 'cfg_assm', 'cfg_hcst', 'cfg_lens2', '-struct', 'corrval_hcst', 'grid', '-v7.3');
corrval_hcst.grid=grid;
corrval_hcst.data=data_hcst_em;
save(matfilename, '-struct', 'corrval_hcst', '-v7.3');

clear data_assm data_lens2 data_assm_em data_lens2_em data_hcst data_hcst_em corrval_hcst corrval_assm corrval_lens2 grid
end
disp('postprocessing complete');


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
        case 'TREFHT'
%             obsname_simple='HadCRUT5';
%             obsname_simple='ERA5';
            obsname_simple='SMYLE';
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
        case 'TREFHT'
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
        case 'TREFHT'
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

