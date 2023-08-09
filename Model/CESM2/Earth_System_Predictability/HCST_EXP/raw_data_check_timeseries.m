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






cfg.var ='photoC_TOT_zint';

cfg.obs_iyears=1970:1974;
cfg.obs_name=f_obs_name(cfg.var);
cfg.obs_fname_mid=f_obs_name_mid(cfg.var);
cfg.obs_varname=f_obs_varname(cfg.var);
cfg.comp=Func_0025_CESM2_cmpname_var(cfg.var);
cfg.obs_fname_module=f_obs_fname_module(cfg.comp);

dirs.hcstroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/HCST_EXP/archive/', cfg.comp, '/', cfg.var];
dirs.obsroot=['/Volumes/kyy_raid/kimyy/Observation/', cfg.obs_name, '/monthly_reg_', cfg.obs_fname_module(2:4)];
dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', cfg.var];
dirs.lens2root=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/LENS2/archive_analysis/', cfg.comp, '/', cfg.var];
dirs.assmroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/ASSM_EXP/archive_analysis/', cfg.comp, '/', cfg.var];


tmp.gridname = [dirs.hcstroot, tmp.fs, '../grid.nc'];
    tmp.maskname = '/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/grids/ocn/RECCAP2_region_masks_all_v20210412_POP2_grid.nc';

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

tmp.varname=cfg.var;

cfg.casename_m=['ens_all'];
cfg.gnm='f09_g17';

dirs.datadir= [dirs.hcstroot, filesep, cfg.casename_m, filesep];
dirs.assmdir= [dirs.assmroot, filesep, cfg.casename_m, filesep];
dirs.lens2dir= [dirs.lens2root, filesep, 'ens_smbb', filesep];

cfg.iyears=1970;
iyear=1970;
tmp.iyear_str=num2str(iyear, '%04i');
tmp.iind=iyear-min(cfg.iyears)+1;
cfg.casename=[cfg.casename_m, '_i', tmp.iyear_str];
for mon=1:60
    tmp.mon=mon;
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
end


tmp.abc=tmp.sdata_assm(:,:,1);
tmp.abc(tmp.abc>10e35)=NaN;
pcolor(tmp.abc'); shading flat; colorbar; 
lonx=195;
laty=360;
plot(squeeze(tmp.sdata(lonx,laty,:)), 'linewidth', 2)
hold on
% plot(squeeze(tmp.sdata_obs(lonx,laty,:)), 'linewidth', 2)
plot(squeeze(tmp.sdata_assm(lonx,laty,:)), 'linewidth', 2)
plot(squeeze(tmp.sdata_lens2(lonx,laty,:)), 'linewidth', 2)
plot(squeeze(tmp.sdata(lonx,laty,:))-squeeze(tmp.sdata_lens2(lonx,laty,:)), 'linewidth', 2)

% legend({'hcst', 'obs', 'assm', 'lens2'});
legend({'hcst', 'assm', 'lens2', 'hc-le'});
hold off
set(gca,'fontsize', 20)

for i=1:12
    tmp.cor_s(i)=corr(squeeze(tmp.sdata_assm(lonx, laty, (0:12:48)+i)), squeeze((tmp.sdata(lonx, laty, (0:12:48)+i)-tmp.sdata_lens2(lonx, laty, (0:12:48)+i))));
end
bar(1:12,tmp.cor_s)

for i=1:12
    tmp.cor_s(i)=corr(squeeze(tmp.sdata_assm(lonx, laty, (0:12:48)+i)), squeeze(tmp.sdata_lens2(lonx, laty, (0:12:48)+i)));
end
bar(1:12,tmp.cor_s)

for i=1:12
    tmp.cor_s(i)=corr(squeeze(tmp.sdata_assm(lonx, laty, (0:12:48)+i)), squeeze(tmp.sdata(lonx, laty, (0:12:48)+i)));
end
bar(1:12,tmp.cor_s)



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
