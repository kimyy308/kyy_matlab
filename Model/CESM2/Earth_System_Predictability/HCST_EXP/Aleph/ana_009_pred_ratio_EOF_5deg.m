% %  Created 08-Jan-2022 by Yong-Yub Kim

clc; clear all; close all;

%% set path
tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

%% model configuration
dirs.root='/mnt/lustre/proj/earth.system.predictability/regrid_5deg';
dirs.assm_root=[dirs.root, '/ASSM_EXP_5deg/archive'];
dirs.hcst_root=[dirs.root, '/HCST_EXP_5deg/archive'];
dirs.lens2_root=[dirs.root, '/LENS2_5deg/archive'];

dirs.saveroot=[dirs.root, '/statistics'];


cfg.iyears=1960:2016;
% cfg.iyears=1965:2020;

cfg.months=1:12;
cfg.scenname='HIST';
cfg.gridname='f09_g17';
cfg.proj_year=5;


cfg.vars={'TS', 'SST', 'PRECT', 'TWS', 'GPP', 'FAREA_BURNED', 'photoC_TOT_zint_100m'};
% cfg.vars={'photoC_TOT_zint_100m'};
cfg.vars={'COL_FIRE_CLOSS', 'SOILWATER_10CM', 'TLAI', 'SSH', 'PSL', 'AEROD_v'};
cfg.vars={'SSH', 'PSL', 'AEROD_v'};
cfg.vars={'TS'};
for vari=1:length(cfg.vars)
    cfg.var=cfg.vars{vari};
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

%% set grid
tmp.gridname=['/mnt/lustre/proj/earth.system.predictability/regrid_5deg/', ...
    'ASSM_EXP_5deg/archive/b.e21.BHISTsmbb.f09_g17.assm.en4.2_ba-10p1/lnd/proc/tseries/month_1/', ...
    'b.e21.BHISTsmbb.f09_g17.assm.en4.2_ba-10p1.r72x36.clm2.h0.TWS.196101-196512.nc'];
grid.lon=ncread(tmp.gridname, 'lon');
grid.lat=ncread(tmp.gridname, 'lat');
[grid.tlat, grid.tlong]=meshgrid(grid.lat, grid.lon);
grid.nlon=length(grid.lon);
grid.nlat=length(grid.lat);

grid.latdist(1:72)=stdist([0 5], [0 0], wgs84Ellipsoid("m"));
for li=1:36
    grid.londist(li)=stdist([-92.5+li*5 -92.5+li*5], [0 5], wgs84Ellipsoid("km"));
end 
grid.area=grid.latdist' * grid.londist;

%% get ASSM members
cfg_assm.list = dir( [dirs.assm_root, '/', '*BHISTsmbb*' ]); % just get ensemble members
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
data_assm.([cfg.var])=NaN(cfg_assm.len_mem, grid.nlon, grid.nlat, cfg.len_t_y*12); % initialization

for mi= 1:cfg_assm.len_mem
    ty=1;
    tmp.member=cfg_assm.members{mi};
    disp(['assm, ', tmp.member]);
    while ty<=cfg.len_t_y+4
        tmp.year=cfg.iyears(ty);
        tmp.scen = f_scen(tmp.year);
        tmp.casenm = ['b.e21.', tmp.scen, '.f09_g17.assm.',tmp.member];
        tmp.fdir=[dirs.assm_root, '/', tmp.casenm, '/' cfg.comp, '/proc/tseries/month_1'];
        tmp.flist=dir( [tmp.fdir, '/','*.',cfg.var, '.*'] );
        for kk = 1: length (tmp.flist)
            tmp.fname_in = tmp.flist(kk).name;
            tmp.fname_split = strsplit(tmp.fname_in, {'.'} );
            tmp.fname_period = tmp.fname_split{12};
            fyear_str   = strsplit( tmp.fname_period, '-' );
            fyear_start = str2num( fyear_str{1}(1:4) );
            fyear_end   = str2num( fyear_str{2}(1:4) );
            if( tmp.year >= fyear_start && tmp.year <= fyear_end )
%                 disp(tmp.fname_period);
                flag_file_in = true;            break;
            end
        end
        tmp.fname=[tmp.fdir, '/', tmp.casenm, '.r72x36', cfg.obs_fname_module, cfg.var, '.', tmp.fname_period, '.nc'];
        ind_start=(tmp.year - fyear_start)*12+1;
        if fyear_end<=max(cfg.iyears)
            ind_count=(fyear_end-tmp.year+1)*12;
        else
            ind_count=(max(cfg.iyears)-tmp.year+1)*12;
        end
        
        ncid = netcdf.open(tmp.fname, 'NOWRITE');
        tmpvarid = netcdf.inqVarID(ncid, cfg.var);
        [tmp.fvarname,tmp.xtype,tmp.dimids,tmp.natts]= netcdf.inqVar(ncid,tmpvarid);
        if length(tmp.dimids)>3
             tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 ind_start-1], [grid.nlon grid.nlat cfg.vlayer_cnt ind_count]));
%                      tmp.dd=tmp.dd.*grid.mask_ocn;
             tmp.dd(abs(tmp.dd)>1e30)=NaN;
             ddd=mean(tmp.dd,3,'omitnan'); % depth mean
             data_assm.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = ddd; %depth averaged value
        else
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

% ensemble mean
data_assm.ensmean=squeeze(mean(data_assm.([cfg.var]),1));
data_assm.ensmean_ym=squeeze(mean(data_assm.([cfg.var, '_ym']),1));
data_assm.ensstd=squeeze(std(data_assm.([cfg.var]),1));
data_assm.ensstd_ym=squeeze(std(data_assm.([cfg.var, '_ym']),1));

%% read LENS2 members
data_lens2.([cfg.var])=NaN(cfg_lens2.len_mem, grid.nlon, grid.nlat, cfg.len_t_y*12); % initialization

for mi= 1:cfg_lens2.len_mem
    ty=1;
    tmp.member=cfg_lens2.members{mi};
    disp(['lens2, ', tmp.member]);
    while ty<=cfg.len_t_y+4
        tmp.year=cfg.iyears(ty);
        tmp.scen = f_scen(tmp.year);
        tmp.casenm = ['b.e21.', tmp.scen, '.f09_g17.',tmp.member];
        tmp.fdir=[dirs.lens2_root, '/', tmp.casenm, '/' cfg.comp, '/proc/tseries/month_1'];
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
        tmp.fname=[tmp.fdir, '/', tmp.casenm, '.r72x36', cfg.obs_fname_module, cfg.var, '.', tmp.fname_period, '.nc'];
        ind_start=(tmp.year - fyear_start)*12+1;
        if fyear_end<=max(cfg.iyears)
            ind_count=(fyear_end-tmp.year+1)*12;
        else
            ind_count=(max(cfg.iyears)-tmp.year+1)*12;
        end
        
        ncid = netcdf.open(tmp.fname, 'NOWRITE');
        tmpvarid = netcdf.inqVarID(ncid, cfg.var);
        [tmp.fvarname,tmp.xtype,tmp.dimids,tmp.natts]= netcdf.inqVar(ncid,tmpvarid);
        if length(tmp.dimids)>3
             tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 ind_start-1], [grid.nlon grid.nlat cfg.vlayer_cnt ind_count]));
%                      tmp.dd=tmp.dd.*grid.mask_ocn;
             tmp.dd(abs(tmp.dd)>1e30)=NaN;
             ddd=mean(tmp.dd,3,'omitnan'); % depth mean
             data_lens2.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = ddd; %depth averaged value
        else
            tmp.dd=  netcdf.getVar(ncid,tmpvarid, [0 0 ind_start-1], [grid.nlon grid.nlat ind_count]);
            tmp.dd(abs(tmp.dd)>1e30)=NaN;
            data_lens2.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd;
        end
        netcdf.close(ncid);

        ty=ty+ind_count/12;
    end
end
% yearly mean
tmp.reshp=reshape(data_lens2.([cfg.var]), [cfg_lens2.len_mem, grid.nlon, grid.nlat, 12, cfg.len_t_y]);
data_lens2.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,4));

% ensemble mean
data_lens2.ensmean=squeeze(mean(data_lens2.([cfg.var]),1));
data_lens2.ensmean_ym=squeeze(mean(data_lens2.([cfg.var, '_ym']),1));
data_lens2.ensstd=squeeze(std(data_lens2.([cfg.var]),1));
data_lens2.ensstd_ym=squeeze(std(data_lens2.([cfg.var, '_ym']),1));

% for lm=1:60 %1:5
%     tmp.lm_str=['lm',num2str(ly, '%02i')];
    %% read HCST members
    data_hcst.([cfg.var])=NaN(cfg_hcst.len_mem, grid.nlon, grid.nlat, cfg.len_t_y*12); % initialization
    
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
           
%             ind_start=(ly - 1)*12+1;
%             ind_count=12;
            ind_start=1;
            ind_count=60;
            
            ncid = netcdf.open(tmp.fname, 'NOWRITE');
            tmpvarid = netcdf.inqVarID(ncid, cfg.var);
            [tmp.fvarname,tmp.xtype,tmp.dimids,tmp.natts]= netcdf.inqVar(ncid,tmpvarid);
            if length(tmp.dimids)>3
                 tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 ind_start-1], [grid.nlon grid.nlat cfg.vlayer_cnt ind_count]));
    %                      tmp.dd=tmp.dd.*grid.mask_ocn;
                 tmp.dd(abs(tmp.dd)>1e30)=NaN;
                 ddd=mean(tmp.dd,3,'omitnan'); % depth mean
                 data_hcst.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = ddd; %depth averaged value
            else
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
    
    % ensemble mean
    data_hcst.ensmean=squeeze(mean(data_hcst.([cfg.var]),1));
    data_hcst.ensmean_ym=squeeze(mean(data_hcst.([cfg.var, '_ym']),1));
    data_hcst.ensstd=squeeze(std(data_hcst.([cfg.var]),1));
    data_hcst.ensstd_ym=squeeze(std(data_hcst.([cfg.var, '_ym']),1));    
   

% end






%% save matfile
matfilename=[dirs.saveroot, '/corr/', 'corr_5deg_all_',cfg.var,'.mat'];
mkdir([dirs.saveroot, '/corr']);
save(matfilename, 'cfg_assm', 'cfg_hcst', 'cfg_lens2', 'corrval', 'grid');


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

