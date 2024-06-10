% %  Created 19-Feb-2024 by Yong-Yub Kim

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
dirs.saveroot=['/mnt/lustre/proj/kimyy/Model/CESM2/ESP', '/statistics'];


dirs.deg5_root='/mnt/lustre/proj/earth.system.predictability/regrid_5deg';
dirs.assm_deg5_root=[dirs.deg5_root, '/ASSM_EXP_5deg/archive'];
dirs.lens2_deg5_root=[dirs.deg5_root, '/LENS2_5deg/archive'];



cfg.iyears=1960:2020;
cfg.months=1:12;
cfg.scenname='HIST';
cfg.gridname='f09_g17';
cfg.proj_year=5;


cfg.vars={'TS', 'SST', 'PRECT', 'TWS', 'GPP', 'FAREA_BURNED', 'photoC_TOT_zint_100m', ...
    'COL_FIRE_CLOSS', 'SOILWATER_10CM', 'TLAI', 'SSH', 'PSL', 'AEROD_v'};
% cfg.vars={'photoC_TOT_zint_100m'};
% cfg.vars={'COL_FIRE_CLOSS', 'SOILWATER_10CM', 'TLAI', 'SSH', 'PSL', 'AEROD_v'};
% cfg.vars={'SSH', 'PSL', 'AEROD_v'};
% cfg.vars={'TS'};
% cfg.vars={'TLAI'};
cfg.vars={'DpCO2_ALT_CO2'};
cfg.vars={'DIC', 'FG_CO2', 'DIC_ALT_CO2', 'DpCO2_ALT_CO2', 'TEMP', 'SALT', 'SSH', 'DpCO2'};
% cfg.vars={'FG_CO2', 'DpCO2_ALT_CO2', 'DpCO2', 'SSH', 'DIC_ALT_CO2',  'TEMP', 'SALT'};
cfg.vars={'DIC_ALT_CO2',  'DIC', 'TEMP', 'SALT'};
cfg.vars={'SALT', 'TEMP'};

% cfg.vars={'FG_CO2', 'DpCO2_ALT_CO2', 'DpCO2', 'SSH'};
% cfg.vars={'FG_ALT_CO2'}; % nonexist
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
dirs.grid_root=['/mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive_transfer/', cfg.comp];

%% set grid
 tmp.gridname = [dirs.grid_root, tmp.fs, 'grid.nc'];
    tmp.maskname = '/mnt/lustre/proj/kimyy/Model/CESM2/ESP/grids/ocn/RECCAP2_region_masks_all_v20210412_POP2_grid.nc';
    
    % grid.region_mask=ncread(tmp.gridname, 'REGION_MASK'); 
    % grid.ocean_mask=NaN(size(grid.region_mask));
    % grid.ocean_mask(grid.region_mask>0)=1;
    % grid.tarea = ncread(tmp.gridname, 'TAREA');
    
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



%% read climate indices
       %% Initialize variables.
        filename = '/mnt/lustre/proj/kimyy/Observation/AMO/amon.sm.long.data.txt';
        startRow = 2;
        endRow = 169;

        formatSpec = '%5f%9f%9f%9f%9f%9f%9f%9f%9f%9f%9f%9f%f%[^\n\r]';
        fileID = fopen(filename,'r');
        
        dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        
        fclose(fileID);
                
        amon = [dataArray{1:end-1}];
        clearvars filename startRow endRow formatSpec fileID dataArray ans;
        
        amon(amon==-99.9900000000000)=NaN;
        clim_year=amon(:,1);
        clim_index=amon(:,2:13);

        cfg.clim_ys = min(cfg.iyears);
        cfg.clim_ye = max(cfg.iyears);
   
        clim_year2=clim_year(find(clim_year==cfg.clim_ys): find(clim_year==cfg.clim_ye));
        clim_index2=clim_index(find(clim_year==cfg.clim_ys): find(clim_year==cfg.clim_ye),:);

        clim_index_ym=mean(clim_index2,2);
        clim_index_jfm=mean(clim_index2(:,1:3),2);
        for titi=1:cfg.len_t_y-1
            clim_index_djf(titi)=clim_index2(titi,12)+clim_index2(titi+1,1)+clim_index2(titi+1,2);
        end
        clim_index_djf(cfg.len_t_y)=NaN;

%% get ASSM members
cfg_assm.list = dir( [dirs.assm_deg5_root, '/', '*BHISTsmbb*' ]); 
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

%% get LENS2 members
cfg_lens2.list = dir( [dirs.lens2_deg5_root, '/', '*BHISTsmbb*' ]); 
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
%             ncid = netcdf.open(tmp.fname, 'NOWRITE');
%             tmpvarid = netcdf.inqVarID(ncid, cfg.var);
%             tmp.data=netcdf.getVar(ncid,tmpvarid);
%             tmp.dimids=size(tmp.data);
%             netcdf.close(ncid);
        end
        if length(tmp.dimids)>3
%              tic;
%              ncid = netcdf.open(tmp.fname, 'NOWRITE');
%              tmpvarid=netcdf.inqVarID(ncid,cfg.var);
%              tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 cfg.vlayer_1st-1 ind_start-1], [grid.nlon grid.nlat cfg.vlayer_cnt ind_count]));
%              tmp.dd(abs(tmp.dd)>1e30)=NaN;
%              data_assm.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd; %depth averaged value
%              toc;
            
%              tic;
             tmp.fname2='/proj/kimyy/tmp/test_mat.nc';
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
%              toc;
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
data_assm.([cfg.var, '_jfm'])=squeeze(mean(tmp.reshp(:,:,:,1:3,:),4));
data_assm.([cfg.var, '_djf'])=NaN(size(data_assm.([cfg.var, '_jfm'])));
for titi=1:cfg.len_t_y-1
    data_assm.([cfg.var, '_djf'])(:,:,:,titi)=(tmp.reshp(:,:,:,12,titi)+tmp.reshp(:,:,:,1,titi+1)+tmp.reshp(:,:,:,2,titi+1))/3;
end

data_assm_em.([cfg.var, '_ym'])=squeeze(mean(data_assm.([cfg.var, '_ym']),1));
data_assm_em.([cfg.var, '_jfm'])=squeeze(mean(data_assm.([cfg.var, '_jfm']),1));
data_assm_em.([cfg.var, '_djf'])=squeeze(mean(data_assm.([cfg.var, '_djf']),1));
toc;

% % % %% read LENS2 members
% % % tic;
% % % data_lens2.([cfg.var])=NaN(cfg_lens2.len_mem, grid.nlon, grid.nlat, cfg.len_t_y*12); % initialization
% % % 
% % % for mi= 1:cfg_lens2.len_mem
% % %     ty=1;
% % %     tmp.member=cfg_lens2.members{mi};
% % %     disp(['lens2, ', tmp.member]);
% % %     while ty<=cfg.len_t_y
% % %         tmp.year=cfg.iyears(ty);
% % %         tmp.scen = f_scen(tmp.year);
% % %         tmp.casenm = ['b.e21.', tmp.scen, '.f09_g17.',tmp.member];
% % %         tmp.fdir=[dirs.lens2_root, '/', tmp.casenm, '/' cfg.comp, '/proc/tseries/month_1'];
% % %         tmp.flist=dir( [tmp.fdir, '/','*.',cfg.var, '.*'] );
% % %         for kk = 1: length (tmp.flist)
% % %             tmp.fname_in = tmp.flist(kk).name;
% % %             tmp.fname_split = strsplit(tmp.fname_in, {'.'} );
% % %             tmp.fname_period = tmp.fname_split{10};
% % %             fyear_str   = strsplit( tmp.fname_period, '-' );
% % %             fyear_start = str2num( fyear_str{1}(1:4) );
% % %             fyear_end   = str2num( fyear_str{2}(1:4) );
% % %             if( tmp.year >= fyear_start && tmp.year <= fyear_end )
% % % %                 disp(tmp.fname_period);
% % %                 flag_file_in = true;            break;
% % %             end
% % %         end
% % %         tmp.fname=[tmp.fdir, '/', tmp.casenm, cfg.obs_fname_module, cfg.var, '.', tmp.fname_period, '.nc'];
% % %         ind_start=(tmp.year - fyear_start)*12+1;
% % %         if fyear_end<=max(cfg.iyears)
% % %             ind_count=(fyear_end-tmp.year+1)*12;
% % %         else
% % %             ind_count=(max(cfg.iyears)-tmp.year+1)*12;
% % %         end
% % %         
% % %         
% % %         if length(tmp.dimids)>3
% % %              tmp.fname2='/proj/kimyy/tmp/test_mat.nc';
% % %              system(['cdo -O -w ', cfg.vm_str, ' -sellevel,', cfg.vl_z1_str,'/',cfg.vl_z2_str, ' ', tmp.fname, ' ', tmp.fname2]);
% % %              try
% % %                  ncid = netcdf.open(tmp.fname2, 'NOWRITE');
% % %                  tmpvarid = netcdf.inqVarID(ncid, cfg.var);
% % %                  tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 ind_start-1], [grid.nlon grid.nlat ind_count]));
% % %                  tmp.dd(abs(tmp.dd)>1e30)=NaN;
% % %                  data_lens2.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd; %depth averaged value
% % %              catch
% % %                  pause(60);
% % %                  system(['rm -f ', tmp.fname2]);
% % %                  system(['cdo -O -w ', cfg.vm_str, ' -sellevel,', cfg.vl_z1_str,'/',cfg.vl_z2_str, ' ', tmp.fname, ' ', tmp.fname2]);
% % %                  ncid = netcdf.open(tmp.fname2, 'NOWRITE');
% % %                  tmpvarid = netcdf.inqVarID(ncid, cfg.var);
% % %                  tmp.dd=squeeze(netcdf.getVar(ncid,tmpvarid, [0 0 ind_start-1], [grid.nlon grid.nlat ind_count]));
% % %                  tmp.dd(abs(tmp.dd)>1e30)=NaN;
% % %                  data_lens2.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd; %depth averaged value
% % %              end
% % % 
% % %         else
% % %             ncid = netcdf.open(tmp.fname, 'NOWRITE');
% % %             tmpvarid = netcdf.inqVarID(ncid, cfg.var);
% % %             tmp.dd=  netcdf.getVar(ncid,tmpvarid, [0 0 ind_start-1], [grid.nlon grid.nlat ind_count]);
% % %             tmp.dd(abs(tmp.dd)>1e30)=NaN;
% % %             data_lens2.([cfg.var])(mi,:,:,(ty-1)*12+1:(ty-1)*12+ind_count) = tmp.dd;
% % %         end
% % %         netcdf.close(ncid);
% % % 
% % %         ty=ty+ind_count/12;
% % %     end
% % % end
% % % tmp.reshp=reshape(data_lens2.([cfg.var]), [cfg_lens2.len_mem, grid.nlon, grid.nlat, 12, cfg.len_t_y]);
% % % data_lens2.([cfg.var, '_ym'])=squeeze(mean(tmp.reshp,4));
% % % data_lens2.([cfg.var, '_jfm'])=squeeze(mean(tmp.reshp(:,:,:,1:3,:),4));
% % % data_lens2.([cfg.var, '_djf'])=NaN(size(data_lens2.([cfg.var, '_jfm'])));
% % % for titi=1:cfg.len_t_y-1
% % %     data_lens2.([cfg.var, '_djf'])(:,:,:,titi)=(tmp.reshp(:,:,:,12,titi)+tmp.reshp(:,:,:,1,titi+1)+tmp.reshp(:,:,:,2,titi+1))/3;
% % % end
% % % 
% % % data_lens2_em.([cfg.var, '_ym'])=squeeze(mean(data_lens2.([cfg.var, '_ym']),1));
% % % data_lens2_em.([cfg.var, '_jfm'])=squeeze(mean(data_lens2.([cfg.var, '_jfm']),1));
% % % data_lens2_em.([cfg.var, '_djf'])=squeeze(mean(data_lens2.([cfg.var, '_djf']),1));
% % % toc;


%% corr, OBS <-> ASSM
tic;
% for mi2= 1:cfg_assm.len_mem
%     tmp.member2=cfg_assm.members{mi2};
%     corrval.obs_assm.corr_member{1, mi2}=[cfg.obs_name, ' x ', tmp.member2];
%     disp(['corr, ', corrval.obs_assm.corr_member{1, mi2}]);
%     for loni=1:grid.nlon
%         for lati=1:grid.nlat
%             [tmp.corr, tmp.corr_p]=corrcoef(data_obs.([cfg.var,'_ym'])(loni,lati,:), ...
%                 data_assm.([cfg.var,'_ym'])(mi2,loni,lati,:), 'Rows', 'complete');
%             corrval.obs_assm.val(mi2, loni,lati)=tmp.corr(1,2);
%             corrval.obs_assm.p(mi2, loni,lati)=tmp.corr_p(1,2);
%         end
%     end
% end
% corrval.obs_assm.val_mean=squeeze(mean(corrval.obs_assm.val,1));
% corrval.obs_assm.val_median=squeeze(median(corrval.obs_assm.val,1));

% ensmean based corr
for loni=1:grid.nlon
    for lati=1:grid.nlat
        [tmp.corr, tmp.corr_p]=corrcoef(clim_index_ym, ...
            data_assm_em.([cfg.var, '_ym'])(loni,lati,:), 'Rows', 'complete');
        corrval.obs_assm_em_ym.val(loni,lati)=tmp.corr(1,2);
        corrval.obs_assm_em_ym.p(loni,lati)=tmp.corr_p(1,2);
    end
end

for loni=1:grid.nlon
    for lati=1:grid.nlat
        [tmp.corr, tmp.corr_p]=corrcoef(clim_index_jfm, ...
            data_assm_em.([cfg.var, '_jfm'])(loni,lati,:), 'Rows', 'complete');
        corrval.obs_assm_em_jfm.val(loni,lati)=tmp.corr(1,2);
        corrval.obs_assm_em_jfm.p(loni,lati)=tmp.corr_p(1,2);
    end
end

for loni=1:grid.nlon
    for lati=1:grid.nlat
        [tmp.corr, tmp.corr_p]=corrcoef(clim_index_djf, ...
            data_assm_em.([cfg.var, '_djf'])(loni,lati,:), 'Rows', 'complete');
        corrval.obs_assm_em_djf.val(loni,lati)=tmp.corr(1,2);
        corrval.obs_assm_em_djf.p(loni,lati)=tmp.corr_p(1,2);
    end
end
toc;

% pcolor(grid.lon, grid.lat, squeeze(corrval.obs_assm_em.val-corrval.obs_assm.val_mean)'); shading flat; colorbar;
% pcolor(grid.lon, grid.lat,squeeze(corrval.obs_assm_em.val-corrval.obs_assm.val_median)'); shading flat; colorbar;
% % caxis([-0.3 0.7])
% pcolor(grid.lon, grid.lat,squeeze(corrval.obs_assm.val_mean)'); shading flat; colorbar;
% pcolor(grid.lon, grid.lat,squeeze(corrval.obs_assm_em.val)'); shading flat; colorbar;




% % % %% corr, OBS <-> LENS2
% % % tic;
% % % % for mi2= 1:cfg_lens2.len_mem
% % % %     tmp.member2=cfg_lens2.members{mi2};
% % % %     corrval.obs_lens2.corr_member{1, mi2}=[cfg.obs_name, ' x ', tmp.member2];
% % % %     disp(['corr, ', corrval.obs_lens2.corr_member{1, mi2}]);
% % % %     for loni=1:grid.nlon
% % % %         for lati=1:grid.nlat
% % % %             [tmp.corr, tmp.corr_p]=corrcoef(data_obs.([cfg.var,'_ym'])(loni,lati,:), ...
% % % %                 data_lens2.([cfg.var,'_ym'])(mi2,loni,lati,:), 'Rows', 'complete');
% % % %             corrval.obs_lens2.val(mi2, loni,lati)=tmp.corr(1,2);
% % % %             corrval.obs_lens2.p(mi2, loni,lati)=tmp.corr_p(1,2);
% % % %         end
% % % %     end
% % % % end
% % % % corrval.obs_lens2.val_mean=squeeze(mean(corrval.obs_lens2.val,1));
% % % % corrval.obs_lens2.val_median=squeeze(median(corrval.obs_lens2.val,1));
% % % 
% % % % ensmean based corr
% % % for loni=1:grid.nlon
% % %     for lati=1:grid.nlat
% % %         [tmp.corr, tmp.corr_p]=corrcoef(clim_index_ym, ...
% % %             data_lens2_em.([cfg.var, '_ym'])(loni,lati,:), 'Rows', 'complete');
% % %         corrval.obs_lens2_em_ym.val(loni,lati)=tmp.corr(1,2);
% % %         corrval.obs_lens2_em_ym.p(loni,lati)=tmp.corr_p(1,2);
% % %     end
% % % end
% % % 
% % % for loni=1:grid.nlon
% % %     for lati=1:grid.nlat
% % %         [tmp.corr, tmp.corr_p]=corrcoef(clim_index_jfm, ...
% % %             data_lens2_em.([cfg.var, '_jfm'])(loni,lati,:), 'Rows', 'complete');
% % %         corrval.obs_lens2_em_jfm.val(loni,lati)=tmp.corr(1,2);
% % %         corrval.obs_lens2_em_jfm.p(loni,lati)=tmp.corr_p(1,2);
% % %     end
% % % end
% % % 
% % % for loni=1:grid.nlon
% % %     for lati=1:grid.nlat
% % %         [tmp.corr, tmp.corr_p]=corrcoef(clim_index_djf, ...
% % %             data_lens2_em.([cfg.var, '_djf'])(loni,lati,:), 'Rows', 'complete');
% % %         corrval.obs_lens2_em_djf.val(loni,lati)=tmp.corr(1,2);
% % %         corrval.obs_lens2_em_djf.p(loni,lati)=tmp.corr_p(1,2);
% % %     end
% % % end
% % % toc;

% pcolor(grid.lon, grid.lat, squeeze(corrval.obs_lens2_em.val-corrval.obs_lens2.val_mean)'); shading flat; colorbar;
% pcolor(grid.lon, grid.lat,squeeze(corrval.obs_lens2_em.val-corrval.obs_lens2.val_median)'); shading flat; colorbar;
% % caxis([-0.3 0.7])
% pcolor(grid.lon, grid.lat,squeeze(corrval.obs_lens2.val_mean)'); shading flat; colorbar;
% pcolor(grid.lon, grid.lat,squeeze(corrval.obs_lens2_em.val)'); shading flat; colorbar;


%% save matfile
matfilename=[dirs.saveroot, '/corr/', 'corr_all_',cfg.var,'_v',num2str(cfg.vlayer_1st), '_v', num2str(cfg.vlayer_cnt), '.mat'];
mkdir([dirs.saveroot, '/corr']);
save(matfilename, 'cfg_assm', 'cfg_lens2', 'corrval', 'grid', 'data_assm_em', 'clim_index_ym', 'clim_index_jfm', 'clim_index_djf');

clear data_assm data_lens2 data_assm_em data_lens2_em

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

