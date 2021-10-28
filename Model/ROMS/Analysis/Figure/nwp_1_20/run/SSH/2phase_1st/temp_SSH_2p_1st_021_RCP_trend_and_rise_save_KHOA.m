% %         calculation of RCP SSH trends (2006-2100) and rise (1986-2005 ~ 2081-2100; raw data)
% %  Updated 08-Oct-2021 by Yong-Yub Kim

close all; clear all; clc;

%  % % configuration of RCM
RCM_info.histname={'test53', 'test54', 'test55', 'test56'};
RCM_info.rcp26name={'test61', 'test62', 'test63', 'test64'};
RCM_info.rcp45name={'test57', 'test58', 'test59', 'test60'};
RCM_info.rcp85name={'test65', 'test66', 'test67', 'test68'};

RCM_info.abbs = {'RCM-IPSL-L', 'RCM-IPSL-M', 'RCM-Nor', 'RCM-MPI'};
RCM_info.model = 'nwp_1_20';
% % RCM_info.dataroot = '/data1/RCM/CMIP6/';
% RCM_info.dataroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
%     'ROMS', filesep, 'nwp_1_20', filesep];    % D:\Data\Model\ROMS\nwp_1_20
% RCM_info.saveroot = '/data1/RCM/CMIP6/';
RCM_info.dataroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep];   % D:\Data\Model\ROMS\nwp_1_20
RCM_info.phase = 'run';
% RCM_info.region = {'NWP', 'AKP4'};
RCM_info.region = {'NWP'};
RCM_info.histyears = 1976:2005;
RCM_info.rcpyears = 2006:2100;
RCM_info.months = 1:12;
RCM_info.scennames = {'rcp26', 'rcp45', 'rcp85'};
RCM_grid.dl = 1/20;

RCM_info.n_hist=length(RCM_info.histname);
RCM_info.n_rcp26=length(RCM_info.rcp26name);
RCM_info.n_rcp45=length(RCM_info.rcp45name);
RCM_info.n_rcp85=length(RCM_info.rcp85name);
RCM_info.n_scen=length(RCM_info.scennames);

% % %  read historical data
for i=1:RCM_info.n_hist
    tmp.testname = RCM_info.histname{i};
    tmp.histfilename = [RCM_info.dataroot, tmp.testname, filesep, RCM_info.phase, ...
        filesep, tmp.testname, '_NWP_ssh_analysis_', ...
        num2str(min(RCM_info.histyears)), '_', num2str(max(RCM_info.histyears)), '.nc'];
    tmp.ncinfo=ncinfo(tmp.histfilename);
    tmp_data.raw_ssh_86_05=ncread(tmp.histfilename, 'raw_ssh', [1 1 121], [inf inf inf]);
    RCM_data.histdata_05(:,:,i)=mean(tmp_data.raw_ssh_86_05(:,:,end-11:end),3);
    RCM_data.histdata_mean_86_05(:,:,i)=mean(tmp_data.raw_ssh_86_05, 3);
    
    RCM_data.histtrend(:,:,i)=ncread(tmp.histfilename, 'trend_filtered', [1 1], [inf inf]);
    RCM_grid.lon_rho = ncread(tmp.histfilename, 'lon_rho');
    RCM_grid.lat_rho = ncread(tmp.histfilename, 'lat_rho');
    RCM_grid.len_lon = size(RCM_grid.lon_rho, 1);
    RCM_grid.len_lat = size(RCM_grid.lat_rho, 2);
end

% % %  read RCP data (rcp26, rcp85)
if (RCM_info.n_rcp26 == RCM_info.n_rcp85)
    for i=1:RCM_info.n_rcp26
        for j=1:RCM_info.n_scen
            tmp.scenname=RCM_info.scennames{j};
            tmp.testname = RCM_info.([tmp.scenname,'name']){i};
            tmp.rcpfilename = [RCM_info.dataroot, tmp.testname, filesep, RCM_info.phase, ...
                filesep, tmp.testname, '_NWP_ssh_trend_', ...
                num2str(min(RCM_info.rcpyears)), '_', num2str(max(RCM_info.rcpyears)), '.nc'];
            
            if (strcmp(tmp.scenname, 'rcp45')~=1)
                tmp.ncinfo=ncinfo(tmp.rcpfilename);
                tmp_data.raw_ssh_06=ncread(tmp.histfilename, 'raw_ssh', [1 1 1], [inf inf 12]);
                RCM_data.([tmp.scenname,'data_06'])(:,:,i)=mean(tmp_data.raw_ssh_06(:,:,:),3);
                
                tmp_data.raw_ssh_81_00=ncread(tmp.rcpfilename, 'raw_ssh', [1 1 1140-239], [inf inf inf]);
                RCM_data.([tmp.scenname,'data_mean_81_00'])(:,:,i)=mean(tmp_data.raw_ssh_81_00, 3);
            end
            RCM_data.([tmp.scenname,'trend'])(:,:,i)=ncread(tmp.rcpfilename, 'trend_filtered', [1 1], [inf inf]);
        end
    end
end

% % %  read RCP data (rcp45)
for i=1:RCM_info.n_rcp45
    tmp.scenname='rcp45';
    tmp.testname = RCM_info.([tmp.scenname,'name']){i};
    
    for j= 1: length(RCM_info.months)
        tmp.rcpfilename = [RCM_info.dataroot, 'backup_surf_monthly', filesep, tmp.testname, filesep, RCM_info.phase, ...
            filesep, 'zeta', filesep, '2006', filesep, 'pck_', tmp.testname, '_zeta_monthly_', ...
            '2006', '_', num2str(RCM_info.months(j), '%02i'), '.nc'];
        tmp.ncinfo=ncinfo(tmp.rcpfilename);
        tmp_data.raw_ssh_06(:,:,j)=ncread(tmp.rcpfilename, 'zeta');
    end
    RCM_data.([tmp.scenname,'data_06'])(:,:,i)=mean(tmp_data.raw_ssh_06(:,:,:),3);
    
    for y=2081:2100
        for j= 1: length(RCM_info.months)
            ind=(y-2081)*12+j;
            tmp.rcpfilename = [RCM_info.dataroot, 'backup_surf_monthly', filesep, tmp.testname, filesep, RCM_info.phase, ...
                filesep, 'zeta', filesep, num2str(y), filesep, 'pck_', tmp.testname, '_zeta_monthly_', ...
                num2str(y), '_', num2str(RCM_info.months(j), '%02i'), '.nc'];
            tmp.ncinfo=ncinfo(tmp.rcpfilename);
            tmp_data.raw_ssh_81_00(:,:,ind)=ncread(tmp.rcpfilename, 'zeta');
        end
    end
    RCM_data.([tmp.scenname,'data_mean_81_00'])(:,:,i)=mean(tmp_data.raw_ssh_81_00, 3);
end


% %  get ensemble mean of historical trend
tmp.histtrend_ens=mean(RCM_data.histtrend,3);
[tmp.histtrend_ens_mean, error_status] = ...
    Func_0011_get_area_weighted_mean(tmp.histtrend_ens, RCM_grid.lon_rho, RCM_grid.lat_rho);
rcp45_correction_val=tmp.histtrend_ens_mean/1000.0;

% %  get ensemble mean of historical sea-level in 2005
RCM_data.histdata_05_ens=mean(RCM_data.histdata_05, 3);
[RCM_data.histdata_05_ens_mean, error_status] = ...
    Func_0011_get_area_weighted_mean(RCM_data.histdata_05_ens, RCM_grid.lon_rho, RCM_grid.lat_rho);

% %  get ensemble mean of rcp45 sea-level in 2006
RCM_data.rcp45data_06_ens=mean(RCM_data.rcp45data_06, 3);
[RCM_data.rcp45data_06_ens_mean, error_status] = ...
    Func_0011_get_area_weighted_mean(RCM_data.rcp45data_06_ens, RCM_grid.lon_rho, RCM_grid.lat_rho);

% %  correction of difference of reference value under rcp45
RCM_data.rcp45data_mean_81_00= ...
    RCM_data.rcp45data_mean_81_00 - RCM_data.rcp45data_06_ens_mean ...
    + RCM_data.histdata_05_ens_mean + rcp45_correction_val;

% %  get ensemble mean of 1986-2005 sea-level under historical periods
RCM_data.histens_mean_86_05 = mean(RCM_data.histdata_mean_86_05, 3);

% %  get ensemble mean of 2081-2100 sea-level under RCP scenarios
for j=1:RCM_info.n_scen
    tmp.scenname=RCM_info.scennames{j};
    RCM_data.([tmp.scenname,'ens_mean_81_00']) = mean(RCM_data.([tmp.scenname,'data_mean_81_00']), 3);
    RCM_findata.([tmp.scenname,'ens_trend']) = mean(RCM_data.([tmp.scenname,'trend']), 3);
end

pcolor(RCM_data.rcp26trend(:,:,3)'); shading flat; colorbar;


% %  get sea-level rise in 1986-2005 ~ 2081-2100 under RCP scenarios
for j=1:RCM_info.n_scen
    tmp.scenname=RCM_info.scennames{j};
    RCM_findata.([tmp.scenname,'ens_mean_slr']) = (RCM_data.([tmp.scenname,'ens_mean_81_00']) - RCM_data.histens_mean_86_05).*100; % m -> cm
end

% pcolor(RCM_data.histdata_05(:,:,4)'); shading flat; colorbar;
% histfilename =

for j=1:RCM_info.n_scen
    tmp.scenname=RCM_info.scennames{j};
    
% % %     make slr file
    tmp_nc.slrfilename = [RCM_info.dataroot, 'RCP_ens_slr', filesep, ...
                'RCM_ensemble_SLR_', tmp.scenname,  '.nc'];
    tmp_nc.slrncid = netcdf.create(tmp_nc.slrfilename,'NETCDF4');
    
    tmp_nc.lon_dimid = netcdf.defDim(tmp_nc.slrncid, 'xi_rho', RCM_grid.len_lon);
    tmp_nc.lat_dimid = netcdf.defDim(tmp_nc.slrncid,'eta_rho', RCM_grid.len_lat);
    
    netcdf.putAtt(tmp_nc.slrncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', ['Northwest Pacific 1/20 ensemble sea-level trend file']);
    netcdf.putAtt(tmp_nc.slrncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'source', ['ROMS RCM NWP 1/20 ', tmp.scenname, ' results']);
    netcdf.putAtt(tmp_nc.slrncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'downscaled GCMs', 'IPSL-CM5A-LR, IPSL_CM5A-MR, NorESM1-M, MPI-ESM-LR');
    netcdf.putAtt(tmp_nc.slrncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'author', ['Created by Y.Y. Kim and Y.K. Cho from MEPL, SNU']);
    netcdf.putAtt(tmp_nc.slrncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'contact', ['Y.Y. Kim - kimyy308@snu.ac.kr']);
    netcdf.putAtt(tmp_nc.slrncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'date', date);
                
    tmp_nc.lon_varid=netcdf.defVar(tmp_nc.slrncid, 'lon_rho', 'NC_DOUBLE', [tmp_nc.lon_dimid tmp_nc.lat_dimid]);
    netcdf.putAtt(tmp_nc.slrncid,tmp_nc.lon_varid,'long_name','longitude at rho grid');
    netcdf.putAtt(tmp_nc.slrncid,tmp_nc.lon_varid,'units','degree_east');
    
    tmp_nc.lat_varid=netcdf.defVar(tmp_nc.slrncid, 'lat_rho', 'NC_DOUBLE', [tmp_nc.lon_dimid tmp_nc.lat_dimid]);
    netcdf.putAtt(tmp_nc.slrncid,tmp_nc.lat_varid,'long_name','latitude at rho grid');
    netcdf.putAtt(tmp_nc.slrncid,tmp_nc.lat_varid,'units','degree_north');
    
    tmp_nc.slrvarid=netcdf.defVar(tmp_nc.slrncid, 'ensemble_slr', 'NC_FLOAT', [tmp_nc.lon_dimid tmp_nc.lat_dimid]);
    netcdf.putAtt(tmp_nc.slrncid,tmp_nc.slrvarid,'long_name','ensemble sea level rise of RCMs ensemble');
    netcdf.putAtt(tmp_nc.slrncid,tmp_nc.slrvarid,'period','mean(2081~2100) - mean(1986~2005)');
    netcdf.putAtt(tmp_nc.slrncid,tmp_nc.slrvarid,'units','cm');

    netcdf.endDef(tmp_nc.slrncid);
    
    netcdf.putVar(tmp_nc.slrncid, tmp_nc.lon_varid, [0 0], [RCM_grid.len_lon RCM_grid.len_lat], RCM_grid.lon_rho);
    netcdf.putVar(tmp_nc.slrncid, tmp_nc.lat_varid, [0 0], [RCM_grid.len_lon RCM_grid.len_lat], RCM_grid.lat_rho);

    netcdf.putVar(tmp_nc.slrncid, tmp_nc.slrvarid, [0 0], [RCM_grid.len_lon RCM_grid.len_lat], ...
        RCM_findata.([tmp.scenname,'ens_mean_slr']));

    netcdf.close(tmp_nc.slrncid);
    
% % %     make trend file
    tmp_nc.trendfilename = [RCM_info.dataroot, 'RCP_ens_slr', filesep, ...
                'RCM_ensemble_trend_', tmp.scenname,  '.nc'];
    tmp_nc.trendncid = netcdf.create(tmp_nc.trendfilename,'NETCDF4');
    
    tmp_nc.lon_dimid = netcdf.defDim(tmp_nc.trendncid, 'xi_rho', RCM_grid.len_lon);
    tmp_nc.lat_dimid = netcdf.defDim(tmp_nc.trendncid,'eta_rho', RCM_grid.len_lat);
    
    netcdf.putAtt(tmp_nc.trendncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', ['Northwest Pacific 1/20 ensemble sea-level trend file']);
    netcdf.putAtt(tmp_nc.trendncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'source', ['ROMS RCM NWP 1/20 ', tmp.scenname, ' results']);
    netcdf.putAtt(tmp_nc.trendncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'downscaled GCMs', 'IPSL-CM5A-LR, IPSL_CM5A-MR, NorESM1-M, MPI-ESM-LR');
    netcdf.putAtt(tmp_nc.trendncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'author', ['Created by Y.Y. Kim and Y.K. Cho from MEPL, SNU']);
    netcdf.putAtt(tmp_nc.trendncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'contact', ['Y.Y. Kim - kimyy308@snu.ac.kr']);
    netcdf.putAtt(tmp_nc.trendncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'date', date);
                
    tmp_nc.lon_varid=netcdf.defVar(tmp_nc.trendncid, 'lon_rho', 'NC_DOUBLE', [tmp_nc.lon_dimid tmp_nc.lat_dimid]);
    netcdf.putAtt(tmp_nc.trendncid,tmp_nc.lon_varid,'long_name','longitude at ROMS rho grid');
    netcdf.putAtt(tmp_nc.trendncid,tmp_nc.lon_varid,'units','degree_east');
    
    tmp_nc.lat_varid=netcdf.defVar(tmp_nc.trendncid, 'lat_rho', 'NC_DOUBLE', [tmp_nc.lon_dimid tmp_nc.lat_dimid]);
    netcdf.putAtt(tmp_nc.trendncid,tmp_nc.lat_varid,'long_name','latitude at ROMS rho grid');
    netcdf.putAtt(tmp_nc.trendncid,tmp_nc.lat_varid,'units','degree_north');
    
    tmp_nc.trendvarid=netcdf.defVar(tmp_nc.trendncid, 'ensemble_trend', 'NC_FLOAT', [tmp_nc.lon_dimid tmp_nc.lat_dimid]);
    netcdf.putAtt(tmp_nc.trendncid,tmp_nc.trendvarid,'long_name','ensemble sea-level trend of RCMs ensemble');
    netcdf.putAtt(tmp_nc.trendncid,tmp_nc.trendvarid,'period','2006 ~ 2100');
    netcdf.putAtt(tmp_nc.trendncid,tmp_nc.trendvarid,'units','mm/yr');

    netcdf.endDef(tmp_nc.trendncid);
    
    netcdf.putVar(tmp_nc.trendncid, tmp_nc.lon_varid, [0 0], [RCM_grid.len_lon RCM_grid.len_lat], RCM_grid.lon_rho);
    netcdf.putVar(tmp_nc.trendncid, tmp_nc.lat_varid, [0 0], [RCM_grid.len_lon RCM_grid.len_lat], RCM_grid.lat_rho);

    netcdf.putVar(tmp_nc.trendncid, tmp_nc.trendvarid, [0 0], [RCM_grid.len_lon RCM_grid.len_lat], ...
        RCM_findata.([tmp.scenname,'ens_trend']));

    netcdf.close(tmp_nc.trendncid);
end

