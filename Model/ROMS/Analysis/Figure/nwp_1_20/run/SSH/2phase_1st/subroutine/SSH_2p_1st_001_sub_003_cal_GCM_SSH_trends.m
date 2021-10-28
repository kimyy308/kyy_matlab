% %         calculation of SSH trends (raw data)      
% %  Updated 07-Oct-2021 by Yong-Yub Kim, 


% % %         sea level analysis
fprintf('\n')
disp('subroutine SSH_2p_1st_001_sub_003_cal_GCM_SSH_trends') 

if (exist(GCM_info.matname , 'file') == 2 || flags.fig_switch(1)~=2)   
    load(RCM_info.matname); 
    load(GCM_info.matname);
    load(CMEMS_info.matname);
    load(GCM_info.matname_interped);
end
GCM_info.matname_trends = [RCM_info.savedir,GCM_info.testname,'_',GCM_info.regionname, '_GCM_ssh_trend_', ...
    num2str(min(GCM_info.years),'%04i'),'_',num2str(max(GCM_info.years),'%04i'),'.mat'];
GCM_info.filename_trends = strcat(RCM_info.savedir, GCM_info.testname,'_',GCM_info.regionname, 'GCM_ssh_trend_', ...
    num2str(min(GCM_info.years),'%04i'),'_',num2str(max(GCM_info.years),'%04i'), '.nc');

GCM_info.matname_trends_interped = [RCM_info.savedir,GCM_info.testname,'_',GCM_info.regionname, '_GCM_ssh_interped_trend_', ...
    num2str(min(GCM_info.years),'%04i'),'_',num2str(max(GCM_info.years),'%04i'),'.mat'];

GCM_data_trend.conf_level=0.95;
GCM_data_interped_trend.conf_level=GCM_data_trend.conf_level;

% % %         trend

if (exist(GCM_info.matname_trends , 'file') ~= 2 || flags.fig_tmp==2)   
   
    
% % %     disp('Trend calculation start ')
% % % % % %   trend  time check
% % %     tmp.lap_time_j=tic;
% % %     tmp.elapsed=toc(tmp.lap_time_j);
% % %     tmp.totlap= (GCM_grid.size_lon*GCM_grid.size_lat);
% % %     tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-1+1)/1);
% % % 
% % %     GCM_data_trend.trend(1:GCM_grid.size_lon,1:GCM_grid.size_lat)=NaN;
% % %     GCM_data_trend.trend_lower_bounds(1:GCM_grid.size_lon,1:GCM_grid.size_lat)=NaN;
% % %     GCM_data_trend.trend_upper_bounds(1:GCM_grid.size_lon,1:GCM_grid.size_lat)=NaN;
% % %     GCM_data_trend.trend_rmse(1:GCM_grid.size_lon,1:GCM_grid.size_lat)=NaN;
% % %     GCM_data_trend.trend_rsquare(1:GCM_grid.size_lon,1:GCM_grid.size_lat)=NaN;
% % %     GCM_data_trend.trend_adjrsquare(1:GCM_grid.size_lon,1:GCM_grid.size_lat)=NaN;
% % %     for i=1:GCM_grid.size_lon
% % %         for j=1:GCM_grid.size_lat
% % % %             [p, S]=polyfit(GCM_time.trendtime, squeeze(GCM_data.all(i,j,:))',1);
% % %             tmp.elapsed=toc(tmp.lap_time_j);
% % %             tmp.templap= i*GCM_grid.size_lat+j;
% % %             fprintf(repmat('\b',1,sum(tmp.nchar)))  % remove printed time
% % %             tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-tmp.templap+1)/tmp.templap);
% % %             
% % %             if isfinite(sum(GCM_data.all(i,j,:)))
% % %                 [tmp_fitobject, tmp_gof] = fit(GCM_time.trendtime', squeeze(GCM_data.all(i,j,:)).*1000.0, 'poly1'); %% m/y -> mm/y
% % %                 GCM_data_trend.tmp_conf_bounds=confint(tmp_fitobject, GCM_data_trend.conf_level);
% % %                 GCM_data_trend.trend_lower_bounds(i,j)=GCM_data_trend.tmp_conf_bounds(1,1);
% % %                 GCM_data_trend.trend_upper_bounds(i,j)=GCM_data_trend.tmp_conf_bounds(2,1);
% % %                 GCM_data_trend.trend(i,j)=tmp_fitobject.p1;
% % %                 GCM_data_trend.trend_rmse(i,j)=tmp_gof.rmse;
% % %                 GCM_data_trend.trend_rsquare(i,j)=tmp_gof.rsquare;
% % %                 GCM_data_trend.trend_adjrsquare(i,j)=tmp_gof.adjrsquare;
% % %             end
% % %         end
% % %     end
% % %     fprintf('\n')
% % %    disp('trend complete') 
   

% % %         yearly trend
   disp('yearly trend calculation start ')

   % % %   yearly trend  time check
    tmp.lap_time_j=tic;
    tmp.elapsed=toc(tmp.lap_time_j);
    tmp.totlap= (GCM_grid.size_lon);
    tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-1+1)/1);
    
    GCM_data_trend.yearly_trend(1:GCM_grid.size_lon,1:GCM_grid.size_lat)=NaN;
    GCM_data_trend.yearly_trend_lower_bounds(1:GCM_grid.size_lon,1:GCM_grid.size_lat)=NaN;
    GCM_data_trend.yearly_trend_upper_bounds(1:GCM_grid.size_lon,1:GCM_grid.size_lat)=NaN;
    GCM_data_trend.yearly_trend_rmse(1:GCM_grid.size_lon,1:GCM_grid.size_lat)=NaN;
    GCM_data_trend.yearly_trend_rsquare(1:GCM_grid.size_lon,1:GCM_grid.size_lat)=NaN;
    GCM_data_trend.yearly_trend_adjrsquare(1:GCM_grid.size_lon,1:GCM_grid.size_lat)=NaN;
    
    for i=1:GCM_grid.size_lon
        tmp.elapsed=toc(tmp.lap_time_j);
        tmp.templap= i-1;
        fprintf(repmat('\b',1,sum(tmp.nchar)))  % remove printed time
        tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-tmp.templap+1)/tmp.templap);

        for j=1:GCM_grid.size_lat
             if isfinite(sum(GCM_data.yearly_mean(i,j,:)))
                [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', squeeze(GCM_data.yearly_mean(i,j,:)).*1000.0, 'poly1'); %% m/y -> mm/y
                GCM_data_trend.tmp_conf_bounds=confint(tmp_fitobject, GCM_data_trend.conf_level);
                GCM_data_trend.yearly_trend_lower_bounds(i,j)=GCM_data_trend.tmp_conf_bounds(1,1);
                GCM_data_trend.yearly_trend_upper_bounds(i,j)=GCM_data_trend.tmp_conf_bounds(2,1);
                GCM_data_trend.yearly_trend(i,j)=tmp_fitobject.p1;
                GCM_data_trend.yearly_trend_rmse(i,j)=tmp_gof.rmse;
                GCM_data_trend.yearly_trend_rsquare(i,j)=tmp_gof.rsquare;
                GCM_data_trend.yearly_trend_adjrsquare(i,j)=tmp_gof.adjrsquare;
            end
        end
    end
%     trend_filtered = trend_filtered * 1000.0; %% m/y -> mm/y
    fprintf('\n')
    disp('yearly trend complete') 

    save(GCM_info.matname_trends, 'GCM_data_trend', '-v7.3');
    
    
% % %         yearly trend (interped)
   disp('yearly trend (interped) calculation start  ')

   % % %   yearly trend  time check
    tmp.lap_time_j=tic;
    tmp.elapsed=toc(tmp.lap_time_j);
    tmp.totlap= (CMEMS_grid.size_lon);
    tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-1+1)/1);
    
    GCM_data_interped_trend.yearly_trend(1:CMEMS_grid.size_lon,1:CMEMS_grid.size_lat)=NaN;
    GCM_data_interped_trend.yearly_trend_lower_bounds(1:CMEMS_grid.size_lon,1:CMEMS_grid.size_lat)=NaN;
    GCM_data_interped_trend.yearly_trend_upper_bounds(1:CMEMS_grid.size_lon,1:CMEMS_grid.size_lat)=NaN;
    GCM_data_interped_trend.yearly_trend_rmse(1:CMEMS_grid.size_lon,1:CMEMS_grid.size_lat)=NaN;
    GCM_data_interped_trend.yearly_trend_rsquare(1:CMEMS_grid.size_lon,1:CMEMS_grid.size_lat)=NaN;
    GCM_data_interped_trend.yearly_trend_adjrsquare(1:CMEMS_grid.size_lon,1:CMEMS_grid.size_lat)=NaN;
    
    for i=1:CMEMS_grid.size_lon
        tmp.elapsed=toc(tmp.lap_time_j);
        tmp.templap= i-1;
        fprintf(repmat('\b',1,sum(tmp.nchar)))  % remove printed time
        tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-tmp.templap+1)/tmp.templap);

        for j=1:CMEMS_grid.size_lat
             if isfinite(sum(GCM_data_interped.yearly_mean(i,j,:)))
                [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', squeeze(GCM_data_interped.yearly_mean(i,j,:)).*1000.0, 'poly1'); %% m/y -> mm/y
                GCM_data_interped_trend.tmp_conf_bounds=confint(tmp_fitobject, GCM_data_interped_trend.conf_level);
                GCM_data_interped_trend.yearly_trend_lower_bounds(i,j)=GCM_data_interped_trend.tmp_conf_bounds(1,1);
                GCM_data_interped_trend.yearly_trend_upper_bounds(i,j)=GCM_data_interped_trend.tmp_conf_bounds(2,1);
                GCM_data_interped_trend.yearly_trend(i,j)=tmp_fitobject.p1;
                GCM_data_interped_trend.yearly_trend_rmse(i,j)=tmp_gof.rmse;
                GCM_data_interped_trend.yearly_trend_rsquare(i,j)=tmp_gof.rsquare;
                GCM_data_interped_trend.yearly_trend_adjrsquare(i,j)=tmp_gof.adjrsquare;
            end
        end
    end
%     trend_filtered = trend_filtered * 1000.0; %% m/y -> mm/y
    fprintf('\n')
    disp('yearly trend complete') 

    save(GCM_info.matname_trends_interped, 'GCM_data_interped_trend', '-v7.3');
    
    
    
% % % % % %         make ncfile
% % %     ncid = netcdf.create(ncoutfilename,'NETCDF4');
% % % 
% % %     lon_dimid = netcdf.defDim(ncid, 'lon', GCM_grid.size_lon);
% % %     lat_dimid = netcdf.defDim(ncid,'lat',GCM_grid.size_lat);
% % %     time_dimid = netcdf.defDim(ncid, 'time', 0);
% % %     clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);
% % % 
% % %     netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
% % %         'type', ['NWP 1/20 _ ', testname, 'model, cmems monthly SSH analysis file']);
% % %     netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
% % %         'title', [' monthly SSH analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
% % %     netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
% % %         'source', [' ROMS NWP 1/20 data from _ ',testname ]);
% % %     netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
% % %         'author', 'Created by Y.Y.Kim');
% % %     netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
% % %         'date', date);
% % % 
% % %     timevarid=netcdf.defVar(ncid, 'time', 'NC_DOUBLE', time_dimid);
% % %     netcdf.putAtt(ncid,timevarid,'long_name','time');
% % %     netcdf.putAtt(ncid,timevarid,'units','days since 1900-12-31 00:00:00');
% % %     netcdf.putAtt(ncid,timevarid,'calendar','gregorian');
% % %     
% % %     yearly_timevarid=netcdf.defVar(ncid, 'yearly_time', 'NC_DOUBLE', yearly_time_dimid);
% % %     netcdf.putAtt(ncid,yearly_timevarid,'long_name','yearly_time');
% % %     netcdf.putAtt(ncid,yearly_timevarid,'units','days since 1900-12-31 00:00:00');
% % %     netcdf.putAtt(ncid,yearly_timevarid,'calendar','gregorian');
% % % 
% % %     lonvarid=netcdf.defVar(ncid, 'lon', 'NC_DOUBLE', [lon_dimid lat_dimid]);
% % %     netcdf.putAtt(ncid,lonvarid,'long_name','lon_model');
% % %     netcdf.putAtt(ncid,lonvarid,'units','degree_east');
% % % 
% % %     latvarid=netcdf.defVar(ncid, 'lat', 'NC_DOUBLE', [lon_dimid lat_dimid]);
% % %     netcdf.putAtt(ncid,latvarid,'long_name','lat_model');
% % %     netcdf.putAtt(ncid,latvarid,'units','degree_north');
% % % 
% % %     raw_sshvarid=netcdf.defVar(ncid, 'raw_ssh', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
% % %     netcdf.putAtt(ncid,raw_sshvarid,'long_name','raw_ssh');
% % %     netcdf.putAtt(ncid,raw_sshvarid,'units','m');
% % % 
% % %     ssh_filteredvarid=netcdf.defVar(ncid, 'yearly_ssh', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
% % %     netcdf.putAtt(ncid,ssh_filteredvarid,'long_name','yearly_ssh');
% % %     netcdf.putAtt(ncid,ssh_filteredvarid,'units','m');
% % % 
% % %     trendvarid=netcdf.defVar(ncid, 'yearly_trend', 'NC_FLOAT', [lon_dimid lat_dimid]);
% % %     netcdf.putAtt(ncid,trendvarid,'long_name','yearly_trend');
% % %     netcdf.putAtt(ncid,trendvarid,'units','mm/year');
% % % 
% % %     clim_sshvarid=netcdf.defVar(ncid, 'clim_ssh', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
% % %     netcdf.putAtt(ncid,clim_sshvarid,'long_name','clim_ssh');
% % %     netcdf.putAtt(ncid,clim_sshvarid,'units','m');
% % % 
% % %     netcdf.endDef(ncid);
% % % 
% % %     netcdf.putVar(ncid, timevarid, 0, length(ROMS_time.ftime), ROMS_time.ftime);
% % %     netcdf.putVar(ncid, yearly_timevarid, 0, length(ROMS_time.yearlytime), ROMS_time.yearlytime);    
% % %     netcdf.putVar(ncid, lonvarid, [0 0], [GCM_grid.size_lon GCM_grid.size_lat], ROMS_grid.lon);
% % %     netcdf.putVar(ncid, latvarid, [0 0], [GCM_grid.size_lon GCM_grid.size_lat], ROMS_grid.lat);
% % %     netcdf.putVar(ncid, raw_sshvarid, [0 0 0], [GCM_grid.size_lon GCM_grid.size_lat length(ROMS_time.ftime)], ROMS_data.all);
% % %     netcdf.putVar(ncid, trendvarid, [0 0], [GCM_grid.size_lon, GCM_grid.size_lat], trend);
% % %     netcdf.putVar(ncid, clim_sshvarid, [0 0 0], [GCM_grid.size_lon, GCM_grid.size_lat length(climtime)], comb_spatial_meanmodel);
% % %     netcdf.close(ncid);


end