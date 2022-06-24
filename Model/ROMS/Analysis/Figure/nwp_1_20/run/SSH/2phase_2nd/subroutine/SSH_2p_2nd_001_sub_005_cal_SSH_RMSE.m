% %         calculation of SSH RMSE (raw data)      
% %  Updated 06-Oct-2021 by Yong-Yub Kim,

% % %         sea level analysis
fprintf('\n')
disp('subroutine SSH_2p_1st_001_sub_002_cal_SSH_RMSE') 

if (exist(RCM_info.matname , 'file') == 2 || flags.fig_switch(1)~=2)   
    load(RCM_info.matname); 
    load(GCM_info.matname);
    load(CMEMS_info.matname);
end
RCM_info.matname_RMSE = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_RCM_ssh_RMSE_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
RCM_info.filename_RMSE = strcat(RCM_info.savedir, RCM_info.testname,'_',RCM_info.regionname, 'RCM_ssh_trend_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '.nc');



1985-2014 AKP4 mean



RCM_data_trend.conf_level=0.95;

% % %         trend

if (exist(RCM_info.matname_RMSE , 'file') ~= 2 || flags.fig_switch(2)==2)   
   
    
% % %     disp('Trend calculation start ')
% % % % % %   trend  time check
% % %     tmp.lap_time_j=tic;
% % %     tmp.elapsed=toc(tmp.lap_time_j);
% % %     tmp.totlap= (RCM_grid.size_lon_rho*RCM_grid.size_lat_rho);
% % %     tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-1+1)/1);
% % % 
% % %     RCM_data_trend.trend(1:RCM_grid.size_lon_rho,1:RCM_grid.size_lat_rho)=NaN;
% % %     RCM_data_trend.trend_lower_bounds(1:RCM_grid.size_lon_rho,1:RCM_grid.size_lat_rho)=NaN;
% % %     RCM_data_trend.trend_upper_bounds(1:RCM_grid.size_lon_rho,1:RCM_grid.size_lat_rho)=NaN;
% % %     RCM_data_trend.trend_rmse(1:RCM_grid.size_lon_rho,1:RCM_grid.size_lat_rho)=NaN;
% % %     RCM_data_trend.trend_rsquare(1:RCM_grid.size_lon_rho,1:RCM_grid.size_lat_rho)=NaN;
% % %     RCM_data_trend.trend_adjrsquare(1:RCM_grid.size_lon_rho,1:RCM_grid.size_lat_rho)=NaN;
% % %     for i=1:RCM_grid.size_lon_rho
% % %         for j=1:RCM_grid.size_lat_rho
% % % %             [p, S]=polyfit(RCM_time.trendtime, squeeze(RCM_data.all(i,j,:))',1);
% % %             tmp.elapsed=toc(tmp.lap_time_j);
% % %             tmp.templap= i*RCM_grid.size_lat_rho+j;
% % %             fprintf(repmat('\b',1,sum(tmp.nchar)))  % remove printed time
% % %             tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-tmp.templap+1)/tmp.templap);
% % %             
% % %             if isfinite(sum(RCM_data.all(i,j,:)))
% % %                 [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime', squeeze(RCM_data.all(i,j,:)).*1000.0, 'poly1'); %% m/y -> mm/y
% % %                 RCM_data_trend.tmp_conf_bounds=confint(tmp_fitobject, RCM_data_trend.conf_level);
% % %                 RCM_data_trend.trend_lower_bounds(i,j)=RCM_data_trend.tmp_conf_bounds(1,1);
% % %                 RCM_data_trend.trend_upper_bounds(i,j)=RCM_data_trend.tmp_conf_bounds(2,1);
% % %                 RCM_data_trend.trend(i,j)=tmp_fitobject.p1;
% % %                 RCM_data_trend.trend_rmse(i,j)=tmp_gof.rmse;
% % %                 RCM_data_trend.trend_rsquare(i,j)=tmp_gof.rsquare;
% % %                 RCM_data_trend.trend_adjrsquare(i,j)=tmp_gof.adjrsquare;
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
    tmp.totlap= (RCM_grid.size_lon_rho);
    tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-1+1)/1);
    
    RCM_data_trend.yearly_trend(1:RCM_grid.size_lon_rho,1:RCM_grid.size_lat_rho)=NaN;
    RCM_data_trend.yearly_trend_lower_bounds(1:RCM_grid.size_lon_rho,1:RCM_grid.size_lat_rho)=NaN;
    RCM_data_trend.yearly_trend_upper_bounds(1:RCM_grid.size_lon_rho,1:RCM_grid.size_lat_rho)=NaN;
    RCM_data_trend.yearly_trend_rmse(1:RCM_grid.size_lon_rho,1:RCM_grid.size_lat_rho)=NaN;
    RCM_data_trend.yearly_trend_rsquare(1:RCM_grid.size_lon_rho,1:RCM_grid.size_lat_rho)=NaN;
    RCM_data_trend.yearly_trend_adjrsquare(1:RCM_grid.size_lon_rho,1:RCM_grid.size_lat_rho)=NaN;
    
    for i=1:RCM_grid.size_lon_rho
        tmp.elapsed=toc(tmp.lap_time_j);
        tmp.templap= i;
        fprintf(repmat('\b',1,sum(tmp.nchar)))  % remove printed time
        tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed/tmp.templap*(tmp.totlap-tmp.templap+1));

        for j=1:RCM_grid.size_lat_rho
             if isfinite(sum(RCM_data.yearly_mean(i,j,:)))
                [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', squeeze(RCM_data.yearly_mean(i,j,:)).*1000.0, 'poly1'); %% m/y -> mm/y
                RCM_data_trend.tmp_conf_bounds=confint(tmp_fitobject, RCM_data_trend.conf_level);
                RCM_data_trend.yearly_trend_lower_bounds(i,j)=RCM_data_trend.tmp_conf_bounds(1,1);
                RCM_data_trend.yearly_trend_upper_bounds(i,j)=RCM_data_trend.tmp_conf_bounds(2,1);
                RCM_data_trend.yearly_trend(i,j)=tmp_fitobject.p1;
                RCM_data_trend.yearly_trend_rmse(i,j)=tmp_gof.rmse;
                RCM_data_trend.yearly_trend_rsquare(i,j)=tmp_gof.rsquare;
                RCM_data_trend.yearly_trend_adjrsquare(i,j)=tmp_gof.adjrsquare;
            end
        end
    end
%     trend_filtered = trend_filtered * 1000.0; %% m/y -> mm/y
    fprintf('\n')
    disp('yearly trend complete') 

    save(RCM_info.matname_RMSE, 'RCM_data_trend', '-v7.3');

% % % % % %         make ncfile
% % %     ncid = netcdf.create(ncoutfilename,'NETCDF4');
% % % 
% % %     lon_dimid = netcdf.defDim(ncid, 'lon', RCM_grid.size_lon_rho);
% % %     lat_dimid = netcdf.defDim(ncid,'lat',RCM_grid.size_lat_rho);
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
% % %     lon_rhovarid=netcdf.defVar(ncid, 'lon_rho', 'NC_DOUBLE', [lon_dimid lat_dimid]);
% % %     netcdf.putAtt(ncid,lon_rhovarid,'long_name','lon_model');
% % %     netcdf.putAtt(ncid,lon_rhovarid,'units','degree_east');
% % % 
% % %     lat_rhovarid=netcdf.defVar(ncid, 'lat_rho', 'NC_DOUBLE', [lon_dimid lat_dimid]);
% % %     netcdf.putAtt(ncid,lat_rhovarid,'long_name','lat_model');
% % %     netcdf.putAtt(ncid,lat_rhovarid,'units','degree_north');
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
% % %     netcdf.putVar(ncid, lon_rhovarid, [0 0], [RCM_grid.size_lon_rho RCM_grid.size_lat_rho], ROMS_grid.lon_rho);
% % %     netcdf.putVar(ncid, lat_rhovarid, [0 0], [RCM_grid.size_lon_rho RCM_grid.size_lat_rho], ROMS_grid.lat_rho);
% % %     netcdf.putVar(ncid, raw_sshvarid, [0 0 0], [RCM_grid.size_lon_rho RCM_grid.size_lat_rho length(ROMS_time.ftime)], ROMS_data.all);
% % %     netcdf.putVar(ncid, trendvarid, [0 0], [RCM_grid.size_lon_rho, RCM_grid.size_lat_rho], trend);
% % %     netcdf.putVar(ncid, clim_sshvarid, [0 0 0], [RCM_grid.size_lon_rho, RCM_grid.size_lat_rho length(climtime)], comb_spatial_meanmodel);
% % %     netcdf.close(ncid);


end