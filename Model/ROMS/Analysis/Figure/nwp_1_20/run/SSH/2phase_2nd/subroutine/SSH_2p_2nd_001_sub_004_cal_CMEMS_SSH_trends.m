% %         calculation of SSH trends (raw data)      
% %  Updated 07-Oct-2021 by Yong-Yub Kim, 


% % %         sea level analysis
fprintf('\n')
disp('subroutine SSH_2p_1st_001_sub_004_cal_CMEMS_SSH_trends') 

if (exist(CMEMS_info.matname , 'file') == 2 || flags.fig_switch(1)~=2)   
    load(RCM_info.matname);
    load(CMEMS_info.matname);
end
CMEMS_info.matname_trends = [RCM_info.savedir,RCM_info.regionname, '_CMEMS_ssh_trend_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
CMEMS_info.filename_trends = strcat(RCM_info.savedir,RCM_info.regionname, 'CMEMS_ssh_trend_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '.nc');

CMEMS_data_trend.conf_level=0.95;

% % %         trend

if (exist(CMEMS_info.matname_trends , 'file') ~= 2 || flags.fig_tmp==2)   
% % %         yearly trend
   disp('yearly trend calculation start ')

   % % %   yearly trend  time check
    tmp.lap_time_j=tic;
    tmp.elapsed=toc(tmp.lap_time_j);
    tmp.totlap= (CMEMS_grid.size_lon);
    tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-1+1)/1);
    
    CMEMS_data_trend.yearly_trend(1:CMEMS_grid.size_lon,1:CMEMS_grid.size_lat)=NaN;
    CMEMS_data_trend.yearly_trend_lower_bounds(1:CMEMS_grid.size_lon,1:CMEMS_grid.size_lat)=NaN;
    CMEMS_data_trend.yearly_trend_upper_bounds(1:CMEMS_grid.size_lon,1:CMEMS_grid.size_lat)=NaN;
    CMEMS_data_trend.yearly_trend_rmse(1:CMEMS_grid.size_lon,1:CMEMS_grid.size_lat)=NaN;
    CMEMS_data_trend.yearly_trend_rsquare(1:CMEMS_grid.size_lon,1:CMEMS_grid.size_lat)=NaN;
    CMEMS_data_trend.yearly_trend_adjrsquare(1:CMEMS_grid.size_lon,1:CMEMS_grid.size_lat)=NaN;
    
    for i=1:CMEMS_grid.size_lon
        tmp.elapsed=toc(tmp.lap_time_j);
        tmp.templap= i-1;
        fprintf(repmat('\b',1,sum(tmp.nchar)))  % remove printed time
        tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed/tmp.templap*(tmp.totlap-tmp.templap+1));

        for j=1:CMEMS_grid.size_lat
             if isfinite(sum(CMEMS_data.sla_yearly_mean(i,j,:)))
                [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', squeeze(CMEMS_data.sla_yearly_mean(i,j,:)).*1000.0, 'poly1'); %% m/y -> mm/y
                CMEMS_data_trend.tmp_conf_bounds=confint(tmp_fitobject, CMEMS_data_trend.conf_level);
                CMEMS_data_trend.yearly_trend_lower_bounds(i,j)=CMEMS_data_trend.tmp_conf_bounds(1,1);
                CMEMS_data_trend.yearly_trend_upper_bounds(i,j)=CMEMS_data_trend.tmp_conf_bounds(2,1);
                CMEMS_data_trend.yearly_trend(i,j)=tmp_fitobject.p1;
                CMEMS_data_trend.yearly_trend_rmse(i,j)=tmp_gof.rmse;
                CMEMS_data_trend.yearly_trend_rsquare(i,j)=tmp_gof.rsquare;
                CMEMS_data_trend.yearly_trend_adjrsquare(i,j)=tmp_gof.adjrsquare;
            end
        end
    end
%     trend_filtered = trend_filtered * 1000.0; %% m/y -> mm/y
    fprintf('\n')
    disp('yearly trend complete') 

    save(CMEMS_info.matname_trends, 'CMEMS_data_trend', '-v7.3');
end