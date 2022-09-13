% %         get model data and cmems data      
% %  Updated 13-Jan-2022 by Yong-Yub Kim, create

disp('subroutine v_2p_2nd_003_sub_001_get_v_monthly') 
disp([RCM_info.testname, ', ', RCM_info.regionname, ', ', tmp.variable, ', ', ...
    num2str(min(RCM_info.years),'%04i'),' ~ ',num2str(max(RCM_info.years),'%04i')])

%% parameters, filename, tic set
tmp.lap_time_j=tic;
run(tmp.param_script);
tmp.ind_for=1;
RCM_info.matname = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_RCM_', tmp.variable, '_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '_', ...
    RCM_info.season, '.mat'];
GCM_info.matname = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_GCM_', tmp.variable, '_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '_', ...
    RCM_info.season, '.mat'];

[RCM_info.months, tmp.error_status] = Func_0019_get_month_from_season(RCM_info.season);

tmp.elapsed=toc(tmp.lap_time_j);
tmp.totlap= (length(RCM_info.years)*length(RCM_info.months));
tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-1+1)/1);
switch RCM_info.regionname
    case 'EKWC2'
        tmp.EKWC_lon_ind_max = 70; % for EKWC2
    otherwise
%         tmp.EKWC_lon_ind_max = 70;
end

%% get data loop (year, month)
if (exist(RCM_info.matname , 'file') ~= 2 || flags.fig_tmp==2)
    RCM_info.size_t=length(RCM_info.years)*length(RCM_info.months);
    RCM_info.size_years=length(RCM_info.years);
    RCM_info.size_months=length(RCM_info.months);
    for yearij = 1:length(RCM_info.years)
        tmp.year = RCM_info.years(yearij);
%         [RCM_info.days, tmp.error_status] = Func_0020_get_day_from_season(RCM_info.months, year2day(tmp.year)-365);
        for monthij = 1:length(RCM_info.months)
            
%% print present time
            tmp.elapsed=toc(tmp.lap_time_j);
            tmp.templap= (yearij-1)*length(RCM_info.months)+monthij-1;
            fprintf(repmat('\b',1,sum(tmp.nchar)))  % remove printed time
            tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-tmp.templap+1)/tmp.templap);
            tmp.month = RCM_info.months(monthij);
            
%% set data filename (raw grid, ES)
            switch(RCM_info.testname)
                otherwise
                    RCM_info.filename = strcat(RCM_info.filedir, tmp.fs, 'run', tmp.fs, ...
                        tmp.variable, tmp.fs, num2str(tmp.year,'%04i'), tmp.fs, ...
                    'NWP_pck_', RCM_info.testname, '_', tmp.variable, '_monthly_', num2str(tmp.year,'%04i'), ...
                    '_', num2str(tmp.month,'%02i'), '.nc');
            end
            
%% get RCM grids
            if (yearij == 1 && monthij == 1)
%                 RCM_grid.filename = RCM_info.filename;
%                 RCM_grid.('lon_rho')= ncread(RCM_info.filename, 'lon_rho');
%                 RCM_grid.('lat_rho')= ncread(RCM_info.filename, 'lat_rho');
%                 RCM_grid.('depth')= ncread(RCM_info.filename, 'depth');
                [RCM_grid.ind_w, RCM_grid.ind_e, RCM_grid.ind_s, RCM_grid.ind_n] = ...
                    findind_Y(1/20, RCM_grid.domain(1:4), RCM_grid.lon_rho, RCM_grid.lat_rho);
                tmp.lonlats={'lon', 'lat'};
                for gridi=1:length(tmp.lonlats)
                    RCM_grid.(['cut_', tmp.lonlats{gridi}, '_rho']) = ...
                        RCM_grid.([tmp.lonlats{gridi}, '_rho'])(RCM_grid.ind_w(1):RCM_grid.ind_e(1), ...
                        RCM_grid.ind_s(1):RCM_grid.ind_n(1));
                end
                
                RCM_grid.size_lon_rho = RCM_grid.ind_e-RCM_grid.ind_w+1;
                RCM_grid.size_lat_rho = RCM_grid.ind_n-RCM_grid.ind_s+1;
                RCM_grid.size_lon_u = RCM_grid.ind_e-RCM_grid.ind_w;
                RCM_grid.size_lat_u = RCM_grid.ind_n-RCM_grid.ind_s+1;
                RCM_grid.size_lon_v = RCM_grid.ind_e-RCM_grid.ind_w+1;
                RCM_grid.size_lat_v = RCM_grid.ind_n-RCM_grid.ind_s;
                
                RCM_grid.mask_region = double(inpolygon(RCM_grid.cut_lon_rho, RCM_grid.cut_lat_rho, ...
                    RCM_grid.refpolygon(:,1), RCM_grid.refpolygon(:,2)));
                RCM_grid.mask_region(RCM_grid.mask_region==0)=NaN;
                RCM_info.data_info = ncinfo(RCM_info.filename, param.varname);
                RCM_grid.mask_ocean= ncread(RCM_info.filename, tmp.variable, [RCM_grid.ind_w RCM_grid.ind_s, 1, 1], ...
                    [RCM_grid.size_lon_rho RCM_grid.size_lat_rho, 1, 1]);
                RCM_grid.mask_ocean(isfinite(RCM_grid.mask_ocean))=1;
                RCM_grid.mask_ocean = RCM_grid.mask_ocean .*  RCM_grid.mask_region;
                RCM_grid.mask_land = RCM_grid.mask_ocean;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                RCM_grid.mask_land(isfinite(RCM_grid.mask_ocean))=NaN;
                RCM_grid.mask_land(isnan(RCM_grid.mask_ocean))=1;
            end
            
%% get RCM data
            RCM_data.tmp_data = ncread(RCM_info.filename, tmp.variable, [RCM_grid.ind_w RCM_grid.ind_s, 1, 1], ...
                    [RCM_grid.size_lon_rho RCM_grid.size_lat_rho, 1, 1]);
            RCM_data.tmp_data=RCM_data.tmp_data.*RCM_grid.mask_region; % [x y z] .* [x y]
            RCM_grid.size_depth = size(RCM_data.tmp_data,3);
            if (isfield(RCM_data, 'all') ~= 1)
                RCM_data.all=NaN(RCM_grid.size_lon_rho, RCM_grid.size_lat_rho, RCM_info.size_t);
%                 RCM_data.clim_mean=(zeros([RCM_grid.size_lon_rho,RCM_grid.size_lat_rho,RCM_info.size_months]));
                RCM_data.yearly_mean=(zeros([RCM_grid.size_lon_rho,RCM_grid.size_lat_rho,RCM_info.size_years]));
                RCM_data.maxval=NaN(RCM_grid.size_lat_rho, RCM_info.size_t);
                RCM_data.maxind=NaN(RCM_grid.size_lat_rho, RCM_info.size_t);
                RCM_data.yearly_maxval=NaN(RCM_grid.size_lat_rho, RCM_info.size_years);
                RCM_data.yearly_maxind=NaN(RCM_grid.size_lat_rho, RCM_info.size_years);
                
%                 RCM_data.all=NaN(RCM_grid.size_lon_rho, RCM_grid.size_lat_rho, RCM_grid.size_depth, RCM_info.size_t);
%                 RCM_data.clim_mean=(zeros([RCM_grid.size_lon_rho,RCM_grid.size_lat_rho,RCM_grid.size_depth,RCM_info.size_months]));
%                 RCM_data.yearly_mean=(zeros([RCM_grid.size_lon_rho,RCM_grid.size_lat_rho,RCM_grid.size_depth,RCM_info.size_years]));
%                 RCM_data.maxval=NaN(RCM_grid.size_lat_rho, RCM_info.size_t);
%                 RCM_data.maxind=NaN(RCM_grid.size_lat_rho, RCM_info.size_t);
%                 RCM_data.yearly_maxval=NaN(RCM_grid.size_lat_rho, RCM_info.size_years);
%                 RCM_data.yearly_maxind=NaN(RCM_grid.size_lat_rho, RCM_info.size_years);
             
            end
            RCM_data.all(:,:,tmp.ind_for) = single(RCM_data.tmp_data);
            RCM_data.yearly_mean(:,:,yearij)=RCM_data.yearly_mean(:,:,yearij) + ...
                RCM_data.tmp_data / double(RCM_info.size_months);
            
%             RCM_data.all(:,:,:,tmp.ind_for) = single(RCM_data.tmp_data);
%             RCM_data.clim_mean(:,:,:,dayij)=RCM_data.clim_mean(:,:,:,dayij) + ...
%                 RCM_data.tmp_data / double(RCM_info.size_years);
%             RCM_data.yearly_mean(:,:,:,yearij)=RCM_data.yearly_mean(:,:,:,yearij) + ...
%                 RCM_data.tmp_data / double(RCM_info.size_months);
            RCM_data.tmp_surface_data=RCM_data.tmp_data(:,:,1);
            switch RCM_info.regionname
                case 'EKWC2'
                    [RCM_data.maxval(:,tmp.ind_for), RCM_data.maxind(:,tmp.ind_for)]=max(RCM_data.tmp_surface_data(1:tmp.EKWC_lon_ind_max,:),[],1); 
                otherwise
                    [RCM_data.maxval(:,tmp.ind_for), RCM_data.maxind(:,tmp.ind_for)]=max(RCM_data.tmp_surface_data,[],1);
            end
            
            tmp.ind_for = tmp.ind_for + 1;
            
            tmptmp=mean(mean(RCM_data.yearly_mean,1,'omitnan'), 2,'omitnan');
            plot(squeeze(tmptmp))
        end % end for day
%         RCM_data.tmp_yearly_surface_data=RCM_data.yearly_mean(:,:,1,yearij);
%         switch RCM_info.regionname
%             case 'EKWC2'
%                 [RCM_data.yearly_maxval(:,yearij), RCM_data.yearly_maxind(:,yearij)]=max(RCM_data.tmp_yearly_surface_data(1:tmp.EKWC_lon_ind_max,:),[],1);                
%             otherwise
%                 [RCM_data.yearly_maxval(:,yearij), RCM_data.yearly_maxind(:,yearij)]=max(RCM_data.tmp_yearly_surface_data,[],1);
%         end

    end % end for year

    save(RCM_info.matname, 'RCM_info', 'RCM_grid', 'RCM_data', '-v7.3');
    
end % end if (exist(RCM_info.matname , 'file') ~= 2 || flags.fig_tmp==2)  
disp(' ')