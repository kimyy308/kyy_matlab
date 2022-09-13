% %  get model wind data (input)
% %  Updated 17-Jan-2021 by Yong-Yub Kim

RCM_info_wind=RCM_info;
RCM_grid_wind=RCM_grid;
disp('subroutine SSH_2p_2nd_001_sub_002_get_wind') 
disp([RCM_info_wind.testname, ', ', RCM_info_wind.regionname, ', ', num2str(min(RCM_info_wind.years),'%04i'),' ~ ',num2str(max(RCM_info_wind.years),'%04i')])

%% parameters, filename, tic set
tmp.lap_time_j=tic;
run(tmp.param_script);
RCM_info_wind.matname = [RCM_info_wind.windsavedir,RCM_info_wind.testname,'_',RCM_info_wind.regionname, '_RCM_data_wind_', ...
    num2str(min(RCM_info_wind.years),'%04i'),'_',num2str(max(RCM_info_wind.years),'%04i'), '_', ...
    RCM_info.season, '.mat'];

tmp.elapsed=toc(tmp.lap_time_j);
tmp.totlap= (length(RCM_info_wind.years)*length(RCM_info_wind.months));
tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-1+1)/1);

%% get data loop (year, month)
if (exist(RCM_info_wind.matname , 'file') ~= 2 || flags.fig_tmp==2)
    RCM_info_wind.size_t=length(RCM_info_wind.years)*length(RCM_info_wind.months);
    RCM_info_wind.size_years=length(RCM_info_wind.years);
    RCM_info_wind.size_months=length(RCM_info_wind.months);
    tmp.wind_vars={'Uwind', 'Vwind'};
    for windij = 1:length(tmp.wind_vars)
        tmp.ind_for=1;
        tmp.variable_wind=tmp.wind_vars{windij};
        for yearij = 1:length(RCM_info_wind.years)
            for monthij = 1:length(RCM_info_wind.months)
    %% print present time
                tmp.elapsed=toc(tmp.lap_time_j);
                tmp.templap= (yearij-1)*length(RCM_info_wind.months)+monthij-1;
                fprintf(repmat('\b',1,sum(tmp.nchar)))  % remove printed time
                tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-tmp.templap+1)/tmp.templap);

                tmp.year = RCM_info_wind.years(yearij);
                tmp.month = RCM_info_wind.months(monthij);

    %% set data filename                        
                switch(RCM_info_wind.testname)
                    otherwise
                        RCM_info_wind.filename = strcat(RCM_info_wind.atmfiledir, tmp.fs, tmp.variable_wind, tmp.fs, ...
                            num2str(tmp.year,'%04i'), tmp.fs, 'NWP_pck_', RCM_info_wind.testname, ...
                            '_', tmp.variable_wind, '_monthly_', num2str(tmp.year,'%04i'), ...
                            '_', num2str(tmp.month,'%02i'), '.nc');
                end

    %% get RCM grids
                if (windij ==1 && yearij == 1 && monthij == 1)
                    RCM_grid_wind.filename_lon_rho = [RCM_info_wind.atmroot, 'NWP_pck_ocean_lon_rho_NWP.nc'];
                    RCM_grid_wind.filename_lon_u = [RCM_info_wind.atmroot, 'NWP_pck_ocean_lon_u_NWP.nc'];
                    RCM_grid_wind.filename_lon_v = [RCM_info_wind.atmroot, 'NWP_pck_ocean_lon_v_NWP.nc'];
                    RCM_grid_wind.filename_lat_rho = [RCM_info_wind.atmroot, 'NWP_pck_ocean_lat_rho_NWP.nc'];
                    RCM_grid_wind.filename_lat_u = [RCM_info_wind.atmroot, 'NWP_pck_ocean_lat_u_NWP.nc'];
                    RCM_grid_wind.filename_lat_v = [RCM_info_wind.atmroot, 'NWP_pck_ocean_lat_v_NWP.nc'];

                    RCM_grid_wind.lon_rho_whole = ncread(RCM_grid_wind.filename_lon_rho,'lon_rho',[1 1],[inf,1]);
                    RCM_grid_wind.lat_rho_whole = ncread(RCM_grid_wind.filename_lat_rho,'lat_rho',[1 1],[1,inf]);
                    [RCM_grid_wind.ind_w, RCM_grid_wind.ind_e, RCM_grid_wind.ind_s, RCM_grid_wind.ind_n] = ...
                        Func_0012_findind_Y(RCM_grid_wind.dl, RCM_grid_wind.domain, RCM_grid_wind.lon_rho_whole, RCM_grid_wind.lat_rho_whole);
                    RCM_grid_wind.size_lon_rho = RCM_grid_wind.ind_e-RCM_grid_wind.ind_w+1;
                    RCM_grid_wind.size_lat_rho = RCM_grid_wind.ind_n-RCM_grid_wind.ind_s+1;
                    RCM_grid_wind.size_lon_u = RCM_grid_wind.ind_e-RCM_grid_wind.ind_w;
                    RCM_grid_wind.size_lat_u = RCM_grid_wind.ind_n-RCM_grid_wind.ind_s+1;
                    RCM_grid_wind.size_lon_v = RCM_grid_wind.ind_e-RCM_grid_wind.ind_w+1;
                    RCM_grid_wind.size_lat_v = RCM_grid_wind.ind_n-RCM_grid_wind.ind_s;

                    RCM_grid_wind.cut_lon_rho = ncread(RCM_grid_wind.filename_lon_rho,'lon_rho', [RCM_grid_wind.ind_w RCM_grid_wind.ind_s], ...
                        [RCM_grid_wind.size_lon_rho RCM_grid_wind.size_lat_rho]);
                    RCM_grid_wind.cut_lat_rho = ncread(RCM_grid_wind.filename_lat_rho,'lat_rho', [RCM_grid_wind.ind_w RCM_grid_wind.ind_s], ...
                        [RCM_grid_wind.size_lon_rho RCM_grid_wind.size_lat_rho]);
                    RCM_grid_wind.cut_lon_u = ncread(RCM_grid_wind.filename_lon_u,'lon_u', [RCM_grid_wind.ind_w RCM_grid_wind.ind_s], ...
                        [RCM_grid_wind.size_lon_u RCM_grid_wind.size_lat_u]);
                    RCM_grid_wind.cut_lat_u = ncread(RCM_grid_wind.filename_lat_u,'lat_u', [RCM_grid_wind.ind_w RCM_grid_wind.ind_s], ...
                        [RCM_grid_wind.size_lon_u RCM_grid_wind.size_lat_u]);
                    RCM_grid_wind.cut_lon_v = ncread(RCM_grid_wind.filename_lon_v,'lon_v', [RCM_grid_wind.ind_w RCM_grid_wind.ind_s], ...
                        [RCM_grid_wind.size_lon_v RCM_grid_wind.size_lat_v]);
                    RCM_grid_wind.cut_lat_v = ncread(RCM_grid_wind.filename_lat_v,'lat_v', [RCM_grid_wind.ind_w RCM_grid_wind.ind_s], ...
                        [RCM_grid_wind.size_lon_v RCM_grid_wind.size_lat_v]);

                    RCM_grid_wind.mask_region = double(inpolygon(RCM_grid_wind.cut_lon_rho, RCM_grid_wind.cut_lat_rho, ...
                        RCM_grid_wind.refpolygon(:,1),RCM_grid_wind.refpolygon(:,2)));
                    RCM_grid_wind.mask_region(RCM_grid_wind.mask_region==0)=NaN;
                    RCM_info_wind.data_info = ncinfo(RCM_info_wind.filename, tmp.variable_wind);
                    RCM_grid_wind.mask_ocean= ncread(RCM_info_wind.filename, tmp.variable_wind, [RCM_grid_wind.ind_w RCM_grid_wind.ind_s, 1], ...
                        [RCM_grid_wind.size_lon_rho RCM_grid_wind.size_lat_rho, 1]);
                    RCM_grid_wind.mask_ocean(isfinite(RCM_grid_wind.mask_ocean))=1;
                    RCM_grid_wind.mask_ocean = RCM_grid_wind.mask_ocean .*  RCM_grid_wind.mask_region;
                    RCM_grid_wind.mask_land = RCM_grid_wind.mask_ocean;
                    RCM_grid_wind.mask_land(isfinite(RCM_grid_wind.mask_ocean))=NaN;
                    RCM_grid_wind.mask_land(isnan(RCM_grid_wind.mask_ocean))=1;
                end

    % %             get RCM data
                RCM_data_wind.tmp_data = ncread(RCM_info_wind.filename, tmp.variable_wind, [RCM_grid_wind.ind_w RCM_grid_wind.ind_s, 1], ...
                        [RCM_grid_wind.size_lon_rho RCM_grid_wind.size_lat_rho, 1]);
                RCM_data_wind.tmp_data=RCM_data_wind.tmp_data.*RCM_grid_wind.mask_region;

                if (yearij == 1 && monthij == 1)
                    RCM_data.(['all_', tmp.variable_wind])=NaN(RCM_grid_wind.size_lon_rho, RCM_grid_wind.size_lat_rho, RCM_info_wind.size_t);
                    RCM_data_wind.(['clim_mean_', tmp.variable_wind])=(zeros([RCM_grid_wind.size_lon_rho,RCM_grid_wind.size_lat_rho,RCM_info_wind.size_months]));
                    RCM_data_wind.(['yearly_mean_', tmp.variable_wind])=(zeros([RCM_grid_wind.size_lon_rho,RCM_grid_wind.size_lat_rho,RCM_info_wind.size_years]));
                end
                RCM_data_wind.(['all_', tmp.variable_wind])(:,:,tmp.ind_for) = single(RCM_data_wind.tmp_data);
                RCM_data_wind.(['clim_mean_', tmp.variable_wind])(:,:,monthij)=RCM_data_wind.(['clim_mean_', tmp.variable_wind])(:,:,monthij) + ...
                    RCM_data_wind.tmp_data / double(length(RCM_info_wind.years));
                RCM_data_wind.(['yearly_mean_', tmp.variable_wind])(:,:,yearij)=RCM_data_wind.(['yearly_mean_', tmp.variable_wind])(:,:,yearij) + ...
                    RCM_data_wind.tmp_data / double(length(RCM_info_wind.months));

                tmp.ind_for = tmp.ind_for + 1;


            end
        end
    end
    save(RCM_info_wind.matname, 'RCM_info_wind', 'RCM_grid_wind', 'RCM_data_wind', '-v7.3');
end
disp(' ')