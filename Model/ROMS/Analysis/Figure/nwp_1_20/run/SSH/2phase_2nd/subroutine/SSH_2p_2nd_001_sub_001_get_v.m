% %         get model data and cmems data      
% %  Updated 13-Jan-2022 by Yong-Yub Kim, create

disp('subroutine v_2p_2nd_001_sub_001_get_v') 
disp([RCM_info.testname, ', ', RCM_info.regionname, ', ', num2str(min(RCM_info.years),'%04i'),' ~ ',num2str(max(RCM_info.years),'%04i')])
%% parameters, filename, tic set
tmp.lap_time_j=tic;
run(tmp.param_script);
tmp.ind_for=1;
RCM_info.matname = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_RCM_', tmp.variable, '_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '_', ...
    RCM_info.season, '.mat'];
% RCM_info.matname_interped = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_RCM_cmems_interped_v_', ...
%     num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
GCM_info.matname = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_GCM_', tmp.variable, '_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '_', ...
    RCM_info.season, '.mat'];
% GCM_info.matname_interped = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_GCM_cmems_interped_v_', ...
%     num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
% CMEMS_info.matname = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_CMEMS_v_', ...
%     num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];

tmp.elapsed=toc(tmp.lap_time_j);
tmp.totlap= (length(RCM_info.years)*length(RCM_info.months));
tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-1+1)/1);
tmp.EKWC_lon_ind_max = 70;

%% get data loop (year, month)
if (exist(RCM_info.matname , 'file') ~= 2 || flags.fig_tmp==2)
    RCM_info.size_t=length(RCM_info.years)*length(RCM_info.months);
    RCM_info.size_years=length(RCM_info.years);
    RCM_info.size_months=length(RCM_info.months);
    for yearij = 1:length(RCM_info.years)
        for monthij = 1:length(RCM_info.months)
            
%% print present time
            tmp.elapsed=toc(tmp.lap_time_j);
            tmp.templap= (yearij-1)*length(RCM_info.months)+monthij-1;
            fprintf(repmat('\b',1,sum(tmp.nchar)))  % remove printed time
            tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-tmp.templap+1)/tmp.templap);

            tmp.year = RCM_info.years(yearij);
            tmp.month = RCM_info.months(monthij);
            
%% set data filename            
            switch(RCM_info.testname)
                otherwise
                    RCM_info.filename = strcat(RCM_info.filedir, num2str(tmp.year,'%04i'), tmp.fs, ...
                    RCM_info.testname, '_monthly_std_', num2str(tmp.year,'%04i'), ...
                    '_', num2str(tmp.month,'%02i'), '.nc');
            end
            
%% get RCM grids
            if (yearij == 1 && monthij == 1)
%                 RCM_grid.filename = RCM_info.filename;
                RCM_grid.('lon_rho')= ncread(RCM_info.filename, 'lon_rho');
                RCM_grid.('lat_rho')= ncread(RCM_info.filename, 'lat_rho');
                RCM_grid.('depth')= ncread(RCM_info.filename, 'depth');
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
                    [RCM_grid.size_lon_rho RCM_grid.size_lat_rho, inf, 1]);
            RCM_data.tmp_data=RCM_data.tmp_data.*RCM_grid.mask_region; % [x y z] .* [x y]
            RCM_grid.size_depth = size(RCM_data.tmp_data,3);
            if (isfield(RCM_data, 'all') ~= 1)
                RCM_data.all=NaN(RCM_grid.size_lon_rho, RCM_grid.size_lat_rho, RCM_grid.size_depth, RCM_info.size_t);
                RCM_data.clim_mean=(zeros([RCM_grid.size_lon_rho,RCM_grid.size_lat_rho,RCM_grid.size_depth,RCM_info.size_months]));
                RCM_data.yearly_mean=(zeros([RCM_grid.size_lon_rho,RCM_grid.size_lat_rho,RCM_grid.size_depth,RCM_info.size_years]));
                RCM_data.maxval=NaN(RCM_grid.size_lat_rho, RCM_info.size_t);
                RCM_data.maxind=NaN(RCM_grid.size_lat_rho, RCM_info.size_t);
                RCM_data.yearly_maxval=NaN(RCM_grid.size_lat_rho, RCM_info.size_years);
                RCM_data.yearly_maxind=NaN(RCM_grid.size_lat_rho, RCM_info.size_years);
                
%                 RCM_data.all_v=NaN(RCM_grid.size_lon_rho, RCM_grid.size_lat_rho, RCM_grid.size_depth, RCM_info.size_t);
%                 RCM_data.clim_mean_v=(zeros([RCM_grid.size_lon_rho,RCM_grid.size_lat_rho,RCM_grid.size_depth,RCM_info.size_months]));
%                 RCM_data.yearly_mean_v=(zeros([RCM_grid.size_lon_rho,RCM_grid.size_lat_rho,RCM_grid.size_depth,RCM_info.size_years]));
%                 RCM_data.maxval_v=NaN(RCM_grid.size_lat_rho, RCM_info.size_t);
%                 RCM_data.maxind_v=NaN(RCM_grid.size_lat_rho, RCM_info.size_t);
%                 RCM_data.yearly_maxval_v=NaN(RCM_grid.size_lat_rho, RCM_info.size_years);
%                 RCM_data.yearly_maxind_v=NaN(RCM_grid.size_lat_rho, RCM_info.size_years);
                
%                 RCM_data.all_u=NaN(RCM_grid.size_lon_rho, RCM_grid.size_lat_rho, RCM_grid.size_depth, RCM_info.size_t);
%                 RCM_data.clim_mean_u=(zeros([RCM_grid.size_lon_rho,RCM_grid.size_lat_rho,RCM_grid.size_depth,RCM_info.size_months]));
%                 RCM_data.yearly_mean_u=(zeros([RCM_grid.size_lon_rho,RCM_grid.size_lat_rho,RCM_grid.size_depth,RCM_info.size_years]));
%                 RCM_data.maxval_u=NaN(RCM_grid.size_lat_rho, RCM_info.size_t);
%                 RCM_data.maxind_u=NaN(RCM_grid.size_lat_rho, RCM_info.size_t);
%                 RCM_data.yearly_maxval_u=NaN(RCM_grid.size_lat_rho, RCM_info.size_years);
%                 RCM_data.yearly_maxind_u=NaN(RCM_grid.size_lat_rho, RCM_info.size_years);
                
                RCM_data.all_transport_w=NaN(RCM_info.size_t,1);
                RCM_data.yearly_mean_transport_w=zeros(RCM_info.size_years,1);
                RCM_data.clim_mean_transport_w=zeros(RCM_info.size_months,1);
                RCM_data.all_transport_e=NaN(RCM_info.size_t,1);
                RCM_data.yearly_mean_transport_e=zeros(RCM_info.size_years,1);
                RCM_data.clim_mean_transport_e=zeros(RCM_info.size_months,1);
            end

            RCM_data.all(:,:,:,tmp.ind_for) = single(RCM_data.tmp_data);
            RCM_data.clim_mean(:,:,:,monthij)=RCM_data.clim_mean(:,:,:,monthij) + ...
                RCM_data.tmp_data / double(RCM_info.size_years);
            RCM_data.yearly_mean(:,:,:,yearij)=RCM_data.yearly_mean(:,:,:,yearij) + ...
                RCM_data.tmp_data / double(RCM_info.size_months);
            RCM_data.tmp_surface_data=RCM_data.tmp_data(:,:,1);
            switch RCM_info.regionname
                case 'EKWC2'
                    [RCM_data.maxval(:,tmp.ind_for), RCM_data.maxind(:,tmp.ind_for)]=max(RCM_data.tmp_surface_data(1:tmp.EKWC_lon_ind_max,:),[],1);                
                otherwise
                    [RCM_data.maxval(:,tmp.ind_for), RCM_data.maxind(:,tmp.ind_for)]=max(RCM_data.tmp_surface_data,[],1);
            end
            
%% get RCM transport across the western Korea Strait [128.75  35.16;   129.35  34.58]
            RCM_info.transfilename = [RCM_info.transfiledir, 'nwp_1_20_monthly_', RCM_info.testname, '_', num2str(tmp.year, '%04i'), '.txt'];

            tmp.startRow = 2;
            tmp.formatSpec = '%8f%9f%9f%9f%9f%9f%9f%9f%s%[^\n\r]';
            tmp.fileID = fopen(RCM_info.transfilename,'r');
            tmp.dataArray = textscan(tmp.fileID, tmp.formatSpec, 'Delimiter', '', ...
                'WhiteSpace', '', 'HeaderLines' , tmp.startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            tmp.dataArray{9} = strtrim(tmp.dataArray{9});
            fclose(tmp.fileID);
            tmp.transport_w = tmp.dataArray{:, 7};
            RCM_data.all_transport_w(tmp.ind_for)=tmp.transport_w(tmp.month);   
            RCM_data.clim_mean_transport_w(monthij)=RCM_data.clim_mean_transport_w(monthij) + ...
                tmp.transport_w(tmp.month) / double(RCM_info.size_years);
            RCM_data.yearly_mean_transport_w(yearij)=RCM_data.yearly_mean_transport_w(yearij) + ...
                tmp.transport_w(tmp.month) /  double(RCM_info.size_months);
%% get RCM transport across the eastern Korea Strait [129.35  34.58;   130.40  33.53]            
            tmp.transport_e = tmp.dataArray{:, 8};
            RCM_data.all_transport_e(tmp.ind_for)=tmp.transport_e(tmp.month);   
            RCM_data.clim_mean_transport_e(monthij)=RCM_data.clim_mean_transport_e(monthij) + ...
                tmp.transport_e(tmp.month) / double(RCM_info.size_years);
            RCM_data.yearly_mean_transport_e(yearij)=RCM_data.yearly_mean_transport_e(yearij) + ...
                tmp.transport_e(tmp.month) /  double(RCM_info.size_months);
            
%% get GCM grid
            GCM_info.filename = strcat(GCM_info.filedir, ...
                    tmp.variable_GCM, '_Omon_', GCM_info.scenario, '_', GCM_info.testname, '_', ...
                    num2str(tmp.year,'%04i'), '.nc');
            switch(GCM_info.testname)
                case{'CNRM-ESM2-1', 'CNRM-CM6-1-HR'}
                    GCM_grid.lonname='lon';
                    GCM_grid.latname='lat';
                    GCM_grid.xname = 'x';
                    GCM_grid.yname = 'y';
                case{'EC-Earth3-Veg', 'ACCESS-CM2', 'CMCC-ESM2'}
                    GCM_grid.lonname='longitude';
                    GCM_grid.latname='latitude';
                    GCM_grid.xname = 'i';
                    GCM_grid.yname = 'j';
            end
            if (yearij == 1 && monthij == 1)
                GCM_grid.filename_lon = [GCM_info.filedir, ...
                    'lon', '_Omon_', GCM_info.scenario, '_', GCM_info.testname, '.nc'];
                GCM_grid.filename_lat = [GCM_info.filedir, ...
                    'lat', '_Omon_', GCM_info.scenario, '_', GCM_info.testname, '.nc'];
                GCM_grid.domain = RCM_grid.domain;
                GCM_grid.refpolygon=RCM_grid.refpolygon;
                GCM_grid.lon_whole = ncread(GCM_grid.filename_lon,GCM_grid.lonname,[1 1],[inf,1]);
                GCM_grid.lat_whole = ncread(GCM_grid.filename_lat,GCM_grid.latname,[1 1],[1,inf]);
                [GCM_grid.ind_w, GCM_grid.ind_e, GCM_grid.ind_s, GCM_grid.ind_n] = ...
                    Func_0012_findind_Y(GCM_grid.dl, GCM_grid.domain, GCM_grid.lon_whole, GCM_grid.lat_whole);
                GCM_grid.size_lon = GCM_grid.ind_e-GCM_grid.ind_w+1;
                GCM_grid.size_lat = GCM_grid.ind_n-GCM_grid.ind_s+1;
                
                GCM_grid.lon= ncread(GCM_grid.filename_lon,GCM_grid.lonname, [GCM_grid.ind_w GCM_grid.ind_s], ...
                    [GCM_grid.size_lon GCM_grid.size_lat]);
                GCM_grid.lat = ncread(GCM_grid.filename_lat,GCM_grid.latname, [GCM_grid.ind_w GCM_grid.ind_s], ...
                    [GCM_grid.size_lon GCM_grid.size_lat]);
                
                GCM_grid.mask_region = double(inpolygon(GCM_grid.lon, GCM_grid.lat, ...
                    GCM_grid.refpolygon(:,1),GCM_grid.refpolygon(:,2)));
                GCM_grid.mask_region(GCM_grid.mask_region==0)=NaN;
                GCM_info.data_info = ncinfo(GCM_info.filename, tmp.variable_GCM);
                GCM_grid.mask_ocean= ncread(GCM_info.filename, tmp.variable_GCM, [GCM_grid.ind_w GCM_grid.ind_s, 1, 1], ...
                    [GCM_grid.size_lon GCM_grid.size_lat, inf, 1]);
                GCM_grid.mask_ocean(isfinite(GCM_grid.mask_ocean))=1;
                GCM_grid.mask_ocean=GCM_grid.mask_ocean .*  GCM_grid.mask_region;
                GCM_grid.mask_land = GCM_grid.mask_ocean;
                GCM_grid.mask_land(isfinite(GCM_grid.mask_ocean))=NaN;
                GCM_grid.mask_land(isnan(GCM_grid.mask_ocean))=1;
            end    
                
%%             get GCM data
            GCM_data.tmp_data = ncread(GCM_info.filename, tmp.variable_GCM, [GCM_grid.ind_w GCM_grid.ind_s, 1, monthij], ...
                    [GCM_grid.size_lon GCM_grid.size_lat, inf, 1]);
            GCM_data.tmp_data=GCM_data.tmp_data.*GCM_grid.mask_region;
            GCM_grid.size_depth = size(GCM_data.tmp_data,3);

            if (isfield(GCM_data, 'all') ~= 1)
                GCM_data.all=NaN(GCM_grid.size_lon, GCM_grid.size_lat, GCM_grid.size_depth, RCM_info.size_t);
                GCM_data.clim_mean=(zeros([GCM_grid.size_lon,GCM_grid.size_lat,GCM_grid.size_depth,RCM_info.size_months]));
                GCM_data.yearly_mean=(zeros([GCM_grid.size_lon,GCM_grid.size_lat,GCM_grid.size_depth,RCM_info.size_years]));         
            else
                GCM_data.all(:,:,:,tmp.ind_for) = single(GCM_data.tmp_data);
                GCM_data.clim_mean(:,:,:,monthij)=GCM_data.clim_mean(:,:,:,monthij) + ...
                    GCM_data.tmp_data / double(length(GCM_info.years));
                GCM_data.yearly_mean(:,:,:,yearij)=GCM_data.yearly_mean(:,:,:,yearij) + ...
                    GCM_data.tmp_data / double(length(GCM_info.months));
            end      
    
            
            tmp.ind_for = tmp.ind_for + 1;
            
            
        end % end for month
        RCM_data.tmp_yearly_surface_data=RCM_data.yearly_mean(:,:,1,yearij);
        switch RCM_info.regionname
            case 'EKWC2'
                [RCM_data.yearly_maxval(:,yearij), RCM_data.yearly_maxind(:,yearij)]=max(RCM_data.tmp_yearly_surface_data(1:tmp.EKWC_lon_ind_max,:),[],1);                
            otherwise
                [RCM_data.yearly_maxval(:,yearij), RCM_data.yearly_maxind(:,yearij)]=max(RCM_data.tmp_yearly_surface_data,[],1);
        end

    end % end for year

     
    save(RCM_info.matname, 'RCM_info', 'RCM_grid', 'RCM_data', '-v7.3');
%     save(RCM_info.matname_interped, 'RCM_info', 'RCM_grid', 'RCM_data_interped', 'CMEMS_grid', '-v7.3');
    save(GCM_info.matname, 'GCM_info', 'GCM_grid', 'GCM_data', '-v7.3');
%     save(GCM_info.matname_interped, 'GCM_info', 'GCM_grid', 'GCM_data_interped', 'CMEMS_grid', '-v7.3');
end % end if (exist(RCM_info.matname , 'file') ~= 2 || flags.fig_tmp==2)  
disp(' ')