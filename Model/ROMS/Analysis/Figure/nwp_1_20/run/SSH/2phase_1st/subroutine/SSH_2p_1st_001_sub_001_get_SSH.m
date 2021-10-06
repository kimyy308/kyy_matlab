% %         get model data and cmems data      
% %  Updated 05-Jul-2021 by Yong-Yub Kim, modified to Func_0012


disp('subroutine SSH_2p_1st_001_sub_001_get_SSH') 
disp([RCM_info.testname, ', ', RCM_info.regionname, ', ', num2str(min(RCM_info.years),'%04i'),' ~ ',num2str(max(RCM_info.years),'%04i')])

% % tmp.lap_time_j=tic;
% % jlen=10;
% % for j=1:jlen
% % tmp.elapsed=toc(tmp.lap_time_j);
% % tmp.totlap= (length(RCM_info.years)*length(RCM_info.months));
% % tmp.templap= yearij*length(RCM_info.months)+monthij;
% % tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-tmp.templap+1)/tmp.templap);
% % fprintf(repmat('\b',1,sum(tmp.nchar)))  % remove printed time
% % end

tmp.lap_time_j=tic;

run(tmp.param_script);
tmp.ind_for=1;
RCM_info.matname = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_RCM_ssh_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
RCM_info.matname_interped = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_RCM_cmems_interped_ssh_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
GCM_info.matname = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_GCM_ssh_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
GCM_info.matname_interped = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_GCM_cmems_interped_ssh_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
CMEMS_info.matname = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_CMEMS_ssh_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];

tmp.elapsed=toc(tmp.lap_time_j);
tmp.totlap= (length(RCM_info.years)*length(RCM_info.months));
tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-1+1)/1);


if (exist(RCM_info.matname , 'file') ~= 2 || flags.fig_switch(1)==2)   
    for yearij = 1:length(RCM_info.years)
        for monthij = 1:length(RCM_info.months)
            tmp.elapsed=toc(tmp.lap_time_j);
            tmp.templap= yearij*length(RCM_info.months)+monthij;
            fprintf(repmat('\b',1,sum(tmp.nchar)))  % remove printed time
            tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-tmp.templap+1)/tmp.templap);
            %             disp([num2str(yearij), 'y_',num2str(monthij),'m'])


            tmp.year = RCM_info.years(yearij);
            tmp.month = RCM_info.months(monthij);
            % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
            RCM_info.filename = strcat(RCM_info.filedir, num2str(tmp.year,'%04i'), tmp.fs, ...
                    'pck_', RCM_info.testname, '_', tmp.variable, '_monthly_', num2str(tmp.year,'%04i'), ...
                    '_', num2str(tmp.month,'%02i'), '.nc');
% %             get RCM_grid
            if (yearij == 1 && monthij == 1)
                RCM_grid.filename_lon_rho = [RCM_info.dataroot, 'NWP_pck_ocean_lon_rho_NWP.nc'];
                RCM_grid.filename_lon_u = [RCM_info.dataroot, 'NWP_pck_ocean_lon_u_NWP.nc'];
                RCM_grid.filename_lon_v = [RCM_info.dataroot, 'NWP_pck_ocean_lon_v_NWP.nc'];
                RCM_grid.filename_lat_rho = [RCM_info.dataroot, 'NWP_pck_ocean_lat_rho_NWP.nc'];
                RCM_grid.filename_lat_u = [RCM_info.dataroot, 'NWP_pck_ocean_lat_u_NWP.nc'];
                RCM_grid.filename_lat_v = [RCM_info.dataroot, 'NWP_pck_ocean_lat_v_NWP.nc'];
                                
                RCM_grid.lon_rho_whole = ncread(RCM_grid.filename_lon_rho,'lon_rho',[1 1],[inf,1]);
                RCM_grid.lat_rho_whole = ncread(RCM_grid.filename_lat_rho,'lat_rho',[1 1],[1,inf]);
                [RCM_grid.ind_w, RCM_grid.ind_e, RCM_grid.ind_s, RCM_grid.ind_n] = ...
                    Func_0012_findind_Y(RCM_grid.dl, RCM_grid.domain, RCM_grid.lon_rho_whole, RCM_grid.lat_rho_whole);
                RCM_grid.size_lon_rho = RCM_grid.ind_e-RCM_grid.ind_w+1;
                RCM_grid.size_lat_rho = RCM_grid.ind_n-RCM_grid.ind_s+1;
                RCM_grid.size_lon_u = RCM_grid.ind_e-RCM_grid.ind_w;
                RCM_grid.size_lat_u = RCM_grid.ind_n-RCM_grid.ind_s+1;
                RCM_grid.size_lon_v = RCM_grid.ind_e-RCM_grid.ind_w+1;
                RCM_grid.size_lat_v = RCM_grid.ind_n-RCM_grid.ind_s;
                
                RCM_grid.lon_rho = ncread(RCM_grid.filename_lon_rho,'lon_rho', [RCM_grid.ind_w RCM_grid.ind_s], ...
                    [RCM_grid.size_lon_rho RCM_grid.size_lat_rho]);
                RCM_grid.lat_rho = ncread(RCM_grid.filename_lat_rho,'lat_rho', [RCM_grid.ind_w RCM_grid.ind_s], ...
                    [RCM_grid.size_lon_rho RCM_grid.size_lat_rho]);
                RCM_grid.lon_u = ncread(RCM_grid.filename_lon_u,'lon_u', [RCM_grid.ind_w RCM_grid.ind_s], ...
                    [RCM_grid.size_lon_u RCM_grid.size_lat_u]);
                RCM_grid.lat_u = ncread(RCM_grid.filename_lat_u,'lat_u', [RCM_grid.ind_w RCM_grid.ind_s], ...
                    [RCM_grid.size_lon_u RCM_grid.size_lat_u]);
                RCM_grid.lon_v = ncread(RCM_grid.filename_lon_v,'lon_v', [RCM_grid.ind_w RCM_grid.ind_s], ...
                    [RCM_grid.size_lon_v RCM_grid.size_lat_v]);
                RCM_grid.lat_v = ncread(RCM_grid.filename_lat_v,'lat_v', [RCM_grid.ind_w RCM_grid.ind_s], ...
                    [RCM_grid.size_lon_v RCM_grid.size_lat_v]);
                
                RCM_grid.mask_region = double(inpolygon(RCM_grid.lon_rho, RCM_grid.lat_rho, ...
                    RCM_grid.refpolygon(:,1),RCM_grid.refpolygon(:,2)));
                RCM_grid.mask_region(RCM_grid.mask_region==0)=NaN;
                RCM_info.data_info = ncinfo(RCM_info.filename, param.varname);
                RCM_grid.mask_ocean= ncread(RCM_info.filename, tmp.variable, [RCM_grid.ind_w RCM_grid.ind_s, 1], ...
                    [RCM_grid.size_lon_rho RCM_grid.size_lat_rho, 1]);
                RCM_grid.mask_ocean(isfinite(RCM_grid.mask_ocean))=1;
                RCM_grid.mask_ocean = RCM_grid.mask_ocean .*  RCM_grid.mask_region;
                RCM_grid.mask_land = RCM_grid.mask_ocean;
                RCM_grid.mask_land(isfinite(RCM_grid.mask_ocean))=NaN;
                RCM_grid.mask_land(isnan(RCM_grid.mask_ocean))=1;
            end
            
% %             get RCM data
            RCM_data.tmp_data = ncread(RCM_info.filename, tmp.variable, [RCM_grid.ind_w RCM_grid.ind_s, 1], ...
                    [RCM_grid.size_lon_rho RCM_grid.size_lat_rho, 1]);
            RCM_data.tmp_data=RCM_data.tmp_data.*RCM_grid.mask_region;

            if (yearij == 1 && monthij == 1)
                RCM_data.clim_mean=(zeros([RCM_grid.size_lon_rho,RCM_grid.size_lat_rho,12]));
                RCM_data.annual_mean=(zeros([RCM_grid.size_lon_rho,RCM_grid.size_lat_rho,length(RCM_info.years)]));
            end
            RCM_data.all(:,:,tmp.ind_for) = single(RCM_data.tmp_data);
            RCM_data.clim_mean(:,:,monthij)=RCM_data.clim_mean(:,:,monthij) + ...
                RCM_data.tmp_data / double(length(RCM_info.years));
            RCM_data.annual_mean(:,:,yearij)=RCM_data.annual_mean(:,:,yearij) + ...
                RCM_data.tmp_data / double(length(RCM_info.months));

% % % % % get GCM grid
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
                    GCM_grid.lonname, '_Omon_', GCM_info.scenario, '_', GCM_info.testname, '.nc'];
                GCM_grid.filename_lat = [GCM_info.filedir, ...
                    GCM_grid.latname, '_Omon_', GCM_info.scenario, '_', GCM_info.testname, '.nc'];
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
                GCM_grid.mask_ocean= ncread(GCM_info.filename, tmp.variable_GCM, [GCM_grid.ind_w GCM_grid.ind_s, 1], ...
                    [GCM_grid.size_lon GCM_grid.size_lat, 1]);
                GCM_grid.mask_ocean(isfinite(GCM_grid.mask_ocean))=1;
                GCM_grid.mask_ocean=GCM_grid.mask_ocean .*  GCM_grid.mask_region;
                GCM_grid.mask_land = GCM_grid.mask_ocean;
                GCM_grid.mask_land(isfinite(GCM_grid.mask_ocean))=NaN;
                GCM_grid.mask_land(isnan(GCM_grid.mask_ocean))=1;
            end    
                
% %             get GCM data
            GCM_data.tmp_data = ncread(GCM_info.filename, tmp.variable_GCM, [GCM_grid.ind_w GCM_grid.ind_s, monthij], ...
                    [GCM_grid.size_lon GCM_grid.size_lat, 1]);
            GCM_data.tmp_data=GCM_data.tmp_data.*GCM_grid.mask_region;

            if (yearij == 1 && monthij == 1)
                GCM_data.clim_mean=(zeros([GCM_grid.size_lon,GCM_grid.size_lat,12]));
                GCM_data.annual_mean=(zeros([GCM_grid.size_lon,GCM_grid.size_lat,length(GCM_info.years)]));                
            end
            GCM_data.all(:,:,tmp.ind_for) = single(GCM_data.tmp_data);
            GCM_data.clim_mean(:,:,monthij)=GCM_data.clim_mean(:,:,monthij) + ...
                GCM_data.tmp_data / double(length(GCM_info.years));
            GCM_data.annual_mean(:,:,yearij)=GCM_data.annual_mean(:,:,yearij) + ...
                GCM_data.tmp_data / double(length(GCM_info.months));    
%         pcolor(GCM_grid.lon', GCM_grid.lat', GCM_data.tmp_data'); shading flat; colorbar; 
                

% %             get CMEMS data
            tmp.variable_CMEMS_daily=[tmp.variable_CMEMS, '_daily'];
            tmp.variable_CMEMS_all=[tmp.variable_CMEMS, '_all'];
            tmp.variable_CMEMS_clim_mean=[tmp.variable_CMEMS, '_clim_mean'];
            tmp.variable_CMEMS_annual_mean=[tmp.variable_CMEMS, '_annual_mean'];

            CMEMS_info.filename = strcat(CMEMS_info.filedir,'cmems_nwp_ssh_', num2str(tmp.year,'%04i'),'.nc');
            if (yearij == 1 && monthij == 1)
                
                CMEMS_info.file_info=ncinfo(CMEMS_info.filename);
                CMEMS_info.lonname='longitude';
                CMEMS_info.latname='latitude';
                tmp.variable_CMEMS='sla';
                CMEMS_info.lon_info=ncinfo(CMEMS_info.filename,CMEMS_info.lonname);
                CMEMS_info.lat_info=ncinfo(CMEMS_info.filename,CMEMS_info.latname);
                CMEMS_grid.domain = RCM_grid.domain;
                CMEMS_grid.lon = ncread(CMEMS_info.filename,CMEMS_info.lonname,1,CMEMS_info.lon_info.Dimensions.Length);
                CMEMS_grid.lat = ncread(CMEMS_info.filename,CMEMS_info.latname,1,CMEMS_info.lat_info.Dimensions.Length);
                [CMEMS_grid.lat2_whole, CMEMS_grid.lon2_whole]= meshgrid(CMEMS_grid.lat, CMEMS_grid.lon);
                [CMEMS_grid.ind_w, CMEMS_grid.ind_e, CMEMS_grid.ind_s, CMEMS_grid.ind_n] = ...
                    Func_0012_findind_Y(CMEMS_grid.dl, CMEMS_grid.domain, CMEMS_grid.lon2_whole, CMEMS_grid.lat2_whole);
                
                CMEMS_grid.lon2 = CMEMS_grid.lon2_whole(CMEMS_grid.ind_w:CMEMS_grid.ind_e, CMEMS_grid.ind_s:CMEMS_grid.ind_n);
                CMEMS_grid.lat2 = CMEMS_grid.lat2_whole(CMEMS_grid.ind_w:CMEMS_grid.ind_e, CMEMS_grid.ind_s:CMEMS_grid.ind_n);

                CMEMS_grid.size_lon=CMEMS_grid.ind_e-CMEMS_grid.ind_w+1;
                CMEMS_grid.size_lat=CMEMS_grid.ind_n-CMEMS_grid.ind_s+1;
                CMEMS_grid.refpolygon=RCM_grid.refpolygon;
                CMEMS_data.adt_clim_mean=(zeros([CMEMS_grid.size_lon, CMEMS_grid.size_lat,12]));
                CMEMS_data.(tmp.variable_CMEMS_clim_mean)=(zeros([CMEMS_grid.size_lon, CMEMS_grid.size_lat,12]));
                CMEMS_data.adt_annual_mean=(zeros([CMEMS_grid.size_lon, CMEMS_grid.size_lat,length(RCM_info.years)]));
                CMEMS_data.(tmp.variable_CMEMS_annual_mean)=(zeros([CMEMS_grid.size_lon, CMEMS_grid.size_lat,length(RCM_info.years)]));

                RCM_data_interped.clim_mean=(zeros([CMEMS_grid.size_lon, CMEMS_grid.size_lat,12]));
                GCM_data_interped.clim_mean=(zeros([CMEMS_grid.size_lon, CMEMS_grid.size_lat,12]));
                RCM_data_interped.annual_mean=(zeros([CMEMS_grid.size_lon, CMEMS_grid.size_lat,length(RCM_info.years)]));
                GCM_data_interped.annual_mean=(zeros([CMEMS_grid.size_lon, CMEMS_grid.size_lat,length(RCM_info.years)]));
                
                CMEMS_grid.mask_region = double(inpolygon(CMEMS_grid.lon2, CMEMS_grid.lat2, ...
                    CMEMS_grid.refpolygon(:,1),CMEMS_grid.refpolygon(:,2)));
                CMEMS_grid.mask_region(CMEMS_grid.mask_region==0)=NaN;
                CMEMS_grid.mask_ocean= ncread(CMEMS_info.filename, tmp.variable_CMEMS, [CMEMS_grid.ind_w CMEMS_grid.ind_s, 1], ...
                    [CMEMS_grid.size_lon CMEMS_grid.size_lat, 1]);
                CMEMS_grid.mask_ocean(isfinite(CMEMS_grid.mask_ocean))=1;
                CMEMS_grid.mask_ocean=CMEMS_grid.mask_ocean .* CMEMS_grid.mask_region;
                CMEMS_grid.mask_land = CMEMS_grid.mask_ocean;
                CMEMS_grid.mask_land(isfinite(CMEMS_grid.mask_ocean))=NaN;
                CMEMS_grid.mask_land(isnan(CMEMS_grid.mask_ocean))=1;
            end
            
            if tmp.month==1
                CMEMS_data.adt_daily=ncread(CMEMS_info.filename,'adt',[CMEMS_grid.ind_w CMEMS_grid.ind_s 1], [CMEMS_grid.size_lon CMEMS_grid.size_lat 31]);
                CMEMS_data.(tmp.variable_CMEMS_daily)=ncread(CMEMS_info.filename,tmp.variable_CMEMS,[CMEMS_grid.ind_w CMEMS_grid.ind_s 1], [CMEMS_grid.size_lon CMEMS_grid.size_lat 31]);
            else
                CMEMS_data.adt_daily=ncread(CMEMS_info.filename,'adt', ...
                    [CMEMS_grid.ind_w CMEMS_grid.ind_s sum(eomday(tmp.year,1:tmp.month-1))+1], ...
                    [CMEMS_grid.size_lon CMEMS_grid.size_lat eomday(tmp.year,tmp.month)]);
                CMEMS_data.(tmp.variable_CMEMS_daily)=ncread(CMEMS_info.filename,tmp.variable_CMEMS, ...
                    [CMEMS_grid.ind_w CMEMS_grid.ind_s sum(eomday(tmp.year,1:tmp.month-1))+1], ...
                    [CMEMS_grid.size_lon CMEMS_grid.size_lat eomday(tmp.year,tmp.month)]);
            end
            CMEMS_data.adt = mean(CMEMS_data.adt_daily,3,'omitnan');
            CMEMS_data.adt = CMEMS_data.adt .* CMEMS_grid.mask_region;
            CMEMS_data.(tmp.variable_CMEMS) = mean(CMEMS_data.(tmp.variable_CMEMS_daily), 3, 'omitnan');
            CMEMS_data.(tmp.variable_CMEMS)=CMEMS_data.(tmp.variable_CMEMS) .* CMEMS_grid.mask_region;

%             cmems_data = ncread(CMEMS_info.filename,varname,[cmems_lon_min(1) cmems_lat_min(1) tmp.month], [cmems_lon_max(1)-cmems_lon_min(1) cmems_lat_max(1)-cmems_lat_min(1) 1]);
%             cmems_data(cmems_data<-1000)=NaN;
%             cmems_data(cmems_data>1000)=NaN;

            CMEMS_data.adt_all(:,:,tmp.ind_for) = CMEMS_data.adt;
            CMEMS_data.(tmp.variable_CMEMS_all)(:,:,tmp.ind_for) = CMEMS_data.(tmp.variable_CMEMS);
            
            CMEMS_data.adt_clim_mean(:,:,monthij)= ...
                CMEMS_data.adt_clim_mean(:,:,monthij) + CMEMS_data.adt / double(length(RCM_info.years));
            CMEMS_data.(tmp.variable_CMEMS_clim_mean)(:,:,monthij) = ...
                CMEMS_data.(tmp.variable_CMEMS_clim_mean)(:,:,monthij) + CMEMS_data.(tmp.variable_CMEMS)/double(length(RCM_info.years));
            CMEMS_data.adt_annual_mean(:,:,yearij)=CMEMS_data.adt_annual_mean(:,:,yearij) + ...
                CMEMS_data.adt / double(length(RCM_info.months));
            CMEMS_data.(tmp.variable_CMEMS_annual_mean)(:,:,yearij)=CMEMS_data.(tmp.variable_CMEMS_annual_mean)(:,:,yearij) + ...
                CMEMS_data.(tmp.variable_CMEMS) / double(length(RCM_info.months));
            
            RCM_data_interped.tmp_data = ...
                griddata(double(RCM_grid.lon_rho), double(RCM_grid.lat_rho), RCM_data.tmp_data, ...
                            double(CMEMS_grid.lon2),double(CMEMS_grid.lat2)); 
            GCM_data_interped.tmp_data = ...
                griddata(double(GCM_grid.lon), double(GCM_grid.lat), GCM_data.tmp_data, ...
                            double(CMEMS_grid.lon2),double(CMEMS_grid.lat2)); 
            
            RCM_data_interped.all(:,:,tmp.ind_for) = RCM_data_interped.tmp_data;
            GCM_data_interped.all(:,:,tmp.ind_for) = GCM_data_interped.tmp_data;

            RCM_data_interped.clim_mean(:,:,monthij) = ...
                RCM_data_interped.clim_mean(:,:,monthij) + RCM_data_interped.tmp_data/double(length(RCM_info.years));
            GCM_data_interped.clim_mean(:,:,monthij) = ...
                GCM_data_interped.clim_mean(:,:,monthij) + GCM_data_interped.tmp_data/double(length(RCM_info.years));
            
            RCM_data_interped.annual_mean(:,:,yearij) = ...
                RCM_data_interped.annual_mean(:,:,yearij) + RCM_data_interped.tmp_data/double(length(RCM_info.months));
            GCM_data_interped.annual_mean(:,:,yearij) = ...
                GCM_data_interped.annual_mean(:,:,yearij) + GCM_data_interped.tmp_data/double(length(RCM_info.months));
            
            
            if (yearij == 1 && monthij == 1)
                RCM_grid.mask_ocean_interped = NaN(size(RCM_data_interped.tmp_data));
                RCM_grid.mask_ocean_interped(isfinite(RCM_data_interped.tmp_data))=1;
                RCM_grid.mask_land_interped = NaN(size(RCM_grid.mask_ocean_interped));
                RCM_grid.mask_land_interped(isnan(RCM_grid.mask_ocean_interped))=1;
                GCM_grid.mask_ocean_interped = NaN(size(GCM_data_interped.tmp_data));
                GCM_grid.mask_ocean_interped(isfinite(GCM_data_interped.tmp_data))=1;
                GCM_grid.mask_land_interped = NaN(size(GCM_grid.mask_ocean_interped));
                GCM_grid.mask_land_interped(isnan(GCM_grid.mask_ocean_interped))=1;
            end
            
            tmp.ind_for = tmp.ind_for + 1;
        end
    end
%     save(matname, 'RCM_grid.mask_region', 'len_lon_model', 'len_lat_model', 'comb_spatial_meanmodel', 'comb_data', ...
%         'cmems_lonsize_cut', 'cmems_latsize_cut', 'comb_spatial_meancmems', 'comb_spatial_meaninterped', ...
%         'mask_cmems', 'comb_cmems_data', 'comb_cmems_adt', 'comb_interped_data', 'comb_spatial_meaninterped', ...
%         'cmems_lon2', 'cmems_lat2', 'lon', 'lat', '-v7.3');
    save(RCM_info.matname, 'RCM_info', 'RCM_grid', 'RCM_data', '-v7.3');
    save(RCM_info.matname_interped, 'RCM_info', 'RCM_grid', 'RCM_data_interped', 'CMEMS_grid', '-v7.3');
    save(GCM_info.matname, 'GCM_info', 'GCM_grid', 'GCM_data', '-v7.3');
    save(GCM_info.matname_interped, 'GCM_info', 'GCM_grid', 'GCM_data_interped', 'CMEMS_grid', '-v7.3');
    save(CMEMS_info.matname, 'CMEMS_info', 'CMEMS_grid', 'CMEMS_data', '-v7.3');
end
