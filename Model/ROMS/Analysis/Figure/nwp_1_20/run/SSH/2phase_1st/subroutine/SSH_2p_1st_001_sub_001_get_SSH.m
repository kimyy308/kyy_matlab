% %         get model data and cmems data      
% %  Updated 05-Jul-2021 by Yong-Yub Kim, modified to Func_0012

run(tmp.param_script);
tmp.ind_for=1;
matname = [tmp.savedir,tmp.testname,'_',tmp.regionname,'model_cmems_ssh_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];

if (exist(matname , 'file') ~= 2 || flags.fig_switch==2)   
    for yearij = 1:length(RCM_info.years)
        for monthij = 1:length(RCM_info.months)
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])
            tic;
            tmp.year = RCM_info.years(yearij);
            tmp.month = RCM_info.months(monthij);
            % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
            filename = strcat(tmp.filedir, num2str(tmp.year,'%04i'), tmp.fs, ...
                    'pck_', tmp.testname, '_monthly_', num2str(tmp.year,'%04i'), ...
                    '_', num2str(tmp.month,'%02i'), '.nc');
% %             get RCM_RCM_grid
            if (exist('RCM_RCM_grid.lon_rho')==0)
                RCM_grid.lon_rho_whole = ncread(filename,'lon_rho',[1 1],[inf,1]);
                RCM_grid.lat_rho_whole = ncread(filename,'lat_rho',[1 1],[1,inf]);
                [RCM_grid.ind_w, RCM_grid.ind_e, RCM_grid.ind_s, RCM_grid.ind_n] = ...
                    Func_0012_findind_Y(RCM_grid.dl, RCM_grid.domain, RCM_grid.lon_rho_whole, RCM_grid.lat_rho_whole);
                RCM_grid.size_lon = RCM_grid.ind_e-RCM_grid.ind_w+1;
                RCM_grid.size_lat = RCM_grid.ind_n-RCM_grid.ind_s+1;
                RCM_grid.lon_rho = ncread(filename,'lon_rho', [RCM_grid.ind_w RCM_grid.ind_s], ...
                    [RCM_grid.size_lon RCM_grid.size_lat]);
                RCM_grid.lat_rho = ncread(filename,'lat_rho', [RCM_grid.ind_w RCM_grid.ind_s], ...
                    [RCM_grid.size_lon RCM_grid.size_lat]);

                RCM_grid.mask_region = double(inpolygon(RCM_grid.lon_rho, RCM_grid.lat_rho, ...
                    RCM_grid.refpolygon(:,1),RCM_grid.refpolygon(:,2)));
                RCM_grid.mask_region(RCM_grid.mask_region==0)=NaN;
                RCM_info.data_info = ncinfo(filename, varname);
                RCM_grid.mask_ocean= ncread(filename, 'zeta', [RCM_grid.ind_w RCM_grid.ind_s, 1], ...
                    [RCM_grid.size_lon RCM_grid.size_lat, 1]);
                RCM_grid.mask_ocean(isfinite(RCM_grid.mask_ocean))=1;
                RCM_grid.mask_land = RCM_grid.mask_ocean;
                RCM_grid.mask_land(isfinite(RCM_grid.mask_ocean))=NaN;
                RCM_grid.mask_land(isnan(RCM_grid.mask_ocean))=1;
            end
            
% %             get RCM data
            RCM_data.tmp_data = ncread(filename, 'zeta', [RCM_grid.ind_w RCM_grid.ind_s, monthij], ...
                    [RCM_grid.size_lon RCM_grid.size_lat, 1]);
            RCM_data.tmp_data=RCM_data.tmp_data.*RCM_grid.mask_region;

            if (exist('comb_spatial_meanmodel')==0)
                RCM_data.clim_mean=(zeros([RCM_grid.size_lon,RCM_grid.size_lat,12]));
            end
            RCM_data.all(:,:,tmp.ind_for) = single(RCM_data.tmp_data);
            RCM_data.clim_mean(:,:,monthij)=RCM_data.clim_mean(:,:,monthij) + ...
                RCM_data.tmp_data / double(length(RCM_info.years));


% % % % % need to GCM handling
ffffffffffffffffffff
% %             get CMEMS data
            CMEMS_info.filename = strcat(tmp.cmemsdir,'cmems_nwp_ssh_', num2str(tmp.year,'%04i'),'.nc');
            if (exist('cmems_lon')==0)
                CMEMS_info.file_info=ncinfo(CMEMS_info.filename);
                CMEMS_info.lon_name='longitude';
                CMEMS_info.lat_name='latitude';
                CMEMS_info.var_name='ssha';
                CMEMS_info.lon_info=ncinfo(CMEMS_info.filename,CMEMS_info.lon_name);
                CMEMS_info.lat_info=ncinfo(CMEMS_info.filename,CMEMS_info.lat_name);
                CMEMS_grid.dl = 1/4;
                CMEMS_grid.domain = RCM_grid.domain;
                CMEMS_grid.lon = ncread(CMEMS_info.filename,CMEMS_info.lon_name,1,cmemsinfo_lon.Dimensions.Length);
                CMEMS_grid.lat = ncread(CMEMS_info.filename,CMEMS_info.lat_name,1,cmemsinfo_lat.Dimensions.Length);
                [CMEMS_grid.lat2_whole, CMEMS_grid.lon2_whole]= meshgrid(CMEMS_grid.lat, CMEMS_grid.lon);
                [CMEMS_grid.ind_w, CMEMS_grid.ind_e, CMEMS_grid.ind_s, CMEMS_grid.ind_n] = ...
                    Func_0012_findind_Y(CMEMS_grid.dl, CMEMS_grid.domain, CMEMS_grid.lon2_whole, CMEMS_grid.lat2_whole);
                
                CMEMS_grid.lon2 = cmems_lon2(CMEMS_grid.ind_w:CMEMS_grid.ind_e, CMEMS_grid.ind_s:CMEMS_grid.ind_n);
                CMEMS_grid.lat2 = cmems_lat2(CMEMS_grid.ind_w:CMEMS_grid.ind_e, CMEMS_grid.ind_s:CMEMS_grid.ind_n);

                CMEMS_grid.size_lon=CMEMS_grid.ind_e-CMEMS_grid.ind_w+1;
                CMEMS_grid.size_lat=CMEMS_grid.ind_n-CMEMS_grid.ind_s+1;
                CMEMS_grid.refpolygon=RCM_grid.refpolygon;
%                 len_lon=cmems_lonsize_cut;
%                 len_lat=cmems_latsize_cut;
%                     comb_spatial_meanrms=(zeros([length(cmems_lon),length(cmems_lat),12]));
%                     comb_spatial_meanbias=(zeros([length(cmems_lon),length(cmems_lat),12]));
                CMEMS_data.clim_mean=(zeros([CMEMS_grid.size_lon, CMEMS_grid.size_lat,12]));
                RCM_data.clim_mean_interped=(zeros([CMEMS_grid.size_lon, CMEMS_grid.size_lat,12]));
%                     comb_spatial_meanmodel=(zeros([length(cmems_lon),length(cmems_lat),12]));

                CMEMS_grid.mask_region = double(inpolygon(CMEMS_grid.lon2, CMEMS_grid.lat2, ...
                    CMEMS_grid.refpolygon(:,1),CMEMS_grid.refpolygon(:,2)));
                CMEMS_grid.mask_region(CMEMS_grid.mask_region==0)=NaN;
                CMEMS_grid.mask_ocean= ncread(filename, 'zeta', [CMEMS_grid.ind_w CMEMS_grid.ind_s, 1], ...
                    [CMEMS_grid.size_lon CMEMS_grid.size_lat, 1]);
                CMEMS_grid.mask_ocean(isfinite(CMEMS_grid.mask_ocean))=1;
                CMEMS_grid.mask_land = CMEMS_grid.mask_ocean;
                CMEMS_grid.mask_land(isfinite(CMEMS_grid.mask_ocean))=NaN;
                CMEMS_grid.mask_land(isnan(CMEMS_grid.mask_ocean))=1;
            end
1111111111111111111
            if tmp.month==1
                cmems_daily_adt=ncread(CMEMS_info.filename,'adt',[cmems_lon_min(1) cmems_lat_min(1) 1], [cmems_lon_max(1)-cmems_lon_min(1)+1 cmems_lat_max(1)-cmems_lat_min(1)+1 31]);
                cmems_daily_data=ncread(CMEMS_info.filename,'sla',[cmems_lon_min(1) cmems_lat_min(1) 1], [cmems_lon_max(1)-cmems_lon_min(1)+1 cmems_lat_max(1)-cmems_lat_min(1)+1 31]);
            else
                cmems_daily_adt=ncread(CMEMS_info.filename,'adt',[cmems_lon_min(1) cmems_lat_min(1) sum(eomday(tmp.year,1:tmp.month-1))+1], [cmems_lon_max(1)-cmems_lon_min(1)+1 cmems_lat_max(1)-cmems_lat_min(1)+1 eomday(tmp.year,tmp.month)]);
                cmems_daily_data=ncread(CMEMS_info.filename,'sla',[cmems_lon_min(1) cmems_lat_min(1) sum(eomday(tmp.year,1:tmp.month-1))+1], [cmems_lon_max(1)-cmems_lon_min(1)+1 cmems_lat_max(1)-cmems_lat_min(1)+1 eomday(tmp.year,tmp.month)]);
            end
            cmems_adt = mean(cmems_daily_adt,3,'omitnan');
            cmems_data = mean(cmems_daily_data,3,'omitnan');
%             cmems_data = ncread(CMEMS_info.filename,varname,[cmems_lon_min(1) cmems_lat_min(1) tmp.month], [cmems_lon_max(1)-cmems_lon_min(1) cmems_lat_max(1)-cmems_lat_min(1) 1]);
%             cmems_data(cmems_data<-1000)=NaN;
%             cmems_data(cmems_data>1000)=NaN;
            cmems_data=cmems_data.*mask_cmems;

            comb_cmems_adt(:,:,ind) = cmems_adt;
            comb_cmems_data(:,:,ind) = cmems_data;
            comb_spatial_meancmems(:,:,monthij)=comb_spatial_meancmems(:,:,monthij)+cmems_adt/double(length(RCM_info.years));        %                 comb_spatial_meancmems(:,:,monthij)=comb_spatial_meancmems(:,:,monthij)+cmems_adt/double(length(RCM_info.years));

            interped_data = griddata(double(lon), double(lat), data,double(cmems_lon2),double(cmems_lat2));   

            comb_interped_data(:,:,ind) = interped_data;
            comb_spatial_meaninterped(:,:,monthij)=comb_spatial_meaninterped(:,:,monthij)+interped_data/double(length(RCM_info.years));

            tmp.ind_for = tmp.ind_for + 1;
            toc;
        end
    end
    save(matname, 'RCM_grid.mask_region', 'len_lon_model', 'len_lat_model', 'comb_spatial_meanmodel', 'comb_data', ...
        'cmems_lonsize_cut', 'cmems_latsize_cut', 'comb_spatial_meancmems', 'comb_spatial_meaninterped', ...
        'mask_cmems', 'comb_cmems_data', 'comb_cmems_adt', 'comb_interped_data', 'comb_spatial_meaninterped', ...
        'cmems_lon2', 'cmems_lat2', 'lon', 'lat', '-v7.3');
else
    load(matname);
end
flags.fig_switch=0;
