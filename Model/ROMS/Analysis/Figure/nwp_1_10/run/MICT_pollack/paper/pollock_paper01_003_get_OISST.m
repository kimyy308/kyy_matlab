close all; clear all;  clc;   
% all_region ={'NWP','ES', 'SS', 'YS', 'ECS'}
% all_region ={'EKB','EK','WB'}
% all_region ={'pollock_egg', 'pollock_egg2'}
% all_region ={'ES', 'pollock_egg'}
all_region ={'pollock_egg3'}

warning off;
for regionind=1:length(all_region)
    clearvars '*' -except regionind all_region

    % % % 
    % % % Read Model SST
    % % % interp
    % % % get sq_er
    % % % get BIAS
    system_name=computer;
    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        dropboxpath='C:\Users\KYY\Dropbox';
        addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
        addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
        addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
        addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
    elseif (strcmp(system_name,'GLNXA64'))
        dropboxpath='/home/kimyy/Dropbox';
        addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
        addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
        addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    end


    shadlev = [0 35];
    sq_er_shadlev = [0 4];
    bias_shadlev = [-4 4];
    conlev  = 0:5:35;
    
    % for snu_desktop
%     testname='avg_ens_10km_mean_monthly_'   % % need to change
    testname='test06'   % % need to change
    inputyear = [1983:2018]; % % put year which you want to plot [year year ...]
    inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
    
    variable = 'SST';
    varname ='temp';  %% reference variable -> temperature
    run('nwp_polygon_point.m');
    regionname=all_region{regionind}
    [error_status, refpolygon, lonlat] = Func_0007_get_polygon_data_from_regionname(regionname);


    switch(testname)
        case('avg_ens_10km_mean_monthly_') %% seo's Reanalysis data
            if (strcmp(system_name,'PCWIN64'))
                % % for windows
                figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd_year\figure\avg_ens_10km_mean\',regionname,'\'); % % where figure files will be saved
                param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_sq_er.m'
        %         filedir = strcat('E:\Data\Observation\OISST\monthly\'); % % where data files are
        %         avhrrdir='E:\Data\Observation\OISST\monthly\';
                filedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\'); % % where data files are
                avhrrdir='E:\Data\Observation\OISST\monthly_kimyy\';
            elseif (strcmp(system_name,'GLNXA64'))
            end
            dl=1/10;
        otherwise
            if (strcmp(system_name,'PCWIN64'))
                % % for windows
%                 figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd_year\figure\nwp_1_10\',testname,'\',regionname,'\'); % % where figure files will be saved
                param_script =['C:\users\user/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_10/run/fig_param/fig_param_kyy_', regionname, '.m'];
        %         filedir = strcat('E:\Data\Observation\OISST\monthly\'); % % where data files are
        %         avhrrdir='E:\Data\Observation\OISST\monthly\';
                filedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA'); % % where data files are
                avhrrdir='Z:\내 드라이브\Data\Observation\OISST\monthly_kimyy\';
            elseif (strcmp(system_name,'GLNXA64'))
            end
            dl=1/10;
    end
    

    run(param_script);
    ind=1;
    for yearij = 1:length(inputyear)
        for monthij = 1:length(inputmonth)
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])
            tic;
            tempyear = inputyear(yearij);
            tempmonth = inputmonth(monthij);
            % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc

            switch(testname)
                case('avg_ens_10km_mean_monthly_') %% seo's Reanalysis data
                    filename = strcat(filedir, ...
                    testname, num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');
                    grdname = strcat(filedir, 'roms_grid_final.nc');
                otherwise
                    filename = strcat(filedir, '\', num2str(tempyear,'%04i'),'\', ...
                    testname, '_monthly_',num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');
                    grdname = filename;
            end
            % read model data
            if (exist('lon')==0)
                modelinfo=ncinfo(grdname);
                lon = ncread(grdname,'lon_rho',[1 1],[modelinfo.Dimensions(5).Length,1]);
                lat = ncread(grdname,'lat_rho',[1 1],[1,modelinfo.Dimensions(6).Length]);

                lon_west = abs(lon - (lonlat(1)-1));
                min_lon_west=min(lon_west);
                lon_east = abs(lon - (lonlat(2)+1));
                min_lon_east=min(lon_east);
                lat_south = abs(lat - (lonlat(3)-1));
                min_lat_south=min(lat_south);
                lat_north = abs(lat - (lonlat(4)+1));
                min_lat_north=min(lat_north);

                lon_min = find(lon_west == min_lon_west);
                lon_max = find(lon_east == min_lon_east);
                lat_min = find(lat_south == min_lat_south);
                lat_max = find(lat_north == min_lat_north);

                lon = ncread(grdname,'lon_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                lat = ncread(grdname,'lat_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);

    %             polygon_ind=NaN(size(refpolygon,1),2);
    %             for i=1:size(refpolygon,1)
    %                 [polygon_ind(i,1), trash_ind, polygon_ind(i,2), trash_ind]=findind_Y(dl, [refpolygon(i,1),refpolygon(i,2)],lon',lat');
    %             end
    %             mask_model = inpolygon(lon,lat,polygon_ind(:,1),polygon_ind(:,2));
                switch(regionname)
                    case('NWP') %% North western Pacific
                        mask_model(1:size(lon,1),1:size(lon,2))=1;
                    otherwise
                        mask_model = double(inpolygon(lon,lat,refpolygon(:,1),refpolygon(:,2)));
                        mask_model(mask_model==0)=NaN;
                end
            end

            data_info = ncinfo(filename, varname);  %% [lon lat depth time] -> [1601 1201 33 1]
            if (length(data_info.Dimensions)==4)
                data = ncread(filename,varname,[lon_min(1) lat_min(1) 40 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);
            else
                data = ncread(filename,varname,[lon_min(1) lat_min(1) 40], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
            end
            data(data<-100)=NaN;
            data(data>100)=NaN;
            data(data==0)=NaN;
            
            data=data.*mask_model;
            
            
            % read OISST DATA
            
%             if (tempyear<2011)
%                 avhrrfilename = strcat(avhrrdir,'avhrr_monthly', num2str(tempyear,'%04i'), ...
%                     '_',num2str(tempmonth,'%02i'), '.nc');
%             else
%                 avhrrfilename = strcat(avhrrdir,'avhrr-monthly_', num2str(tempyear,'%04i'), ...
%                     num2str(tempmonth,'%02i'), '.nc');
%             end
            avhrrfilename = strcat(avhrrdir,'avhrr_only_monthly_v2_', num2str(tempyear,'%04i'), '.nc');
                
%             if (exist('avhrr_lon')==0)
                avhrrinfo=ncinfo(avhrrfilename);
%                 if (tempyear<2011)
%                     avhrr_lon = ncread(avhrrfilename,'long',[1 1],[avhrrinfo.Dimensions(2).Length,1]);
%                     avhrr_lat = ncread(avhrrfilename,'lat',[1 1],[1,avhrrinfo.Dimensions(1).Length]);
%                 else
%                     avhrr_lon = ncread(avhrrfilename,'lon',[1],[avhrrinfo.Dimensions(2).Length]);
%                     avhrr_lat = ncread(avhrrfilename,'lat',[1],[avhrrinfo.Dimensions(1).Length]);
%                 end
                avhrr_lon = ncread(avhrrfilename,'lon',[1],[avhrrinfo.Dimensions(2).Length]);
                avhrr_lat = ncread(avhrrfilename,'lat',[1],[avhrrinfo.Dimensions(1).Length]);
                
                avhrr_lon_west = abs(avhrr_lon - (lonlat(1)-1));
                min_avhrr_lon_west=min(avhrr_lon_west);
                avhrr_lon_east = abs(avhrr_lon - (lonlat(2)+1));
                min_avhrr_lon_east=min(avhrr_lon_east);
                avhrr_lat_south = abs(avhrr_lat - (lonlat(3)-1));
                min_avhrr_lat_south=min(avhrr_lat_south);
                avhrr_lat_north = abs(avhrr_lat - (lonlat(4)+1));
                min_avhrr_lat_north=min(avhrr_lat_north);

                avhrr_lon_min = find(avhrr_lon_west == min_avhrr_lon_west);
                avhrr_lon_max = find(avhrr_lon_east == min_avhrr_lon_east);
                avhrr_lat_min = find(avhrr_lat_south == min_avhrr_lat_south);
                avhrr_lat_max = find(avhrr_lat_north == min_avhrr_lat_north);

        %         ncinfo('E:\Data\Observation\OISST\monthly\avhrr_monthly1983_11.nc');
% %                 if (tempyear<2011)
% %                     avhrr_lon = ncread(avhrrfilename,'long', [avhrr_lon_min(1) avhrr_lat_min(1)], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1)]);
% %                     avhrr_lat = ncread(avhrrfilename,'lat', [avhrr_lon_min(1) avhrr_lat_min(1)], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1)]);
% % 
% %                 else
% %                     avhrr_lon = ncread(avhrrfilename,'lon', [avhrr_lon_min(1)], [avhrr_lon_max(1)-avhrr_lon_min(1)]);
% %                     avhrr_lat = ncread(avhrrfilename,'lat', [avhrr_lat_min(1)], [avhrr_lat_max(1)-avhrr_lat_min(1)]);
% %                 end
                avhrr_lon = ncread(avhrrfilename,'lon', [avhrr_lon_min(1)], [avhrr_lon_max(1)-avhrr_lon_min(1)]);
                avhrr_lat = ncread(avhrrfilename,'lat', [avhrr_lat_min(1)], [avhrr_lat_max(1)-avhrr_lat_min(1)]);
                
                if (exist('comb_spatial_meansq_er')==0)
                    comb_spatial_meansq_er=(zeros([length(avhrr_lon),length(avhrr_lat),12]));
                    comb_spatial_meanbias=(zeros([length(avhrr_lon),length(avhrr_lat),12]));
                    comb_spatial_meanavhrr=(zeros([length(avhrr_lon),length(avhrr_lat),12]));
                    comb_spatial_meanmodel=(zeros([length(avhrr_lon),length(avhrr_lat),12]));
                end
                switch(regionname)
                    case('NWP') %% North western Pacific
                        mask_avhrr(1:size(avhrr_lon,1),1:size(avhrr_lon,2))=1;  %% before 2011
                    otherwise
%                         if (tempyear<2011)
%                             avhrr_lat2=avhrr_lat;
%                             avhrr_lon2=avhrr_lon;
%                         else
%                             [avhrr_lat2 avhrr_lon2]=meshgrid(avhrr_lat, avhrr_lon);
%                         end
                        [avhrr_lat2 avhrr_lon2]=meshgrid(avhrr_lat, avhrr_lon);
                        mask_avhrr = double(inpolygon(avhrr_lon2,avhrr_lat2,refpolygon(:,1),refpolygon(:,2)));
                        mask_avhrr(mask_avhrr==0)=NaN;
                end
%             end
            
%             if (tempyear<2011)
%                 len_lon = length(avhrr_lon(:,1));
%                 len_lat = length(avhrr_lat(1,:));
%                 avhrr_data = ncread(avhrrfilename,varname,[avhrr_lon_min(1) avhrr_lat_min(1)], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1)]);
%             else
%                 len_lon = length(avhrr_lon(:));
%                 len_lat = length(avhrr_lat(:));
%                 avhrr_data = ncread(avhrrfilename,'sst',[avhrr_lon_min(1) avhrr_lat_min(1) 1], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1) 1]);
%             end
            len_lon = length(avhrr_lon(:));
            len_lat = length(avhrr_lat(:));
            len_lon_model = size(data,1);
            len_lat_model = size(data,2);
            avhrr_data = ncread(avhrrfilename,'temp',[avhrr_lon_min(1) avhrr_lat_min(1) tempmonth], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1) 1]);            
            avhrr_err = ncread(avhrrfilename,'err',[avhrr_lon_min(1) avhrr_lat_min(1) tempmonth], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1) 1]);            
            
            avhrr_data(avhrr_data<-9)=NaN;
            avhrr_data(avhrr_data>1000)=NaN;
            avhrr_data=avhrr_data.*mask_avhrr;
            
            avhrr_err(avhrr_err<-100)=NaN;
            avhrr_err(avhrr_err>1000)=NaN;
            avhrr_err=avhrr_err.*mask_avhrr;
            
            interped_data = griddata(double(lon), double(lat), double(data), double(avhrr_lon), double(avhrr_lat)')';   
            bias = interped_data-avhrr_data;  
%             rms = sqrt((avhrr_data-interped_data).^2);  
            sq_er = (avhrr_data-interped_data).^2;            
            meanbias = mean(mean(bias,'omitnan'),'omitnan');
            meansq_er = mean(mean(sq_er,'omitnan'),'omitnan');
            
            comb_interped_data(:,:,ind) = interped_data;
            comb_avhrr_data(:,:,ind) = avhrr_data;
            comb_avhrr_err(:,:,ind) =avhrr_err;
            comb_sq_er_data(:,:,ind) = sq_er;
            comb_bias_data(:,:,ind) = bias;
            
            comb_meansq_er(yearij,monthij)=meansq_er;
            comb_meanbias(yearij,monthij)=meanbias;
            comb_spatial_meansq_er(:,:,monthij)=comb_spatial_meansq_er(:,:,monthij)+sq_er/double(length(inputyear));
            comb_spatial_meanbias(:,:,monthij)=comb_spatial_meanbias(:,:,monthij)+bias/double(length(inputyear));
            comb_spatial_meanavhrr(:,:,monthij)=comb_spatial_meanavhrr(:,:,monthij)+avhrr_data/double(length(inputyear));
            comb_spatial_meanmodel(:,:,monthij)=comb_spatial_meanmodel(:,:,monthij)+interped_data/double(length(inputyear));
        
            comb_data(:,:,ind) = data;
%     pcolor(avhrr_data'-273.15); shading flat; colorbar        
            comb_interped_data(:,:,ind) = interped_data;
%     pcolor(interped_data'); shading flat; colorbar        
            ind = ind + 1;
            toc;
        end
    end
    comb_avhrr_clim_divided=reshape(comb_avhrr_data,[len_lon, len_lat, 12, length(inputyear)]);
    comb_interped_clim_divided=reshape(comb_interped_data,[len_lon, len_lat, 12, length(inputyear)]);
    avhrr_land = ncread(avhrrfilename,'temp',[avhrr_lon_min(1) avhrr_lat_min(1) tempmonth], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1) 1]);
    avhrr_land(isnan(avhrr_land))=5000;
    avhrr_land(avhrr_land<5000)=NaN;
    avhrr_land(isfinite(avhrr_land))=1;
    
    model_land = ncread(filename,varname,[lon_min(1) lat_min(1) 40 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);
    model_land(isnan(model_land))=5000;
    model_land(model_land<5000)=NaN;
    model_land(isfinite(model_land))=1;
    
%     pcolor(avhrr_land'); shading flat; colorbar
    %     pcolor(model_land'); shading flat; colorbar

    trendtime=inputyear(1):1/12:inputyear(end)+1-1/12;
    trend(1:len_lon,1:len_lat)=NaN;
    for i=1:len_lon
        for j=1:len_lat
            p=polyfit(trendtime,squeeze(comb_interped_data(i,j,:))',1);
            trend(i,j)=p(1);
        end
    end
    
    avhrr_trend(1:len_lon,1:len_lat)=NaN;
    for i=1:len_lon
        for j=1:len_lat
            p=polyfit(trendtime,squeeze(comb_avhrr_data(i,j,:))',1);
            avhrr_trend(i,j)=p(1);
        end
    end

    for t=1:length(inputyear)
        comb_interped_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_interped_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meanmodel;
        comb_avhrr_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_avhrr_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meanavhrr;
    end
    
    trend_filtered(1:len_lon,1:len_lat)=NaN;
    avhrr_trend_filtered(1:len_lon,1:len_lat)=NaN;
    
    for i=1:len_lon
        for j=1:len_lat
            p=polyfit(trendtime,squeeze(comb_interped_data_filtered(i,j,:))',1);
            trend_filtered(i,j)=p(1);
        end
    end
    
    for i=1:len_lon
        for j=1:len_lat
            p=polyfit(trendtime,squeeze(comb_avhrr_data_filtered(i,j,:))',1);
            avhrr_trend_filtered(i,j)=p(1);
        end
    end
    
    mean_trend=mean(mean(trend,'omitnan'),'omitnan');
    mean_trend_filtered=mean(mean(trend_filtered,'omitnan'),'omitnan');
    mean_avhrr_trend=mean(mean(avhrr_trend,'omitnan'),'omitnan');
    mean_avhrr_trend_filtered=mean(mean(avhrr_trend_filtered,'omitnan'),'omitnan');
    
    clim_avhrr_trend_filtered(1:len_lon,1:len_lat,1:12)=NaN;
    for k=1:12
        clim_trendtime(:,k)=inputyear(1)+(k-1)/12 : 1 : inputyear(end)+(k-1)/12;
        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(clim_trendtime(:,k)',squeeze(comb_avhrr_clim_divided(i,j,k,:))',1);
                clim_avhrr_trend_divided(i,j,k)=p(1);
            end
        end
    end
    
    clim_interped_trend_filtered(1:len_lon,1:len_lat,1:12)=NaN;
    for k=1:12
        clim_trendtime(:,k)=inputyear(1)+(k-1)/12 : 1 : inputyear(end)+(k-1)/12;
        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(clim_trendtime(:,k)',squeeze(comb_interped_clim_divided(i,j,k,:))',1);
                clim_interped_trend_divided(i,j,k)=p(1);
            end
        end
    end

    mean_avhrr_trend=mean(mean(avhrr_trend,'omitnan'),'omitnan');
    mean_avhrr_trend_filtered=mean(mean(avhrr_trend_filtered,'omitnan'),'omitnan');
    mean_clim_avhrr_trend_divided=mean(mean(clim_avhrr_trend_divided,'omitnan'),'omitnan');
    mean_clim_interped_trend_divided=mean(mean(clim_interped_trend_divided,'omitnan'),'omitnan');

%     figdir=[figrawdir,'CLIM\'];
%     outfile = strcat(figdir,regionname);
%     if (exist(strcat(figdir) , 'dir') ~= 7)
%         mkdir(strcat(figdir));
%     end 

    save([filedir,'\', regionname,'sst_sq_er_and_bias_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

    ncid = netcdf.create(strcat(filedir, '\', testname,regionname,'_sq_er_clim_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc'),'NETCDF4');

    onedimid = netcdf.defDim(ncid,'one', 1);
    lon_dimid = netcdf.defDim(ncid, 'lon', len_lon);
    lat_dimid = netcdf.defDim(ncid,'lat',len_lat);
    xidimid = netcdf.defDim(ncid, 'xi_rho', len_lon_model);
    etadimid = netcdf.defDim(ncid,'eta_rho',len_lat_model);
    time_dimid = netcdf.defDim(ncid, 'time', 0);
    clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ['NWP 1/10 Reanalysis_ ', testname, 'monthly SST sq_er/BIAS file']);
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', ['NWP 1/10 Reanalysis SST trend (',num2str(inputyear(1)),'-',num2str(inputyear(end)),') ']);
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', [' NWP 1/10 Reanalysis ' ]);
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'author', 'Created by Y.Y.Kim');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'date', date);

    timevarid=netcdf.defVar(ncid, 'time', 'NC_DOUBLE', time_dimid);
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days since 1900-12-31 00:00:00');
    netcdf.putAtt(ncid,timevarid,'calendar','gregorian');

    clim_timevarid=netcdf.defVar(ncid, 'clim_time', 'NC_DOUBLE', clim_time_dimid);
    netcdf.putAtt(ncid,clim_timevarid,'long_name','clim_time');
    netcdf.putAtt(ncid,clim_timevarid,'units','days since 1900-12-31 00:00:00');
    netcdf.putAtt(ncid,clim_timevarid,'calendar','gregorian');

    lonvarid=netcdf.defVar(ncid, 'lon', 'NC_DOUBLE', lon_dimid);
    netcdf.putAtt(ncid,lonvarid,'long_name','longitude');
    netcdf.putAtt(ncid,lonvarid,'units','degree_east');

    latvarid=netcdf.defVar(ncid, 'lat', 'NC_DOUBLE', lat_dimid);
    netcdf.putAtt(ncid,latvarid,'long_name','latitude');
    netcdf.putAtt(ncid,latvarid,'units','degree_north');
    
    lon_rhovarid=netcdf.defVar(ncid, 'lon_rho', 'NC_DOUBLE', [xidimid etadimid]);
    netcdf.putAtt(ncid,lon_rhovarid,'long_name','lon_model');
    netcdf.putAtt(ncid,lon_rhovarid,'units','degree_east');

    lat_rhovarid=netcdf.defVar(ncid, 'lat_rho', 'NC_DOUBLE', [xidimid etadimid]);
    netcdf.putAtt(ncid,lat_rhovarid,'long_name','lat_model');
    netcdf.putAtt(ncid,lat_rhovarid,'units','degree_north');
    
    raw_sstvarid=netcdf.defVar(ncid, 'raw_sst', 'NC_FLOAT', [xidimid etadimid time_dimid]);
    netcdf.putAtt(ncid,raw_sstvarid,'long_name','raw_sst');
    netcdf.putAtt(ncid,raw_sstvarid,'units','Celsius');

    sstvarid=netcdf.defVar(ncid, 'sst', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,sstvarid,'long_name','sst');
    netcdf.putAtt(ncid,sstvarid,'units','Celsius');

    sst_filteredvarid=netcdf.defVar(ncid, 'sst_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,sst_filteredvarid,'long_name','sst_filtered');
    netcdf.putAtt(ncid,sst_filteredvarid,'units','Celsius');
    
    avhrr_sstvarid=netcdf.defVar(ncid, 'avhrr_sst', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,avhrr_sstvarid,'long_name','avhrr_sst');
    netcdf.putAtt(ncid,avhrr_sstvarid,'units','Celsius');
    
    avhrr_errvarid=netcdf.defVar(ncid, 'avhrr_err', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,avhrr_errvarid,'long_name','ahvrr_err');
    netcdf.putAtt(ncid,avhrr_errvarid,'units','Celsius');

    avhrr_sst_filteredvarid=netcdf.defVar(ncid, 'avhrr_sst_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,avhrr_sst_filteredvarid,'long_name','avhrr_sst_filtered');
    netcdf.putAtt(ncid,avhrr_sst_filteredvarid,'units','Celsius');
    
    model_landvarid=netcdf.defVar(ncid, 'model_land', 'NC_FLOAT', [xidimid etadimid]);
    netcdf.putAtt(ncid,model_landvarid,'long_name','model_land');
    
    avhrr_landvarid=netcdf.defVar(ncid, 'avhrr_land', 'NC_FLOAT', [lon_dimid lat_dimid]);
    netcdf.putAtt(ncid,avhrr_landvarid,'long_name','avhrr_land');
    
    trendvarid=netcdf.defVar(ncid, 'trend', 'NC_FLOAT', [lon_dimid lat_dimid]);
    netcdf.putAtt(ncid,trendvarid,'long_name','trend');
    netcdf.putAtt(ncid,trendvarid,'units','Celsius/year');

    trend_filteredvarid=netcdf.defVar(ncid, 'trend_filtered', 'NC_FLOAT', [lon_dimid lat_dimid]);
    netcdf.putAtt(ncid,trend_filteredvarid,'long_name','trend_filtered');
    netcdf.putAtt(ncid,trend_filteredvarid,'units','Celsius/year');
    
    avhrr_trendvarid=netcdf.defVar(ncid, 'avhrr_trend', 'NC_FLOAT', [lon_dimid lat_dimid]);
    netcdf.putAtt(ncid,avhrr_trendvarid,'long_name','avhrr_trend');
    netcdf.putAtt(ncid,avhrr_trendvarid,'units','Celsius/year');

    avhrr_trend_filteredvarid=netcdf.defVar(ncid, 'avhrr_trend_filtered', 'NC_FLOAT', [lon_dimid lat_dimid]);
    netcdf.putAtt(ncid,avhrr_trend_filteredvarid,'long_name','avhrr_trend_filtered');
    netcdf.putAtt(ncid,avhrr_trend_filteredvarid,'units','Celsius/year');
    
    clim_avhrr_trend_dividedvarid=netcdf.defVar(ncid, 'clim_avhrr_trend_divided', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_avhrr_trend_dividedvarid,'long_name','clim_avhrr_trend_divided');
    netcdf.putAtt(ncid,clim_avhrr_trend_dividedvarid,'units','Celsius/year');
    
    clim_interped_trend_dividedvarid=netcdf.defVar(ncid, 'clim_interped_trend_divided', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_interped_trend_dividedvarid,'long_name','clim_interped_trend_divided');
    netcdf.putAtt(ncid,clim_interped_trend_dividedvarid,'units','Celsius/year');
    
    sq_ervarid=netcdf.defVar(ncid, 'sq_er', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,sq_ervarid,'long_name','sq_er');
    netcdf.putAtt(ncid,sq_ervarid,'units','Celsius');

    biasvarid=netcdf.defVar(ncid, 'bias', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,biasvarid,'long_name','bias');
    netcdf.putAtt(ncid,biasvarid,'units','Celsius');

    mean_trendvarid=netcdf.defVar(ncid, 'mean_trend', 'NC_FLOAT', onedimid);
    netcdf.putAtt(ncid,mean_trendvarid,'long_name','mean_trend');
    netcdf.putAtt(ncid,mean_trendvarid,'units','Celsius/year');

    mean_trend_filteredvarid=netcdf.defVar(ncid, 'mean_trend_filtered', 'NC_FLOAT', onedimid);
    netcdf.putAtt(ncid,mean_trend_filteredvarid,'long_name','mean_trend_filtered');
    netcdf.putAtt(ncid,mean_trend_filteredvarid,'units','Celsius/year');
    
    mean_avhrr_trendvarid=netcdf.defVar(ncid, 'mean_avhrr_trend', 'NC_FLOAT', onedimid);
    netcdf.putAtt(ncid,mean_avhrr_trendvarid,'long_name','mean_avhrr_trend');
    netcdf.putAtt(ncid,mean_avhrr_trendvarid,'units','Celsius/year');

    mean_avhrr_trend_filteredvarid=netcdf.defVar(ncid, 'mean_avhrr_trend_filtered', 'NC_FLOAT', onedimid);
    netcdf.putAtt(ncid,mean_avhrr_trend_filteredvarid,'long_name','mean_avhrr_trend_filtered');
    netcdf.putAtt(ncid,mean_avhrr_trend_filteredvarid,'units','Celsius/year');
    
    clim_sstvarid=netcdf.defVar(ncid, 'clim_sst', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_sstvarid,'long_name','clim_sst');
    netcdf.putAtt(ncid,clim_sstvarid,'units','Celsius');
    
    mean_clim_avhrr_trend_dividedvarid=netcdf.defVar(ncid, 'mean_clim_avhrr_trend_divided', 'NC_FLOAT', clim_time_dimid);
    netcdf.putAtt(ncid,mean_clim_avhrr_trend_dividedvarid,'long_name','mean_clim_avhrr_trend_divided');
    netcdf.putAtt(ncid,mean_clim_avhrr_trend_dividedvarid,'units','Celsius/year');
    
    mean_clim_interped_trend_dividedvarid=netcdf.defVar(ncid, 'mean_clim_interped_trend_divided', 'NC_FLOAT', clim_time_dimid);
    netcdf.putAtt(ncid,mean_clim_interped_trend_dividedvarid,'long_name','mean_clim_interped_trend_divided');
    netcdf.putAtt(ncid,mean_clim_interped_trend_dividedvarid,'units','Celsius/year');

    clim_avhrrvarid=netcdf.defVar(ncid, 'clim_avhrr_sst', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_avhrrvarid,'long_name','clim_avhrr_sst');
    netcdf.putAtt(ncid,clim_avhrrvarid,'units','Celsius');
    
    clim_sq_ervarid=netcdf.defVar(ncid, 'clim_sq_er', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_sq_ervarid,'long_name','clim_sq_er');
    netcdf.putAtt(ncid,clim_sq_ervarid,'units','Celsius');

    clim_biasvarid=netcdf.defVar(ncid, 'clim_bias', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_biasvarid,'long_name','clim_bias');
    netcdf.putAtt(ncid,clim_biasvarid,'units','Celsius');

    netcdf.endDef(ncid);

    tind=1;
    for yearij = 1:length(inputyear)
        for month=1:12 
            tempyear = inputyear(yearij);
            ftime(tind) = datenum(tempyear,month,15) - datenum(1900,12,31);
            tind=tind+1;
        end
    end
    for month=1:12 
            tempyear = inputyear(yearij);
            climtime(month) = datenum(1950,month,15) - datenum(1900,12,31);
    end

    netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
    netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
%     if (tempyear<2011)
%         netcdf.putVar(ncid, lonvarid, 0, len_lon, avhrr_lon(:,1));
%         netcdf.putVar(ncid, latvarid, 0, len_lat, avhrr_lat(1,:));
%     else
        netcdf.putVar(ncid, lonvarid, 0, len_lon, avhrr_lon(:));
        netcdf.putVar(ncid, latvarid, 0, len_lat, avhrr_lat(:));
%     end
%     netcdf.putVar(ncid, avhrr_sstvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_avhrr_data);
%     netcdf.putVar(ncid, avhrr_sst_filteredvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_avhrr_data_filtered);
%     netcdf.putVar(ncid, clim_avhrrvarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meanavhrr);
%     netcdf.putVar(ncid, avhrr_trendvarid, [0 0], [len_lon len_lat], avhrr_trend);
%     netcdf.putVar(ncid, avhrr_trend_filteredvarid, [0 0], [len_lon len_lat], avhrr_trend_filtered);
%     netcdf.putVar(ncid, mean_avhrr_trendvarid, [0], [1], mean_avhrr_trend);
%     netcdf.putVar(ncid, mean_avhrr_trend_filteredvarid, [0], [1], mean_avhrr_trend_filtered);

    netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
    netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
%     netcdf.putVar(ncid, lonvarid, 0, len_lon, avhrr_lon(:,1));
%     netcdf.putVar(ncid, latvarid, 0, len_lat, avhrr_lat(1,:));
    netcdf.putVar(ncid, lon_rhovarid, [0 0], [len_lon_model len_lat_model], lon);
    netcdf.putVar(ncid, lat_rhovarid, [0 0], [len_lon_model len_lat_model], lat);
    netcdf.putVar(ncid, raw_sstvarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_data);
    netcdf.putVar(ncid, sstvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_interped_data);
    netcdf.putVar(ncid, sst_filteredvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_interped_data_filtered);
    netcdf.putVar(ncid, avhrr_sstvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_avhrr_data);
    netcdf.putVar(ncid, avhrr_errvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_avhrr_err);
    netcdf.putVar(ncid, avhrr_sst_filteredvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_avhrr_data_filtered);
    netcdf.putVar(ncid, sq_ervarid, [0 0 0], [len_lon len_lat length(ftime)], comb_sq_er_data);
    netcdf.putVar(ncid, biasvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_bias_data);
    netcdf.putVar(ncid, clim_sstvarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meanmodel);
    netcdf.putVar(ncid, clim_avhrrvarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meanavhrr);
    netcdf.putVar(ncid, clim_sq_ervarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meansq_er);
    netcdf.putVar(ncid, clim_biasvarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meanbias);
    netcdf.putVar(ncid, avhrr_landvarid, [0 0], [len_lon len_lat], avhrr_land);
    netcdf.putVar(ncid, model_landvarid, [0 0], [len_lon_model len_lat_model], model_land);
    netcdf.putVar(ncid, trendvarid, [0 0], [len_lon len_lat], trend);
    netcdf.putVar(ncid, trend_filteredvarid, [0 0], [len_lon len_lat], trend_filtered);
    netcdf.putVar(ncid, avhrr_trendvarid, [0 0], [len_lon len_lat], avhrr_trend);
    netcdf.putVar(ncid, avhrr_trend_filteredvarid, [0 0], [len_lon len_lat], avhrr_trend_filtered);
    netcdf.putVar(ncid, mean_trendvarid, [0], [1], mean_trend);
    netcdf.putVar(ncid, mean_trend_filteredvarid, [0], [1], mean_trend_filtered);
    netcdf.putVar(ncid, mean_avhrr_trendvarid, [0], [1], mean_avhrr_trend);
    netcdf.putVar(ncid, mean_avhrr_trend_filteredvarid, [0], [1], mean_avhrr_trend_filtered);
    netcdf.putVar(ncid, clim_avhrr_trend_dividedvarid, [0 0 0], [len_lon len_lat length(climtime)], clim_avhrr_trend_divided);
    netcdf.putVar(ncid, clim_interped_trend_dividedvarid, [0 0 0], [len_lon len_lat length(climtime)], clim_interped_trend_divided);
    netcdf.putVar(ncid, mean_clim_avhrr_trend_dividedvarid, [0], [length(climtime)], mean_clim_avhrr_trend_divided);
    netcdf.putVar(ncid, mean_clim_interped_trend_dividedvarid, [0], [length(climtime)], mean_clim_interped_trend_divided);

    netcdf.close(ncid);
end