close all; clear all;  clc;   
% %  get SSS RMS, SSS BIAS data and plot(climatological value time series)

% all_region ={'ES','SS', 'YS'}
% all_region ={'ES', 'SS', 'YS', 'ECS'}

% all_testname = {'ens03','test53','test54','test55','test56'};
all_testname = {'test49', 'test52'};

% all_region ={'AKP', 'NWP', 'ES', 'SS', 'YS'};
all_region ={'NWP', 'AKP2'};

for testnameind=1:length(all_testname)
    for regionind=1:length(all_region)
        clearvars '*' -except regionind testnameind all_region all_testname

        % % % 
        % % % Read Model SSS
        % % % interp
        % % % get RMS
        % % % get BIAS
        system_name=computer;
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            dropboxpath='C:\Users\KYY\Dropbox';
            addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
            addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
            addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
            addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
        elseif (strcmp(system_name,'GLNXA64'))
            dropboxpath='/home/kimyy/Dropbox';
            addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
            addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
            addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
            addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
        end


        shadlev = [0 35];
        rms_shadlev = [0 5];
        bias_shadlev = [-5 5];
        conlev  = 0:5:35;
        dl=1/20;
        % for snu_desktop
        testname=all_testname{testnameind} 
        inputyear = [1980:2008]; % % put year which you want to plot [year year ...]
        inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]


        varname ='salt';  %% reference variable -> temperature
        run('nwp_polygon_point.m');
        regionname=all_region{regionind}
        switch(regionname)
            case('NWP') %% North western Pacific
                lonlat = [115, 164, 15, 52];  %% whole data area
                refpolygon(1,1)=lonlat(1);
                refpolygon(2,1)=lonlat(2);
                refpolygon(1,2)=lonlat(3);
                refpolygon(2,2)=lonlat(4);
            case('ES') %% East Sea
                refpolygon=espolygon;
            case('SS') %% South Sea
                refpolygon=sspolygon;
            case('YS') %% Yellow Sea
                refpolygon=yspolygon;
            case('ECS') %% East China Sea
                refpolygon=ecspolygon;
            case('AKP') %% Around Korea Peninsula
                refpolygon=akppolygon;
            case('AKP2') %% Around Korea Peninsula
                refpolygon=akp2polygon;
            case('EKB') %% Around Korea Peninsula
                refpolygon=ekbpolygon;
            otherwise
                ('?')
        end
        lonlat(1)=min(refpolygon(:,1));
        lonlat(2)=max(refpolygon(:,1));
        lonlat(3)=min(refpolygon(:,2));
        lonlat(4)=max(refpolygon(:,2));
        % % % for EKB
        % regionname='EKB';
        % lonlat = [127, 129.5, 38, 40.5];

        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
            param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
            filedir = strcat('E:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            avhrrdir='E:\Data\Observation\OISST\monthly\';
        elseif (strcmp(system_name,'GLNXA64'))
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
                filename = strcat(filedir, num2str(tempyear,'%04i'), '\', ...
                        testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');
                % read model data
                if (exist('lon')==0)
                    modelinfo=ncinfo(filename);
                    lon = ncread(filename,'lon_rho',[1 1],[modelinfo.Dimensions(5).Length,1]);
                    lat = ncread(filename,'lat_rho',[1 1],[1,modelinfo.Dimensions(6).Length]);

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

                    lon = ncread(filename,'lon_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                    lat = ncread(filename,'lat_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);

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

                data = ncread(filename,varname,[lon_min(1) lat_min(1) 40 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);
                data=data.*mask_model;

                % read OISST DATA
                avhrrfilename = strcat(avhrrdir,'avhrr_monthly', num2str(1982,'%04i'), ...
                        '_',num2str(tempmonth,'%02i'), '.nc');
                if (exist('avhrr_lon')==0)
                    avhrrinfo=ncinfo(avhrrfilename);
                    avhrr_lon = ncread(avhrrfilename,'long',[1 1],[avhrrinfo.Dimensions(2).Length,1]);
                    avhrr_lat = ncread(avhrrfilename,'lat',[1 1],[1,avhrrinfo.Dimensions(1).Length]);

                    avhrr_lon_west = abs(avhrr_lon - (lonlat(1)));
                    min_avhrr_lon_west=min(avhrr_lon_west);
                    avhrr_lon_east = abs(avhrr_lon - (lonlat(2)));
                    min_avhrr_lon_east=min(avhrr_lon_east);
                    avhrr_lat_south = abs(avhrr_lat - (lonlat(3)));
                    min_avhrr_lat_south=min(avhrr_lat_south);
                    avhrr_lat_north = abs(avhrr_lat - (lonlat(4)));
                    min_avhrr_lat_north=min(avhrr_lat_north);

                    avhrr_lon_min = find(avhrr_lon_west == min_avhrr_lon_west);
                    avhrr_lon_max = find(avhrr_lon_east == min_avhrr_lon_east);
                    avhrr_lat_min = find(avhrr_lat_south == min_avhrr_lat_south);
                    avhrr_lat_max = find(avhrr_lat_north == min_avhrr_lat_north);

            %         ncinfo('E:\Data\Observation\OISST\monthly\avhrr_monthly1983_11.nc');

                    avhrr_lon = ncread(avhrrfilename,'long', [avhrr_lon_min(1) avhrr_lat_min(1)], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1)]);
                    avhrr_lat = ncread(avhrrfilename,'lat', [avhrr_lon_min(1) avhrr_lat_min(1)], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1)]);

                    comb_spatial_meanrms=(zeros([size(avhrr_lon),12]));
                    comb_spatial_meanbias=(zeros([size(avhrr_lon),12]));
                    comb_spatial_meanavhrr=(zeros([size(avhrr_lon),12]));
                    comb_spatial_meanmodel=(zeros([size(avhrr_lon),12]));

                    switch(regionname)
                        case('NWP') %% North western Pacific
                            mask_avhrr(1:size(avhrr_lon,1),1:size(avhrr_lon,2))=1;
                            mask_avhrr(mask_avhrr==0)=NaN;
                        otherwise
                            mask_avhrr = double(inpolygon(avhrr_lon,avhrr_lat,refpolygon(:,1),refpolygon(:,2)));
                            mask_avhrr(mask_avhrr==0)=NaN;
                    end
                end
                len_lon = length(avhrr_lon(:,1));
                len_lat = length(avhrr_lat(1,:));
                len_lon_model = size(data,1);
                len_lat_model = size(data,2);

                interped_data = griddata(double(lon), double(lat), data,double(avhrr_lon),double(avhrr_lat));   
                
                comb_interped_data(:,:,ind) = interped_data;
                comb_spatial_meanmodel(:,:,monthij)=comb_spatial_meanmodel(:,:,monthij)+interped_data/double(length(inputyear));

                comb_data(:,:,ind) = data;
                comb_interped_data(:,:,ind) = interped_data;

                ind = ind + 1;
                toc;
            end
        end

        trendtime=inputyear(1):1/12:inputyear(end)+1-1/12;
        trend(1:len_lon,1:len_lat)=NaN;
        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(trendtime,squeeze(comb_interped_data(i,j,:))',1);
                trend(i,j)=p(1);
            end
        end

        
        for t=1:length(inputyear)
            comb_interped_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_interped_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meanmodel;
        end

        trend_filtered(1:len_lon,1:len_lat)=NaN;
        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(trendtime,squeeze(comb_interped_data_filtered(i,j,:))',1);
                trend_filtered(i,j)=p(1);
            end
        end


        mean_trend=mean(mean(trend,'omitnan'),'omitnan');
        mean_trend_filtered=mean(mean(trend_filtered,'omitnan'),'omitnan');

        figdir=[figrawdir,'CLIM\'];
        outfile = strcat(figdir,regionname);
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        
        save([filedir,regionname,'sss_model_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

        ncid = netcdf.create(strcat(filedir, testname, '_',regionname,'_clim_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc'),'NETCDF4');


        onedimid = netcdf.defDim(ncid,'one', 1);
        lon_dimid = netcdf.defDim(ncid, 'lon', len_lon);
        lat_dimid = netcdf.defDim(ncid,'lat',len_lat);
        xidimid = netcdf.defDim(ncid, 'xi_rho', len_lon_model);
        etadimid = netcdf.defDim(ncid,'eta_rho',len_lat_model);
        time_dimid = netcdf.defDim(ncid, 'time', 0);
        clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
            'type', ['NWP 1/20 _ ', testname, 'monthly SSS file']);
        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
            'title', [' monthly SSS  (', num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ']);
        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
            'source', [' ROMS NWP 1/20 data from _ ',testname ]);
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

        raw_sssvarid=netcdf.defVar(ncid, 'raw_sss', 'NC_FLOAT', [xidimid etadimid time_dimid]);
        netcdf.putAtt(ncid,raw_sssvarid,'long_name','raw_sss');
        netcdf.putAtt(ncid,raw_sssvarid,'units','Celsius');

        sssvarid=netcdf.defVar(ncid, 'sss', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,sssvarid,'long_name','sss');
        netcdf.putAtt(ncid,sssvarid,'units','Celsius');

        sss_filteredvarid=netcdf.defVar(ncid, 'sss_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,sss_filteredvarid,'long_name','sss_filtered');
        netcdf.putAtt(ncid,sss_filteredvarid,'units','Celsius');

        trendvarid=netcdf.defVar(ncid, 'trend', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,trendvarid,'long_name','trend');
        netcdf.putAtt(ncid,trendvarid,'units','Celsius/year');

        trend_filteredvarid=netcdf.defVar(ncid, 'trend_filtered', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,trend_filteredvarid,'long_name','trend_filtered');
        netcdf.putAtt(ncid,trend_filteredvarid,'units','Celsius/year');

        mean_trendvarid=netcdf.defVar(ncid, 'mean_trend', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_trendvarid,'long_name','mean_trend');
        netcdf.putAtt(ncid,mean_trendvarid,'units','Celsius/year');

        mean_trend_filteredvarid=netcdf.defVar(ncid, 'mean_trend_filtered', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_trend_filteredvarid,'long_name','mean_trend_filtered');
        netcdf.putAtt(ncid,mean_trend_filteredvarid,'units','Celsius/year');

        clim_sssvarid=netcdf.defVar(ncid, 'clim_sss', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_sssvarid,'long_name','clim_sss');
        netcdf.putAtt(ncid,clim_sssvarid,'units','Celsius');

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
        netcdf.putVar(ncid, lonvarid, 0, len_lon, avhrr_lon(:,1));
        netcdf.putVar(ncid, latvarid, 0, len_lat, avhrr_lat(1,:));
        netcdf.putVar(ncid, lon_rhovarid, [0 0], [len_lon_model len_lat_model], lon);
        netcdf.putVar(ncid, lat_rhovarid, [0 0], [len_lon_model len_lat_model], lat);
        netcdf.putVar(ncid, raw_sssvarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_data);
        netcdf.putVar(ncid, sssvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_interped_data);
        netcdf.putVar(ncid, sss_filteredvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_interped_data_filtered);
        netcdf.putVar(ncid, clim_sssvarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meanmodel);
        netcdf.putVar(ncid, trendvarid, [0 0], [len_lon len_lat], trend);
        netcdf.putVar(ncid, trend_filteredvarid, [0 0], [len_lon len_lat], trend_filtered);
        netcdf.putVar(ncid, mean_trendvarid, [0], [1], mean_trend);
        netcdf.putVar(ncid, mean_trend_filteredvarid, [0], [1], mean_trend_filtered);

        netcdf.close(ncid);
    end
end

% SSH_4th_mid_report3