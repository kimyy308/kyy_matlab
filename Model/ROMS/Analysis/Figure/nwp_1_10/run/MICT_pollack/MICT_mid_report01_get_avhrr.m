close all; clear all;  clc;   
% all_region ={'NWP','ES', 'SS', 'YS', 'ECS'}
all_region ={'NWP', 'AKP2'}

for regionind=1:length(all_region)
    clearvars '*' -except regionind all_region

    % % % 
    % % % Read Model SST
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
    elseif (strcmp(system_name,'GLNXA64'))
        dropboxpath='/home/kimyy/Dropbox';
        addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
        addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
        addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    end


    shadlev = [0 35];
    rms_shadlev = [0 4];
    bias_shadlev = [-4 4];
    conlev  = 0:5:35;
    dl=1/20;
    % for snu_desktop
    testname='test06'   % % need to change
    inputyear = [1989:2017]; % % put year which you want to plot [year year ...]
    inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
%     inputmonth = [1 2 3]; % % put month which you want to plot [month month ...]


    
    varname ='temp';  %% reference variable -> temperature
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
        case('EKB') %% North western Pacific
            lonlat = [127, 131, 37, 42];  %% East Korea Bay
            refpolygon(1,1)=lonlat(1);
            refpolygon(2,1)=lonlat(2);
            refpolygon(1,2)=lonlat(3);
            refpolygon(2,2)=lonlat(4);
        case('AKP2')
            refpolygon=akp2polygon;
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
        figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd_year\figure\nwp_1_10\',testname,'\DA\',regionname,'\'); % % where figure files will be saved
        param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
%         filedir = strcat('E:\Data\Observation\OISST\monthly\'); % % where data files are
%         avhrrdir='E:\Data\Observation\OISST\monthly\';
        filedir = strcat('E:\Data\Observation\OISST\monthly_kimyy\'); % % where data files are
        avhrrdir='E:\Data\Observation\OISST\monthly_kimyy\';
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
                
                if (exist('comb_spatial_meanrms')==0)
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
                avhrr_data = ncread(avhrrfilename,'temp',[avhrr_lon_min(1) avhrr_lat_min(1) monthij], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1) 1]);

            
            
            avhrr_data(avhrr_data<-9)=NaN;
            avhrr_data(avhrr_data>1000)=NaN;
            avhrr_data=avhrr_data.*mask_avhrr;
            
            comb_avhrr_data(:,:,ind) = avhrr_data;
            comb_spatial_meanavhrr(:,:,monthij)=comb_spatial_meanavhrr(:,:,monthij)+avhrr_data/double(length(inputyear));

            ind = ind + 1;
            toc;
        end
    end
    comb_avhrr_clim_divided=reshape(comb_avhrr_data,[len_lon, len_lat, length(inputmonth), length(inputyear)]);
    
    trendtime=inputyear(1):1/12:inputyear(end)+1-1/12;
    avhrr_trend(1:len_lon,1:len_lat)=NaN;
    for i=1:len_lon
        for j=1:len_lat
            p=polyfit(trendtime,squeeze(comb_avhrr_data(i,j,:))',1);
            avhrr_trend(i,j)=p(1);
        end
    end

    for t=1:length(inputyear)
        comb_avhrr_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_avhrr_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meanavhrr;
    end

    avhrr_trend_filtered(1:len_lon,1:len_lat)=NaN;

    for i=1:len_lon
        for j=1:len_lat
            p=polyfit(trendtime,squeeze(comb_avhrr_data_filtered(i,j,:))',1);
            avhrr_trend_filtered(i,j)=p(1);
        end
    end
    
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

    mean_avhrr_trend=mean(mean(avhrr_trend,'omitnan'),'omitnan');
    mean_avhrr_trend_filtered=mean(mean(avhrr_trend_filtered,'omitnan'),'omitnan');
    mean_clim_avhrr_trend_divided=mean(mean(clim_avhrr_trend_divided,'omitnan'),'omitnan');
    figdir=[figrawdir,'CLIM\'];
    outfile = strcat(figdir,regionname);
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 

    save([filedir,regionname,'sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

    ncid = netcdf.create(strcat(filedir, testname, '_',regionname,'_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc'),'NETCDF4');

    onedimid = netcdf.defDim(ncid,'one', 1);
    lon_dimid = netcdf.defDim(ncid, 'lon', len_lon);
    lat_dimid = netcdf.defDim(ncid,'lat',len_lat);
    time_dimid = netcdf.defDim(ncid, 'time', 0);
    clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ['AVHRR(OISST) _ ', testname, 'monthly SST trend file']);
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', [' AVHRR SST trend (',num2str(inputyear(1)),'-',num2str(inputyear(end)),') ']);
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', [' AVHRR(OISST) satellite data ' ]);
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

    avhrr_sstvarid=netcdf.defVar(ncid, 'avhrr_sst', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,avhrr_sstvarid,'long_name','avhrr_sst');
    netcdf.putAtt(ncid,avhrr_sstvarid,'units','Celsius');

    avhrr_sst_filteredvarid=netcdf.defVar(ncid, 'avhrr_sst_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,avhrr_sst_filteredvarid,'long_name','avhrr_sst_filtered');
    netcdf.putAtt(ncid,avhrr_sst_filteredvarid,'units','Celsius');

    avhrr_trendvarid=netcdf.defVar(ncid, 'avhrr_trend', 'NC_FLOAT', [lon_dimid lat_dimid]);
    netcdf.putAtt(ncid,avhrr_trendvarid,'long_name','avhrr_trend');
    netcdf.putAtt(ncid,avhrr_trendvarid,'units','Celsius/year');

    avhrr_trend_filteredvarid=netcdf.defVar(ncid, 'avhrr_trend_filtered', 'NC_FLOAT', [lon_dimid lat_dimid]);
    netcdf.putAtt(ncid,avhrr_trend_filteredvarid,'long_name','avhrr_trend_filtered');
    netcdf.putAtt(ncid,avhrr_trend_filteredvarid,'units','Celsius/year');
    
    clim_avhrr_trend_dividedvarid=netcdf.defVar(ncid, 'clim_avhrr_trend_divided', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_avhrr_trend_dividedvarid,'long_name','clim_avhrr_trend_divided');
    netcdf.putAtt(ncid,clim_avhrr_trend_dividedvarid,'units','Celsius/year');

    mean_avhrr_trendvarid=netcdf.defVar(ncid, 'mean_avhrr_trend', 'NC_FLOAT', onedimid);
    netcdf.putAtt(ncid,mean_avhrr_trendvarid,'long_name','mean_avhrr_trend');
    netcdf.putAtt(ncid,mean_avhrr_trendvarid,'units','Celsius/year');

    mean_avhrr_trend_filteredvarid=netcdf.defVar(ncid, 'mean_avhrr_trend_filtered', 'NC_FLOAT', onedimid);
    netcdf.putAtt(ncid,mean_avhrr_trend_filteredvarid,'long_name','mean_avhrr_trend_filtered');
    netcdf.putAtt(ncid,mean_avhrr_trend_filteredvarid,'units','Celsius/year');
    
    mean_clim_avhrr_trend_dividedvarid=netcdf.defVar(ncid, 'mean_clim_avhrr_trend_divided', 'NC_FLOAT', clim_time_dimid);
    netcdf.putAtt(ncid,mean_clim_avhrr_trend_dividedvarid,'long_name','mean_clim_avhrr_trend_divided');
    netcdf.putAtt(ncid,mean_clim_avhrr_trend_dividedvarid,'units','Celsius/year');

    clim_avhrrvarid=netcdf.defVar(ncid, 'clim_avhrr_sst', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_avhrrvarid,'long_name','clim_avhrr_sst');
    netcdf.putAtt(ncid,clim_avhrrvarid,'units','Celsius');

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
    if (tempyear<2011)
        netcdf.putVar(ncid, lonvarid, 0, len_lon, avhrr_lon(:,1));
        netcdf.putVar(ncid, latvarid, 0, len_lat, avhrr_lat(1,:));
    else
        netcdf.putVar(ncid, lonvarid, 0, len_lon, avhrr_lon(:));
        netcdf.putVar(ncid, latvarid, 0, len_lat, avhrr_lat(:));
    end
    netcdf.putVar(ncid, avhrr_sstvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_avhrr_data);
    netcdf.putVar(ncid, avhrr_sst_filteredvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_avhrr_data_filtered);
    netcdf.putVar(ncid, clim_avhrrvarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meanavhrr);
    netcdf.putVar(ncid, avhrr_trendvarid, [0 0], [len_lon len_lat], avhrr_trend);
    netcdf.putVar(ncid, avhrr_trend_filteredvarid, [0 0], [len_lon len_lat], avhrr_trend_filtered);
    netcdf.putVar(ncid, clim_avhrr_trend_dividedvarid, [0 0 0], [len_lon len_lat length(climtime)], clim_avhrr_trend_divided);
    netcdf.putVar(ncid, mean_avhrr_trendvarid, [0], [1], mean_avhrr_trend);
    netcdf.putVar(ncid, mean_avhrr_trend_filteredvarid, [0], [1], mean_avhrr_trend_filtered);
    netcdf.putVar(ncid, mean_clim_avhrr_trend_dividedvarid, [0], [length(climtime)], mean_clim_avhrr_trend_divided);

    netcdf.close(ncid);
end
