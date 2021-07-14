close all; clear all;  clc;   
% all_region ={'NWP','ES', 'SS', 'YS', 'ECS'}
all_region ={'NWP3'}

for regionind=1:length(all_region)
    clearvars '*' -except regionind all_region

    % % % 
    % % % Read Model ssh
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
%     testname='test49'   % % need to change
    inputyear = [1994:2017]; % % put year which you want to plot [year year ...]
    inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
    

    
    varname ='SOSSHEIG';  %% reference variable -> temperature
    run('nwp_polygon_point.m');
    regionname=all_region{regionind}
    switch(regionname)
        case('NWP') %% North western Pacific
            lonlat = [115, 164, 15, 52];  %% whole data area
            refpolygon(1,1)=lonlat(1);
            refpolygon(2,1)=lonlat(2);
            refpolygon(1,2)=lonlat(3);
            refpolygon(2,2)=lonlat(4);
        case('NWP2') %% North western Pacific
            lonlat = [115, 145, 25, 52];  %% whole data area
            refpolygon(1,1)=lonlat(1);
            refpolygon(2,1)=lonlat(2);
            refpolygon(1,2)=lonlat(3);
            refpolygon(2,2)=lonlat(4);
        case('NWP3') %% North western Pacific
            lonlat = [114, 166, 14, 53];  %% whole data area
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
        figrawdir =strcat('D:\OneDrive - ������б�\MEPL\project\MICT_pollack\3rd year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
        param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
        filedir = strcat('E:\Data\Observation\OIssh\monthly\'); % % where data files are
        ORAS5dir='E:\Data\Observation\OIssh\monthly\';
    elseif (strcmp(system_name,'GLNXA64'))
        figrawdir =strcat('/data1/stlee/ext_hdd/ORAS5/',regionname,'/'); % % where figure files will be saved
        param_script ='/home/kimyy/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run/fig_param/fig_param_kyy_EKB_RMS.m'
        filedir = strcat('/data2/kimyy/Reanalysis/ORAS5/'); % % where data files are
        ORAS5dir='/data2/kimyy/Reanalysis/ORAS5/';
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

            % read OIssh DATA
            ORAS5filename = strcat(ORAS5dir,'ORAS5_SSH_', num2str(min(inputyear),'%04i'),'_', num2str(max(inputyear), '%04i'), '.nc');
            if (exist('ORAS5_lon')==0)
                ORAS5info=ncinfo(ORAS5filename);
%                 ORAS5_lon = ncread(ORAS5filename, 'LON111_171',[1],[ORAS5info.Dimensions(3).Length]);
%                 ORAS5_lat = ncread(ORAS5filename, 'LAT101_150',[1],[ORAS5info.Dimensions(2).Length]);
                
                ORAS5_lon = ncread(ORAS5filename, 'LON111_171',[1],[inf]);
                ORAS5_lat = ncread(ORAS5filename, 'LAT101_150',[1],[inf]);
                
                ORAS5_lon_west = abs(ORAS5_lon - (lonlat(1)));
                min_ORAS5_lon_west=min(ORAS5_lon_west);
                ORAS5_lon_east = abs(ORAS5_lon - (lonlat(2)));
                min_ORAS5_lon_east=min(ORAS5_lon_east);
                ORAS5_lat_south = abs(ORAS5_lat - (lonlat(3)));
                min_ORAS5_lat_south=min(ORAS5_lat_south);
                ORAS5_lat_north = abs(ORAS5_lat - (lonlat(4)));
                min_ORAS5_lat_north=min(ORAS5_lat_north);

                ORAS5_lon_min = find(ORAS5_lon_west == min_ORAS5_lon_west);
                ORAS5_lon_max = find(ORAS5_lon_east == min_ORAS5_lon_east);
                ORAS5_lat_min = find(ORAS5_lat_south == min_ORAS5_lat_south);
                ORAS5_lat_max = find(ORAS5_lat_north == min_ORAS5_lat_north);

        %         ncinfo('E:\Data\Observation\OIssh\monthly\ORAS5_monthly1983_11.nc');

                ORAS5_lon = ncread(ORAS5filename,'LON111_171', [ORAS5_lon_min(1)], [ORAS5_lon_max(1)-ORAS5_lon_min(1)]);
                ORAS5_lat = ncread(ORAS5filename,'LAT101_150', [ORAS5_lat_min(1)], [ORAS5_lat_max(1)-ORAS5_lat_min(1)]);

                comb_spatial_meanrms=(zeros([length(ORAS5_lon),length(ORAS5_lat),12]));
                comb_spatial_meanbias=(zeros([length(ORAS5_lon),length(ORAS5_lat),12]));
                comb_spatial_meanORAS5=(zeros([length(ORAS5_lon),length(ORAS5_lat),12]));
                comb_spatial_meanmodel=(zeros([length(ORAS5_lon),length(ORAS5_lat),12]));
                
                [ORAS5_lat2 ORAS5_lon2]=meshgrid(ORAS5_lat, ORAS5_lon);
                switch(regionname)
                    case('NWP') %% North western Pacific
                        mask_ORAS5(1:size(ORAS5_lon,1),1:size(ORAS5_lon,2))=1;
                    otherwise
                        mask_ORAS5 = double(inpolygon(ORAS5_lon2,ORAS5_lat2,refpolygon(:,1),refpolygon(:,2)));
                        mask_ORAS5(mask_ORAS5==0)=NaN;
                end
            end
            len_lon = length(ORAS5_lon(:));
            len_lat = length(ORAS5_lat(:));

            ORAS5_data = ncread(ORAS5filename,varname, ...
                [ORAS5_lon_min(1) ORAS5_lat_min(1) (yearij-1)*12+tempmonth], ...
                [ORAS5_lon_max(1)-ORAS5_lon_min(1) ORAS5_lat_max(1)-ORAS5_lat_min(1) 1]);
            ORAS5_data(ORAS5_data<-1000)=NaN;
            ORAS5_data(ORAS5_data>1000)=NaN;
            ORAS5_data=ORAS5_data.*mask_ORAS5;

            comb_ORAS5_data(:,:,ind) = ORAS5_data;
            comb_spatial_meanORAS5(:,:,monthij)=comb_spatial_meanORAS5(:,:,monthij)+ORAS5_data/double(length(inputyear));

            ind = ind + 1;
            toc;
        end
    end
    comb_ORAS5_clim_divided=reshape(comb_ORAS5_data,[len_lon, len_lat, 12, length(inputyear)]);
    
    trendtime=inputyear(1):1/12:inputyear(end)+1-1/12;
    ORAS5_trend(1:len_lon,1:len_lat)=NaN;
    for i=1:len_lon
        for j=1:len_lat
            p=polyfit(trendtime,squeeze(comb_ORAS5_data(i,j,:))',1);
            ORAS5_trend(i,j)=p(1) * 1000.0 ;
        end
    end

    for t=1:length(inputyear)
        comb_ORAS5_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_ORAS5_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meanORAS5;
    end

    ORAS5_trend_filtered(1:len_lon,1:len_lat)=NaN;

    for i=1:len_lon
        for j=1:len_lat
            p=polyfit(trendtime,squeeze(comb_ORAS5_data_filtered(i,j,:))',1);
            ORAS5_trend_filtered(i,j)=p(1) * 1000.0 ;
        end
    end
    
    clim_ORAS5_trend_filtered(1:len_lon,1:len_lat,1:12)=NaN;
    for k=1:12
        clim_trendtime(:,k)=inputyear(1)+(k-1)/12 : 1 : inputyear(end)+(k-1)/12;
        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(clim_trendtime(:,k)',squeeze(comb_ORAS5_clim_divided(i,j,k,:))',1);
                clim_ORAS5_trend_divided(i,j,k)=p(1) * 1000.0 ;
            end
        end
    end

    mean_ORAS5_trend=mean(mean(ORAS5_trend,'omitnan'),'omitnan');
    mean_ORAS5_trend_filtered=mean(mean(ORAS5_trend_filtered,'omitnan'),'omitnan');
    mean_clim_ORAS5_trend_divided=mean(mean(clim_ORAS5_trend_divided,'omitnan'),'omitnan');
    figdir=[figrawdir,'CLIM\'];
    outfile = strcat(figdir,regionname);
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 

    save([filedir,regionname,'ORAS5_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

    ncid = netcdf.create(strcat(filedir,regionname,'ORAS5_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc'),'NETCDF4');

    onedimid = netcdf.defDim(ncid,'one', 1);
    lon_dimid = netcdf.defDim(ncid, 'lon', len_lon);
    lat_dimid = netcdf.defDim(ncid,'lat',len_lat);
    time_dimid = netcdf.defDim(ncid, 'time', 0);
    clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ['ORAS5 _ ', 'monthly ssh trend file']);
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', [' ORAS5 ssh trend (',num2str(inputyear(1)),'-',num2str(inputyear(end)),') ']);
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', [' ORAS5(OIssh) satellite data ' ]);
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

    ORAS5_sshvarid=netcdf.defVar(ncid, 'ORAS5_ssh', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,ORAS5_sshvarid,'long_name','ORAS5_ssh');
    netcdf.putAtt(ncid,ORAS5_sshvarid,'units','mm ');

    ORAS5_ssh_filteredvarid=netcdf.defVar(ncid, 'ORAS5_ssh_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,ORAS5_ssh_filteredvarid,'long_name','ORAS5_ssh_filtered');
    netcdf.putAtt(ncid,ORAS5_ssh_filteredvarid,'units','mm ');

    ORAS5_trendvarid=netcdf.defVar(ncid, 'ORAS5_trend', 'NC_FLOAT', [lon_dimid lat_dimid]);
    netcdf.putAtt(ncid,ORAS5_trendvarid,'long_name','ORAS5_trend');
    netcdf.putAtt(ncid,ORAS5_trendvarid,'units','mm /year');

    ORAS5_trend_filteredvarid=netcdf.defVar(ncid, 'ORAS5_trend_filtered', 'NC_FLOAT', [lon_dimid lat_dimid]);
    netcdf.putAtt(ncid,ORAS5_trend_filteredvarid,'long_name','ORAS5_trend_filtered');
    netcdf.putAtt(ncid,ORAS5_trend_filteredvarid,'units','mm /year');
    
    clim_ORAS5_trend_dividedvarid=netcdf.defVar(ncid, 'clim_ORAS5_trend_divided', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_ORAS5_trend_dividedvarid,'long_name','clim_ORAS5_trend_divided');
    netcdf.putAtt(ncid,clim_ORAS5_trend_dividedvarid,'units','mm /year');

    mean_ORAS5_trendvarid=netcdf.defVar(ncid, 'mean_ORAS5_trend', 'NC_FLOAT', onedimid);
    netcdf.putAtt(ncid,mean_ORAS5_trendvarid,'long_name','mean_ORAS5_trend');
    netcdf.putAtt(ncid,mean_ORAS5_trendvarid,'units','mm /year');

    mean_ORAS5_trend_filteredvarid=netcdf.defVar(ncid, 'mean_ORAS5_trend_filtered', 'NC_FLOAT', onedimid);
    netcdf.putAtt(ncid,mean_ORAS5_trend_filteredvarid,'long_name','mean_ORAS5_trend_filtered');
    netcdf.putAtt(ncid,mean_ORAS5_trend_filteredvarid,'units','mm /year');
    
    mean_clim_ORAS5_trend_dividedvarid=netcdf.defVar(ncid, 'mean_clim_ORAS5_trend_divided', 'NC_FLOAT', clim_time_dimid);
    netcdf.putAtt(ncid,mean_clim_ORAS5_trend_dividedvarid,'long_name','mean_clim_ORAS5_trend_divided');
    netcdf.putAtt(ncid,mean_clim_ORAS5_trend_dividedvarid,'units','mm /year');

    clim_ORAS5varid=netcdf.defVar(ncid, 'clim_ORAS5_ssh', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_ORAS5varid,'long_name','clim_ORAS5_ssh');
    netcdf.putAtt(ncid,clim_ORAS5varid,'units','mm ');

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
    netcdf.putVar(ncid, lonvarid, 0, len_lon, ORAS5_lon(:));
    netcdf.putVar(ncid, latvarid, 0, len_lat, ORAS5_lat(:));
    netcdf.putVar(ncid, ORAS5_sshvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_ORAS5_data);
    netcdf.putVar(ncid, ORAS5_ssh_filteredvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_ORAS5_data_filtered);
    netcdf.putVar(ncid, clim_ORAS5varid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meanORAS5);
    netcdf.putVar(ncid, ORAS5_trendvarid, [0 0], [len_lon len_lat], ORAS5_trend);
    netcdf.putVar(ncid, ORAS5_trend_filteredvarid, [0 0], [len_lon len_lat], ORAS5_trend_filtered);
    netcdf.putVar(ncid, clim_ORAS5_trend_dividedvarid, [0 0 0], [len_lon len_lat length(climtime)], clim_ORAS5_trend_divided);
    netcdf.putVar(ncid, mean_ORAS5_trendvarid, [0], [1], mean_ORAS5_trend);
    netcdf.putVar(ncid, mean_ORAS5_trend_filteredvarid, [0], [1], mean_ORAS5_trend_filtered);
    netcdf.putVar(ncid, mean_clim_ORAS5_trend_dividedvarid, [0], [length(climtime)], mean_clim_ORAS5_trend_divided);

    netcdf.close(ncid);
end