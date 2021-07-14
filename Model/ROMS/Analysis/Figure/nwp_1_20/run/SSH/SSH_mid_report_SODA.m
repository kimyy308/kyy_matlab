close all; clear all;  clc;   
% all_region ={'NWP','ES', 'SS', 'YS', 'ECS'}
all_region ={'NWP3'}
% all_region ={'GLO'}

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
    testname='test49'   % % need to change
    inputyear = [1994:2014]; % % put year which you want to plot [year year ...]
    inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
    

    
    varname ='ssh';  %% reference variable -> temperature
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
        case('GLO')
            lonlat = [-360, 360, -180, 180];  %% East Korea Bay
            refpolygon(1,1)=lonlat(1);
            refpolygon(2,1)=lonlat(2);
            refpolygon(1,2)=lonlat(3);
            refpolygon(2,2)=lonlat(4);
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
        sodadir='E:\Data\Observation\OIssh\monthly\';
    elseif (strcmp(system_name,'GLNXA64'))
        figrawdir =strcat('/data2/kimyy/Reanalysis/SODA/SODA_3_4_2/monthly/figures/',regionname,'/'); % % where figure files will be saved
        param_script ='/home/kimyy/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run/fig_param/fig_param_kyy_EKB_RMS.m'
        filedir = strcat('/data2/kimyy/ext_hdd/Reanalysis/SODA/SODA_3_4_2/monthly/'); % % where data files are
        sodadir='/data2/kimyy/ext_hdd/Reanalysis/SODA/SODA_3_4_2/monthly/';
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

            % read ssh DATA
            sodafilename = strcat(sodadir,'soda3.4.2_mn_ocean_reg_', num2str(tempyear,'%04i'),'.nc');
            if (exist('soda_lon')==0)
                sodainfo=ncinfo(sodafilename);
                soda_lon = ncread(sodafilename,'xt_ocean',[1],[sodainfo.Dimensions(1).Length]);
                soda_lat = ncread(sodafilename,'yt_ocean',[1],[sodainfo.Dimensions(2).Length]);

                soda_lon_west = abs(soda_lon - (lonlat(1)));
                min_soda_lon_west=min(soda_lon_west);
                soda_lon_east = abs(soda_lon - (lonlat(2)));
                min_soda_lon_east=min(soda_lon_east);
                soda_lat_south = abs(soda_lat - (lonlat(3)));
                min_soda_lat_south=min(soda_lat_south);
                soda_lat_north = abs(soda_lat - (lonlat(4)));
                min_soda_lat_north=min(soda_lat_north);

                soda_lon_min = find(soda_lon_west == min_soda_lon_west);
                soda_lon_max = find(soda_lon_east == min_soda_lon_east);
                soda_lat_min = find(soda_lat_south == min_soda_lat_south);
                soda_lat_max = find(soda_lat_north == min_soda_lat_north);

        %         ncinfo('E:\Data\Observation\OIssh\monthly\soda_monthly1983_11.nc');

                soda_lon = ncread(sodafilename,'xt_ocean', [soda_lon_min(1)], [soda_lon_max(1)-soda_lon_min(1)]);
                soda_lat = ncread(sodafilename,'yt_ocean', [soda_lat_min(1)], [soda_lat_max(1)-soda_lat_min(1)]);

                comb_spatial_meanrms=(zeros([length(soda_lon),length(soda_lat),12]));
                comb_spatial_meanbias=(zeros([length(soda_lon),length(soda_lat),12]));
                comb_spatial_meansoda=(zeros([length(soda_lon),length(soda_lat),12]));
                comb_spatial_meanmodel=(zeros([length(soda_lon),length(soda_lat),12]));
                
                [soda_lat2 soda_lon2]=meshgrid(soda_lat, soda_lon);
                switch(regionname)
                    case('NWP') %% North western Pacific
                        mask_soda(1:size(soda_lon,1),1:size(soda_lon,2))=1;
                    otherwise
                        mask_soda = double(inpolygon(soda_lon2,soda_lat2,refpolygon(:,1),refpolygon(:,2)));
                        mask_soda(mask_soda==0)=NaN;
                end
            end
            len_lon = length(soda_lon(:));
            len_lat = length(soda_lat(:));

            soda_data = ncread(sodafilename,varname,[soda_lon_min(1) soda_lat_min(1)  tempmonth], [soda_lon_max(1)-soda_lon_min(1) soda_lat_max(1)-soda_lat_min(1) 1]);
            soda_data(soda_data<-1000)=NaN;
            soda_data(soda_data>1000)=NaN;
            soda_data=soda_data.*mask_soda;

            comb_soda_data(:,:,ind) = soda_data;
            comb_spatial_meansoda(:,:,monthij)=comb_spatial_meansoda(:,:,monthij)+soda_data/double(length(inputyear));

            ind = ind + 1;
            toc;
        end
    end
    comb_soda_clim_divided=reshape(comb_soda_data,[len_lon, len_lat, 12, length(inputyear)]);
    
    trendtime=inputyear(1):1/12:inputyear(end)+1-1/12;
    soda_trend(1:len_lon,1:len_lat)=NaN;
    for i=1:len_lon
        for j=1:len_lat
            p=polyfit(trendtime,squeeze(comb_soda_data(i,j,:))',1);
            soda_trend(i,j)=p(1) * 1000.0 ;
        end
    end

    for t=1:length(inputyear)
        comb_soda_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_soda_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meansoda;
    end

    soda_trend_filtered(1:len_lon,1:len_lat)=NaN;

    for i=1:len_lon
        for j=1:len_lat
            p=polyfit(trendtime,squeeze(comb_soda_data_filtered(i,j,:))',1);
            soda_trend_filtered(i,j)=p(1) * 1000.0 ;
        end
    end
    
    clim_soda_trend_filtered(1:len_lon,1:len_lat,1:12)=NaN;
    for k=1:12
        clim_trendtime(:,k)=inputyear(1)+(k-1)/12 : 1 : inputyear(end)+(k-1)/12;
        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(clim_trendtime(:,k)',squeeze(comb_soda_clim_divided(i,j,k,:))',1);
                clim_soda_trend_divided(i,j,k)=p(1) * 1000.0 ;
            end
        end
    end

    mean_soda_trend=mean(mean(soda_trend,'omitnan'),'omitnan');
    mean_soda_trend_filtered=mean(mean(soda_trend_filtered,'omitnan'),'omitnan');
    mean_clim_soda_trend_divided=mean(mean(clim_soda_trend_divided,'omitnan'),'omitnan');
    figdir=[figrawdir,'CLIM\'];
    outfile = strcat(figdir,regionname);
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 

    save([filedir,regionname,'soda_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

    ncid = netcdf.create(strcat(filedir,regionname,'soda_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc'),'NETCDF4');

    onedimid = netcdf.defDim(ncid,'one', 1);
    lon_dimid = netcdf.defDim(ncid, 'lon', len_lon);
    lat_dimid = netcdf.defDim(ncid,'lat',len_lat);
    time_dimid = netcdf.defDim(ncid, 'time', 0);
    clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ['soda(OIssh) _ ', testname, 'monthly ssh trend file']);
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', [' soda ssh trend (',num2str(inputyear(1)),'-',num2str(inputyear(end)),') ']);
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', [' soda(OIssh) satellite data ' ]);
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

    soda_sshvarid=netcdf.defVar(ncid, 'soda_ssh', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,soda_sshvarid,'long_name','soda_ssh');
    netcdf.putAtt(ncid,soda_sshvarid,'units','mm ');

    soda_ssh_filteredvarid=netcdf.defVar(ncid, 'soda_ssh_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,soda_ssh_filteredvarid,'long_name','soda_ssh_filtered');
    netcdf.putAtt(ncid,soda_ssh_filteredvarid,'units','mm ');

    soda_trendvarid=netcdf.defVar(ncid, 'soda_trend', 'NC_FLOAT', [lon_dimid lat_dimid]);
    netcdf.putAtt(ncid,soda_trendvarid,'long_name','soda_trend');
    netcdf.putAtt(ncid,soda_trendvarid,'units','mm /year');

    soda_trend_filteredvarid=netcdf.defVar(ncid, 'soda_trend_filtered', 'NC_FLOAT', [lon_dimid lat_dimid]);
    netcdf.putAtt(ncid,soda_trend_filteredvarid,'long_name','soda_trend_filtered');
    netcdf.putAtt(ncid,soda_trend_filteredvarid,'units','mm /year');
    
    clim_soda_trend_dividedvarid=netcdf.defVar(ncid, 'clim_soda_trend_divided', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_soda_trend_dividedvarid,'long_name','clim_soda_trend_divided');
    netcdf.putAtt(ncid,clim_soda_trend_dividedvarid,'units','mm /year');

    mean_soda_trendvarid=netcdf.defVar(ncid, 'mean_soda_trend', 'NC_FLOAT', onedimid);
    netcdf.putAtt(ncid,mean_soda_trendvarid,'long_name','mean_soda_trend');
    netcdf.putAtt(ncid,mean_soda_trendvarid,'units','mm /year');

    mean_soda_trend_filteredvarid=netcdf.defVar(ncid, 'mean_soda_trend_filtered', 'NC_FLOAT', onedimid);
    netcdf.putAtt(ncid,mean_soda_trend_filteredvarid,'long_name','mean_soda_trend_filtered');
    netcdf.putAtt(ncid,mean_soda_trend_filteredvarid,'units','mm /year');
    
    mean_clim_soda_trend_dividedvarid=netcdf.defVar(ncid, 'mean_clim_soda_trend_divided', 'NC_FLOAT', clim_time_dimid);
    netcdf.putAtt(ncid,mean_clim_soda_trend_dividedvarid,'long_name','mean_clim_soda_trend_divided');
    netcdf.putAtt(ncid,mean_clim_soda_trend_dividedvarid,'units','mm /year');

    clim_sodavarid=netcdf.defVar(ncid, 'clim_soda_ssh', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_sodavarid,'long_name','clim_soda_ssh');
    netcdf.putAtt(ncid,clim_sodavarid,'units','mm ');

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
    netcdf.putVar(ncid, lonvarid, 0, len_lon, soda_lon(:));
    netcdf.putVar(ncid, latvarid, 0, len_lat, soda_lat(:));
    netcdf.putVar(ncid, soda_sshvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_soda_data);
    netcdf.putVar(ncid, soda_ssh_filteredvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_soda_data_filtered);
    netcdf.putVar(ncid, clim_sodavarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meansoda);
    netcdf.putVar(ncid, soda_trendvarid, [0 0], [len_lon len_lat], soda_trend);
    netcdf.putVar(ncid, soda_trend_filteredvarid, [0 0], [len_lon len_lat], soda_trend_filtered);
    netcdf.putVar(ncid, clim_soda_trend_dividedvarid, [0 0 0], [len_lon len_lat length(climtime)], clim_soda_trend_divided);
    netcdf.putVar(ncid, mean_soda_trendvarid, [0], [1], mean_soda_trend);
    netcdf.putVar(ncid, mean_soda_trend_filteredvarid, [0], [1], mean_soda_trend_filtered);
    netcdf.putVar(ncid, mean_clim_soda_trend_dividedvarid, [0], [length(climtime)], mean_clim_soda_trend_divided);

    netcdf.close(ncid);
end
