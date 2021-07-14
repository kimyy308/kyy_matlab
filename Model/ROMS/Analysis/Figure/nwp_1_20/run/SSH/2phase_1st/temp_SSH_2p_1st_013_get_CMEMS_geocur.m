close all; clear all;  clc;   
% all_region ={'NWP','ES', 'SS', 'YS', 'ECS'}
% all_region ={'NWP', 'AKP2'}
all_region ={'NWP'}

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
        dropboxpath='C:\Users\user\Dropbox';
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
    inputyear = [1993:2014]; % % put year which you want to plot [year year ...]
    inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

    
    varname ='sla';  %% reference variable -> temperature
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
        case('AKP2') %% Around Korean Peninsula (except kuroshio)
            refpolygon=akp2polygon;
        case('AKP3') %% Around Korea Peninsula
            refpolygon=akp3polygon;
        case('AKP4') %% Around Korea Peninsula
            refpolygon=akp4polygon;
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
        figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\6th_year\figure\CMEMS\',regionname,'\'); % % where figure files will be saved
        param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
        filedir = strcat('Z:\내 드라이브\Data\Observation\CMEMS\'); % % where data files are
        cmemsdir='Z:\내 드라이브\Data\Observation\CMEMS\';
    elseif (strcmp(system_name,'GLNXA64'))
        figrawdir =strcat('/data2/kimyy/Reanalysis/SODA/SODA_3_4_2/monthly/figures/',regionname,'/'); % % where figure files will be saved
        param_script ='/home/kimyy/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run/fig_param/fig_param_kyy_EKB_RMS.m'
        filedir = strcat('/data2/kimyy/Reanalysis/SODA/SODA_3_4_2/monthly/'); % % where data files are
        cmemsdir='/data2/kimyy/Reanalysis/SODA/SODA_3_4_2/monthly/';
    end
    variable = 'SSH';
    run(param_script);
    ind=1;
    for yearij = 1:length(inputyear)
%         inputmonth = 1:yeardays(inputyear(yearij));
        for monthij = 1:length(inputmonth)
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])
            tic;
            tempyear = inputyear(yearij);
            tempmonth = inputmonth(monthij);
            % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc

            % read CMEMS DATA
            cmemsfilename = strcat(cmemsdir,'cmems_nwp_ssh_', num2str(tempyear,'%04i'),'.nc');
            if (exist('cmems_lon')==0)
                cmemsinfo=ncinfo(cmemsfilename);
                cmems_lon = ncread(cmemsfilename,'longitude',[1],[cmemsinfo.Dimensions(3).Length]);
                cmems_lat = ncread(cmemsfilename,'latitude',[1],[cmemsinfo.Dimensions(2).Length]);

                cmems_lon_west = abs(cmems_lon - (lonlat(1)));
                min_cmems_lon_west=min(cmems_lon_west);
                cmems_lon_east = abs(cmems_lon - (lonlat(2)));
                min_cmems_lon_east=min(cmems_lon_east);
                cmems_lat_south = abs(cmems_lat - (lonlat(3)));
                min_cmems_lat_south=min(cmems_lat_south);
                cmems_lat_north = abs(cmems_lat - (lonlat(4)));
                min_cmems_lat_north=min(cmems_lat_north);

                cmems_lon_min = find(cmems_lon_west == min_cmems_lon_west);
                cmems_lon_max = find(cmems_lon_east == min_cmems_lon_east);
                cmems_lat_min = find(cmems_lat_south == min_cmems_lat_south);
                cmems_lat_max = find(cmems_lat_north == min_cmems_lat_north);

        %         ncinfo('E:\Data\Observation\OIssh\monthly\cmems_monthly1983_11.nc');

                cmems_lon = ncread(cmemsfilename,'longitude', [cmems_lon_min(1)], [cmems_lon_max(1)-cmems_lon_min(1)]);
                cmems_lat = ncread(cmemsfilename,'latitude', [cmems_lat_min(1)], [cmems_lat_max(1)-cmems_lat_min(1)]);

                comb_spatial_mean_u=(zeros([length(cmems_lon),length(cmems_lat),12]));
                comb_spatial_mean_v=(zeros([length(cmems_lon),length(cmems_lat),12]));
                
                [cmems_lat2 cmems_lon2]=meshgrid(cmems_lat, cmems_lon);
                switch(regionname)
                    case('NWP') %% North western Pacific
                        mask_cmems(1:size(cmems_lon,1),1:size(cmems_lon,2))=1;
                    otherwise
                        mask_cmems = double(inpolygon(cmems_lon2,cmems_lat2,refpolygon(:,1),refpolygon(:,2)));
                        mask_cmems(mask_cmems==0)=NaN;
                end
            end
            len_lon = length(cmems_lon(:));
            len_lat = length(cmems_lat(:));
            if tempmonth==1
                cmems_daily_u=ncread(cmemsfilename,'ugos',[cmems_lon_min(1) cmems_lat_min(1) 1], [cmems_lon_max(1)-cmems_lon_min(1) cmems_lat_max(1)-cmems_lat_min(1) 31]);
                cmems_daily_v=ncread(cmemsfilename,'vgos',[cmems_lon_min(1) cmems_lat_min(1) 1], [cmems_lon_max(1)-cmems_lon_min(1) cmems_lat_max(1)-cmems_lat_min(1) 31]);
            else
                cmems_daily_u=ncread(cmemsfilename,'ugos',[cmems_lon_min(1) cmems_lat_min(1) sum(eomday(tempyear,1:tempmonth-1))+1], [cmems_lon_max(1)-cmems_lon_min(1) cmems_lat_max(1)-cmems_lat_min(1) eomday(tempyear,tempmonth)]);
                cmems_daily_v=ncread(cmemsfilename,'vgos',[cmems_lon_min(1) cmems_lat_min(1) sum(eomday(tempyear,1:tempmonth-1))+1], [cmems_lon_max(1)-cmems_lon_min(1) cmems_lat_max(1)-cmems_lat_min(1) eomday(tempyear,tempmonth)]);
            end
            cmems_u = nanmean(cmems_daily_u,3);
            cmems_v = nanmean(cmems_daily_v,3);
            cmems_u(cmems_u<-1000)=NaN;
            cmems_u(cmems_u>1000)=NaN;
            cmems_u=cmems_u.*mask_cmems;
            cmems_v(cmems_v<-1000)=NaN;
            cmems_v(cmems_v>1000)=NaN;
            cmems_v=cmems_v.*mask_cmems;
            
            comb_cmems_u(:,:,ind) = cmems_u;
            comb_cmems_v(:,:,ind) = cmems_v;

            ind = ind + 1;
            toc;
        end
    end
    mean_u=mean(comb_cmems_u,3,'omitnan');
    mean_v=mean(comb_cmems_v,3,'omitnan');

    figdir=[figrawdir,'CLIM\'];
    outfile = strcat(figdir,regionname);
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 

    save([filedir,regionname,'cmems_gos_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

    ncid = netcdf.create(strcat(filedir,regionname,'cmems_gos_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc'),'NETCDF4');

    onedimid = netcdf.defDim(ncid,'one', 1);
    lon_dimid = netcdf.defDim(ncid, 'lon', len_lon);
    lat_dimid = netcdf.defDim(ncid,'lat',len_lat);
    time_dimid = netcdf.defDim(ncid, 'time', 0);
    clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ['cmems(ssh) _ ', 'monthly ssh trend file']);
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', [' cmems ssh trend (',num2str(inputyear(1)),'-',num2str(inputyear(end)),') ']);
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', [' cmems(OIssh) satellite data ' ]);
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
       
    cmems_uvarid=netcdf.defVar(ncid, 'cmems_u', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,cmems_uvarid,'long_name','cmems_u');
    netcdf.putAtt(ncid,cmems_uvarid,'units','m ');
    netcdf.putAtt(ncid,cmems_uvarid,'source','monthly Absolute geostrophic velocity : zonal component');
    
    cmems_vvarid=netcdf.defVar(ncid, 'cmems_v', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,cmems_vvarid,'long_name','cmems_v');
    netcdf.putAtt(ncid,cmems_vvarid,'units','m ');
    netcdf.putAtt(ncid,cmems_vvarid,'source','monthly Absolute geostrophic velocity : meridional component');
    
    mean_cmems_uvarid=netcdf.defVar(ncid, 'mean_cmems_u', 'NC_FLOAT', [lon_dimid lat_dimid]);
    netcdf.putAtt(ncid,mean_cmems_uvarid,'long_name','mean_cmems_u');
    netcdf.putAtt(ncid,mean_cmems_uvarid,'units','m ');
    netcdf.putAtt(ncid,mean_cmems_uvarid,'source','monthly Absolute geostrophic velocity : zonal component');
    
    mean_cmems_vvarid=netcdf.defVar(ncid, 'mean_cmems_v', 'NC_FLOAT', [lon_dimid lat_dimid]);
    netcdf.putAtt(ncid,mean_cmems_vvarid,'long_name','mean_cmems_v');
    netcdf.putAtt(ncid,mean_cmems_vvarid,'units','m ');
    netcdf.putAtt(ncid,mean_cmems_vvarid,'source','monthly Absolute geostrophic velocity : meridional component');
    
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
    netcdf.putVar(ncid, lonvarid, 0, len_lon, cmems_lon(:));
    netcdf.putVar(ncid, latvarid, 0, len_lat, cmems_lat(:));
    netcdf.putVar(ncid, mean_cmems_uvarid, [0 0], [len_lon len_lat], mean_u);
    netcdf.putVar(ncid, mean_cmems_vvarid, [0 0], [len_lon len_lat], mean_v);
    netcdf.putVar(ncid, cmems_uvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_cmems_u);
    netcdf.putVar(ncid, cmems_vvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_cmems_v);

    netcdf.close(ncid);
end