close all; clear all;  clc;   
% all_region ={'NWP','AKP2','ES', 'SS', 'YS', 'ECS'}
all_region ={'NWP'}
% all_region ={'BOH'}

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
    inputyear = [1980:2005]; % % put year which you want to plot [year year ...]
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
        case('NES') %% Northern East Sea
            refpolygon=nespolygon;
        case('SES') %% Southern East Sea
            refpolygon=sespolygon;
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
        case('BOH') %% Around Korea Peninsula
            refpolygon=bohpolygon;
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
        figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\5th_year\figure\recon\',regionname,'\'); % % where figure files will be saved
%         param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
%         filedir = strcat('E:\Data\Observation\CMEMS\'); % % where data files are
        filedir_out = 'E:\Data\Reanalysis\CSEOF_reconstructed_SSHA\';
        recondir='E:\Data\Reanalysis\CSEOF_reconstructed_SSHA\';
    elseif (strcmp(system_name,'GLNXA64'))
    end

%     run(param_script);
    ind=1;
    for yearij = 1:length(inputyear)
        for monthij = 1:length(inputmonth)
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])
            tic;
            tempyear = inputyear(yearij);
            tempmonth = inputmonth(monthij);
            % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc

            % read CMEMS DATA
            reconfilename = strcat(recondir,'CCAR_recon_sea_level_monthly_merged', '.nc');
            if (exist('recon_lon')==0)
                reconinfo=ncinfo(reconfilename);
                recon_lonname='lon';
                recon_latname='lat';
                recon_varname='ssha';
                reconinfo_lon=ncinfo(reconfilename,recon_lonname);
                reconinfo_lat=ncinfo(reconfilename,recon_latname);
                recon_lon = ncread(reconfilename,recon_lonname,1,reconinfo_lon.Dimensions.Length);
                recon_lat = ncread(reconfilename,recon_latname,1,reconinfo_lat.Dimensions.Length);
                [recon_lat2, recon_lon2]= meshgrid(recon_lat, recon_lon);
                [recon_lon_min, recon_lon_max, recon_lat_min, recon_lat_max] = findind_Y(1/20, lonlat(1:4), recon_lon2, recon_lat2);
                recon_lon2 = recon_lon2(recon_lon_min(1):recon_lon_max(1), recon_lat_min(1):recon_lat_max(1));
                recon_lat2 = recon_lat2(recon_lon_min(1):recon_lon_max(1), recon_lat_min(1):recon_lat_max(1));
                
                recon_lonsize_cut=recon_lon_max(1)-recon_lon_min(1)+1;
                recon_latsize_cut=recon_lat_max(1)-recon_lat_min(1)+1;
                comb_spatial_meanrms=(zeros([recon_lonsize_cut,recon_latsize_cut,12]));
                comb_spatial_meanbias=(zeros([recon_lonsize_cut,recon_latsize_cut,12]));
                comb_spatial_meanrecon=(zeros([recon_lonsize_cut,recon_latsize_cut,12]));
                comb_spatial_meanmodel=(zeros([recon_lonsize_cut,recon_latsize_cut,12]));
                
                switch(regionname)
                    case('NWP') %% North western Pacific
                        mask_recon(1:recon_lonsize_cut,1:recon_latsize_cut)=1;
                    otherwise
                        mask_recon = double(inpolygon(recon_lon2,recon_lat2,refpolygon(:,1),refpolygon(:,2)));
                        mask_recon(mask_recon==0)=NaN;
                end
            end

            recon_data_ind=(tempyear-1950)*12+tempmonth-1;
            recon_data = ncread(reconfilename,'ssha',[recon_lon_min(1) recon_lat_min(1) recon_data_ind],  ...
                [recon_lonsize_cut recon_latsize_cut 1]) * 0.01;  %% get data (cm) and change (cm) to (m)
            recon_data=recon_data.*mask_recon;
            recon_data(recon_data<-1000)=NaN;
            recon_data(recon_data>1000)=NaN;
            recon_data=recon_data.*mask_recon;
            
            comb_recon_data(:,:,ind) = recon_data;

            comb_spatial_meanrecon(:,:,monthij)=comb_spatial_meanrecon(:,:,monthij)+recon_data/double(length(inputyear));

            ind = ind + 1;
            toc;
        end
    end

    for loni=1:size(recon_data,1)
        for lati=1:size(recon_data,2)
            recon_sla_var(loni,lati)=var(comb_recon_data(loni,lati,:));
        end
    end
                
    comb_recon_clim_divided=reshape(comb_recon_data,[recon_lonsize_cut, recon_latsize_cut, 12, length(inputyear)]);
    
    trendtime=inputyear(1):1/12:inputyear(end)+1-1/12;
    recon_trend(1:recon_lonsize_cut,1:recon_latsize_cut)=NaN;
    for i=1:recon_lonsize_cut
        for j=1:recon_latsize_cut
            p=polyfit(trendtime,squeeze(comb_recon_data(i,j,:))',1);
            recon_trend(i,j)=p(1) * 1000.0 ;
        end
    end

    for t=1:length(inputyear)
        comb_recon_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_recon_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meanrecon;
    end

    recon_trend_filtered(1:recon_lonsize_cut,1:recon_latsize_cut)=NaN;

    for i=1:recon_lonsize_cut
        for j=1:recon_latsize_cut
            p=polyfit(trendtime,squeeze(comb_recon_data_filtered(i,j,:))',1);
            recon_trend_filtered(i,j)=p(1) * 1000.0 ;
        end
    end
    
    clim_recon_trend_divided(1:recon_lonsize_cut,1:recon_latsize_cut,1:12)=NaN;
    for k=1:12
        clim_trendtime(:,k)=inputyear(1)+(k-1)/12 : 1 : inputyear(end)+(k-1)/12;
        for i=1:recon_lonsize_cut
            for j=1:recon_latsize_cut
                p=polyfit(clim_trendtime(:,k)',squeeze(comb_recon_clim_divided(i,j,k,:))',1);
                clim_recon_trend_divided(i,j,k)=p(1) * 1000.0 ;
            end
        end
    end

    mean_recon_trend=mean(mean(recon_trend,'omitnan'),'omitnan');
    mean_recon_trend_filtered=mean(mean(recon_trend_filtered,'omitnan'),'omitnan');
    mean_clim_recon_trend_divided=mean(mean(clim_recon_trend_divided,'omitnan'),'omitnan');

    save([filedir_out,regionname,'recon_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

    ncid = netcdf.create(strcat(filedir_out,regionname,'recon_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc'),'NETCDF4');

    onedimid = netcdf.defDim(ncid,'one', 1);
    lon_dimid = netcdf.defDim(ncid, 'lon', recon_lonsize_cut);
    lat_dimid = netcdf.defDim(ncid,'lat',recon_latsize_cut);
    time_dimid = netcdf.defDim(ncid, 'time', 0);
    clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ['recon(ssh) _ ', 'monthly ssh trend file']);
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', [' recon ssh trend (',num2str(inputyear(1)),'-',num2str(inputyear(end)),') ']);
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', [' recon(OIssh) satellite data ' ]);
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

    recon_slavarid=netcdf.defVar(ncid, 'recon_sla', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,recon_slavarid,'long_name','recon_sla');
    netcdf.putAtt(ncid,recon_slavarid,'units','m ');
    
    recon_adtvarid=netcdf.defVar(ncid, 'recon_adt', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,recon_adtvarid,'long_name','recon_adt');
    netcdf.putAtt(ncid,recon_adtvarid,'units','m ');
    netcdf.putAtt(ncid,recon_adtvarid,'source','mdt from AVISO + sla');
    
    recon_errvarid=netcdf.defVar(ncid, 'recon_err', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,recon_errvarid,'long_name','recon_err');
    netcdf.putAtt(ncid,recon_errvarid,'units','m ');
    netcdf.putAtt(ncid,recon_errvarid,'source','Formal mapping error');
    
    recon_ssh_filteredvarid=netcdf.defVar(ncid, 'recon_ssh_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
    netcdf.putAtt(ncid,recon_ssh_filteredvarid,'long_name','recon_ssh_filtered');
    netcdf.putAtt(ncid,recon_ssh_filteredvarid,'units','mm ');
    
    recon_sla_varvarid=netcdf.defVar(ncid, 'recon_sla_var', 'NC_FLOAT', [lon_dimid lat_dimid]);
    netcdf.putAtt(ncid,recon_sla_varvarid,'long_name','recon_sla_var');
    netcdf.putAtt(ncid,recon_sla_varvarid,'units','mm /year');
    
    recon_trendvarid=netcdf.defVar(ncid, 'recon_trend', 'NC_FLOAT', [lon_dimid lat_dimid]);
    netcdf.putAtt(ncid,recon_trendvarid,'long_name','recon_trend');
    netcdf.putAtt(ncid,recon_trendvarid,'units','mm /year');

    recon_trend_filteredvarid=netcdf.defVar(ncid, 'recon_trend_filtered', 'NC_FLOAT', [lon_dimid lat_dimid]);
    netcdf.putAtt(ncid,recon_trend_filteredvarid,'long_name','recon_trend_filtered');
    netcdf.putAtt(ncid,recon_trend_filteredvarid,'units','mm /year');
    
    clim_recon_trend_dividedvarid=netcdf.defVar(ncid, 'clim_recon_trend_divided', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_recon_trend_dividedvarid,'long_name','clim_recon_trend_divided');
    netcdf.putAtt(ncid,clim_recon_trend_dividedvarid,'units','mm /year');

    mean_recon_trendvarid=netcdf.defVar(ncid, 'mean_recon_trend', 'NC_FLOAT', onedimid);
    netcdf.putAtt(ncid,mean_recon_trendvarid,'long_name','mean_recon_trend');
    netcdf.putAtt(ncid,mean_recon_trendvarid,'units','mm /year');

    mean_recon_trend_filteredvarid=netcdf.defVar(ncid, 'mean_recon_trend_filtered', 'NC_FLOAT', onedimid);
    netcdf.putAtt(ncid,mean_recon_trend_filteredvarid,'long_name','mean_recon_trend_filtered');
    netcdf.putAtt(ncid,mean_recon_trend_filteredvarid,'units','mm /year');
    
    mean_clim_recon_trend_dividedvarid=netcdf.defVar(ncid, 'mean_clim_recon_trend_divided', 'NC_FLOAT', clim_time_dimid);
    netcdf.putAtt(ncid,mean_clim_recon_trend_dividedvarid,'long_name','mean_clim_recon_trend_divided');
    netcdf.putAtt(ncid,mean_clim_recon_trend_dividedvarid,'units','mm /year');

    clim_reconvarid=netcdf.defVar(ncid, 'clim_recon_ssh', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_reconvarid,'long_name','clim_recon_ssh');
    netcdf.putAtt(ncid,clim_reconvarid,'units','mm ');

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
    netcdf.putVar(ncid, lonvarid, 0, recon_lonsize_cut, recon_lon(:));
    netcdf.putVar(ncid, latvarid, 0, recon_latsize_cut, recon_lat(:));
    netcdf.putVar(ncid, recon_slavarid, [0 0 0], [recon_lonsize_cut recon_latsize_cut length(ftime)], comb_recon_data);
    netcdf.putVar(ncid, recon_adtvarid, [0 0 0], [recon_lonsize_cut recon_latsize_cut length(ftime)], comb_recon_adt);
    netcdf.putVar(ncid, recon_errvarid, [0 0 0], [recon_lonsize_cut recon_latsize_cut length(ftime)], comb_recon_err);    
    netcdf.putVar(ncid, recon_ssh_filteredvarid, [0 0 0], [recon_lonsize_cut recon_latsize_cut length(ftime)], comb_recon_data_filtered);
    netcdf.putVar(ncid, clim_reconvarid, [0 0 0], [recon_lonsize_cut recon_latsize_cut length(climtime)], comb_spatial_meanrecon);
    netcdf.putVar(ncid, recon_sla_varvarid, [0 0], [recon_lonsize_cut recon_latsize_cut], recon_sla_var);
    netcdf.putVar(ncid, recon_trendvarid, [0 0], [recon_lonsize_cut recon_latsize_cut], recon_trend);
    netcdf.putVar(ncid, recon_trend_filteredvarid, [0 0], [recon_lonsize_cut recon_latsize_cut], recon_trend_filtered);
    netcdf.putVar(ncid, clim_recon_trend_dividedvarid, [0 0 0], [recon_lonsize_cut recon_latsize_cut length(climtime)], clim_recon_trend_divided);
    netcdf.putVar(ncid, mean_recon_trendvarid, [0], [1], mean_recon_trend);
    netcdf.putVar(ncid, mean_recon_trend_filteredvarid, [0], [1], mean_recon_trend_filtered);
    netcdf.putVar(ncid, mean_clim_recon_trend_dividedvarid, [0], [length(climtime)], mean_clim_recon_trend_divided);

    netcdf.close(ncid);
end

