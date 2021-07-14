close all; clear all;  clc;   
% all_region ={'NWP','ES', 'SS', 'YS', 'ECS'}
% all_region ={'YS', 'SS'}

% all_var ={'temp', 'salt', 'u', 'v', 'w', 'zeta', 'swrad', 'shflux', 'Uwind', 'Vwind'};
all_var ={'temp','salt','zeta'}

for varind=1:length(all_var)
    clearvars '*' -except varind all_var
    
    % % % 
    % % % Read Model variable at surface, bottom, 
    % % % boundary(north, south, east, west)

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
    rms_shadlev = [0 4];
    bias_shadlev = [-4 4];
    conlev  = 0:5:35;
    dl=1/20;
    % for snu_desktop
    testname='test03'   % % need to change
    inputyear = [1980:2015]; % % put year which you want to plot [year year ...]
    inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
    varname =all_var{varind};  %% reference variable -> temperature
    run('nwp_polygon_point.m');
%     regionname=all_region{regionind}
    regionname='NWP'
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
        otherwise
            ('?')
    end
    lonlat(1)=min(refpolygon(:,1));
    lonlat(2)=max(refpolygon(:,1));
    lonlat(3)=min(refpolygon(:,2));
    lonlat(4)=max(refpolygon(:,2));
    
    switch(varname)
        case('zeta') %% North western Pacific
            ndim=3;
            data_units = 'm';
        case('temp') %% North western Pacific
            ndim=4;
            data_units = 'Celsius'
        otherwise
            ndim=4;
            data_units = '?'
    end
    
    
    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd_year\figure\nwp_1_10\',testname,'\',regionname,'\'); % % where figure files will be saved
        param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
        filedir = strcat('E:\Data\Model\ROMS\nwp_1_10\', testname, '\run\'); % % where data files are
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
                lat = ncread(filename,'lat_rho',[1 1],[1,modelinfo.Dimensions(4).Length]);
                s_rho = ncread(filename,'s_rho');

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
                n_z = modelinfo.Dimensions(2).Length;
                
                if (ndim==3)
                    comb_spatial_mean_surf=(zeros([size(lon,1), size(lon,2), 12]));
                    comb_spatial_mean_bot=(zeros([size(lon,1), size(lon,2), 12]));
                    comb_spatial_mean_north=(zeros([size(lon,1), 12]));
                    comb_spatial_mean_south=(zeros([size(lon,1), 12]));
                    comb_spatial_mean_east=(zeros([size(lat,2), 12]));
                    comb_spatial_mean_west=(zeros([size(lat,2), 12]));
                elseif (ndim==4)
                    comb_spatial_mean_surf=(zeros([size(lon,1), size(lon,2), 12]));
                    comb_spatial_mean_bot=(zeros([size(lon,1), size(lon,2), 12]));
                    comb_spatial_mean_north=(zeros([size(lon,1), n_z, 12]));
                    comb_spatial_mean_south=(zeros([size(lon,1), n_z, 12]));
                    comb_spatial_mean_east=(zeros([size(lat,2), n_z, 12]));
                    comb_spatial_mean_west=(zeros([size(lat,2), n_z, 12]));
                end
            end

            data_info = ncinfo(filename, varname);  %% [lon lat depth time] -> [1601 1201 33 1]
            
            
            if (ndim==3)
                surf_data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]); %% NWP : [1:980, 1:920, 1]
                surf_data=surf_data.*mask_model;
                bot_data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]); %% NWP : [1:980, 1:920, 1]
                bot_data=bot_data.*mask_model;
                north_data = squeeze(ncread(filename,varname,[lon_min(1) lat_max(1) 1], [lon_max(1)-lon_min(1)+1 1 1])); %% NWP : [1:980, 920, 1]
                south_data = squeeze(ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 1 1])); %% NWP : [1:980, 1, 1]
                east_data = squeeze(ncread(filename,varname,[lon_max(1) lat_min(1) 1],  [1 lat_max(1)-lat_min(1)+1 1])); %% NWP : [980, 1:920, 1]
                west_data = squeeze(ncread(filename,varname,[lon_min(1) lat_min(1) 1],  [1 lat_max(1)-lat_min(1)+1 1])); %% NWP : [1, 1:920, 1]
            elseif (ndim==4)
                surf_data = ncread(filename,varname,[lon_min(1) lat_min(1) 40 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]); %% NWP : [1:980, 1:920, 40, 1]
                surf_data=surf_data.*mask_model;
                bot_data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]); %% NWP : [1:980, 1:920, 1, 1]
                bot_data=bot_data.*mask_model;
                north_data = squeeze(ncread(filename,varname,[lon_min(1) lat_max(1) 1 1], [lon_max(1)-lon_min(1)+1 1 n_z 1])); %% NWP : [1:980, 920, 1:40, 1]
                south_data = squeeze(ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 1 n_z 1])); %% NWP : [1:980, 1, 1:40, 1]
                east_data = squeeze(ncread(filename,varname,[lon_max(1) lat_min(1) 1 1],  [1 lat_max(1)-lat_min(1)+1 n_z 1])); %% NWP : [980, 1:920, 1:40, 1]
                west_data = squeeze(ncread(filename,varname,[lon_min(1) lat_min(1) 1 1],  [1 lat_max(1)-lat_min(1)+1 n_z 1])); %% NWP : [1, 1:920, 1:40, 1]
            end
            
            len_lon_model = size(surf_data,1);
            len_lat_model = size(surf_data,2);

%             meanbias = mean(mean(bias,'omitnan'),'omitnan');
%             meanrms = mean(mean(rms,'omitnan'),'omitnan');
            
            comb_surf_data(:,:,ind)=surf_data;
            comb_bot_data(:,:,ind)=bot_data;
            
            comb_spatial_mean_surf(:,:,monthij)=comb_spatial_mean_surf(:,:,monthij)+surf_data/double(length(inputyear));  
            comb_spatial_mean_bot(:,:,monthij)=comb_spatial_mean_bot(:,:,monthij)+bot_data/double(length(inputyear)); 
            if (ndim==3)
                comb_north_data(:,ind)=north_data;
                comb_south_data(:,ind)=south_data;
                comb_east_data(:,ind)=east_data;
                comb_west_data(:,ind)=west_data;
                comb_spatial_mean_north(:,monthij)=comb_spatial_mean_north(:,monthij)+north_data/double(length(inputyear));  
                comb_spatial_mean_south(:,monthij)=comb_spatial_mean_south(:,monthij)+south_data/double(length(inputyear));  
                comb_spatial_mean_east(:,monthij)=comb_spatial_mean_east(:,monthij)+east_data'/double(length(inputyear));  
                comb_spatial_mean_west(:,monthij)=comb_spatial_mean_west(:,monthij)+west_data'/double(length(inputyear));
            elseif (ndim==4)
                comb_north_data(:,:,ind)=north_data;
                comb_south_data(:,:,ind)=south_data;
                comb_east_data(:,:,ind)=east_data;
                comb_west_data(:,:,ind)=west_data;
                comb_spatial_mean_north(:,:,monthij)=comb_spatial_mean_north(:,:,monthij)+north_data/double(length(inputyear));  
                comb_spatial_mean_south(:,:,monthij)=comb_spatial_mean_south(:,:,monthij)+south_data/double(length(inputyear));  
                comb_spatial_mean_east(:,:,monthij)=comb_spatial_mean_east(:,:,monthij)+east_data/double(length(inputyear));  
                comb_spatial_mean_west(:,:,monthij)=comb_spatial_mean_west(:,:,monthij)+west_data/double(length(inputyear));
            end
              
            hold off
            close all;
            ind = ind + 1;
            toc;
        end
    end

%     trendtime=inputyear(1):1/12:inputyear(end)+1-1/12;
%     trend(1:len_lon,1:len_lat)=NaN;
%     for i=1:len_lon
%         for j=1:len_lat
%             p=polyfit(trendtime,squeeze(comb_interped_data(i,j,:))',1);
%             trend(i,j)=p(1);
%         end
%     end
% 
%     avhrr_trend(1:len_lon,1:len_lat)=NaN;
%     for i=1:len_lon
%         for j=1:len_lat
%             p=polyfit(trendtime,squeeze(comb_avhrr_data(i,j,:))',1);
%             avhrr_trend(i,j)=p(1);
%         end
%     end
% 
% 
%     for t=1:length(inputyear)
%         comb_interped_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_interped_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meanmodel;
%         comb_avhrr_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_avhrr_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meanavhrr;
%     end
% 
% 
%     trend_filtered(1:len_lon,1:len_lat)=NaN;
%     avhrr_trend_filtered(1:len_lon,1:len_lat)=NaN;
%     for i=1:len_lon
%         for j=1:len_lat
%             p=polyfit(trendtime,squeeze(comb_interped_data_filtered(i,j,:))',1);
%             trend_filtered(i,j)=p(1);
%         end
%     end
% 
%     for i=1:len_lon
%         for j=1:len_lat
%             p=polyfit(trendtime,squeeze(comb_avhrr_data_filtered(i,j,:))',1);
%             avhrr_trend_filtered(i,j)=p(1);
%         end
%     end
% 
%     mean_trend=mean(mean(trend,'omitnan'),'omitnan');
%     mean_trend_filtered=mean(mean(trend_filtered,'omitnan'),'omitnan');
%     mean_avhrr_trend=mean(mean(avhrr_trend,'omitnan'),'omitnan');
%     mean_avhrr_trend_filtered=mean(mean(avhrr_trend_filtered,'omitnan'),'omitnan');
% 
% 

    figdir=[figrawdir,'CLIM\'];
    outfile = strcat(figdir,regionname);
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 

%     rmsplot=plot(mean(comb_meanrms,1),'k')
%     jpgname=strcat(outfile, '_', testname,'_',regionname,'_climrms', '.jpg'); %% ~_year_month.jpg
%     xlabel('month')
%     ylabel('RMS(^oC)')
%     title(['RMS, climatology(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
%     ylim([0 4])
%     set(rmsplot,'LineWidth',2);
%     set(gca,'FontSize',15);
%     grid on
%     saveas(gcf,jpgname,'jpg');
%     grid off
% 
%     biasplot=plot(mean(comb_meanbias,1) ,'k')
%     jpgname=strcat(outfile, '_', testname,'_',regionname, '_climbias', '.jpg'); %% ~_year_month.jpg
%     xlabel('month')
%     ylabel('bias(^o)')
%     title(['BIAS, climatology(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
%     ylim([-4 4])
%     set(biasplot,'LineWidth',2);
%     grid on
%     saveas(gcf,jpgname,'jpg');
%     grid off

    save([filedir,regionname,'_',varname,'_all_time_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

    ncid = netcdf.create(strcat(filedir, testname, '_',regionname,'_',varname,'_all_time_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc'),'NETCDF4');

    onedimid = netcdf.defDim(ncid,'one', 1);
    xidimid = netcdf.defDim(ncid, 'xi_rho', len_lon_model);
    etadimid = netcdf.defDim(ncid,'eta_rho',len_lat_model);
    s_rhodimid = netcdf.defDim(ncid, 's_rho',n_z);
    time_dimid = netcdf.defDim(ncid, 'time', 0);
    clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ['NWP 1/20 _ ', testname, 'all time monthly ',varname, ' file']);
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', [' all time monthly ', varname, ' (' , num2str(inputyear(1)),'-',num2str(inputyear(end)),')']);
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
    
    xivarid=netcdf.defVar(ncid, 'xi_rho', 'NC_DOUBLE', [xidimid]);
    netcdf.putAtt(ncid,xivarid,'long_name','xi_rho');
    
    etavarid=netcdf.defVar(ncid, 'eta_rho', 'NC_DOUBLE', [etadimid]);
    netcdf.putAtt(ncid,etavarid,'long_name','eta_rho');
    
    s_rhovarid=netcdf.defVar(ncid, 's_rho', 'NC_DOUBLE', [s_rhodimid]);
    netcdf.putAtt(ncid,s_rhovarid,'long_name','s_rho');
    
    lon_rhovarid=netcdf.defVar(ncid, 'lon_rho', 'NC_DOUBLE', [xidimid etadimid]);
    netcdf.putAtt(ncid,lon_rhovarid,'long_name','lon_model');
    netcdf.putAtt(ncid,lon_rhovarid,'units','degree_east');

    lat_rhovarid=netcdf.defVar(ncid, 'lat_rho', 'NC_DOUBLE', [xidimid etadimid]);
    netcdf.putAtt(ncid,lat_rhovarid,'long_name','lat_model');
    netcdf.putAtt(ncid,lat_rhovarid,'units','degree_north');

    surf_datavarid=netcdf.defVar(ncid, 'surf_data', 'NC_FLOAT', [xidimid etadimid time_dimid]);
    netcdf.putAtt(ncid,surf_datavarid,'long_name','surf_data');
    netcdf.putAtt(ncid,surf_datavarid,'units', data_units);
    
    bot_datavarid=netcdf.defVar(ncid, 'bot_data', 'NC_FLOAT', [xidimid etadimid time_dimid]);
    netcdf.putAtt(ncid,bot_datavarid,'long_name','bot_data');
    netcdf.putAtt(ncid,bot_datavarid,'units', data_units);
    
    if (ndim==3)
        north_datavarid=netcdf.defVar(ncid, 'north_data', 'NC_FLOAT', [xidimid time_dimid]);
        netcdf.putAtt(ncid,north_datavarid,'long_name','north_data');
        netcdf.putAtt(ncid,north_datavarid,'units', data_units);

        south_datavarid=netcdf.defVar(ncid, 'south_data', 'NC_FLOAT', [xidimid time_dimid]);
        netcdf.putAtt(ncid,south_datavarid,'long_name','south_data');
        netcdf.putAtt(ncid,south_datavarid,'units', data_units);

        east_datavarid=netcdf.defVar(ncid, 'east_data', 'NC_FLOAT', [etadimid time_dimid]);
        netcdf.putAtt(ncid,east_datavarid,'long_name','east_data');
        netcdf.putAtt(ncid,east_datavarid,'units', data_units);

        west_datavarid=netcdf.defVar(ncid, 'west_data', 'NC_FLOAT', [etadimid time_dimid]);
        netcdf.putAtt(ncid,west_datavarid,'long_name','west_data');
        netcdf.putAtt(ncid,west_datavarid,'units', data_units);
        
        clim_northvarid=netcdf.defVar(ncid, 'clim_north', 'NC_FLOAT', [xidimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_northvarid,'long_name','clim_north');
        netcdf.putAtt(ncid,clim_northvarid,'units', data_units);

        clim_southvarid=netcdf.defVar(ncid, 'clim_south', 'NC_FLOAT', [xidimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_southvarid,'long_name','clim_south');
        netcdf.putAtt(ncid,clim_southvarid,'units', data_units);

        clim_eastvarid=netcdf.defVar(ncid, 'clim_east', 'NC_FLOAT', [etadimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_eastvarid,'long_name','clim_east');
        netcdf.putAtt(ncid,clim_eastvarid,'units', data_units);

        clim_westvarid=netcdf.defVar(ncid, 'clim_west', 'NC_FLOAT', [etadimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_westvarid,'long_name','clim_west');
        netcdf.putAtt(ncid,clim_westvarid,'units', data_units);
        
    elseif (ndim==4)
        north_datavarid=netcdf.defVar(ncid, 'north_data', 'NC_FLOAT', [xidimid s_rhodimid time_dimid]);
        netcdf.putAtt(ncid,north_datavarid,'long_name','north_data');
        netcdf.putAtt(ncid,north_datavarid,'units',data_units);

        south_datavarid=netcdf.defVar(ncid, 'south_data', 'NC_FLOAT', [xidimid s_rhodimid time_dimid]);
        netcdf.putAtt(ncid,south_datavarid,'long_name','south_data');
        netcdf.putAtt(ncid,south_datavarid,'units',data_units);

        east_datavarid=netcdf.defVar(ncid, 'east_data', 'NC_FLOAT', [etadimid s_rhodimid time_dimid]);
        netcdf.putAtt(ncid,east_datavarid,'long_name','east_data');
        netcdf.putAtt(ncid,east_datavarid,'units',data_units);

        west_datavarid=netcdf.defVar(ncid, 'west_data', 'NC_FLOAT', [etadimid s_rhodimid time_dimid]);
        netcdf.putAtt(ncid,west_datavarid,'long_name','west_data');
        netcdf.putAtt(ncid,west_datavarid,'units',data_units);

        clim_northvarid=netcdf.defVar(ncid, 'clim_north', 'NC_FLOAT', [xidimid s_rhodimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_northvarid,'long_name','clim_north');
        netcdf.putAtt(ncid,clim_northvarid,'units',data_units);

        clim_southvarid=netcdf.defVar(ncid, 'clim_south', 'NC_FLOAT', [xidimid s_rhodimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_southvarid,'long_name','clim_south');
        netcdf.putAtt(ncid,clim_southvarid,'units',data_units);

        clim_eastvarid=netcdf.defVar(ncid, 'clim_east', 'NC_FLOAT', [etadimid s_rhodimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_eastvarid,'long_name','clim_east');
        netcdf.putAtt(ncid,clim_eastvarid,'units',data_units);

        clim_westvarid=netcdf.defVar(ncid, 'clim_west', 'NC_FLOAT', [etadimid s_rhodimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_westvarid,'long_name','clim_west');
        netcdf.putAtt(ncid,clim_westvarid,'units',data_units);
    end

    clim_surfvarid=netcdf.defVar(ncid, 'clim_surf', 'NC_FLOAT', [xidimid etadimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_surfvarid,'long_name','clim_surf');
    netcdf.putAtt(ncid,clim_surfvarid,'units',data_units);
    
    clim_botvarid=netcdf.defVar(ncid, 'clim_bot', 'NC_FLOAT', [xidimid etadimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_botvarid,'long_name','clim_bot');
    netcdf.putAtt(ncid,clim_botvarid,'units',data_units);
    
    netcdf.endDef(ncid);
    if (ndim==3)
    elseif (ndim==4)
    end
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
    netcdf.putVar(ncid, xivarid, 0, len_lon_model, 1:len_lon_model);
    netcdf.putVar(ncid, etavarid, 0, len_lat_model, 1:len_lat_model);

    netcdf.putVar(ncid, s_rhovarid, 0, n_z, s_rho);
    netcdf.putVar(ncid, lon_rhovarid, [0 0], [len_lon_model len_lat_model], lon);
    netcdf.putVar(ncid, lat_rhovarid, [0 0], [len_lon_model len_lat_model], lat);
    netcdf.putVar(ncid, surf_datavarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_surf_data);
    netcdf.putVar(ncid, bot_datavarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_bot_data);
    
    if (ndim==3)
        netcdf.putVar(ncid, north_datavarid, [0 0], [len_lon_model length(ftime)], comb_north_data);
        netcdf.putVar(ncid, south_datavarid, [0 0], [len_lon_model length(ftime)], comb_south_data);
        netcdf.putVar(ncid, east_datavarid, [0 0], [len_lat_model length(ftime)], comb_east_data);
        netcdf.putVar(ncid, west_datavarid, [0 0], [len_lat_model length(ftime)], comb_west_data);
        netcdf.putVar(ncid, clim_northvarid, [0 0], [len_lon_model   length(climtime)], comb_spatial_mean_north);
        netcdf.putVar(ncid, clim_southvarid, [0 0], [len_lon_model   length(climtime)], comb_spatial_mean_south);
        netcdf.putVar(ncid, clim_eastvarid, [0 0], [len_lat_model   length(climtime)], comb_spatial_mean_east);
        netcdf.putVar(ncid, clim_westvarid, [0 0], [len_lat_model   length(climtime)], comb_spatial_mean_west);
    elseif (ndim==4)
        netcdf.putVar(ncid, north_datavarid, [0 0 0], [len_lon_model n_z length(ftime)], comb_north_data);
        netcdf.putVar(ncid, south_datavarid, [0 0 0], [len_lon_model n_z length(ftime)], comb_south_data);
        netcdf.putVar(ncid, east_datavarid, [0 0 0], [len_lat_model n_z length(ftime)], comb_east_data);
        netcdf.putVar(ncid, west_datavarid, [0 0 0], [len_lat_model n_z length(ftime)], comb_west_data);
        netcdf.putVar(ncid, clim_northvarid, [0 0 0], [len_lon_model n_z length(climtime)], comb_spatial_mean_north);
        netcdf.putVar(ncid, clim_southvarid, [0 0 0], [len_lon_model n_z length(climtime)], comb_spatial_mean_south);
        netcdf.putVar(ncid, clim_eastvarid, [0 0 0], [len_lat_model n_z length(climtime)], comb_spatial_mean_east);
        netcdf.putVar(ncid, clim_westvarid, [0 0 0], [len_lat_model n_z length(climtime)], comb_spatial_mean_west);
    end

    netcdf.putVar(ncid, clim_surfvarid, [0 0 0], [len_lon_model len_lat_model length(climtime)], comb_spatial_mean_surf);
    netcdf.putVar(ncid, clim_botvarid, [0 0 0], [len_lon_model len_lat_model length(climtime)], comb_spatial_mean_bot);
    
    netcdf.close(ncid);
end