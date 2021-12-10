close all; clear all;  clc;   
% %  get SST RMS, SST BIAS data and plot(climatological value time series)
% %  Updated 29-Nov-2021 by Yong-Yub Kim 


RCM_info.name={'test2117', 'test2118', 'test2119', 'test2120', 'test2121'};
% RCM_info.name={'test2127', 'test2128', 'test2129', 'test2130', 'test2131'};

RCM_info.abbs = {'RCM-CNE', 'RCM-ECV', 'RCM-ACC', 'RCM-CNH', 'RCM-CMC'};

% RCM_info.region ={'TEST'};

RCM_info.model = 'nwp_1_20';
% RCM_info.dataroot = '/data1/RCM/CMIP6/';
RCM_info.dataroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'backup_surf', filesep];
% RCM_info.saveroot = '/data1/RCM/CMIP6/';
RCM_info.saveroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'backup_surf', filesep];
RCM_info.phase = 'run';
% RCM_info.region = {'NWP', 'AKP4'};
RCM_info.region = {'AKP4'};
RCM_info.years = 1985:2014;
% RCM_info.years = 2015:2050;
RCM_info.months = 1:12;
RCM_grid.dl = 1/20;

% % % configuration of GCM
GCM_info.name={'CNRM-ESM2-1', 'EC-Earth3-Veg', 'ACCESS-CM2', 'CNRM-CM6-1-HR', 'CMCC-ESM2'};
GCM_info.abbs = {'GCM-CNE', 'GCM-ECV', 'GCM-ACC', 'GCM-CNH', 'GCM-CMC'};
% GCM_info.name={'CNRM-ESM2-1'};
% GCM_info.abbs = {  'GCM-CNRM'};
GCM_info.model = GCM_info.name;
GCM_info.dataroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'CMIP6', filesep, 'NWP', filesep];
GCM_info.saveroot = GCM_info.dataroot;
GCM_info.phase = RCM_info.phase;
GCM_info.region = RCM_info.region;
GCM_info.years = RCM_info.years;
GCM_info.months = RCM_info.months;
GCM_grid.dl = 1/2;

% % % configuration of OISST
OISST_info.filedir = ['D:', filesep, 'Data', filesep, 'Observation', ...
    filesep, 'OISST', filesep, 'monthly_kimyy', filesep];
OISST_grid.dl = 1/4;



for testnameind=1:length(RCM_info.name)
    for regionind=1:length(RCM_info.region)
        clearvars '*' -except regionind testnameind flags flg_lev param ...
            RCM_info RCM_grid GCM_info GCM_grid OISST_info OISST_grid
        
        tmp.variable = 'temp';
        tmp.variable_GCM = 'thetao';
        tmp.variable_OISST = 'temp';
        tmp.fs = filesep; % file separator win = '\', linux = '/'
        
        
        % %     set dropbox path
        if (strcmp(computer,'PCWIN64'))
            tmp.dropboxpath = 'C:\Users\User\Dropbox';
        else
            tmp.dropboxpath = '/home/kimyy/Dropbox';
        end
        addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));


        
        % %     set temporary variables (testname, regionname, filesep, ...)
        RCM_info.testname = RCM_info.name{testnameind};
        RCM_info.regionname = RCM_info.region{regionind};
        RCM_info.abb = RCM_info.abbs{testnameind};
        [RCM_info.scenario, tmp.error_status] = Func_0013_RCM_CMIP6_scenname(RCM_info.testname);
        
        GCM_info.testname = GCM_info.name{testnameind};
        GCM_info.regionname = RCM_info.regionname;
        GCM_info.abb = GCM_info.abbs{testnameind};
        GCM_info.scenario = RCM_info.scenario;

        [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
        addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Model' ...
            tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', tmp.fs, 'nwp_1_20', tmp.fs ...
            'run', tmp.fs, 'SSH', tmp.fs, '2phase_1st', tmp.fs, 'subroutine']));
        
        
        [RCM_grid.refpolygon, RCM_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(RCM_info.regionname);

        tmp.param_script = [tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', ...
            tmp.fs, 'Model', tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', ...
            tmp.fs, 'nwp_1_20', tmp.fs, 'run', tmp.fs, 'fig_param', tmp.fs, ...
            'fig_param2_kyy_', RCM_info.regionname, '.m'];
        RCM_info.filedir = [RCM_info.dataroot, RCM_info.testname, tmp.fs, 'run', tmp.fs, ...
            tmp.variable, tmp.fs];
        RCM_info.savedir = [RCM_info.dataroot, RCM_info.testname, tmp.fs, 'run', tmp.fs, ...
            tmp.variable, tmp.fs];
        GCM_info.filedir = [GCM_info.dataroot, tmp.variable_GCM, tmp.fs, GCM_info.scenario, tmp.fs, ...
            'Omon', tmp.fs, GCM_info.testname, tmp.fs];
        
% %         
% % 
% %         if (strcmp(system_name,'PCWIN64'))
% %             % % for windows
% %             figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
% %             param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
% %             filedir = strcat('H:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
% %             avhrrdir='E:\Data\Observation\OISST\monthly\';
% %         elseif (strcmp(system_name,'GLNXA64'))
% %         end

% % %         flag configuration (process)
        for folding=1:1
            flags.fig_name{1}='get model (RCM, GCM) SST and OISST data';

            for flagi=1:length(flags.fig_name)
                flags.fig_switch(flagi)=0;
            end
            flags.fig_switch(1)=2;
            flags.fig_switch(2)=0;
            flags.fig_switch(3)=0;
            flags.fig_switch(4)=0;
            flags.fig_switch(5)=0;
        end
        
        if flags.fig_switch(1) > 0
            flags.fig_tmp = flags.fig_switch(1);
            SSH_2p_1st_011_sub_001_get_SST;
        end

        if (exist(RCM_info.matname , 'file') == 2 || flags.fig_switch(1)~=2)   
            load(RCM_info.matname_interped); 
            load(GCM_info.matname_interped);
            load(OISST_info.matname);
            
            [tmp.m_value, tmp.error_status] = ...
                Func_0011_get_area_weighted_mean(RCM_data_interped.rms, OISST_grid.lon2, OISST_grid.lat2);
            disp(['RCM SST RMS: ', num2str(tmp.m_value), ])
            [tmp.m_value, tmp.error_status] = ...
                Func_0011_get_area_weighted_mean(GCM_data_interped.rms, OISST_grid.lon2, OISST_grid.lat2);
            disp(['GCM SST RMS: ', num2str(tmp.m_value), ])
            [tmp.m_value, tmp.error_status] = ...
                Func_0011_get_area_weighted_mean(RCM_data_interped.yearly_rms, OISST_grid.lon2, OISST_grid.lat2);
            disp(['RCM SST yearly RMS: ', num2str(tmp.m_value), ])
            [tmp.m_value, tmp.error_status] = ...
                Func_0011_get_area_weighted_mean(GCM_data_interped.yearly_rms, OISST_grid.lon2, OISST_grid.lat2);
            disp(['GCM SST yearly RMS: ', num2str(tmp.m_value), ])
        end

% %         run(param_script);
% %         ind=1;
% %         for yearij = 1:length(inputyear)
% %             for monthij = 1:length(inputmonth)
% %                 disp([num2str(yearij), 'y_',num2str(monthij),'m'])
% %                 tic;
% %                 tempyear = inputyear(yearij);
% %                 tempmonth = inputmonth(monthij);
% %                 % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
% %                 filename = strcat(filedir, num2str(tempyear,'%04i'), '\', ...
% %                         testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');
% %                 % read model data
% %                 if (exist('lon')==0)
% %                     modelinfo=ncinfo(filename);
% %                     lon = ncread(filename,'lon_rho',[1 1],[modelinfo.Dimensions(5).Length,1]);
% %                     lat = ncread(filename,'lat_rho',[1 1],[1,modelinfo.Dimensions(6).Length]);
% % 
% %                     lon_west = abs(lon - (lonlat(1)-1));
% %                     min_lon_west=min(lon_west);
% %                     lon_east = abs(lon - (lonlat(2)+1));
% %                     min_lon_east=min(lon_east);
% %                     lat_south = abs(lat - (lonlat(3)-1));
% %                     min_lat_south=min(lat_south);
% %                     lat_north = abs(lat - (lonlat(4)+1));
% %                     min_lat_north=min(lat_north);
% % 
% %                     lon_min = find(lon_west == min_lon_west);
% %                     lon_max = find(lon_east == min_lon_east);
% %                     lat_min = find(lat_south == min_lat_south);
% %                     lat_max = find(lat_north == min_lat_north);
% % 
% %                     lon = ncread(filename,'lon_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
% %                     lat = ncread(filename,'lat_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
% % 
% %         %             polygon_ind=NaN(size(refpolygon,1),2);
% %         %             for i=1:size(refpolygon,1)
% %         %                 [polygon_ind(i,1), trash_ind, polygon_ind(i,2), trash_ind]=findind_Y(dl, [refpolygon(i,1),refpolygon(i,2)],lon',lat');
% %         %             end
% %         %             mask_model = inpolygon(lon,lat,polygon_ind(:,1),polygon_ind(:,2));
% %                     switch(regionname)
% %                         case('NWP') %% North western Pacific
% %                             mask_model(1:size(lon,1),1:size(lon,2))=1;
% %                         otherwise
% %                             mask_model = double(inpolygon(lon,lat,refpolygon(:,1),refpolygon(:,2)));
% %                             mask_model(mask_model==0)=NaN;
% %                     end
% %                 end
% % 
% %                 data_info = ncinfo(filename, varname);  %% [lon lat depth time] -> [1601 1201 33 1]
% % 
% %                 data = ncread(filename,varname,[lon_min(1) lat_min(1) 40 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);
% %                 data=data.*mask_model;
% % 
% %                 % read OISST DATA
% %                 avhrrfilename = strcat(avhrrdir,'avhrr_monthly', num2str(tempyear,'%04i'), ...
% %                         '_',num2str(tempmonth,'%02i'), '.nc');
% %                 if (exist('avhrr_lon')==0)
% %                     avhrrinfo=ncinfo(avhrrfilename);
% %                     avhrr_lon = ncread(avhrrfilename,'long',[1 1],[avhrrinfo.Dimensions(2).Length,1]);
% %                     avhrr_lat = ncread(avhrrfilename,'lat',[1 1],[1,avhrrinfo.Dimensions(1).Length]);
% % 
% %                     avhrr_lon_west = abs(avhrr_lon - (lonlat(1)));
% %                     min_avhrr_lon_west=min(avhrr_lon_west);
% %                     avhrr_lon_east = abs(avhrr_lon - (lonlat(2)));
% %                     min_avhrr_lon_east=min(avhrr_lon_east);
% %                     avhrr_lat_south = abs(avhrr_lat - (lonlat(3)));
% %                     min_avhrr_lat_south=min(avhrr_lat_south);
% %                     avhrr_lat_north = abs(avhrr_lat - (lonlat(4)));
% %                     min_avhrr_lat_north=min(avhrr_lat_north);
% % 
% %                     avhrr_lon_min = find(avhrr_lon_west == min_avhrr_lon_west);
% %                     avhrr_lon_max = find(avhrr_lon_east == min_avhrr_lon_east);
% %                     avhrr_lat_min = find(avhrr_lat_south == min_avhrr_lat_south);
% %                     avhrr_lat_max = find(avhrr_lat_north == min_avhrr_lat_north);
% % 
% %             %         ncinfo('E:\Data\Observation\OISST\monthly\avhrr_monthly1983_11.nc');
% % 
% %                     avhrr_lon = ncread(avhrrfilename,'long', [avhrr_lon_min(1) avhrr_lat_min(1)], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1)]);
% %                     avhrr_lat = ncread(avhrrfilename,'lat', [avhrr_lon_min(1) avhrr_lat_min(1)], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1)]);
% % 
% %                     comb_spatial_meanrms=(zeros([size(avhrr_lon),12]));
% %                     comb_spatial_meanbias=(zeros([size(avhrr_lon),12]));
% %                     comb_spatial_meanavhrr=(zeros([size(avhrr_lon),12]));
% %                     comb_spatial_meanmodel=(zeros([size(avhrr_lon),12]));
% % 
% %                     switch(regionname)
% %                         case('NWP') %% North western Pacific
% %                             mask_avhrr(1:size(avhrr_lon,1),1:size(avhrr_lon,2))=1;
% %                         otherwise
% %                             mask_avhrr = double(inpolygon(avhrr_lon,avhrr_lat,refpolygon(:,1),refpolygon(:,2)));
% %                             mask_avhrr(mask_avhrr==0)=NaN;
% %                     end
% %                 end
% %                 len_lon = length(avhrr_lon(:,1));
% %                 len_lat = length(avhrr_lat(1,:));
% %                 len_lon_model = size(data,1);
% %                 len_lat_model = size(data,2);
% % 
% % 
% %                 avhrr_data = ncread(avhrrfilename,varname,[avhrr_lon_min(1) avhrr_lat_min(1)], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1)]);
% %                 avhrr_data=avhrr_data.*mask_avhrr;
% % 
% %                 interped_data = griddata(double(lon), double(lat), data,double(avhrr_lon),double(avhrr_lat));   
% %                 bias = interped_data-avhrr_data;  
% %                 rms = sqrt((interped_data-avhrr_data).^2);  
% % %                 meanbias = mean(mean(bias,'omitnan'),'omitnan');
% % %                 meanrms = mean(mean(rms,'omitnan'),'omitnan');
% %                 meanbias = mean(bias(:),'omitnan');
% %                 meanrms = mean(rms(:),'omitnan');
% % 
% %                 comb_interped_data(:,:,ind) = interped_data;
% %                 comb_avhrr_data(:,:,ind) = avhrr_data;
% %                 comb_rms_data(:,:,ind) = rms;
% %                 comb_bias_data(:,:,ind) = bias;
% % 
% %                 comb_meanrms(yearij,monthij)=meanrms;
% %                 comb_meanbias(yearij,monthij)=meanbias;
% %                 comb_spatial_meanrms(:,:,monthij)=comb_spatial_meanrms(:,:,monthij)+rms/double(length(inputyear));
% %                 comb_spatial_meanbias(:,:,monthij)=comb_spatial_meanbias(:,:,monthij)+bias/double(length(inputyear));
% %                 comb_spatial_meanavhrr(:,:,monthij)=comb_spatial_meanavhrr(:,:,monthij)+avhrr_data/double(length(inputyear));
% %                 comb_spatial_meanmodel(:,:,monthij)=comb_spatial_meanmodel(:,:,monthij)+interped_data/double(length(inputyear));
% % 
% %                 comb_data(:,:,ind) = data;
% %                 comb_interped_data(:,:,ind) = interped_data;
% %                 comb_avhrr_data(:,:,ind) = avhrr_data;
% % 
% %                 hold off
% %                 close all;
% %                 ind = ind + 1;
% %                 toc;
% %             end
% %         end
% % 
% %         trendtime=inputyear(1):1/12:inputyear(end)+1-1/12;
% %         trend(1:len_lon,1:len_lat)=NaN;
% %         for i=1:len_lon
% %             for j=1:len_lat
% %                 p=polyfit(trendtime,squeeze(comb_interped_data(i,j,:))',1);
% %                 trend(i,j)=p(1);
% %             end
% %         end
% % 
% %         avhrr_trend(1:len_lon,1:len_lat)=NaN;
% %         for i=1:len_lon
% %             for j=1:len_lat
% %                 p=polyfit(trendtime,squeeze(comb_avhrr_data(i,j,:))',1);
% %                 avhrr_trend(i,j)=p(1);
% %             end
% %         end
% % 
% % 
% %         for t=1:length(inputyear)
% %             comb_interped_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_interped_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meanmodel;
% %             comb_avhrr_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_avhrr_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meanavhrr;
% %         end
% % 
% %         trend_filtered(1:len_lon,1:len_lat)=NaN;
% %         avhrr_trend_filtered(1:len_lon,1:len_lat)=NaN;
% %         for i=1:len_lon
% %             for j=1:len_lat
% %                 p=polyfit(trendtime,squeeze(comb_interped_data_filtered(i,j,:))',1);
% %                 trend_filtered(i,j)=p(1);
% %             end
% %         end
% % 
% %         for i=1:len_lon
% %             for j=1:len_lat
% %                 p=polyfit(trendtime,squeeze(comb_avhrr_data_filtered(i,j,:))',1);
% %                 avhrr_trend_filtered(i,j)=p(1);
% %             end
% %         end
% % 
% %         mean_trend=mean(mean(trend,'omitnan'),'omitnan');
% %         mean_trend_filtered=mean(mean(trend_filtered,'omitnan'),'omitnan');
% %         mean_avhrr_trend=mean(mean(avhrr_trend,'omitnan'),'omitnan');
% %         mean_avhrr_trend_filtered=mean(mean(avhrr_trend_filtered,'omitnan'),'omitnan');
% % 
% %   
% % 
% % 
% %         figdir=[figrawdir,'CLIM\'];
% %         outfile = strcat(figdir,regionname);
% %         if (exist(strcat(figdir) , 'dir') ~= 7)
% %             mkdir(strcat(figdir));
% %         end 
% % 
% %         rmsplot=plot(mean(comb_meanrms,1),'k')
% %         jpgname=strcat(outfile, '_', testname,'_',regionname,'_climrms_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
% %         xlabel('month')
% %         ylabel('RMS(^oC)')
% %         title(['RMS, climatology(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
% %         ylim([0 4])
% %         set(rmsplot,'LineWidth',2);
% %         set(gca,'FontSize',15);
% %         grid on
% %         saveas(gcf,jpgname,'jpg');
% %         grid off
% % 
% %         biasplot=plot(mean(comb_meanbias,1) ,'k')
% %         jpgname=strcat(outfile, '_', testname,'_',regionname, '_climbias_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
% %         xlabel('month')
% %         ylabel('bias(^o)')
% %         title(['BIAS, climatology(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
% %         ylim([-4 4])
% %         set(biasplot,'LineWidth',2);
% %         grid on
% %         saveas(gcf,jpgname,'jpg');
% %         grid off
% % 
% %         save([filedir,regionname,'sst_rms_and_bias_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
% % 
% %         ncid = netcdf.create(strcat(filedir, testname, '_',regionname,'_rms_clim_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc'),'NETCDF4');
% % 
% % 
% %         onedimid = netcdf.defDim(ncid,'one', 1);
% %         lon_dimid = netcdf.defDim(ncid, 'lon', len_lon);
% %         lat_dimid = netcdf.defDim(ncid,'lat',len_lat);
% %         xidimid = netcdf.defDim(ncid, 'xi_rho', len_lon_model);
% %         etadimid = netcdf.defDim(ncid,'eta_rho',len_lat_model);
% %         time_dimid = netcdf.defDim(ncid, 'time', 0);
% %         clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);
% % 
% %         netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
% %             'type', ['NWP 1/20 _ ', testname, 'monthly SST RMS/BIAS file']);
% %         netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
% %             'title', ' monthly SST RMS/BIAS (1982-2009) ');
% %         netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
% %             'source', [' ROMS NWP 1/20 data from _ ',testname ]);
% %         netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
% %             'author', 'Created by Y.Y.Kim');
% %         netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
% %             'date', date);
% % 
% %         timevarid=netcdf.defVar(ncid, 'time', 'NC_DOUBLE', time_dimid);
% %         netcdf.putAtt(ncid,timevarid,'long_name','time');
% %         netcdf.putAtt(ncid,timevarid,'units','days since 1900-12-31 00:00:00');
% %         netcdf.putAtt(ncid,timevarid,'calendar','gregorian');
% % 
% %         clim_timevarid=netcdf.defVar(ncid, 'clim_time', 'NC_DOUBLE', clim_time_dimid);
% %         netcdf.putAtt(ncid,clim_timevarid,'long_name','clim_time');
% %         netcdf.putAtt(ncid,clim_timevarid,'units','days since 1900-12-31 00:00:00');
% %         netcdf.putAtt(ncid,clim_timevarid,'calendar','gregorian');
% % 
% %         lonvarid=netcdf.defVar(ncid, 'lon', 'NC_DOUBLE', lon_dimid);
% %         netcdf.putAtt(ncid,lonvarid,'long_name','longitude');
% %         netcdf.putAtt(ncid,lonvarid,'units','degree_east');
% % 
% %         latvarid=netcdf.defVar(ncid, 'lat', 'NC_DOUBLE', lat_dimid);
% %         netcdf.putAtt(ncid,latvarid,'long_name','latitude');
% %         netcdf.putAtt(ncid,latvarid,'units','degree_north');
% % 
% %         lon_rhovarid=netcdf.defVar(ncid, 'lon_rho', 'NC_DOUBLE', [xidimid etadimid]);
% %         netcdf.putAtt(ncid,lon_rhovarid,'long_name','lon_model');
% %         netcdf.putAtt(ncid,lon_rhovarid,'units','degree_east');
% % 
% %         lat_rhovarid=netcdf.defVar(ncid, 'lat_rho', 'NC_DOUBLE', [xidimid etadimid]);
% %         netcdf.putAtt(ncid,lat_rhovarid,'long_name','lat_model');
% %         netcdf.putAtt(ncid,lat_rhovarid,'units','degree_north');
% % 
% %         raw_sstvarid=netcdf.defVar(ncid, 'raw_sst', 'NC_FLOAT', [xidimid etadimid time_dimid]);
% %         netcdf.putAtt(ncid,raw_sstvarid,'long_name','raw_sst');
% %         netcdf.putAtt(ncid,raw_sstvarid,'units','Celsius');
% % 
% %         sstvarid=netcdf.defVar(ncid, 'sst', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
% %         netcdf.putAtt(ncid,sstvarid,'long_name','sst');
% %         netcdf.putAtt(ncid,sstvarid,'units','Celsius');
% % 
% %         sst_filteredvarid=netcdf.defVar(ncid, 'sst_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
% %         netcdf.putAtt(ncid,sst_filteredvarid,'long_name','sst_filtered');
% %         netcdf.putAtt(ncid,sst_filteredvarid,'units','Celsius');
% % 
% %         avhrr_sstvarid=netcdf.defVar(ncid, 'avhrr_sst', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
% %         netcdf.putAtt(ncid,avhrr_sstvarid,'long_name','avhrr_sst');
% %         netcdf.putAtt(ncid,avhrr_sstvarid,'units','Celsius');
% % 
% %         avhrr_sst_filteredvarid=netcdf.defVar(ncid, 'avhrr_sst_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
% %         netcdf.putAtt(ncid,avhrr_sst_filteredvarid,'long_name','avhrr_sst_filtered');
% %         netcdf.putAtt(ncid,avhrr_sst_filteredvarid,'units','Celsius');
% % 
% %         trendvarid=netcdf.defVar(ncid, 'trend', 'NC_FLOAT', [lon_dimid lat_dimid]);
% %         netcdf.putAtt(ncid,trendvarid,'long_name','trend');
% %         netcdf.putAtt(ncid,trendvarid,'units','Celsius/year');
% % 
% %         trend_filteredvarid=netcdf.defVar(ncid, 'trend_filtered', 'NC_FLOAT', [lon_dimid lat_dimid]);
% %         netcdf.putAtt(ncid,trend_filteredvarid,'long_name','trend_filtered');
% %         netcdf.putAtt(ncid,trend_filteredvarid,'units','Celsius/year');
% % 
% %         avhrr_trendvarid=netcdf.defVar(ncid, 'avhrr_trend', 'NC_FLOAT', [lon_dimid lat_dimid]);
% %         netcdf.putAtt(ncid,avhrr_trendvarid,'long_name','avhrr_trend');
% %         netcdf.putAtt(ncid,avhrr_trendvarid,'units','Celsius/year');
% % 
% %         avhrr_trend_filteredvarid=netcdf.defVar(ncid, 'avhrr_trend_filtered', 'NC_FLOAT', [lon_dimid lat_dimid]);
% %         netcdf.putAtt(ncid,avhrr_trend_filteredvarid,'long_name','avhrr_trend_filtered');
% %         netcdf.putAtt(ncid,avhrr_trend_filteredvarid,'units','Celsius/year');
% % 
% %         rmsvarid=netcdf.defVar(ncid, 'rms', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
% %         netcdf.putAtt(ncid,rmsvarid,'long_name','rms');
% %         netcdf.putAtt(ncid,rmsvarid,'units','Celsius');
% % 
% %         biasvarid=netcdf.defVar(ncid, 'bias', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
% %         netcdf.putAtt(ncid,biasvarid,'long_name','bias');
% %         netcdf.putAtt(ncid,biasvarid,'units','Celsius');
% % 
% %         mean_trendvarid=netcdf.defVar(ncid, 'mean_trend', 'NC_FLOAT', onedimid);
% %         netcdf.putAtt(ncid,mean_trendvarid,'long_name','mean_trend');
% %         netcdf.putAtt(ncid,mean_trendvarid,'units','Celsius/year');
% % 
% %         mean_trend_filteredvarid=netcdf.defVar(ncid, 'mean_trend_filtered', 'NC_FLOAT', onedimid);
% %         netcdf.putAtt(ncid,mean_trend_filteredvarid,'long_name','mean_trend_filtered');
% %         netcdf.putAtt(ncid,mean_trend_filteredvarid,'units','Celsius/year');
% % 
% %         mean_avhrr_trendvarid=netcdf.defVar(ncid, 'mean_avhrr_trend', 'NC_FLOAT', onedimid);
% %         netcdf.putAtt(ncid,mean_avhrr_trendvarid,'long_name','mean_avhrr_trend');
% %         netcdf.putAtt(ncid,mean_avhrr_trendvarid,'units','Celsius/year');
% % 
% %         mean_avhrr_trend_filteredvarid=netcdf.defVar(ncid, 'mean_avhrr_trend_filtered', 'NC_FLOAT', onedimid);
% %         netcdf.putAtt(ncid,mean_avhrr_trend_filteredvarid,'long_name','mean_avhrr_trend_filtered');
% %         netcdf.putAtt(ncid,mean_avhrr_trend_filteredvarid,'units','Celsius/year');
% % 
% %         clim_sstvarid=netcdf.defVar(ncid, 'clim_sst', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
% %         netcdf.putAtt(ncid,clim_sstvarid,'long_name','clim_sst');
% %         netcdf.putAtt(ncid,clim_sstvarid,'units','Celsius');
% % 
% %         clim_avhrrvarid=netcdf.defVar(ncid, 'clim_avhrr_sst', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
% %         netcdf.putAtt(ncid,clim_avhrrvarid,'long_name','clim_avhrr_sst');
% %         netcdf.putAtt(ncid,clim_avhrrvarid,'units','Celsius');
% % 
% %         clim_rmsvarid=netcdf.defVar(ncid, 'clim_rms', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
% %         netcdf.putAtt(ncid,clim_rmsvarid,'long_name','clim_rms');
% %         netcdf.putAtt(ncid,clim_rmsvarid,'units','Celsius');
% % 
% %         clim_biasvarid=netcdf.defVar(ncid, 'clim_bias', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
% %         netcdf.putAtt(ncid,clim_biasvarid,'long_name','clim_bias');
% %         netcdf.putAtt(ncid,clim_biasvarid,'units','Celsius');
% % 
% %         netcdf.endDef(ncid);
% % 
% %         tind=1;
% %         for yearij = 1:length(inputyear)
% %             for month=1:12 
% %                 tempyear = inputyear(yearij);
% %                 ftime(tind) = datenum(tempyear,month,15) - datenum(1900,12,31);
% %                 tind=tind+1;
% %             end
% %         end
% %         for month=1:12 
% %                 tempyear = inputyear(yearij);
% %                 climtime(month) = datenum(1950,month,15) - datenum(1900,12,31);
% %         end
% % 
% %         netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
% %         netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
% %         netcdf.putVar(ncid, lonvarid, 0, len_lon, avhrr_lon(:,1));
% %         netcdf.putVar(ncid, latvarid, 0, len_lat, avhrr_lat(1,:));
% %         netcdf.putVar(ncid, lon_rhovarid, [0 0], [len_lon_model len_lat_model], lon);
% %         netcdf.putVar(ncid, lat_rhovarid, [0 0], [len_lon_model len_lat_model], lat);
% %         netcdf.putVar(ncid, raw_sstvarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_data);
% %         netcdf.putVar(ncid, sstvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_interped_data);
% %         netcdf.putVar(ncid, sst_filteredvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_interped_data_filtered);
% %         netcdf.putVar(ncid, avhrr_sstvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_avhrr_data);
% %         netcdf.putVar(ncid, avhrr_sst_filteredvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_avhrr_data_filtered);
% %         netcdf.putVar(ncid, rmsvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_rms_data);
% %         netcdf.putVar(ncid, biasvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_bias_data);
% %         netcdf.putVar(ncid, clim_sstvarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meanmodel);
% %         netcdf.putVar(ncid, clim_avhrrvarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meanavhrr);
% %         netcdf.putVar(ncid, clim_rmsvarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meanrms);
% %         netcdf.putVar(ncid, clim_biasvarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meanbias);
% %         netcdf.putVar(ncid, trendvarid, [0 0], [len_lon len_lat], trend);
% %         netcdf.putVar(ncid, trend_filteredvarid, [0 0], [len_lon len_lat], trend_filtered);
% %         netcdf.putVar(ncid, avhrr_trendvarid, [0 0], [len_lon len_lat], avhrr_trend);
% %         netcdf.putVar(ncid, avhrr_trend_filteredvarid, [0 0], [len_lon len_lat], avhrr_trend_filtered);
% %         netcdf.putVar(ncid, mean_trendvarid, [0], [1], mean_trend);
% %         netcdf.putVar(ncid, mean_trend_filteredvarid, [0], [1], mean_trend_filtered);
% %         netcdf.putVar(ncid, mean_avhrr_trendvarid, [0], [1], mean_avhrr_trend);
% %         netcdf.putVar(ncid, mean_avhrr_trend_filteredvarid, [0], [1], mean_avhrr_trend_filtered);
% % 
% %         netcdf.close(ncid);
    end
end

% SSH_4th_mid_report3