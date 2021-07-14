close all; clear all;  clc;
% %  get reconstructed SSH. compare soda and reSSH. save. 
% all_region ={'AKP4'}
% all_region ={'NWP', 'AKP2', 'ES', 'YS', 'NES', 'SES'}
% all_region ={'YSECS', 'ECS2', 'ES', 'YS', 'NES', 'SES'}

% all_testname = {'test57', 'test58', 'test59'};
all_testname = {'MyOcean'};
all_region ={'AKP4'};
for testnameind=1:length(all_testname)
    for regionind=1:length(all_region)
        clearvars '*' -except regionind testnameind all_region all_testname
        % % % 
        % % % Read Model SSH
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

        shadlev = [-2 2];
        % rms_shadlev = [0 4];
        trend_shadlev = [0 4];
        % bias_shadlev = [-4 4];
        conlev  = 0:5:35;
        dl=1/10;

        % for snu_desktop
        testname=all_testname{testnameind} 
        inputyear = [1994:2017]; % % put year which you want to plot [year year ...]
        inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

        varname ='zeta'
        regionname=all_region{regionind}
        run('nwp_polygon_point.m');
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
            case('ECS2') %% East China Sea
                refpolygon=ecs2polygon;
            case('YSECS') %% YS & East China Sea
                refpolygon=ysecspolygon;
            case('AKP') %% Around Korea Peninsula
                refpolygon=akppolygon;
            case('AKP2') %% Around Korea Peninsula
                refpolygon=akp2polygon;
            case('AKP3') %% Around Korea Peninsula
                refpolygon=akp3polygon;
            case('AKP4') %% Around Korea Peninsula
                refpolygon=akp4polygon;
            case('CA') %% Coastal Area around korea peninsula
                refpolygon=capolygon;
            case('EKB') %% Coastal Area around korea peninsula
                refpolygon=ekbpolygon;
            case('BOH') %% Coastal Area around korea peninsula
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

        system_name=computer;
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\4th_year\figure\',testname,'\',regionname,'\'); % % where figure files will be saved
            param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
            filedir = strcat('E:\Data\Reanalysis\', testname, '\'); % % where data files are
            cmemsdir='E:\Data\Observation\CMEMS\';
        elseif (strcmp(system_name,'GLNXA64'))
        end

        run(param_script);
        ind=1;
        readsst=0;
        for yearij = 1:length(inputyear)
            for monthij = 1:length(inputmonth)
                disp([num2str(yearij), 'y_',num2str(monthij),'m'])
                tic;
                tempyear = inputyear(yearij);
                tempmonth = inputmonth(monthij);
                % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
                filename = strcat(filedir, 'NWPMyOcean_ssh_trend_', num2str(inputyear(1),'%04i'), ...
                            '_',num2str(inputyear(end),'%04i'), '.nc');
                if readsst==1
                    sstfilename = strcat(filedir, 'NWPsoda_temp_trend_', num2str(inputyear(1),'%04i'), ...
                                '_',num2str(inputyear(end),'%04i'), '.nc');
                end
                % read soda data
                if (exist('lon')==0)
                    sodainfo=ncinfo(filename);
%                     lon1= ncread(filename, 'lon');
%                     lat1= ncread(filename, 'lat');
%                     [lat, lon]=meshgrid(lat1, lon1);
                    lon= ncread(filename, 'lon');
                    lat= ncread(filename, 'lat');
                    
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

%                     lon = ncread(filename,'lon_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
%                     lat = ncread(filename,'lat_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
%                     
%                     lon=lon(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
%                     lat=lat(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    lon2 = ncread(filename,'lon', [lon_min(1)], [lon_max(1)-lon_min(1)+1]);
                    lat2 = ncread(filename,'lat', [lat_min(1)], [lat_max(1)-lat_min(1)+1]);
                    
                    [lat lon]=meshgrid(lat2,lon2);
                    switch(regionname)
                        case('NWP') %% North western Pacific
                            mask_soda(1:size(lon,1),1:size(lon,2))=1;
                        otherwise
                            mask_soda = double(inpolygon(lon,lat,refpolygon(:,1),refpolygon(:,2)));
                            mask_soda(mask_soda==0)=NaN;
                    end
                end


                data_info = ncinfo(filename, 'MyOcean_ssh');  %% [lon lat depth time] -> [1601 1201 33 1]

                data = ncread(filename,'MyOcean_ssh',[lon_min(1) lat_min(1) (yearij-1)*12+monthij], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                data=data.*mask_soda;
                if readsst==1
                    sst = ncread(sstfilename,'soda_temp',[lon_min(1) lat_min(1) (yearij-1)*12+monthij], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                    sst=sst.*mask_soda;
                end
                
                len_lon_soda = size(data,1);
                len_lat_soda = size(data,2);
                len_lon=len_lon_soda;
                len_lat=len_lat_soda;
                if (exist('comb_spatial_meansoda')==0)
                    comb_spatial_meansoda=(zeros([len_lon_soda,len_lat_soda,12]));
                    if readsst==1
                        comb_spatial_meansoda_sst=(zeros([len_lon_soda,len_lat_soda,12]));
                    end
                end
                
                comb_data(:,:,ind) = single(data);
                if readsst==1
                    comb_sst(:,:,ind) = single(sst);
                end
                comb_spatial_meansoda(:,:,monthij)=comb_spatial_meansoda(:,:,monthij)+data/double(length(inputyear));
                if readsst==1
                    comb_spatial_meansoda_sst(:,:,monthij)=comb_spatial_meansoda_sst(:,:,monthij)+sst/double(length(inputyear));
                end
                
                
% % % %                 % read CMEMS DATA
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

%                     comb_spatial_meanrms=(zeros([length(cmems_lon),length(cmems_lat),12]));
%                     comb_spatial_meanbias=(zeros([length(cmems_lon),length(cmems_lat),12]));
                    comb_spatial_meancmems=(zeros([length(cmems_lon),length(cmems_lat),12]));
                    comb_spatial_meaninterped=(zeros([length(cmems_lon),length(cmems_lat),12]));
                    if readsst==1
                        comb_spatial_meaninterped_sst=(zeros([length(cmems_lon),length(cmems_lat),12]));
                    end
%                     comb_spatial_meansoda=(zeros([length(cmems_lon),length(cmems_lat),12]));

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
                    cmems_daily_adt=ncread(cmemsfilename,'adt',[cmems_lon_min(1) cmems_lat_min(1) 1], [cmems_lon_max(1)-cmems_lon_min(1) cmems_lat_max(1)-cmems_lat_min(1) 31]);
                    cmems_daily_data=ncread(cmemsfilename,'sla',[cmems_lon_min(1) cmems_lat_min(1) 1], [cmems_lon_max(1)-cmems_lon_min(1) cmems_lat_max(1)-cmems_lat_min(1) 31]);
                else
                    cmems_daily_adt=ncread(cmemsfilename,'adt',[cmems_lon_min(1) cmems_lat_min(1) sum(eomday(tempyear,1:tempmonth-1))+1], [cmems_lon_max(1)-cmems_lon_min(1) cmems_lat_max(1)-cmems_lat_min(1) eomday(tempyear,tempmonth)]);
                    cmems_daily_data=ncread(cmemsfilename,'sla',[cmems_lon_min(1) cmems_lat_min(1) sum(eomday(tempyear,1:tempmonth-1))+1], [cmems_lon_max(1)-cmems_lon_min(1) cmems_lat_max(1)-cmems_lat_min(1) eomday(tempyear,tempmonth)]);
                end
                cmems_adt = nanmean(cmems_daily_adt,3);
                cmems_data = nanmean(cmems_daily_data,3);
    %             cmems_data = ncread(cmemsfilename,varname,[cmems_lon_min(1) cmems_lat_min(1) tempmonth], [cmems_lon_max(1)-cmems_lon_min(1) cmems_lat_max(1)-cmems_lat_min(1) 1]);
                cmems_data(cmems_data<-1000)=NaN;
                cmems_data(cmems_data>1000)=NaN;
                cmems_data=cmems_data.*mask_cmems;
                
                comb_cmems_adt(:,:,ind) = cmems_adt;
                comb_cmems_data(:,:,ind) = cmems_data;
                comb_spatial_meancmems(:,:,monthij)=comb_spatial_meancmems(:,:,monthij)+cmems_adt/double(length(inputyear));

                interped_data = griddata(double(lon), double(lat), double(data), double(cmems_lon),double(cmems_lat)')';
                comb_interped_data(:,:,ind) = interped_data;
                comb_spatial_meaninterped(:,:,monthij)=comb_spatial_meaninterped(:,:,monthij)+interped_data/double(length(inputyear));
                
                if readsst==1
                    interped_sst = griddata(double(lon), double(lat), double(sst), double(cmems_lon),double(cmems_lat)')';
                    comb_interped_sst(:,:,ind) = interped_sst;
                    comb_spatial_meaninterped_sst(:,:,monthij)=comb_spatial_meaninterped_sst(:,:,monthij)+interped_sst/double(length(inputyear));
                end
                ind = ind + 1;
                toc;
            end
        end
        
        %%% get sea level anomaly
        cmems_sla = comb_cmems_data;
        for sla_i=1:size(cmems_sla,1)
            for sla_j=1:size(cmems_sla,2)
                cmems_sla_mean(sla_i,sla_j)=mean(cmems_sla(sla_i,sla_j,:),'omitnan');
            end
        end
        interped_ssh=comb_interped_data;
        for sla_i=1:size(cmems_sla,1)
            for sla_j=1:size(cmems_sla,2)
                interped_sla_mean(sla_i,sla_j)=mean(interped_ssh(sla_i,sla_j,:),'omitnan');
                interped_sla(sla_i,sla_j,:)=interped_ssh(sla_i,sla_j,:)-(interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
            end
        end
        interped_sla_divided=reshape(interped_sla,[size(cmems_sla,1), size(cmems_sla,2), 12, length(inputyear)]);
        clim_interped_sla=mean(interped_sla_divided,4,'omitnan');
        for t=1:length(inputyear)
            interped_sla_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(interped_sla(:,:,(t-1)*12+1:(t-1)*12+12)-clim_interped_sla);
        end
        cmems_sla_divided=reshape(cmems_sla,[size(cmems_sla,1), size(cmems_sla,2), 12, length(inputyear)]);
        clim_cmems_sla=mean(cmems_sla_divided,4,'omitnan');
        for t=1:length(inputyear)
            cmems_sla_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(cmems_sla(:,:,(t-1)*12+1:(t-1)*12+12)-clim_cmems_sla);
        end
            
        
% %        2, 3, 5year lowpass filter 
        sample_freq=1;   %sampling frequency
        filt_order=5;  %filtering order
        nq_freq=sample_freq/2;  %Nyquist frequency
        ftype='low';  %filter type 
        
%         nyears=[1,2,3,4,5,6,7,8,9,10];
    nyears=[1,2,3,4,5];
%     nyears=[2];
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            cutoff_freq=1/(12*nyear);
            [filt_coef_b,filt_coef_a]=butter(filt_order,cutoff_freq/nq_freq, ftype); %butterworth filter
            for i=1:len_lon
                for j=1:len_lat
                    eval(['comb_interped_',num2str(nyear),'y_lowpass(i,j,:)=filter(filt_coef_b,filt_coef_a,comb_interped_data(i,j,:)-mean(comb_interped_data(i,j,:)))+mean(comb_interped_data(i,j,:));'])
                    eval(['comb_cmems_',num2str(nyear),'y_lowpass(i,j,:)=filter(filt_coef_b,filt_coef_a,comb_cmems_data(i,j,:)-mean(comb_cmems_data(i,j,:)))+mean(comb_cmems_data(i,j,:));'])
                    eval(['comb_interped_',num2str(nyear),'y_lowpass(i,j,1:6*nyear)=NaN;'])
                    eval(['comb_interped_',num2str(nyear),'y_lowpass(i,j,end-6*nyear+1:end)=NaN;'])
                    eval(['comb_cmems_',num2str(nyear),'y_lowpass(i,j,1:6*nyear)=NaN;'])
                    eval(['comb_cmems_',num2str(nyear),'y_lowpass(i,j,end-6*nyear+1:end)=NaN;'])
                end
            end
        end
    
% % %          sea level anomaly low pass filter
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            cutoff_freq=1/(12*nyear);
            [filt_coef_b,filt_coef_a]=butter(filt_order,cutoff_freq/nq_freq, ftype); %butterworth filter
            for i=1:len_lon
                for j=1:len_lat
                    eval(['comb_interped_sla_',num2str(nyear),'y_lowpass(i,j,:)=filter(filt_coef_b,filt_coef_a,interped_sla(i,j,:)-mean(interped_sla(i,j,:)))+mean(interped_sla(i,j,:));'])
                    eval(['comb_cmems_sla_',num2str(nyear),'y_lowpass(i,j,:)=filter(filt_coef_b,filt_coef_a,cmems_sla(i,j,:)-mean(cmems_sla(i,j,:)))+mean(cmems_sla(i,j,:));'])
                    eval(['comb_interped_sla_',num2str(nyear),'y_lowpass(i,j,1:6*nyear)=NaN;'])
                    eval(['comb_interped_sla_',num2str(nyear),'y_lowpass(i,j,end-6*nyear+1:end)=NaN;'])
                    eval(['comb_cmems_sla_',num2str(nyear),'y_lowpass(i,j,1:6*nyear)=NaN;'])
                    eval(['comb_cmems_sla_',num2str(nyear),'y_lowpass(i,j,end-6*nyear+1:end)=NaN;'])
                end
            end
        end
        
%         for l=1:288
%             tempmsl=interped_sla(:,:,l);
%             msl(l)=mean(tempmsl(:),'omitnan');
%             tempmsl_lp=comb_interped_sla_2y_lowpass(:,:,l);
%             msl_lp(l)=mean(tempmsl_lp(:),'omitnan');
%         end
%         plot(msl,'k');
%         hold on;
%         plot(msl_lp,'r');
%         hold off;
        
        
% % %         trend
        trendtime=inputyear(1):1/12:inputyear(end)+1-1/12;
        trend(1:len_lon_soda,1:len_lat_soda)=NaN;
        for i=1:len_lon_soda
            for j=1:len_lat_soda
                p=polyfit(trendtime,squeeze(comb_data(i,j,:))',1);
                trend(i,j)=p(1);
            end
        end
        trend = trend * 1000.0; %% m/y -> mm/y
       disp('trend complete') 

       
% % %         trend (seasonal filtered)
        for t=1:length(inputyear)
            comb_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(comb_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meansoda);
        end
        
        trend_filtered(1:len_lon_soda,1:len_lat_soda)=NaN;
        for i=1:len_lon_soda
            for j=1:len_lat_soda
                p=polyfit(trendtime,squeeze(comb_data_filtered(i,j,:))',1);
                trend_filtered(i,j)=p(1);
            end
        end
        trend_filtered = trend_filtered * 1000.0; %% m/y -> mm/y
       disp('trend_filtered complete') 
       
% % %        sst trend
        if readsst==1        
            trendtime=inputyear(1):1/12:inputyear(end)+1-1/12;
            sst_trend(1:len_lon_soda,1:len_lat_soda)=NaN;
            for i=1:len_lon_soda
                for j=1:len_lat_soda
                    p=polyfit(trendtime,squeeze(comb_sst(i,j,:))',1);
                    sst_trend(i,j)=p(1);
                end
            end
            sst_trend = sst_trend; %% m/y -> mm/y
            disp('sst trend complete') 
        end
       
% % %        sst trend (seasonal filtered)
        if readsst==1
            for t=1:length(inputyear)
                comb_sst_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(comb_sst(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meansoda_sst);
            end

            sst_trend_filtered(1:len_lon_soda,1:len_lat_soda)=NaN;
            for i=1:len_lon_soda
                for j=1:len_lat_soda
                    p=polyfit(trendtime,squeeze(comb_sst_filtered(i,j,:))',1);
                    sst_trend_filtered(i,j)=p(1);
                end
            end
            
            sst_trend_filtered = sst_trend_filtered; %% m/y -> mm/y
            disp('sst trend_filtered complete') 
        end
% % %         mean trend (raw, seasonal filtered)
        mean_trend=mean(trend(:),'omitnan');
        mean_trend_filtered=mean(trend_filtered(:),'omitnan');
        
% % %         climatological trend 
        comb_spatial_data=reshape(comb_data, [len_lon_soda len_lat_soda 12 length(inputyear)]);
        climtrendtime=inputyear(1):inputyear(end);
        for i=1:len_lon_soda
            for j=1:len_lat_soda
                for k=1:12  % month
                    p=polyfit(climtrendtime,squeeze(comb_spatial_data(i,j,k,:))',1);
                    trend_clim(i,j,k)=p(1);
                end
            end
        end
       disp('climatological trend complete') 

        
% % %         cmems trend 
        comb_cmems_clim_divided=reshape(comb_cmems_data,[len_lon, len_lat, 12, length(inputyear)]);
    
        cmems_trend(1:len_lon,1:len_lat)=NaN;
        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(trendtime,squeeze(comb_cmems_data(i,j,:))',1);
                cmems_trend(i,j)=p(1) * 1000.0 ;
            end
        end
       disp('cmems trend complete') 

% % %         cmems trend (seasonal filtered)
        for t=1:length(inputyear)
            comb_cmems_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_cmems_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meancmems;
        end

        cmems_trend_filtered(1:len_lon,1:len_lat)=NaN;

        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(trendtime,squeeze(comb_cmems_data_filtered(i,j,:))',1);
                cmems_trend_filtered(i,j)=p(1) * 1000.0 ;
            end
        end
       disp('cmems trend_filtered complete') 

        
        
% % %         climatological cmems trend 

        clim_cmems_trend_divided(1:len_lon,1:len_lat,1:12)=NaN;
        for k=1:12
            clim_trendtime(:,k)=inputyear(1)+(k-1)/12 : 1 : inputyear(end)+(k-1)/12;
            for i=1:len_lon
                for j=1:len_lat
                    p=polyfit(clim_trendtime(:,k)',squeeze(comb_cmems_clim_divided(i,j,k,:))',1);
                    clim_cmems_trend_divided(i,j,k)=p(1) * 1000.0 ;
                end
            end
        end
       disp('cmems climatological trend complete') 

        mean_cmems_trend=mean(mean(cmems_trend,'omitnan'),'omitnan');
        mean_cmems_trend_filtered=mean(mean(cmems_trend_filtered,'omitnan'),'omitnan');
        mean_clim_cmems_trend_divided=mean(mean(clim_cmems_trend_divided,'omitnan'),'omitnan');
        
% % %         interped trend 
        comb_interped_clim_divided=reshape(comb_interped_data,[len_lon, len_lat, 12, length(inputyear)]);
    
        interped_trend(1:len_lon,1:len_lat)=NaN;
        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(trendtime,squeeze(comb_interped_data(i,j,:))',1);
                interped_trend(i,j)=p(1) * 1000.0 ;
            end
        end
        disp('interped trend complete') 


% % %         interped trend (seasonal filtered)
        for t=1:length(inputyear)
            comb_interped_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_interped_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meaninterped;
        end

        interped_trend_filtered(1:len_lon,1:len_lat)=NaN;

        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(trendtime,squeeze(comb_interped_data_filtered(i,j,:))',1);
                interped_trend_filtered(i,j)=p(1) * 1000.0 ;
            end
        end
        disp('interped trend_filtered complete') 

% % %         interped sst trend
        if readsst==1
            comb_interped_sst_clim_divided=reshape(comb_interped_sst,[len_lon, len_lat, 12, length(inputyear)]);

            interped_sst_trend(1:len_lon,1:len_lat)=NaN;
            for i=1:len_lon
                for j=1:len_lat
                    p=polyfit(trendtime,squeeze(comb_interped_sst(i,j,:))',1);
                    interped_sst_trend(i,j)=p(1) * 1000.0 ;
                end
            end
            disp('interped sst trend complete') 
        end

% % %         interped sst trend (seasonal filtered)
        if readsst==1
            for t=1:length(inputyear)
                comb_interped_sst_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_interped_sst(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meaninterped_sst;
            end

            interped_sst_trend_filtered(1:len_lon,1:len_lat)=NaN;

            for i=1:len_lon
                for j=1:len_lat
                    p=polyfit(trendtime,squeeze(comb_interped_sst_filtered(i,j,:))',1);
                    interped_sst_trend_filtered(i,j)=p(1) * 1000.0 ;
                end
            end
            disp('interped sst trend_filtered complete') 
        end
        
% % %   correlation coefficient between soda and cmems
        for i=1:len_lon
            for j=1:len_lat
                temp_corr=corrcoef(squeeze(comb_interped_data(i,j,:))',squeeze(comb_cmems_data(i,j,:))');
                corr_interped(i,j)=temp_corr(1,2);
            end
        end
        disp('corr coef complete') 
        
% % %   correlation coefficient between soda ssh and  sst
        if readsst==1
            for i=1:len_lon
                for j=1:len_lat
                    temp_corr=corrcoef(squeeze(comb_interped_data(i,j,:))',squeeze(comb_interped_sst(i,j,:))');
                    corr_interped_sst(i,j)=temp_corr(1,2);
                end
            end
            disp('corr coef_sst complete') 
        end
        
% % %   correlation coefficient between soda_filtered and cmems_filtered 
        for i=1:len_lon
            for j=1:len_lat
                temp_corr=corrcoef(squeeze(comb_interped_data_filtered(i,j,:))',squeeze(comb_cmems_data_filtered(i,j,:))');
                corr_interped_filtered(i,j)=temp_corr(1,2);
            end
        end
        disp('corr coef_filtered complete') 

% % %   correlation coefficient between model_climatology and cmems_climatology 
        for i=1:len_lon
            for j=1:len_lat
                temp_corr=corrcoef(squeeze(comb_spatial_meaninterped(i,j,:))',squeeze(comb_spatial_meancmems(i,j,:))');
                corr_spatial_mean(i,j)=temp_corr(1,2);
            end
        end
        disp('corr coef_spatial_mean complete')      
        
% % %   correlation coefficient between soda_filtered and soda_sst_filtered 
        if readsst==1
            for i=1:len_lon
                for j=1:len_lat
                    temp_corr=corrcoef(squeeze(comb_interped_data_filtered(i,j,:))',squeeze(comb_interped_sst_filtered(i,j,:))');
                    corr_interped_sst_filtered(i,j)=temp_corr(1,2);
                end
            end
            disp('corr coef_sst_filtered complete') 
        end
        
% % %   correlation coefficient between soda_low_passed and cmems_low_passed
%         nyears=[2,3,5];
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            for i=1:len_lon
                for j=1:len_lat
%                     eval(['temp_corr=corrcoef(squeeze(comb_interped_', num2str(nyear), 'y_lowpass(i,j,:)),squeeze(comb_cmems_,num2str(nyear),'y_lowpass(i,j,:)));'])
                    eval(['temp_corr=corrcoef(squeeze(comb_interped_', num2str(nyear), 'y_lowpass(i,j,:)),squeeze(comb_cmems_', num2str(nyear),'y_lowpass(i,j,:)),','''', 'rows','''',',','''','complete','''',');']);
                    eval(['numnan_interped=sum(~isnan(comb_interped_', num2str(nyear), 'y_lowpass(i,j,:)));']);
                    eval(['numnan_cmems=sum(~isnan(comb_cmems_', num2str(nyear), 'y_lowpass(i,j,:)));']);
                    eval(['tlen=length(comb_cmems_', num2str(nyear), 'y_lowpass(i,j,:));']);
                    if (numnan_interped < tlen-nyear*12 || numnan_cmems < tlen-nyear*12 )
                        eval(['corr_interped_',num2str(nyear),'y_lowpass(i,j)=NaN;'])
                    else
                        eval(['corr_interped_',num2str(nyear),'y_lowpass(i,j)=temp_corr(1,2);'])
                    end
% %                     corrcoef(squeeze(comb_interped_8y_lowpass(95,82,:)),squeeze(comb_cmems_8y_lowpass(95,82,:)),'rows','complete')
                end
            end
%             eval(['m_corr_',num2str(nyear),'y_lowpass=mean(corr_interped_',num2str(nyear),'y_lowpass(:),','''','omitnan','''',');']);
        end
disp('corr coef_low_passed complete') 
        
% % %   correlation coefficient between climatological ssh and climatological cmems ssh 
        for i=1:len_lon
            for j=1:len_lat
                for k=1:12
                    temp_corr=corrcoef(squeeze(comb_interped_clim_divided(i,j,k,:))',squeeze(comb_cmems_clim_divided(i,j,k,:))');
                    corr_clim(i,j,k)=temp_corr(1,2);
                end
            end
        end
        disp('corr coef_clim complete') 
        
        
% % %   correlation coefficient between soda_ssh_filtered_detrended and cmems_filtered_detrended
    
        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
        end
        
        for i=1:len_lon
            for j=1:len_lat
                comb_interped_detrended(i,j,:)=comb_interped_data_filtered(i,j,:)-mean(comb_interped_data_filtered(i,j,:));
                ppp=polyfit(xData,squeeze(comb_interped_detrended(i,j,:))',1);
                comb_interped_linear(i,j,:)=xData*ppp(1)+ppp(2);
                comb_interped_detrended(i,j,:)=comb_interped_detrended(i,j,:)-comb_interped_linear(i,j,:);
                comb_cmems_detrended(i,j,:)=comb_cmems_data_filtered(i,j,:)-mean(comb_cmems_data_filtered(i,j,:));
                ppp=polyfit(xData,squeeze(comb_cmems_detrended(i,j,:))',1);
                comb_cmems_linear(i,j,:)=xData*ppp(1)+ppp(2);
                comb_cmems_detrended(i,j,:)=comb_cmems_detrended(i,j,:)-comb_cmems_linear(i,j,:);
            end
        end
        
        for i=1:len_lon
            for j=1:len_lat
                temp_corr=corrcoef(squeeze(comb_interped_detrended(i,j,:))',squeeze(comb_cmems_detrended(i,j,:))');
                corr_interped_detrended(i,j)=temp_corr(1,2);
            end
        end
        disp('corr coef_detrended complete') 

% %        2, 3, 5year detrended data lowpass filter 
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            cutoff_freq=1/(12*nyear);
            [filt_coef_b,filt_coef_a]=butter(filt_order,cutoff_freq/nq_freq, ftype); %butterworth filter
            for i=1:len_lon
                for j=1:len_lat
                    eval(['comb_interped_detrended_',num2str(nyear),'y_lowpass(i,j,:)=filter(filt_coef_b,filt_coef_a,comb_interped_detrended(i,j,:));'])
                    eval(['comb_cmems_detrended_',num2str(nyear),'y_lowpass(i,j,:)=filter(filt_coef_b,filt_coef_a,comb_cmems_detrended(i,j,:));'])
                    eval(['comb_interped_detrended_',num2str(nyear),'y_lowpass(i,j,1:6*nyear)=NaN;'])
                    eval(['comb_interped_detrended_',num2str(nyear),'y_lowpass(i,j,end-6*nyear+1:end)=NaN;'])
                    eval(['comb_cmems_detrended_',num2str(nyear),'y_lowpass(i,j,1:6*nyear)=NaN;'])
                    eval(['comb_cmems_detrended_',num2str(nyear),'y_lowpass(i,j,end-6*nyear+1:end)=NaN;'])
                end
            end
        end
        
        % % %   correlation coefficient between detrended soda_low_passed and detrended cmems_low_passed
%         nyears=[2,3,5];
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            for i=1:len_lon
                for j=1:len_lat
%                     eval(['temp_corr=corrcoef(squeeze(comb_interped_detrended_', num2str(nyear), 'y_lowpass(i,j,:)),squeeze(comb_cmems_detrended_', num2str(nyear), 'y_lowpass(i,j,:)));'])
                    eval(['temp_corr=corrcoef(squeeze(comb_interped_detrended_', num2str(nyear), 'y_lowpass(i,j,:)),squeeze(comb_cmems_detrended_', num2str(nyear), 'y_lowpass(i,j,:)),','''', 'rows','''',',','''','complete','''',');'])
                    eval(['numnan_interped=sum(~isnan(comb_interped_detrended_', num2str(nyear), 'y_lowpass(i,j,:)));']);
                    eval(['numnan_cmems=sum(~isnan(comb_cmems_detrended_', num2str(nyear), 'y_lowpass(i,j,:)));']);
                    eval(['tlen=length(comb_cmems_detrended_', num2str(nyear), 'y_lowpass(i,j,:));']);
                    if (numnan_interped < tlen-nyear*12 || numnan_cmems < tlen-nyear*12 )
                        eval(['corr_interped_detrended_',num2str(nyear),'y_lowpass(i,j)=NaN;'])
                    else
                        eval(['corr_interped_detrended_',num2str(nyear),'y_lowpass(i,j)=temp_corr(1,2);'])
                    end
                end
            end 
%             eval(['m_corr_',num2str(nyear),'y_lowpass=mean(corr_interped_',num2str(nyear),'y_lowpass(:),','''','omitnan','''',');']);
        end
        disp('corr detrended coef_low_passed complete') 
        
        
% %         figure;
% %         pcolor(corr_clim(:,:,2)')
% %         shading flat
% %         colorbar
% %         disp('corr coef_detrended complete') 
        
        

        
%         for l=1:300
%             tempm=squeeze(comb_interped_data_filtered(:,:,l));
%             m_interped(l)=mean(tempm(:),'omitnan');
%             tempm=squeeze(comb_cmems_data_filtered(:,:,l));
%             m_cmems(l)=mean(tempm(:),'omitnan');
%         end
%         corrcoef(m_interped, m_cmems)
%         pcolor(corr_interped')
%         shading flat
%         colorbar;

%         save([filedir,testname,'_',regionname,'ssh_trend_corr_coef_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'], '-v7.3', '-nocompression');
        save([filedir,testname,'_',regionname,'ssh_trend_corr_coef_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'], '-v7.3');
%         load([filedir,testname,'_',regionname,'ssh_trend_corr_coef_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
        
        
        len_lon_cmems=length(cmems_lon);
        len_lat_cmems=length(cmems_lat);
        ncid = netcdf.create(strcat(filedir, testname,'_',regionname, '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc'),'NETCDF4');

        onedimid = netcdf.defDim(ncid,'one', 1);
        lon_dimid = netcdf.defDim(ncid, 'lon', len_lon_soda);
        lat_dimid = netcdf.defDim(ncid,'lat',len_lat_soda);
        xidimid = netcdf.defDim(ncid, 'xi_rho', len_lon_soda);
        etadimid = netcdf.defDim(ncid,'eta_rho',len_lat_soda);
        lon_cmems_dimid = netcdf.defDim(ncid, 'lon_cmems', len_lon_cmems);
        lat_cmems_dimid = netcdf.defDim(ncid,'lat_cmems', len_lat_cmems);
        time_dimid = netcdf.defDim(ncid, 'time', 0);
        clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
            'type', ['NWP 1/20 _ ', testname, 'soda, recon monthly SSH analysis file']);
        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
            'title', [' monthly SSH analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

        lonvarid=netcdf.defVar(ncid, 'lon', 'NC_DOUBLE', [lon_dimid]);
        netcdf.putAtt(ncid,lonvarid,'long_name','longitude');
        netcdf.putAtt(ncid,lonvarid,'units','degree_east');

        latvarid=netcdf.defVar(ncid, 'lat', 'NC_DOUBLE', [lat_dimid]);
        netcdf.putAtt(ncid,latvarid,'long_name','latitude');
        netcdf.putAtt(ncid,latvarid,'units','degree_north');
        
        lon_cmemsvarid=netcdf.defVar(ncid, 'lon_cmems', 'NC_DOUBLE', lon_cmems_dimid);
        netcdf.putAtt(ncid,lon_cmemsvarid,'long_name','longitude');
        netcdf.putAtt(ncid,lon_cmemsvarid,'units','degree_east');

        lat_cmemsvarid=netcdf.defVar(ncid, 'lat_cmems', 'NC_DOUBLE', lat_cmems_dimid);
        netcdf.putAtt(ncid,lat_cmemsvarid,'long_name','latitude');
        netcdf.putAtt(ncid,lat_cmemsvarid,'units','degree_north');
        
        lon_rhovarid=netcdf.defVar(ncid, 'lon_rho', 'NC_DOUBLE', [xidimid etadimid]);
        netcdf.putAtt(ncid,lon_rhovarid,'long_name','lon_soda');
        netcdf.putAtt(ncid,lon_rhovarid,'units','degree_east');

        lat_rhovarid=netcdf.defVar(ncid, 'lat_rho', 'NC_DOUBLE', [xidimid etadimid]);
        netcdf.putAtt(ncid,lat_rhovarid,'long_name','lat_soda');
        netcdf.putAtt(ncid,lat_rhovarid,'units','degree_north');

        raw_sshvarid=netcdf.defVar(ncid, 'raw_ssh', 'NC_FLOAT', [xidimid etadimid time_dimid]);
        netcdf.putAtt(ncid,raw_sshvarid,'long_name','raw_ssh');
        netcdf.putAtt(ncid,raw_sshvarid,'units','m');
        
        raw_sstvarid=netcdf.defVar(ncid, 'raw_sst', 'NC_FLOAT', [xidimid etadimid time_dimid]);
        netcdf.putAtt(ncid,raw_sstvarid,'long_name','raw_sst');
        netcdf.putAtt(ncid,raw_sstvarid,'units','^oC');
        
        ssh_filteredvarid=netcdf.defVar(ncid, 'ssh_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,ssh_filteredvarid,'long_name','ssh_filtered');
        netcdf.putAtt(ncid,ssh_filteredvarid,'units','m');
        if readsst==1
            sst_filteredvarid=netcdf.defVar(ncid, 'sst_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
            netcdf.putAtt(ncid,sst_filteredvarid,'long_name','sst_filtered');
            netcdf.putAtt(ncid,sst_filteredvarid,'units','^oC');
        end
        trendvarid=netcdf.defVar(ncid, 'trend', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,trendvarid,'long_name','trend');
        netcdf.putAtt(ncid,trendvarid,'units','mm/year');

        trend_filteredvarid=netcdf.defVar(ncid, 'trend_filtered', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,trend_filteredvarid,'long_name','trend_filtered');
        netcdf.putAtt(ncid,trend_filteredvarid,'units','mm/year');

        mean_trendvarid=netcdf.defVar(ncid, 'mean_trend', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_trendvarid,'long_name','mean_trend');
        netcdf.putAtt(ncid,mean_trendvarid,'units','mm/year');

        mean_trend_filteredvarid=netcdf.defVar(ncid, 'mean_trend_filtered', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_trend_filteredvarid,'long_name','mean_trend_filtered');
        netcdf.putAtt(ncid,mean_trend_filteredvarid,'units','mm/year');

        clim_sshvarid=netcdf.defVar(ncid, 'clim_ssh', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_sshvarid,'long_name','clim_ssh');
        netcdf.putAtt(ncid,clim_sshvarid,'units','m');
        
        clim_ssh_trendvarid=netcdf.defVar(ncid, 'clim_ssh_trend', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_ssh_trendvarid,'long_name','clim_ssh_trend');
        netcdf.putAtt(ncid,clim_ssh_trendvarid,'units','m');
        
        interped_sshvarid=netcdf.defVar(ncid, 'interped_ssh', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid time_dimid]);
        netcdf.putAtt(ncid,interped_sshvarid,'long_name','interped_ssh');
        netcdf.putAtt(ncid,interped_sshvarid,'units','m ');
        
        interped_slavarid=netcdf.defVar(ncid, 'interped_sla', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid time_dimid]);
        netcdf.putAtt(ncid,interped_slavarid,'long_name','interped_sla');
        netcdf.putAtt(ncid,interped_slavarid,'units','m ');
        
        interped_sla_filteredvarid=netcdf.defVar(ncid, 'interped_sla_filtered', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid time_dimid]);
        netcdf.putAtt(ncid,interped_sla_filteredvarid,'long_name','interped_sla_filtered');
        netcdf.putAtt(ncid,interped_sla_filteredvarid,'units','m ');
        
        cmems_slavarid=netcdf.defVar(ncid, 'cmems_sla', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid time_dimid]);
        netcdf.putAtt(ncid,cmems_slavarid,'long_name','cmems_sla');
        netcdf.putAtt(ncid,cmems_slavarid,'units','m ');

        cmems_adtvarid=netcdf.defVar(ncid, 'cmems_adt', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid time_dimid]);
        netcdf.putAtt(ncid,cmems_adtvarid,'long_name','cmems_adt');
        netcdf.putAtt(ncid,cmems_adtvarid,'units','m ');
        netcdf.putAtt(ncid,cmems_adtvarid,'source','mdt from AVISO + sla');
    
        clim_cmemsvarid=netcdf.defVar(ncid, 'clim_cmems_ssh', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_cmemsvarid,'long_name','clim_cmems_ssh');
        netcdf.putAtt(ncid,clim_cmemsvarid,'units','m ');
    
        cmems_ssh_filteredvarid=netcdf.defVar(ncid, 'cmems_ssh_filtered', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid time_dimid]);
        netcdf.putAtt(ncid,cmems_ssh_filteredvarid,'long_name','cmems_ssh_filtered');
        netcdf.putAtt(ncid,cmems_ssh_filteredvarid,'units','m ');
        
        nc_varname_prefixes={'cmems_', 'interped_', 'interped_sla_', 'cmems_detrended_', 'interped_detrended_'};
        for nc_varnameij=1:length(nc_varname_prefixes)
            nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
            for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            nc_varname=[nc_varname_prefix, num2str(nyear), 'y_lowpass'];
            eval([nc_varname,'varid=netcdf.defVar(ncid,','''', ...
                nc_varname,'''',',', '''', 'NC_FLOAT','''',',', '[lon_cmems_dimid lat_cmems_dimid time_dimid]);'])
            eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'long_name','''', ',', '''', nc_varname,'''',');'])
            eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'units'    ,'''', ',', '''', 'm', '''', ');']);
            end
        end
        
        cmems_trendvarid=netcdf.defVar(ncid, 'cmems_trend', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
        netcdf.putAtt(ncid,cmems_trendvarid,'long_name','cmems_trend');
        netcdf.putAtt(ncid,cmems_trendvarid,'units','mm /year');

        cmems_trend_filteredvarid=netcdf.defVar(ncid, 'cmems_trend_filtered', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
        netcdf.putAtt(ncid,cmems_trend_filteredvarid,'long_name','cmems_trend_filtered');
        netcdf.putAtt(ncid,cmems_trend_filteredvarid,'units','mm /year');

        clim_cmems_trend_dividedvarid=netcdf.defVar(ncid, 'clim_cmems_trend_divided', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_cmems_trend_dividedvarid,'long_name','clim_cmems_trend_divided');
        netcdf.putAtt(ncid,clim_cmems_trend_dividedvarid,'units','mm /year');
        
        interped_trendvarid=netcdf.defVar(ncid, 'interped_trend', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
        netcdf.putAtt(ncid,interped_trendvarid,'long_name','interped_trend');
        netcdf.putAtt(ncid,interped_trendvarid,'units','mm /year');
        if readsst==1
            interped_sst_trendvarid=netcdf.defVar(ncid, 'interped_sst_trend', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
            netcdf.putAtt(ncid,interped_sst_trendvarid,'long_name','interped_sst_trend');
            netcdf.putAtt(ncid,interped_sst_trendvarid,'units','^oC /year');
        end
        interped_trend_filteredvarid=netcdf.defVar(ncid, 'interped_trend_filtered', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
        netcdf.putAtt(ncid,interped_trend_filteredvarid,'long_name','interped_trend_filtered_trend');
        netcdf.putAtt(ncid,interped_trend_filteredvarid,'units','mm /year');
        if readsst==1
            interped_sst_trend_filteredvarid=netcdf.defVar(ncid, 'interped_sst_trend_filtered', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
            netcdf.putAtt(ncid,interped_sst_trend_filteredvarid,'long_name','interped_sst_trend_filtered_trend');
            netcdf.putAtt(ncid,interped_sst_trend_filteredvarid,'units','^oC /year');
        end
        corr_interpedvarid=netcdf.defVar(ncid, 'corr_interped', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
        netcdf.putAtt(ncid,corr_interpedvarid,'long_name','corr_interped');
        netcdf.putAtt(ncid,corr_interpedvarid,'units',' ');
        if readsst==1
            corr_interped_sstvarid=netcdf.defVar(ncid, 'corr_interped_sst', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
            netcdf.putAtt(ncid,corr_interped_sstvarid,'long_name','corr_interped_sst');
            netcdf.putAtt(ncid,corr_interped_sstvarid,'units',' ');
        end
        corr_interped_filteredvarid=netcdf.defVar(ncid, 'corr_interped_filtered', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
        netcdf.putAtt(ncid,corr_interped_filteredvarid,'long_name','corr_interped_filtered');
        netcdf.putAtt(ncid,corr_interped_filteredvarid,'units',' ');
        
        corr_interped_sst_filteredvarid=netcdf.defVar(ncid, 'corr_interped_sst_filtered', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
        netcdf.putAtt(ncid,corr_interped_sst_filteredvarid,'long_name','corr_interped_sst_filtered');
        netcdf.putAtt(ncid,corr_interped_sst_filteredvarid,'units',' ');
        
        corr_spatial_meanvarid=netcdf.defVar(ncid, 'corr_spatial_mean', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
        netcdf.putAtt(ncid,corr_spatial_meanvarid,'long_name','corr_spatial_mean');
        netcdf.putAtt(ncid,corr_spatial_meanvarid,'units',' ');
        
        corr_interped_detrendedvarid=netcdf.defVar(ncid, 'corr_interped_detrended', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
        netcdf.putAtt(ncid,corr_interped_detrendedvarid,'long_name','corr_interped_detrended');
        netcdf.putAtt(ncid,corr_interped_detrendedvarid,'units',' ');
        
        corr_climvarid=netcdf.defVar(ncid, 'corr_clim', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,corr_climvarid,'long_name','corr_clim');
        netcdf.putAtt(ncid,corr_climvarid,'units',' ');
        
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            nc_varname=['corr_interped_',num2str(nyear),'y_lowpass'];
            eval([nc_varname,'varid=netcdf.defVar(ncid,','''', ...
                nc_varname,'''',',', '''', 'NC_FLOAT','''',',', '[lon_cmems_dimid lat_cmems_dimid]);'])
            eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'long_name','''', ',', '''', nc_varname,'''',');'])
            eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'units'    ,'''', ',', '''', ' ', '''', ');']);
        end
         
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            nc_varname=['corr_interped_detrended_',num2str(nyear),'y_lowpass'];
            eval([nc_varname,'varid=netcdf.defVar(ncid,','''', ...
                nc_varname,'''',',', '''', 'NC_FLOAT','''',',', '[lon_cmems_dimid lat_cmems_dimid]);'])
            eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'long_name','''', ',', '''', nc_varname,'''',');'])
            eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'units'    ,'''', ',', '''', ' ', '''', ');']);
         end
        
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
        netcdf.putVar(ncid, lon_rhovarid, [0 0], [len_lon_soda len_lat_soda], lon);
        netcdf.putVar(ncid, lat_rhovarid, [0 0], [len_lon_soda len_lat_soda], lat);
        netcdf.putVar(ncid, raw_sshvarid, [0 0 0], [len_lon_soda len_lat_soda length(ftime)], comb_data);
        netcdf.putVar(ncid, ssh_filteredvarid, [0 0 0], [len_lon_soda len_lat_soda length(ftime)], comb_data_filtered);
        if readsst==1
            netcdf.putVar(ncid, raw_sstvarid, [0 0 0], [len_lon_soda len_lat_soda length(ftime)], comb_sst);
            netcdf.putVar(ncid, sst_filteredvarid, [0 0 0], [len_lon_soda len_lat_soda length(ftime)], comb_sst_filtered);
        end
        netcdf.putVar(ncid, trendvarid, [0 0], [len_lon_soda, len_lat_soda], trend);
        netcdf.putVar(ncid, trend_filteredvarid, [0 0], [len_lon_soda, len_lat_soda], trend_filtered);
        netcdf.putVar(ncid, mean_trendvarid, [0], [1], mean_trend);
        netcdf.putVar(ncid, mean_trend_filteredvarid, [0], [1], mean_trend_filtered);
        netcdf.putVar(ncid, clim_sshvarid, [0 0 0], [len_lon_soda, len_lat_soda length(climtime)], comb_spatial_meansoda);
        netcdf.putVar(ncid, clim_ssh_trendvarid, [0 0 0], [len_lon_soda, len_lat_soda length(climtime)], trend_clim);
        netcdf.putVar(ncid, lon_cmemsvarid, 0, len_lon_cmems, cmems_lon(:));
        netcdf.putVar(ncid, lat_cmemsvarid, 0, len_lat_cmems, cmems_lat(:));
        netcdf.putVar(ncid, interped_sshvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], comb_interped_data);
        netcdf.putVar(ncid, interped_slavarid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], interped_sla);        
        netcdf.putVar(ncid, cmems_slavarid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], comb_cmems_data);
        netcdf.putVar(ncid, cmems_adtvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], comb_cmems_adt);
        netcdf.putVar(ncid, cmems_ssh_filteredvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], comb_cmems_data_filtered);
        netcdf.putVar(ncid, interped_sla_filteredvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], interped_sla_filtered);
        nc_varname_prefixes={'cmems_', 'interped_', 'interped_sla_', 'cmems_detrended_', 'interped_detrended_'};
        for nc_varnameij=1:length(nc_varname_prefixes)
            nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                eval(['netcdf.putVar(ncid,', nc_varname,'varid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], comb_', nc_varname, ');']);
            end
        end
        nc_varname_prefixes={'corr_interped_', 'corr_interped_detrended_'};
        for nc_varnameij=1:length(nc_varname_prefixes)
            nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                nc_varname=[nc_varname_prefix,num2str(nyear),'y_lowpass'];
                eval(['netcdf.putVar(ncid,', nc_varname,'varid, [0 0], [len_lon_cmems len_lat_cmems],', nc_varname, ');']);
            end
        end
        netcdf.putVar(ncid, clim_cmemsvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(climtime)], comb_spatial_meancmems);
        netcdf.putVar(ncid, cmems_trendvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_trend);
        netcdf.putVar(ncid, cmems_trend_filteredvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_trend_filtered);
        netcdf.putVar(ncid, clim_cmems_trend_dividedvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(climtime)], clim_cmems_trend_divided);
        netcdf.putVar(ncid, interped_trendvarid, [0 0], [len_lon_cmems len_lat_cmems], interped_trend);
        
        if readsst==1
            netcdf.putVar(ncid, interped_sst_trendvarid, [0 0], [len_lon_cmems len_lat_cmems], interped_sst_trend);
            netcdf.putVar(ncid, interped_sst_trend_filteredvarid, [0 0], [len_lon_cmems len_lat_cmems], interped_sst_trend_filtered);
            netcdf.putVar(ncid, corr_interped_sstvarid, [0 0], [len_lon_cmems len_lat_cmems], corr_interped_sst);
            netcdf.putVar(ncid, corr_interped_sst_filteredvarid, [0 0], [len_lon_cmems len_lat_cmems], corr_interped_sst_filtered);
        end
        netcdf.putVar(ncid, interped_trend_filteredvarid, [0 0], [len_lon_cmems len_lat_cmems], interped_trend_filtered);
        netcdf.putVar(ncid, corr_interpedvarid, [0 0], [len_lon_cmems len_lat_cmems], corr_interped);
        netcdf.putVar(ncid, corr_interped_filteredvarid, [0 0], [len_lon_cmems len_lat_cmems], corr_interped_filtered);
        netcdf.putVar(ncid, corr_spatial_meanvarid, [0 0], [len_lon_cmems len_lat_cmems], corr_spatial_mean);
        netcdf.putVar(ncid, corr_interped_detrendedvarid, [0 0], [len_lon_cmems len_lat_cmems], corr_interped_detrended);
        netcdf.putVar(ncid, corr_climvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(climtime)], corr_clim);
        netcdf.close(ncid);
    end
end
