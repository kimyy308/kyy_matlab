close all; clear all; clc;

system_name=computer;
if (strcmp(system_name,'PCWIN64'))
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox'; %% SNU_desktop, kyy_laptop
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run']));

elseif (strcmp(system_name,'GLNXA64'))
    dropboxpath='/home/kimyy/Dropbox'; %% DAMO
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run']));
end

warning off;

workdir='/data2/kimyy/Reanalysis/SODA/SODA_3_4_2/monthly/';
inputyear=2018;
inputmonth=1:8;
myoceanvarname={'thetao', 'so', 'uo', 'vo', 'zos'};
sodavarname={'temp', 'salt', 'u', 'v', 'ssh'};
% myoceanvarname={'zos'};
% sodavarname={'ssh'};
section= [115-2, 164+2, 15-2, 52+2];

for numyear=1:length(inputyear)
    tempyear=inputyear(numyear);
    rawsodaname=[workdir,'soda3.4.2_mn_ocean_reg_',num2str(2016,'%04i'),'.nc'];
%     rawsodaname=[workdir,'soda3.4.2_mn_ocean_reg_',num2str(year,'%04i'),'.nc'];
    soda_myocean_name=[workdir,'soda3.4.2_plus_myocean_',num2str(tempyear,'%04i'),'.nc'];
    if (system_name=='GLNXA64')  %% LINUX
            system(['cp ',rawsodaname,' ',soda_myocean_name])  %% make a tide file from copied grd file
    elseif (system_name=='WIN64')  %% WIN
        system(['copy ',rawsodaname,' ',soda_myocean_name])  %% make a tide file from copied grd file
    else
        disp('cannot recognize operation type');
    end
    
    for nummonth=1:length(inputmonth)
        tempmonth=inputmonth(nummonth);
        myocean_ano_name=['/data2/kimyy/Reanalysis/myocean/monthly/diff/',num2str(tempyear,'%04i'), ...
            '/MyOcean_monthly_anomaly_from_2016_',num2str(tempyear,'%04i'),num2str(tempmonth,'%02i'),'.nc'];
        for numvar=1:length(myoceanvarname)
            temp_mo_varname=myoceanvarname{numvar};
            temp_soda_varname=sodavarname{numvar};
            
            lon_soda_global = ncread(soda_myocean_name,'xt_ocean');
            lat_soda_global = ncread(soda_myocean_name,'yt_ocean');
            depth_soda = ncread(soda_myocean_name,'st_ocean');

            dl=mean(diff(lat_soda_global));
            [soda_lon_min, soda_lon_max, soda_lat_min, soda_lat_max] = findind_Y(dl, section(1:4), lon_soda_global, lat_soda_global);
            
            lon_soda_section = lon_soda_global(soda_lon_min:soda_lon_max);
            lat_soda_section = lat_soda_global(soda_lat_min:soda_lat_max);
            
% %             read soda 2016 from newly copied file
            
            if (strcmp(temp_mo_varname,'zos')==1)
                data = ncread(soda_myocean_name, temp_soda_varname,  ...
                [soda_lon_min, soda_lat_min, tempmonth], [soda_lon_max-soda_lon_min+1, soda_lat_max-soda_lat_min+1, 1]);
            else   
            data = ncread(soda_myocean_name, temp_soda_varname,  ...
                [soda_lon_min, soda_lat_min, 1, tempmonth], [soda_lon_max-soda_lon_min+1, soda_lat_max-soda_lat_min+1, length(depth_soda), 1]);
            end
            data(data<-1000)=NaN;
            data(data>1000)=NaN;

% %             get anomdata
            lon_mo_global = ncread(myocean_ano_name,'longitude');
            lat_mo_global = ncread(myocean_ano_name,'latitude');
            depth_mo = ncread(myocean_ano_name,'depth');
            
            
            dl=mean(diff(lat_mo_global));
            [lon_min, lon_max, lat_min, lat_max] = findind_Y(dl, section(1:4), lon_mo_global, lat_mo_global);
            
            lon_mo_section = lon_mo_global(lon_min:lon_max);
            lat_mo_section = lat_mo_global(lat_min:lat_max);
            
            if (strcmp(temp_mo_varname,'zos')==1)
                anomdata = ncread(myocean_ano_name, temp_mo_varname, ...
                    [lon_min, lat_min, 1], [lon_max-lon_min+1, lat_max-lat_min+1, 1]);
            else  
                anomdata = ncread(myocean_ano_name, temp_mo_varname, ...
                    [lon_min, lat_min, 1, 1], [lon_max-lon_min+1, lat_max-lat_min+1, length(depth_mo), 1]);
            end
            anomdata(anomdata<-1000)=NaN;
            anomdata(anomdata>1000)=NaN;
                        
% %             horizontal interpolation of anomdata to sodagrid
            
            if (strcmp(temp_mo_varname,'zos')==1)
                interped_anom_data(:,:)=griddata(double(lon_mo_section),double(lat_mo_section)', ...
                        squeeze(anomdata(:,:))', lon_soda_section, lat_soda_section')';
            else  
                for numdepth=1:length(depth_mo)
                    interped_anom_data(:,:,numdepth)=griddata(double(lon_mo_section),double(lat_mo_section)', ...
                        squeeze(anomdata(:,:,numdepth))', lon_soda_section, lat_soda_section')';
                end
            end
           
% %             vertical interpolation of anomdata to sodagrid
            
            if (strcmp(temp_mo_varname,'zos')==1)
               interped_anom_data(isnan(interped_anom_data)==1)=0;
            else  
                for numx=1:size(interped_anom_data,1)
                    for numy=1:size(interped_anom_data,2)
                        interped_anom_data2(numx,numy,:) = ...
                            interp1(depth_mo,squeeze(interped_anom_data(numx,numy,:)),depth_soda);
                    end
                end
                interped_anom_data2(isnan(interped_anom_data2)==1)=0;
            end
            
          
% %             sum anomdata and soda 2016
            if (strcmp(temp_mo_varname,'zos')==1)
                newdata=data+interped_anom_data;
                ncwrite(soda_myocean_name, temp_soda_varname, newdata, ...
                [soda_lon_min, soda_lat_min, tempmonth]);
            else  
                newdata=data+interped_anom_data2;
                ncwrite(soda_myocean_name, temp_soda_varname, newdata, ...
                [soda_lon_min, soda_lat_min, 1, tempmonth]);
            end    
            clear interped_anom_data 
            clear interped_anom_data2
        end
        xData(nummonth) = datenum([num2str(tempyear),'-',num2str(tempmonth,'%02i'),'-15']);
        reftime=datenum([num2str(1980),'-',num2str(01,'%02i'),'-01']);
        xData(nummonth) = xData(nummonth)-reftime;
        
        
    end
% %     change time data
    ncwrite(soda_myocean_name, 'time', xData);
end