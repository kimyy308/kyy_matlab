clc;clear all; close all;
%set variables
model_name = {'CanESM2','CESM1-BGC','EC-EARTH','HadGEM2-CC','NorESM1-M','MPI-ESM-P','CESM1-WACCM'};

obs_path = '/home/ygkim/ipcc/data/heatcontent/tanom_year'; %tanom_ path

root_path = '/data/ygkim/tmp/';
data_type = 'Omon' ; %Amom Omon

sena_list = {'historical'}; %{'rcp45', 'rcp85'};
vari_list = {'thetao','uo','vo','so'};

for sena = 1:length(sena_list)
for vari = 1:length(vari_list)    
    
    variable = vari_list{vari};
    senario = sena_list{sena};
    data_path = strcat(root_path,variable,'/',senario,'/',data_type,'/');
    
    %Setting variable and path
    %i : model number
    for i = 1:length(model_name)
           %--------STEP1. Load obs grid data--------%
            nc = netcdf(strcat(obs_path,'/','tanom_5555.nc'));
            lat = nc{'lat'}(:);         %unit : degrees_north
            lon = nc{'lon'}(:)+180;     %unit : degrees_east
            clear nc;
            %salinity file
            nc = netcdf(strcat('/home/ygkim/ipcc/data/heatcontent','/climatology/salinity_annual_1deg.nc'));
            depth = nc{'depth'}(:);     %unit : meters
            
            %find NWP area 
             lon_idx = [115.5,164.5]; %lon : 115 ~164
             lat_idx = [15.5,52.5];   %lat : 15 ~ 52

             lon_fin = find(lon>=lon_idx(1) & lon<= lon_idx(2));
             lat_fin = find(lat>=lat_idx(1) & lat <=lat_idx(2));
             lon = lon(lon_fin);
             lat = lat(lat_fin);

            [lat_obs, depth_obs, lon_obs] = meshgrid(lat,depth,lon); 

           %--------STEP2. Setting file path,time range ...------%
            %load thetao folder
            file_thetao =dir(strcat(data_path,variable,'_',data_type,'_',model_name{i},'*'));

            %set time range
            %-------determine range-------%
            ini_yr = 1976 ;ini_mo = 1 ; %initial year,month
            end_yr = 2005 ;end_mo = 12; %end year,month
            %-----------------------------%
            date = file_thetao(1).name;
            date_yr = str2num(date(end-15:end-12));
            date_mo = str2num(date(end-11:end-10));   
            ini_time = (12*(ini_yr-date_yr)+(ini_mo-date_mo)+1); 
            end_time = (12*(end_yr-date_yr)+(end_mo-date_mo)+1);
            num = 0;num2 = 0;

           %-----STEP3. Load Model grid information ------%
            %3.1 load reference nc file
            nc = netcdf(strcat(data_path,file_thetao(1).name));
            depth = nc{'lev'}(:); %depth
            lat = nc{'lat'}(:);   %lattitude
            lon = nc{'lon'}(:);   %longitude
            dat = nc{variable}(1,:,:,:);
            dat(dat>10000) = NaN;
            fin = find(isnan(dat)==1);

            %3.2 convert to lon,lat format
             %longitude : -180~+180 => 0~360 degree
            if min(lon(:)) <0 
                lon = lon+360;
                for ii = 1 :length(lon)
                    if lon(ii)>360
                        lon(ii) =lon(ii)-360;
                    end
                end
            end

            if max(lon(:)) >360
               idx = find(lon>360);
               lon(idx) = lon(idx)-360;
            end


            % make grid lon,lat (if lon,lat are not grid data)
            [m, n] = size(lat);
            if n == 1
                [lon, lat] = meshgrid(lon,lat);
            end

            %3.3 Coverting to wanted 3-D grid format
            [m,n] = size(lat);
            depth_mod = zeros(length(depth),m,n);
            lat_mod = zeros(length(depth),m,n);
            lon_mod = zeros(length(depth),m,n);

            for zz = 1:length(depth)
                depth_mod(zz,:,:) = depth_mod(zz,:,:) + depth(zz);
                lat_mod(zz,:,:) = lat(:,:);
                lon_mod(zz,:,:) = lon(:,:);
            end

            %3.4. convert depth(sigma cooridinate)
            if max(depth) <10
                depth2 = nc{'depth'}(:);
                depth3 = zeros(size(depth_mod));
                for jj = 1: 40
                    depth3(jj,:) = depth2(:)*depth(jj);
                end
                depth_mod = -depth3;
            end
            depth_mod(fin) = [];
            lon_mod(fin) = [];
            lat_mod(fin) = [];
      
            %Limit to NWP area
            fin2 = find(lon_mod<110 | lon_mod>170 | lat_mod <10 | lat_mod>60);
            lon_mod(fin2) = [];
            lat_mod(fin2) = [];
            depth_mod(fin2) = [];
            %STEP4. interpolation process
            if ini_time ==1
                [mm ,nn, qq] = size(lat_obs); %size of obs data
                thetao = zeros((end_time-ini_time+1),mm,nn,qq);
                tim_index = 0;
                % j : file number of same model

                for j = 1:length(file_thetao)
                    nc = netcdf(strcat(data_path,file_thetao(j).name));
                    time =nc{'time'}(:);


                    for k = 1:length(time)
                        tim_index = tim_index+1;
                        if tim_index > end_time-ini_time+1
                            break
                        end
                        data = squeeze(nc{variable}(k,:,:,:));
                        %Nan data process
                           data(data>10000) = NaN;

                           data(fin) = [];data(fin2) = [];
                           data = griddata(lon_mod,lat_mod,depth_mod,data,lon_obs,lat_obs,depth_obs);
                           thetao(tim_index,:,:,:) = data(:,:,:);
                    end
                end
            else
                [mm ,nn, qq] = size(lat_obs); %size of obs data
                thetao = zeros((end_time-ini_time+1),mm,nn,qq);
                tim_index = 0;
                % j : file number of same model

                for j = 1:length(file_thetao)
                    nc = netcdf(strcat(data_path,file_thetao(j).name));
                    time =nc{'time'}(:);


                    for k = 1:length(time)
                        tim_index = tim_index+1;
                        data = squeeze(nc{variable}(k,:,:,:));
                        %Nan data process
                           data(data>10000) = NaN;

                           data(fin) = [];data(fin2) = [];
                           data = griddata(lon_mod,lat_mod,depth_mod,data,lon_obs,lat_obs,depth_obs);
                           thetao(tim_index,:,:,:) = data(:,:,:);
                    end
                end
                thetao(1,:,:,:) = [];
            end
            %------------STEP5. Save preprocess data-----------------%
            savepath = strcat(data_path,model_name{i},'_','preprocess.mat');
            parsave(savepath,thetao);
            disp(i);
    end
end
end
