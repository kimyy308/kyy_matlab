clc;clear all;close all;
%set variables
model_name = {'CCSM4','CESM1-CAM5', 'CSIRO-Mk3-6-0', 'GFDL-ESM2G', 'GFDL-ESM2M', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'MIROC-ESM', 'MIROC5', 'MIROC-ESM-CHEM', 'NorESM1-M', 'NorESM1-ME', 'bcc-csm1-1', 'bcc-csm1-1-m', 'CESM1-WACCM', 'CNRM-CM5', 'CanESM2', 'EC-EARTH', 'FGOALS-g2', 'GISS-E2-R', 'MPI-ESM-LR,', 'MPI-ESM-MR', 'MRI-CGCM3'};

obs_path = '/home/ygkim/ipcc/data/heatcontent/tanom_year'; %tanom_ path

% root_path = '/data/ygkim/tmp/';
root_path = '/data/ygkim/tmp1/';
data_type = 'Omon' ; %Amom Omon
% data_type = 'Amon' ; %Amom Omon

% sena_list = {'rcp45', 'rcp85'};
sena_list = {'historical'};
% vari_list = {'hurs','psl','rsds','uas','vas'};
vari_list = {'thetao','uo','vo','so','zos'};

zzz = 0

for sena = 1:length(sena_list)
for vari = 1:length(vari_list)

    variable = vari_list{vari};
    senario = sena_list{sena};

%     data_path = strcat(root_path,variable,'/',senario,'/',data_type,'/');
    data_path = strcat(root_path,variable,'/',senario,'/',data_type,'/');

    %Setting variable and path
    %i : model number
    zzz = 1
    for i = 1:5
           %--------STEP1. Load obs grid data--------%
            nc = netcdf(strcat(obs_path,'/','tanom_5555.nc'));
            lat = nc{'lat'}(:);         %unit : degrees_north
            lon = nc{'lon'}(:)+180;     %unit : degrees_east
            
            
            %find NWP area 
             lon_idx = [115.5,164.5]; %lon : 115 ~164
             lat_idx = [15.5,52.5];   %lat : 15 ~ 52

             lon_fin = find(lon>=lon_idx(1) & lon<= lon_idx(2));
             lat_fin = find(lat>=lat_idx(1) & lat <=lat_idx(2));
             lon = lon(lon_fin);
             lat = lat(lat_fin);
             
            [lon_obs,lat_obs] = meshgrid(lon,lat); 
            

           %--------STEP2. Setting file path,time range ...------%
            %load zos folder
            file_zos =dir(strcat(data_path,variable,'_',data_type,'_',model_name{i},'*'));
            
            %set time range
            %-------determine range-------%
            ini_yr = 1980 ;ini_mo = 1 ; %initial year,month
            end_yr = 2005 ;end_mo = 12; %end year,month
            %-----------------------------%
            
            date = file_zos(1).name;
            date_yr = str2num(date(end-15:end-12));
            date_mo = str2num(date(end-11:end-10));   
            ini_time = (12*(ini_yr-date_yr)+(ini_mo-date_mo)+1); 
            end_time = (12*(end_yr-date_yr)+(end_mo-date_mo)+1);
            num = 0;num2 = 0;
            
           %-----STEP3. Load Model grid information ------%
            %3.1 load reference nc file
            nc = netcdf(strcat(data_path,file_zos(1).name));

            lat = nc{'lat'}(:);   %lattitude
            lon = nc{'lon'}(:);   %longitude
            dat = nc{variable}(1,:,:);
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
            
            lon_mod = lon;
            lat_mod = lat;

            
            % make grid lon,lat (if lon,lat are not grid data)
            [m, n] = size(lat);
            if n == 1
                [lon_mod, lat_mod] = meshgrid(lon,lat);
            end
            lon_mod(fin) = [];
            lat_mod(fin) = [];
            zzz = 2
            
            %Limit to NWP area
            fin2 = find(lon_mod<110 | lon_mod>160 | lat_mod <10 | lat_mod>60);
            lon_mod(fin2) = [];
            lat_mod(fin2) = [];

            %STEP4. interpolation process
            [mm ,nn] = size(lat_obs); %size of obs data
            zos = zeros((end_time-ini_time+1),mm,nn);
            tim_index = 0;
            % j : file number of same model
            
            for j = 1:length(file_zos)
                nc = netcdf(strcat(data_path,file_zos(j).name));
                time =nc{'time'}(:);
               
        
                for k = 1:length(time)
                    tim_index = tim_index+1;
                    if tim_index > end_time-ini_time+1
                        break
                    end
                    data = squeeze(nc{variable}(k,:,:));
                    %Nan data process
                       data(data>10000) = NaN;
                       
                       data(fin) = [];data(fin2) = [];
                       data = griddata(lon_mod,lat_mod,data,lon_obs,lat_obs);
                       zos(tim_index,:,:,:) = data(:,:,:);
                end
            end
 
            %------------STEP5. Save preprocess data-----------------%
            savepath = strcat(data_path,model_name{i},'_','preprocess_csh.mat');
            parsave(savepath,zos);
            disp(i);
    end
end
end