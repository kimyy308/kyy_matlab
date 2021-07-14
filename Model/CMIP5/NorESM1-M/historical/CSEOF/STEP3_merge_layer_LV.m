% clc;clear all; close all;
% 
close all; clc; clear all;

load('STEP2_5_output.mat')

filepath = [workdir, '/data/']
% filename = 'increment_v1.5_050.nc'
filepath = [workdir, '/data/NorESM1-M_mon/']  
  filename = [varname, '_mon_', modelname, '_', scenname, '_', num2str(tempyear,'%04i'), '.nc']
file = strcat(filepath,filename)
%nc = netcdf(file)

%lat2 = nc{'YT_OCEAN'}(:)
%lon2 = nc{'XT_OCEAN'}(:);
% lat2 = ncread(file,'lon');
% lon2 = ncread(file,'lat');

%----------STEP1. load lon,lat information-------------%
%load atm lon,lat information
lon_atm = lon;
lat_atm = lat;
lon_uv = lon;
lat_uv = lat;
lon2 = lon;
lat2 = lat;

%% ATM component
%tas
%for layer = 1 : 1
%    %set filepath , variables
%    filepath ='/home/kimyy/ygkim/cseof_analysis/exercise_kyy/cseofs/';
%    LV_name = strcat('LV_increment1','.dat');
%    time = 408 ; % time(month)
%    T = 12  ;  %period
%    mode = 10 ; %number of LV mode
%    lon = lon_atm; %lon
%    lat =  lat_atm; %lat
%    mm  = length(lon);
%    nn  = length(lat);
%    
%    %%%%Loading vector
%    %load LV data and reshape
%    LVname = strcat(filepath,LV_name);
%    raw_LV = fopen(LVname);raw_LV = fread(raw_LV,'float');
%    space = length(raw_LV)/mode/T;
%    LV = reshape(raw_LV,nn,mm,T,mode);
%    LV_tas(layer,:,:,:,:) = LV;
%end


%% ocean component

%thetao
T= LVnumber;
mode = cseofmodenum;
%LV_increment1 = zeros(1,length(lat2),length(lon2),T,mode);
LV_increment1 = zeros(1,length(lon2),length(lat2),T,mode);
for layer = 1 : 1
    %set filepath , variables
    filepath =[workdir, '/cseofs/'];
    LV_name = strcat(filepath,'LV_NorESM1-M_tas.dat');
%     LV_name = strcat(filepath,'eof_NorESM1-M_tas.dat');
%     LV_name = strcat('/data1/kimyy/etc/DASK_CSEOF/ygkim/cseof_analysis/exercise_kyy/cseofs/eof_increment1.dat');
    %load mask data
    mask_name = strcat(filepath,'ocean_mask1.txt');
    ocean_mask = importdata(mask_name);
    time = size(data,2) ; % time(month)
    T = LVnumber  ;  %period
%     mode = str2num(mode) ; %number of LV mode
    lon = lon2; %lon
    lat=  lat2; %lat
    mm = length(lon);
    nn = length(lat);
    
    %%%%Loading vector
    %load LV data and reshape
    LVname = strcat(filepath,LV_name)
    raw_LV = fopen(LV_name);
    raw_LV = fread(raw_LV,'float');
    space = length(raw_LV)/mode/T;
    raw_LV = reshape(raw_LV,space,T*mode);
    LV_grid = zeros(nn*mm,T*mode) ;
    LV_grid(ocean_mask,:) = raw_LV(:,:);
    LV = reshape(LV_grid,mm,nn,T,mode);
    LV_increment1(layer,:,:,:,:) = LV;
end

save([workdir, '/cseofs/lv_layer_merge.mat'], '-v7.3');
