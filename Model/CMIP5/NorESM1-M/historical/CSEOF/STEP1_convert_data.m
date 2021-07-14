clc;clear all;close all;
% This is the code to convert acsii file from binary file, 
%fit to cseof analysis
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

%% -----case 1. ATM file (no mask)----------%

%-----------STEP1. set the variables----------------%
%input file path,name
%%filepath = '/home/kimyy/ygkim/cseof_analysis/exercise_kyy/data/';
%%filename = 'increment_v1.5_050.nc';
%%file = strcat(filepath,filename);

%open nc file and set the study area
%%nc = netcdf(file) ;
%%lat = nc{'YT_OCEAN'}(:);lon = nc{'XT_OCEAN'}(:);time = nc{'time'}(:);
%%lat_ini = find(lat<=60, 1, 'last' );lat_end = find(lat>=-60, 1 );
%%lon_ini = find(lon<-50, 1, 'last' ); lon_end = find(lon>=-279, 1 );

%load land_mask file 
%%psl= nc{'psl'}(:,lat_ini:lat_end,lon_ini:lon_end);

%--------STEP2. convert acsii file -------------------%
%extract psl data
%%data =zeros(length(lat2)*length(lon2),length(time));
%%for i = 1:length(time) 
%%    data(:,i) = psl(i,:)';
%%end

%save data
%%save('psl.data','data','-ASCII');

%% ---------case 2. ocean file(ocean mask)-------------%
%example of temperature (40 layer)
%-----------STEP1. set the variables----------------%
%input file path,name
%path in the damo server
tfreq = 1;  %% 1 : monthly, 2 : daily

workdir='/data1/kimyy/etc/NorESM1-M_CSEOF';
varname = ['tas']  % tas, psl, hur, rsds, uo, vo
modelname = 'NorESM1-M';
scenname = 'historical_r1i1p1';
inputyear = 1976:2005;
dl = 0.5;
section = [115 164 15 52];

for nyear = 1:length(inputyear)
  tempyear = inputyear(nyear);
  filepath = [workdir, '/data/NorESM1-M_mon/']
  
%   filename = ['NorESM-M_', varname, '_', num2str(tempyear,'%04i'), '_v2.nc']
  filename = [varname, '_mon_', modelname, '_', scenname, '_', num2str(tempyear,'%04i'), '.nc']

  file = strcat(filepath,filename)
  
  %open nc file and set the study area
  %nc = netcdf.open(file)
  %lat = nc{'YT_OCEAN'}(:)
  %lon = nc{'XT_OCEAN'}(:);
  %lev = nc{'ST_OCEAN'}(:);
  if (nyear == 1)
    lon_glo = ncread(file,'lon');
    lat_glo = ncread(file,'lat');
%    lev = ncread(file,'S');
    xlen_glo=length(lon_glo);
    ylen_glo=length(lat_glo);
%    zlen=length(lev)
    [indw, inde, inds, indn] = findind_Y(dl, section, lon_glo, lat_glo);
    
    lon = ncread(file,'lon',indw, inde-indw+1);
    lat = ncread(file,'lat',inds, indn-inds+1);
    xlen = length(lon);
    ylen = length(lat);
  end 
  % first time is 16-JAN-1979, monthly data.
  %time = nc{'TIME1'}(49:408);
  %time = ncread(file, 'TIME1')
  
  if (tfreq==1)
      tstart=(nyear-1)*12+1;
      tend=nyear*12;
      temptime = ncread(file,'time',1,12);
      ttemp = ncread(file,varname,[indw inds 1], [xlen ylen 12]);
  elseif (tfreq==2)
      tstart=(nyear-1)*365+1;
      tend=nyear*365;
      temptime = ncread(file,'time',1,365);
      ttemp = ncread(file,varname,[indw inds 1], [xlen ylen 365]);
  end
  
%   temptime = ncread(file,'time');
  time(tstart:tend)=temptime;
  tlen=length(time)
  
  temp(:,:,tstart:tend)=ttemp;
end

%--------STEP2. convert acsii file -------------------%
%extract temperature data
    
%    temp=permute(temp, [3 2 1]); %%[t y x]
    %temp = ncread(file,'ODA_GOT_T',[1,1,1,49],[xlen,ylen,zlen,tlen]); 
    
% % %     for ocean mask
    size(temp)
for j = 1 : 1  
    %setting land mask and save mask file
%     temp2 = temp(:,:,1);
%     temp2(:,1,:,:) = temp(:,:,:);
%     size(temp2)
    n=length(lon);
    q=length(lat);
    m=length(time);
%     get salt data for land mask
    saltname = [workdir, '/data/NorESM1-M_mon/so_interp_NorESM1-M_historical_r1i1p1_1976.nc'];
    salt = ncread(saltname,'so');
    lon_salt =ncread(saltname, 'lon');
    lat_salt =ncread(saltname, 'lat');
    salt_interped=griddata(double(lon_salt),double(lat_salt'), double(squeeze(salt(:,:,1,1))'), double(lon), double(lat'))';
    land_idx = zeros(n*q,1);
    idx = find(isnan(salt_interped)==1) ; %land mask: temp==0
    land_idx(idx) = 1 ;
    land_mask = find(land_idx ==1);
    ocean_mask = find(land_idx ==0);
    save(strcat(workdir, '/data/NorESM1-M_mon/ocean_mask',num2str(j),'.txt'),'ocean_mask','-ascii');

    %extract temp data(NaN data remove)
    data =zeros(n*q- length(land_mask),m);
    for i = 1:m 
        temp2=squeeze(temp(:,:,i));
	%temp2(i,j,:,:)=permute(temp(:,:,j,i), [4 3 2 1]); %temp : lon, lat , depth, time
        temp2(land_mask) = [];
        data(:,i) = temp2(:)' ; % [space time]
%	size(temp2)
    end
    
    sizedata=size(data); 
    %% for fast calculation speed, if [M*N], M must be higher than N.
    if (sizedata(1) >= sizedata(2)) 
        data2 = data;
    else
        data2 = data';
    end
    data_len=n*q-length(land_mask)
    size(data2)
    slen = size(data2,1);
    %save data
    tt = strcat(workdir, '/data/','NorESM1-M_',varname,num2str(j),'.data');
    save(tt,'data2','-ASCII');
    disp(['dimension of sampling stations : ',num2str(slen)]);
    disp(['number of sampling points at each station : ',num2str(tlen)]);
end

save('STEP1_output.mat')