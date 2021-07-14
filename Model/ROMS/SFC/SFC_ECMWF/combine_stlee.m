clear all;close all; clc;

warning off;

linux=1; windows=0;

if (windows ==1)
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox';
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
%     addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_ktotalday']));
elseif (linux==1)
    dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_ktotalday']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Roms_tools/Run']));
end
start



grdfile='/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_20/input/roms_grid_combine2_test37.nc';

year_start =2002;
year_end = 2002;
rawdata_dir='/data1/temp/ECMWF_interim/';
output_dir='/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_20/forcing_matlab/';

%% make D
for year= year_start : year_end
    tic;
    var_dir='dewt/'
    rawdata=[rawdata_dir,var_dir,'ECMWF_interim_dewt_',num2str(year),'.nc'];
    nc1=netcdf(rawdata);
    totalday=year2day(year);

    longitude_raw=nc1{'longitude'}(:);
    latitude_raw=nc1{'latitude'}(:);
    time=[0.5:1:totalday-0.5];
      
    Dair=nc1{'d2m'}(:);

    [longitude,latitude]=meshgrid(longitude_raw,latitude_raw);
    section=[115 164 15 52];
    dist=sqrt((longitude-section(1)).^2+(latitude-section(3)).^2); %% get distance from station 1
    min_dist=min(min(dist));
    dist2=sqrt((longitude-section(2)).^2+(latitude-section(4)).^2);  %% get distance from station 2
    min_dist2=min(min(dist2));                
    [x1,y1]=find(dist==min_dist);  %% find closest point from station 1. [lat lon]
    [x2,y2]=find(dist2==min_dist2); %% find closest point from station 2. [lat lon]

    cut_data=Dair(:,x2(1):x1(1),y1(1):y2(1));
    cut_lon=longitude(x2(1):x1(1),y1(1):y2(1));
    cut_lat=latitude(x2(1):x1(1),y1(1):y2(1));
    
    Dsfactor=nc1{'d2m'}.scale_factor(:);
    Daoffset=nc1{'d2m'}.add_offset(:);
    size_D=size(cut_data);
    close(nc1);
    cut_data=cut_data.*Dsfactor+Daoffset-273.15;


    D_daily=ones(totalday,size_D(2),size_D(3))*NaN;

    for i=1:1:totalday;
    %     aa=sum(Dair(1+2*(i-1):2*i,:,:));
        aa=sum(cut_data(1+4*(i-1):4*i,:,:));
    %     bb=sum(Dair((totalday*4+1)+4*(i-1):(totalday*4)+4*i,:,:));
        D_daily(i,:,:)=(aa)/4;
    end

    clear Dair nc1;


    grd=netcdf(grdfile);
    lat_rho=grd{'lat_rho'}(:);
    lon_rho=grd{'lon_rho'}(:);
    close(grd);

    size_inp_D=size(lat_rho);
    inp_D_daily=ones(totalday,size_inp_D(1),size_inp_D(2))*NaN;
%     [Xi,Yi]=meshgrid(cut_lon,cut_lat);

    for ii=1:1:totalday;
        Zi=squeeze(D_daily(ii,:,:));
        Z=griddata(cut_lon,cut_lat,Zi,lon_rho,lat_rho);
        inp_D_daily(ii,:,:)=Z;
    end

    clear D_daily grd Z lat_u lon_u Xi Yi Zi

    fname=[output_dir,'nwp_1_20_',num2str(year),'_Dair_',num2str(totalday),'.nc'];
    toc;
    % fname=['easttest_',num2str(year),'Dair_',num2str(totalday),'.nc'];
    create_roms_forcing_D(fname,inp_D_daily,time,totalday,totalday);
    % clear all;close all;
%     version -date
end


%% make P
close all;clear all;

for year=2012
% grdfile='D:\add2_ini_bry_grd\grid\roms_grid2_ADD_08_2_ep.nc';
% grdfile='d:\NWP_2\grid\roms_grid_combine.nc';
% grdfile='f:\Input\grid_gy_v11_s.nc';
% grdfile='g:\ideal\domain\ideal3_sym_linear_grid.nc';
grdfile='D:\MEPL\project\NWP\make_init\grid\roms_grid_combine2_edited_nwp_ktotalday.nc';

rawdata=['ECMWF_Interim_msl_',num2str(year),'.nc'];
nc1=netcdf(rawdata);
totalday=year2day(year);

longitude=nc1{'longitude'}(:);
latitude=nc1{'latitude'}(:);
time=[0.5:1:totalday-0.5];
Pair=nc1{'msl'}(:);
Psfactor=nc1{'msl'}.scale_factor(:);
Paoffset=nc1{'msl'}.add_offset(:);
size_P=size(Pair);
close(nc1);
Pair=(Pair.*Psfactor+Paoffset)./100;


P_daily=ones(totalday,size_P(2),size_P(3))*NaN;

for i=1:1:totalday;
%     aa=sum(Pair(1+2*(i-1):2*i,:,:));
    aa=sum(Pair(1+4*(i-1):4*i,:,:));
%     bb=sum(Pair((totalday*4+1)+4*(i-1):(totalday*4)+4*i,:,:));
    P_daily(i,:,:)=(aa)/4;
end

clear Pair nc1;



grd=netcdf(grdfile);
lat_rho=grd{'lat_rho'}(:);
lon_rho=grd{'lon_rho'}(:);
close(grd);

size_inp_P=size(lat_rho);
inp_P_daily=ones(totalday,size_inp_P(1),size_inp_P(2))*NaN;
[Xi,Yi]=meshgrid(longitude,latitude);

for ii=1:1:totalday;
    Zi=squeeze(P_daily(ii,:,:));
    Z=griddata(Xi,Yi,Zi,lon_rho,lat_rho);
    inp_P_daily(ii,:,:)=Z;
end

clear P_daily grd Z lat_u lon_ Xi Yi Zi

% inp_P_daily_2007 = inp_P_daily;
% 
% save 2007-Pair.mat inp_P_daily_2007

fname=['nwp_1_20_',num2str(year),'_Pair_',num2str(totalday),'.nc'];

% fname=['easttest_',num2str(year),'Pair_',num2str(totalday),'.nc'];
create_roms_forcing_P(fname,inp_P_daily,time,totalday,totalday)
clear all;close all;
end

%% make swrad
close all;clear all;

for year=2012
% grdfile='D:\add2_ini_bry_grd\grid\roms_grid2_ADD_08_2_ep.nc';
% % grdfile='D:\add_ini_bry_grd\grid\roms_grid_ADD_10_ep.nc';
% % grdfile='D:\add2_ini_bry_grd\grid\east_grid.nc';
% grdfile='d:\NWP_2\grid\roms_grid_combine.nc';
% grdfile='f:\Input\grid_gy_v11_s.nc';
% grdfile='g:\ideal\domain\ideal3_sym_linear_grid.nc';
grdfile='D:\MEPL\project\NWP\make_init\grid\roms_grid_combine2_edited_nwp_ktotalday.nc';

rawdata=['ECMWF_Interim_ssrd_',num2str(year),'.nc'];
nc1=netcdf(rawdata);
totalday=year2day(year);

longitude=nc1{'longitude'}(:);
latitude=nc1{'latitude'}(:);
time=[0.5:1:totalday-0.5];
swrad=nc1{'ssrd'}(:);
SWsfactor=nc1{'ssrd'}.scale_factor(:);
SWaoffset=nc1{'ssrd'}.add_offset(:);
size_SW=size(swrad);
close(nc1);
swrad=(swrad.*SWsfactor+SWaoffset);


SW_daily=ones(totalday,size_SW(2),size_SW(3))*NaN;
% i_c=1;
% for i=1:4:totalday*4;
%     SW_three(1+(i_c-1)*8,:,:)=swrad(i,:,:)/10800/4;  % 12�� ���� 00�ñ��� ���� -00��
%     SW_three(3+(i_c-1)*8,:,:)=swrad(i+1,:,:)/10800/2; % ���ʿ�
%     SW_three(5+(i_c-1)*8,:,:)=swrad(i+2,:,:)/10800/4; % 00�� ���� 12�ñ��� ���� -12��
%     SW_three(7+(i_c-1)*8,:,:)=swrad(i+3,:,:)/10800/2; % ���ʿ�
%     i_c=i_c+1;
% end
% i_c=1;
% for i=1+totalday*4:4:totalday*8;
%     SW_three(2+(i_c-1)*8,:,:)=swrad(i,:,:)/10800/1; % ���ʿ�
%     SW_three(4+(i_c-1)*8,:,:)=swrad(i+1,:,:)/10800/3; % ���ʿ�
%     SW_three(6+(i_c-1)*8,:,:)=swrad(i+2,:,:)/10800/1; % ���ʿ�
%     SW_three(8+(i_c-1)*8,:,:)=swrad(i+3,:,:)/10800/3; % ���ʿ�
%     i_c=i_c+1;
% end


for i=1:1:totalday;
%     aa=sum(swrad(1+4*(i-1):4*i,:,:));
%     bb=sum(swrad((totalday*4+1)+4*(i-1):(totalday*4)+4*i,:,:));
%     SW_daily(i,:,:)=(aa+bb)/8;
%     SW_daily(i,:,:)=(SW_three(1+(i-1)*8,:,:)+SW_three(5+(i-1)*8,:,:))/2;
%     SW_daily(i,:,:)=(swrad(i*2-1,:,:)+swrad(i*2,:,:))/(24*60*60);
    SW_daily(i,:,:)=(swrad(4+(i-1)*8,:,:)+swrad(8+(i-1)*8,:,:))/(24*60*60);
end

clear swrad nc1;



grd=netcdf(grdfile);
lat_rho=grd{'lat_rho'}(:);
lon_rho=grd{'lon_rho'}(:);
close(grd);

size_inp_SW=size(lat_rho);
inp_SW_daily=ones(totalday,size_inp_SW(1),size_inp_SW(2))*NaN;
[Xi,Yi]=meshgrid(longitude,latitude);

for ii=1:1:totalday;
    Zi=squeeze(SW_daily(ii,:,:));
    Z=griddata(Xi,Yi,Zi,lon_rho,lat_rho);
    inp_SW_daily(ii,:,:)=Z;
end

clear SW_daily grd Z lat_u lon_u Xi Yi Zi

% inp_SW_daily_2007 = inp_SW_daily;
% 
% save 2007-SWair.mat inp_SW_daily_2007

fname=['nwp_1_20_',num2str(year),'_SWair_',num2str(totalday),'.nc'];

% fname=['easttest_',num2str(year),'swradd_',num2str(totalday),'.nc'];
create_roms_forcing_SWD(fname,inp_SW_daily,time,totalday,totalday)
clear all;close all;
end

%% make T
close all;clear all;
for year=2012
grdfile='D:\MEPL\project\NWP\make_init\grid\roms_grid_combine2_edited_nwp_ktotalday.nc';

rawdata=['ECMWF_Interim_airT_',num2str(year),'.nc'];
nc1=netcdf(rawdata);
totalday=year2day(year);

longitude=nc1{'longitude'}(:);
latitude=nc1{'latitude'}(:);
tt=nc1{'time'}(:);
% ssrd=nc1{'ssrd'}(:);
time=[0.5:1:totalday-0.5];
Tair=nc1{'t2m'}(:);
Tsfactor=nc1{'t2m'}.scale_factor(:);
Taoffset=nc1{'t2m'}.add_offset(:);
size_T=size(Tair);
close(nc1);
Tair=Tair.*Tsfactor+Taoffset-273.15;


T_daily=ones(totalday,size_T(2),size_T(3))*NaN;

for i=1:1:totalday;
%     aa=sum(Tair(1+2*(i-1):2*i,:,:));
        aa=sum(Tair(1+4*(i-1):4*i,:,:));
%     bb=sum(Tair((totalday*4+1)+4*(i-1):(totalday*4)+4*i,:,:));
%     T_daily(i,:,:)=(aa+bb)/8;
    T_daily(i,:,:)=(aa)/4;
end

clear Tair nc1;

grd=netcdf(grdfile);
lat_rho=grd{'lat_rho'}(:);
lon_rho=grd{'lon_rho'}(:);
close(grd);

size_inp_T=size(lat_rho);
inp_T_daily=ones(totalday,size_inp_T(1),size_inp_T(2))*NaN;
[Xi,Yi]=meshgrid(longitude,latitude);

for ii=1:1:totalday;
    Zi=squeeze(T_daily(ii,:,:));
    Z=griddata(Xi,Yi,Zi,lon_rho,lat_rho);
    inp_T_daily(ii,:,:)=Z;
end

clear T_daily grd Z lat_u lon_u Xi Yi Zi

% inp_T_daily_2010 = inp_T_daily;
% 
% save 2010-Tair.mat inp_T_daily_2010

fname=['nwp_1_20_',num2str(year),'_Tair_',num2str(totalday),'.nc'];
% fname=['easttest_',num2str(year),'Tair_',num2str(totalday),'.nc'];
create_roms_forcing_T(fname,inp_T_daily,time,totalday,totalday)
clear all;close all;
end

%% make u wind
close all;clear all;

for year=2012
grdfile='D:\MEPL\project\NWP\make_init\grid\roms_grid_combine2_edited_nwp_ktotalday.nc';

rawdata=['ECMWF_Interim_u10_',num2str(year),'.nc'];
nc1=netcdf(rawdata);
totalday=year2day(year);
time=[0.5:1:totalday-0.5];
longitude=nc1{'longitude'}(:);
latitude=nc1{'latitude'}(:);
U_wind=nc1{'u10'}(:);
Usfactor=nc1{'u10'}.scale_factor(:);
Uaoffset=nc1{'u10'}.add_offset(:);
size_U=size(U_wind);
close(nc1);
U_wind=U_wind.*Usfactor+Uaoffset;


U_daily=ones(totalday,size_U(2),size_U(3))*NaN;

for i=1:1:totalday;
%     aa=sum(U_wind(1+2*(i-1):2*i,:,:));
    aa=sum(U_wind(1+4*(i-1):4*i,:,:));
%     bb=sum(U_wind((totalday*4+1)+4*(i-1):(totalday*4)+4*i,:,:));
%     U_daily(i,:,:)=(aa+bb)/8;
    U_daily(i,:,:)=aa/4;
end

clear U_wind nc1;



grd=netcdf(grdfile);
lat_rho=grd{'lat_rho'}(:);
lon_rho=grd{'lon_rho'}(:);
close(grd);

size_inp_U=size(lat_rho);
inp_U_daily=ones(totalday,size_inp_U(1),size_inp_U(2))*NaN;
[Xi,Yi]=meshgrid(longitude,latitude);

for ii=1:1:totalday;
    Zi=squeeze(U_daily(ii,:,:));
    Z=griddata(Xi,Yi,Zi,lon_rho,lat_rho);
    inp_U_daily(ii,:,:)=Z;
end

clear U_daily grd Z lat_u lon_u Xi Yi Zi

% inp_U_daily_2011 = inp_U_daily;
% 
% save 2011-Uair.mat inp_U_daily_2011

fname=['nwp_1_20_',num2str(year),'_Uwind_',num2str(totalday),'.nc'];
% fname=['easttest_',num2str(year),'Tair_',num2str(totalday),'.nc'];
create_roms_forcing_U(fname,inp_U_daily,time,totalday,totalday)
clear all;close all
end
%% make V wind
for year=2012
    
grdfile='D:\MEPL\project\NWP\make_init\grid\roms_grid_combine2_edited_nwp_ktotalday.nc';

rawdata=['ECMWF_Interim_v10_',num2str(year),'.nc'];
nc1=netcdf(rawdata);
totalday=year2day(year);
longitude=nc1{'longitude'}(:);
latitude=nc1{'latitude'}(:);
time=[0.5:1:totalday-0.5];
V_wind=nc1{'v10'}(:);
Vsfactor=nc1{'v10'}.scale_factor(:);
Vaoffset=nc1{'v10'}.add_offset(:);
size_V=size(V_wind);
close(nc1);
V_wind=V_wind.*Vsfactor+Vaoffset;



V_daily=ones(totalday,size_V(2),size_V(3))*NaN;

for i=1:1:totalday;
%     aa=sum(V_wind(1+2*(i-1):2*i,:,:));
     aa=sum(V_wind(1+4*(i-1):4*i,:,:));
%     bb=sum(V_wind((totalday*4+1)+4*(i-1):(totalday*4)+4*i,:,:));
%     V_daily(i,:,:)=(aa+bb)/8;
    V_daily(i,:,:)=aa/4;
end

clear V_wind nc1;



grd=netcdf(grdfile);
lat_rho=grd{'lat_rho'}(:);
lon_rho=grd{'lon_rho'}(:);
close(grd);

size_inp_V=size(lat_rho);
inp_V_daily=ones(totalday,size_inp_V(1),size_inp_V(2))*NaN;
[Xi,Yi]=meshgrid(longitude,latitude);

for ii=1:1:totalday;
    Zi=squeeze(V_daily(ii,:,:));
    Z=griddata(Xi,Yi,Zi,lon_rho,lat_rho);
    inp_V_daily(ii,:,:)=Z;
end

clear V_daily grd Z lat_v lon_v Xi Yi Zi

% inp_V_daily_2009 = inp_V_daily;
% 
% save 2009-Vair.mat inp_V_daily_2009

fname=['nwp_1_20_',num2str(year),'_Vwind_',num2str(totalday),'.nc'];
% fname=['easttest_',num2str(year),'Vwind_',num2str(totalday),'.nc'];
create_roms_forcing_V(fname,inp_V_daily,time,totalday,totalday)
clear all;close all;
end
%% make Q
close all;clear all;

for year=2012
    clear t_time
    P='Pair';T='Tair';D='Dair';Q='Qair';foot='.nc';
    
    totalday=year2day(year);strtotalday=num2str(totalday);
%     rawdata=['frc_ecmwf_ideal3_',num2str(year),P,'_',strtotalday,foot];
    rawdata=['nwp_1_20_',num2str(year),'_',P,'_',strtotalday,foot];
%     rawdata=['easttest_',num2str(year),P,'_',strtotalday,foot];
    Pnc=netcdf(rawdata);
%     rawdata=['frc_ecmwf_ideal3_',num2str(year),T,'_',strtotalday,foot];
    rawdata=['nwp_1_20_',num2str(year),'_',T,'_',strtotalday,foot];
%     rawdata=['easttest_',num2str(year),T,'_',strtotalday,foot];
    Tnc=netcdf(rawdata);
%     rawdata=['frc_ecmwf_ideal3_',num2str(year),D,'_',strtotalday,foot];
    rawdata=['nwp_1_20_',num2str(year),'_',D,'_',strtotalday,foot];
%     rawdata=['easttest_',num2str(year),D,'_',strtotalday,foot];
    Dnc=netcdf(rawdata);

    lon=Pnc{'lon_rho'}(:);
    lat=Pnc{'lat_rho'}(:);
    t_time = Pnc{'Pair_time'}(:);
    ttt = length(t_time) ;
        for dd=1:1:ttt;
    
            P_value = Pnc{P}(dd,:,:);
            T_value = Tnc{T}(dd,:,:);
            D_value = Dnc{D}(dd,:,:);

            Qair(dd,:,:) = (qsat_yg(D_value,P_value)./qsat_yg(T_value,P_value)).*100;
            
        end
    close(Pnc);  close(Tnc);  close(Dnc);
    
%     inp_Q_daily_2007 = Qair;
%     
%     save 2007-Qair.mat inp_Q_daily_2007

    fname=['nwp_1_20_',num2str(year),'_Qair_',num2str(totalday),'.nc'];
    
%     fname=['easttest_',num2str(year),'Qair_',strtotalday,foot];
    create_roms_forcing_Qair(fname,Qair,t_time,totalday,totalday)
    clear all;close all;
end

%% combine
% clear all;close all;

