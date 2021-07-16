clc;close all;clear all;
warning off;


%TEST
system_name=computer;
if (strcmp(system_name,'PCWIN64'))
    % % for windows
    dropboxpath='C:/Users/KYY/Dropbox'; %% SNU_desktop
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
%     addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/seo_nwp_1_10/run']));
elseif (strcmp(system_name,'GLNXA64'))    % % for linux
%     dropboxpath='/home01/kimyy/Dropbox'; %% ROMS server
    dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
%     addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/seo_nwp_1_10/run']));
end

run('nwp_polygon_point.m');


filepath = '/data2/kimyy/Observation/OISST/avhrr_only/www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/access/avhrr-only/';
year  = num2str(1982);
month = num2str(02,'%02i');
len_day = 1;

ii=1;
for aa = 1:len_day
%     aa=1;
% filepath = 'F:\ROMS\��м��ڷ�?(�������?)\�ڷᵿȭ\�����ڷ�(OSTIA)\';
name1 = strcat(year,month,'/','avhrr-only-v2.',year,month);
name2 = '.nc';
filename=strcat(filepath,name1,num2str(aa,'%2.2d'),name2);
% filename = strcat(filepath,filename);
lat = ncread(filename,'lat');
lon = ncread(filename,'lon');

sst_s = ncread(filename,'sst');
sst_e = ncread(filename,'err');
if (abs(mean(mean(sst_s,'omitnan'),'omitnan'))>100)
   sst_s=sst_s * 0.01 %% scale factor
   sst_e=sst_e * 0.01
end
% sst_er = ncread(filename,'analysis_error');
lat_s = ncread(filename,'lat');
lon_s = ncread(filename,'lon');

sst_s(sst_s<=-1000)=NaN;
sst_s(sst_s>=1000)=NaN;

lon1 = lon_s+360;
lon_s = [lon_s;lon1];
% sst = sst*scale_factor;

% lon_L = 98; lon_R = 180+(180-76); lat_H = 65; lat_L = -20;
lon_L = 115; lon_R = 162; lat_H = 52; lat_L = 15;

[px py] = find(lon_s >=lon_L & lon_s<=lon_R);
[qx qy] = find(lat_s >=lat_L & lat_s<=lat_H);
[xi, yi]=meshgrid(lon_s(px),lat_s(qx));
xi = double(xi);
yi = double(yi);
% SST = sst_s'-273.15;
SST = sst_s';
ERR = sst_e';

% sst_er = sst_er';
clear sst_s sst_e lon_s 
SST = [SST SST];
ERR = [ERR ERR];
% sst_er = [sst_er sst_er];
% pcolor(SST)

SST_A = SST(qx,px);  %% sampled SST(NWP) from global data
ERR_A = ERR(qx,px);

interval= 8; % 2 degree interval
xi_2d = xi(1:interval:end,1:interval:end);
yi_2d = yi(1:interval:end,1:interval:end);
SST_A_2d = SST_A(1:interval:end,1:interval:end);
ERR_A_2d = ERR_A(1:interval:end,1:interval:end);
interval=4; % 1 degree interval
xi_2d_5 = xi(1:interval:end,1:interval:end);
yi_2d_5 = yi(1:interval:end,1:interval:end);
SST_A_2d_5 = SST_A(1:interval:end,1:interval:end);
ERR_A_2d_5 = ERR_A(1:interval:end,1:interval:end);

% sst_er_a = sst_er_a(1:40:end,1:40:end);
% sst_er_a(find(sst_er_a>5))=nan;
% mean(nanmean(sst_er_a));
clear  SST ERR

% SST_A_2d = SST_A_2d(1:2:end,1:2:end);

xindex_10 = find(xi_2d(1,:)>117 & xi_2d(1,:)<160);
yindex_10 = find(yi_2d(:,1)>17 & yi_2d(:,1)<50);
% xindex_10 = find(xi_2d(1,:)<117 | xi_2d(1,:)>160);
% yindex_10 = find(yi_2d(:,1)<17 | yi_2d(:,1)>50);

% xindex_5 = find(xi_2d_5(1,:)>116 & xi_2d_5(1,:)<162);
% yindex_5 = find(yi_2d_5(:,1)>16 & yi_2d_5(:,1)<49);

mask_es = double(inpolygon(xi_2d_5,yi_2d_5,espolygon(:,1),espolygon(:,2)));
dim1_index_5 = find(mask_es==1);

x_10 = xi_2d(yindex_10,xindex_10);
y_10 = yi_2d(yindex_10,xindex_10);
% x_5 = xi_2d_5(yindex_5,xindex_5);
% y_5 = yi_2d_5(yindex_5,xindex_5);
x_5 = xi_2d_5(dim1_index_5);
y_5 = yi_2d_5(dim1_index_5);

SST_A_2d_10=SST_A_2d;
SST_A_2d(:,:) = nan;
SST_A_2d(yindex_10,xindex_10) =SST_A_2d_10(yindex_10,xindex_10);
% SST_5 = SST_A_2d_5(yindex_5,xindex_5);
SST_5 = SST_A_2d_5(dim1_index_5);
SST_5(y_5>50)=nan;
ERR_A_2d_10=ERR_A_2d;
ERR_A_2d(:,:) = nan;
ERR_A_2d(yindex_10,xindex_10) = ERR_A_2d_10(yindex_10,xindex_10);

% ERR_5 = ERR_A_2d_5(yindex_5,xindex_5);
ERR_5 = ERR_A_2d_5(dim1_index_5);
ERR_5(y_5>50)=nan;

index = 1;
for i = 1:length(SST_A_2d(:,1))
    for j = 1:length(SST_A_2d(1,:))
        
        if (isnan(SST_A_2d(i,j)))
            
        else
            lon_obs(index) = xi_2d(i,j);
            lat_obs(index) = yi_2d(i,j);
            sst_obs(index) = SST_A_2d(i,j);
            obserr(index) = ERR_A_2d(i,j);
            index = index + 1;
        end
    end
end
for i = 1:length(SST_5(:,1))
    for j = 1:length(SST_5(1,:))
        
        if (isnan(SST_5(i,j)))
            
        else
            lon_obs(index) = x_5(i,j);
            lat_obs(index) = y_5(i,j);
            sst_obs(index) = SST_5(i,j);
            obserr(index) = ERR_5(i,j);
            index = index + 1;
        end
    end
end

% % % % 2015 (added point)
% % % nn =[397 398];
% % % lon_obs(nn) = [];
% % % lat_obs(nn) = [];
% % % sst_obs(nn) = [];
% % % 
% % % 
% % % % sst = sst*scale_factor;
% % % % south boundary
% % % nn = [1:32];
% % % % west boundary
% % % nn = [nn 33 68 104 139 174 209];
% % % % land
% % % nn = [nn 244 245 277 278 301:303 322:326 345 346 364 431 440:442];
% % % % land 2
% % % nn = [nn 430 438 439];
% % % % nwp
% % % nn = [nn 529 542 562 575 584 588];
% % % 
% % % lon_obs(nn) = [];
% % % lat_obs(nn) = [];
% % % sst_obs(nn) = [];

% 
% im = [72 137 141 149 109 110 113];
% ix = [301 581 621 1001 1101 1401 2301];
% iy = [701 1081 1201 1321 1501 1501 1501];
% 
%     lon_obs(im) = xi(1,ix);
%     lat_obs(im) = yi(iy,1);
%  for mm = 1:7
%     sst_obs(im(mm)) = SST_A(iy(mm),ix(mm));
%  end
%  
  
% im = 1:16;
% ix1 = 1:200:401;
% ix2 = 1201:200:3601;
% ix = [ix1 ix2];
% iy = 21*ones(1,16);
% 
%     lon_obs(im) = xi(1,ix);
%     lat_obs(im) = yi(iy,1);
%  for mm = 1:7
%     sst_obs(im(mm)) = SST_A(iy(mm),ix(mm));
%  end
%  
%   nn = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 35 54];
% lon_obs(nn) = [];
% lat_obs(nn) = [];
% sst_obs(nn) = [];

% figure
% m_proj('mercator','lon',[98 284],'lat',[-20 65]);
% m_grid('fontsize',15)
% hold on
% m_plot(lon_obs,lat_obs,'ro','Markersize',5)
% m_gshhs_i('color','k')  



num_data = length(lon_obs);
outname = strcat(filepath,'obs_',num2str(aa,'%2.2d'),'.nc')
nccreate(outname,'ixt','Dimensions',{'xt_i',num_data});
nccreate(outname,'rlon','Dimensions',{'xt_i',num_data,'time',1});
nccreate(outname,'rlat','Dimensions',{'xt_i',num_data,'time',1});
nccreate(outname,'rdepth','Dimensions',{'xt_i',num_data,'time',1});
nccreate(outname,'obsdata','Dimensions',{'xt_i',num_data,'time',1});
nccreate(outname,'obserr','Dimensions',{'xt_i',num_data,'time',1});
nccreate(outname,'dindex','Dimensions',{'xt_i',num_data,'time',1});
nccreate(outname,'ndata','Dimensions',{'time',1});
ixt = 1:num_data;
dindex = 2*ones(1,num_data)';
% obserr = 0.67*ones(1,num_data);
ncwrite(outname,'ixt',ixt);
ncwrite(outname,'dindex',dindex);
ncwrite(outname,'ndata',num_data);
ncwrite(outname,'rdepth',zeros(num_data,1));
ncwrite(outname,'rlon',lon_obs');
ncwrite(outname,'rlat',lat_obs');
ncwrite(outname,'obsdata',sst_obs');
ncwrite(outname,'obserr',obserr');

ii=ii+1;


end
% 

% 
figure
m_proj('mercator','lon',[115 162],'lat',[15 52]);
m_grid('fontsize',15)
hold on
m_plot(lon_obs,lat_obs,'ro','Markersize',5)
m_gshhs_i('color','k')  
% sti = 397;
% m_plot(lon_obs(sti),lat_obs(sti),'bo','Markersize',5)