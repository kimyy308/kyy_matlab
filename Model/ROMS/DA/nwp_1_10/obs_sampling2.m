clc;close all;clear all;
warning off;

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
% year  = num2str(1982);
% month = num2str(01,'%02i');
% len_day = 13511;   %% 12412 =34y  %% 13511--> to 2018
len_day = 13511;

% % datenum(1982,1,1)

ii=1;
emptyval=0;
for aa = 12412:7:len_day
%     aa=1;
    initdaynum=datenum(1981,12,31);
    nowdate=datestr(initdaynum+aa,'yyyymmdd')
    year=nowdate(1:4);
    month=nowdate(5:6);
    day=nowdate(7:8);
    % filepath = 'F:\ROMS\��м��ڷ�(�������)\�ڷᵿȭ\�����ڷ�(OSTIA)\';
    name1 = strcat(year,month,'/','avhrr-only-v2.',year,month);
    name2 = '.nc';
    filename=strcat(filepath,name1,num2str(day,'%2.2d'),name2);
    % filename = strcat(filepath,filename);
    lat = ncread(filename,'lat');
    lon = ncread(filename,'lon');

    sst_s = ncread(filename,'sst');
    sst_e = ncread(filename,'err');
    if (abs(mean(mean(sst_s,'omitnan'),'omitnan'))>100)
       sst_s=sst_s * 0.01; %% scale factor
       sst_e=sst_e * 0.01;
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
    xi_2d_es = xi(1:interval:end,1:interval:end);
    yi_2d_es = yi(1:interval:end,1:interval:end);
    SST_A_2d_es = SST_A(1:interval:end,1:interval:end);
    ERR_A_2d_es = ERR_A(1:interval:end,1:interval:end);

    xi_2d_kuro = xi(1:interval:end,1:interval:end);
    yi_2d_kuro = yi(1:interval:end,1:interval:end);
    SST_A_2d_kuro = SST_A(1:interval:end,1:interval:end);
    ERR_A_2d_kuro = ERR_A(1:interval:end,1:interval:end);
    
    interval=2; % 0.5 degree interval
    xi_2d_ekwc = xi(1:interval:end,1:interval:end);
    yi_2d_ekwc = yi(1:interval:end,1:interval:end);
    SST_A_2d_ekwc = SST_A(1:interval:end,1:interval:end);
    ERR_A_2d_ekwc = ERR_A(1:interval:end,1:interval:end);



    % sst_er_a = sst_er_a(1:40:end,1:40:end);
    % sst_er_a(find(sst_er_a>5))=nan;
    % mean(nanmean(sst_er_a));
    clear  SST ERR

    % SST_A_2d = SST_A_2d(1:2:end,1:2:end);

    xindex_10 = find(xi_2d(1,:)>117 & xi_2d(1,:)<160);
    yindex_10 = find(yi_2d(:,1)>17 & yi_2d(:,1)<50);
    % xindex_10 = find(xi_2d(1,:)<117 | xi_2d(1,:)>160);
    % yindex_10 = find(yi_2d(:,1)<17 | yi_2d(:,1)>50);

    % xindex_es = find(xi_2d_es(1,:)>116 & xi_2d_es(1,:)<162);
    % yindex_es = find(yi_2d_es(:,1)>16 & yi_2d_es(:,1)<49);

    mask_es = double(inpolygon(xi_2d_es,yi_2d_es,espolygon(:,1),espolygon(:,2)));
    dim1_index_es = find(mask_es==1);

    mask_kuro = double(inpolygon(xi_2d_kuro,yi_2d_kuro,kuropolygon(:,1),kuropolygon(:,2)));
    dim1_index_kuro = find(mask_kuro==1);
    
    mask_ekwc = double(inpolygon(xi_2d_ekwc,yi_2d_ekwc,ekwcpolygon(:,1),ekwcpolygon(:,2)));
    dim1_index_ekwc = find(mask_ekwc==1);


    x_10 = xi_2d(yindex_10,xindex_10);
    y_10 = yi_2d(yindex_10,xindex_10);
    % x_es = xi_2d_es(yindex_es,xindex_es);
    % y_es = yi_2d_es(yindex_es,xindex_es);
    x_es = xi_2d_es(dim1_index_es);
    y_es = yi_2d_es(dim1_index_es);
    x_kuro = xi_2d_kuro(dim1_index_kuro);
    y_kuro = yi_2d_kuro(dim1_index_kuro);
    x_ekwc = xi_2d_ekwc(dim1_index_ekwc);
    y_ekwc = yi_2d_ekwc(dim1_index_ekwc);



    SST_A_2d_10=SST_A_2d;
    SST_A_2d(:,:) = nan;
    SST_A_2d(yindex_10,xindex_10) =SST_A_2d_10(yindex_10,xindex_10);
    ERR_A_2d_10=ERR_A_2d;
    ERR_A_2d(:,:) = nan;
    ERR_A_2d(yindex_10,xindex_10) = ERR_A_2d_10(yindex_10,xindex_10);

    SST_es = SST_A_2d_es(dim1_index_es);
    SST_es(y_es>50)=nan;
    ERR_es = ERR_A_2d_es(dim1_index_es);
    ERR_es(y_es>50)=nan;

    SST_kuro = SST_A_2d_kuro(dim1_index_kuro);
    SST_kuro(y_kuro>50)=nan;
    ERR_kuro = ERR_A_2d_kuro(dim1_index_kuro);
    ERR_kuro(y_kuro>50)=nan;
    
    SST_ekwc = SST_A_2d_ekwc(dim1_index_ekwc);
    SST_ekwc(y_ekwc>50)=nan;
    ERR_ekwc = ERR_A_2d_ekwc(dim1_index_ekwc);
    ERR_ekwc(y_ekwc>50)=nan;


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
    for i = 1:length(SST_es(:,1))
        for j = 1:length(SST_es(1,:))

            if (isnan(SST_es(i,j)))

            else
                lon_obs(index) = x_es(i,j);
                lat_obs(index) = y_es(i,j);
                sst_obs(index) = SST_es(i,j);
                obserr(index) = ERR_es(i,j);
                index = index + 1;
            end
        end
    end
    for i = 1:length(SST_kuro(:,1))
        for j = 1:length(SST_kuro(1,:))

            if (isnan(SST_kuro(i,j)))

            else
                lon_obs(index) = x_kuro(i,j);
                lat_obs(index) = y_kuro(i,j);
                sst_obs(index) = SST_kuro(i,j);
                obserr(index) = ERR_kuro(i,j);
                index = index + 1;
            end
        end
    end
    for i = 1:length(SST_ekwc(:,1))
        for j = 1:length(SST_ekwc(1,:))

            if (isnan(SST_ekwc(i,j)))

            else
                lon_obs(index) = x_ekwc(i,j);
                lat_obs(index) = y_ekwc(i,j);
                sst_obs(index) = SST_ekwc(i,j);
                obserr(index) = ERR_ekwc(i,j);
                index = index + 1;
            end
        end
    end







%     if(emptyval~=1)
        % % % % 2015 (added point)
        wrong_point = [119.125, 25.125 ; 
                        121.125, 25.125 ; 
                        131.125, 31.125 ;  
                        129.125, 33.125 ;  
                        133.125, 33.125 ;
                        135.125, 34.125 ;
                        129.125, 35.125 ;
                        139.125, 35.125 ;
                        140.125, 35.125 ;
                        136.125, 36.125 ;
                        137.125, 37.125 ;
                        138.125, 37.125 ;
                        129.125, 37.625 ;
                        128.625, 38.125 ;
                        127.625, 39.125 ;
                        125.125, 39.125 ;
                        127.625, 39.625 ;
                        141.125, 41.125 ;
                        140.125, 42.125 ;
                        132.125, 43.125 ;
                        144.125, 44.125 ;
                        145.125, 44.125 ;
                        148.125, 45.125 ;
                        138.125, 46.125 ;
                        143.125, 47.125 ;
                        142.125, 48.125 ;
                        143.125, 49.125 ]
        for search_i=1:length(lon_obs)
            for search_j=1:length(lon_obs)
                if(search_i<length(lon_obs) && search_j<length(lon_obs) && search_i ~= search_j && lon_obs(search_i)==lon_obs(search_j))%%duplication point
                    if(lat_obs(search_i)==lat_obs(search_j))
                        disp(['duplication point : ', num2str(search_j)])   
                        lon_obs(search_j) = [];
                        lat_obs(search_j) = [];
                        sst_obs(search_j) = [];
                        obserr(search_j) = [];
                    end
                end
            end
        end
        for search_i=1:length(lon_obs)
            for search_j=1:length(lon_obs)
                if(search_i<length(lon_obs) && search_j<length(lon_obs) && search_i ~= search_j && lon_obs(search_i)==lon_obs(search_j))%%duplication point
                    if(lat_obs(search_i)==lat_obs(search_j))
                        disp(['duplication point : ', num2str(search_j)])   
                        lon_obs(search_j) = [];
                        lat_obs(search_j) = [];
                        sst_obs(search_j) = [];
                        obserr(search_j) = [];
                    end
                end
            end
        end
        for wrong_i=1:size(wrong_point,1)
            new_ind=find(lon_obs==wrong_point(wrong_i,1) & lat_obs==wrong_point(wrong_i,2));
            if exist('nn')==0
                nn=new_ind;
            else
                nn=[nn new_ind];
            end
        end
        lon_obs(nn) = [];
        lat_obs(nn) = [];
        sst_obs(nn) = [];
        obserr(nn) = [];
%     end
    % nn =[87]; %% near taiwan point index, first wrong point(lowest latitude)
    % lon_obs(nn) = [];
    % lat_obs(nn) = [];
    % sst_obs(nn) = [];


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
    outname = strcat('/data2/kimyy/Observation/OISST/','nwp_1_10_obs_',num2str(aa,'%5.5d'),'.nc')
    if (exist(outname,'file'))
        system(['rm -f ',outname]);
    end
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

%     emptyval=1;
end
% 

gridname=strcat('/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_10/input/test06/spinup/roms_grid_nwp_1_10_test06.nc');
% filename = strcat(filepath,filename);
lat_rho = ncread(gridname,'lat_rho');
lon_rho = ncread(gridname,'lon_rho');
mask_rho = ncread(gridname, 'mask_rho');

% % 
% figure
% % m_proj('mercator','lon',[115 162],'lat',[15 52]);
% % m_grid('fontsize',15)
% hold on
% m_pcolor(lon_rho,lat_rho,mask_rho);
% shading flat;
% m_plot(lon_obs,lat_obs,'ro','Markersize',5)
% % m_gshhs_i('color','k')  
% 
% % sti = 397;
% % m_plot(lon_obs(sti),lat_obs(sti),'bo','Markersize',5)

% 
figure
% m_proj('mercator','lon',[115 162],'lat',[15 52]);
% m_grid('fontsize',15)
hold on
pcolor(lon_rho,lat_rho,mask_rho);
shading flat;
plot(lon_obs,lat_obs,'ro','Markersize',5)
% m_gshhs_i('color','k')  

% sti = 397;
% m_plot(lon_obs(sti),lat_obs(sti),'bo','Markersize',5)
