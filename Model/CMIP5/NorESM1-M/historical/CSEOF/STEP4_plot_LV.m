clc;clear all; close all;
load('/data1/kimyy/etc/NorESM1-M_CSEOF/STEP2_5_output.mat')

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
rehash toolboxcache
%load loading vector file
load([workdir, '/cseofs/lv_layer_merge.mat']);
% lon = lon_atm;
% lat = lat_atm;
mode = 3;
%load cseof pc time series
pct_data = importdata([workdir, '/cseofs/cpct_NorESM1-M_tas.d'])';
pct = pct_data(:)';
pct = reshape(pct,tlen,length(pct)/tlen);
%pct(:,2) = -pct(:,2);
%plot pc time series

for nyear = 1:length(inputyear)
    tempyear = inputyear(nyear);
    for month=1:12
        xData((12*(nyear-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
    end
end

for i = 1:mode
    plot(xData, pct(:,i));
    tt = strcat('TAS -',num2str(i),'mode');
    title(tt,'fontsize',18,'fontweight','bold');  
    datetick('x', 'yy', 'keepticks')
    xlabel('Time(year)','FontSize',20) ;
    set(gca,'FontSize',14);
    saveas(gcf, strcat([workdir, '/figure/',tt,'.png']),'png');
    close all
    i
end

%LV = LV(layer,lon,lat,period,mode)
LV_increment1(LV_increment1==0) = NaN;
%LV_thetao(:,:,:,:,2) = -LV_thetao(:,:,:,:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%% surface plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 hh = 1;
%plot - surface
[lon2, lat2] = meshgrid(lon,lat);
%[lat2, lon2] = meshgrid(lat,lon);

%sst
for i = 1 : mode
    for period = 1:T
    %dat3 = squeeze(LV_increment1(hh,:,:,period,i));
    dat3 = squeeze(LV_increment1(hh,:,:,period,i));
    dat3(dat3==1e20) = NaN;
    m_proj('mercator','lon',[115 164],'lat',[15 52]);
    m_gshhs_l('color','k');
    m_gshhs_l('patch',[.8,.8,.8]);
    m_grid('box','fancy','tickdir','in','linewidth',1);
    hold on;
    %[C,h] = m_contour(lon2,lat2,dat3,'LineWidth',2);
    m_pcolor(lon2,lat2,dat3');
    %clabel(C,h,'FontSize',14,'Color','k','Rotation',0,'fontweight','bold');
    shading interp;
    colorbar;
    colormap jet;
    xlabel(['Longitude (^o E)'],'fontsize',18,'fontweight','bold','fontname','times new roman');
    ylabel(['Latitude(^o N)'],'fontsize',18,'fontweight','bold','fontname','times new roman');
    tt = strcat('SST -',num2str(i),'mode','(period-',num2str(period),')');
    title(tt,'fontsize',18,'fontweight','bold');   %��ٲٱ�
    caxis([-4 4]);
    axis tight;
    saveas(gcf,strcat([workdir, '/figure/surface/',tt,'.png']),'png');
    close all;
    [i period]
    end
end

